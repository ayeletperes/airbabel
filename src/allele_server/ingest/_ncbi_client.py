"""NCBI Entrez (E-utilities) client for the NCBI xref enrichment pass.

Stdlib ``urllib`` only (no new runtime dep; ingest-only). Rails:

* **Disk cache first** — every esearch/efetch body cached under
  ``cache_dir/<sha256(query)>.json``; re-runs are a no-op, bit-for-bit.
* **Graceful** — on network/parse failure, serve the cached body with
  ``stale=True``; never raise for transient errors.
* **Polite** — rate limited (2/s, 8/s with an API key), identifying User-Agent.

BLAST is intentionally omitted (out of scope here). Resolution is name-based
(``esearch[Title]`` -> ``efetch`` -> substring-confirm) plus a direct
``efetch(accession)`` path for alleles that already carry a GenBank accession.
"""

from __future__ import annotations

import hashlib
import json
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal

ENTREZ_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
USER_AGENT = "allele-server/0.1 (contact: aperes7@gmail.com)"

MatchConfidence = Literal["exact_substring", "name_only", "none"]

_KMER = 30
_KMER_THRESHOLD = 0.50


def _kmer_overlap(query: str, record: str, k: int = _KMER) -> float:
    if not query or len(query) < k or not record:
        return 0.0
    qset = {query[i:i + k] for i in range(len(query) - k + 1)}
    if not qset:
        return 0.0
    rset = {record[i:i + k] for i in range(len(record) - k + 1)}
    return len(qset & rset) / len(qset)


@dataclass(frozen=True)
class GBRecord:
    accession: str | None
    sequence: str | None
    definition: str | None
    organism: str | None
    gene_id: str | None
    pubmed_ids: tuple[str, ...] = ()


@dataclass
class ResolveResult:
    accession: str | None
    url: str | None
    gene_id: str | None
    definition: str | None
    confidence: MatchConfidence
    resolved_from_name: str
    stale: bool
    error: str | None
    candidates: list[str] = field(default_factory=list)
    # PMIDs the GenBank record cites (its reference elements). High precision:
    # emitted as *verified* references only when the sequence is confirmed.
    linked_pubmed_ids: list[str] = field(default_factory=list)


def nuccore_url(accession: str) -> str:
    return f"https://www.ncbi.nlm.nih.gov/nuccore/{accession}"


def gene_url(gene_id: str) -> str:
    return f"https://www.ncbi.nlm.nih.gov/gene/{gene_id}"


@dataclass
class NCBIClient:
    cache_dir: Path
    api_key: str | None = None
    rate_limit_per_s: float | None = None
    timeout_s: float = 15.0
    # Injectable HTTP for tests: (url, params) -> body text. Default uses urllib.
    fetcher: Callable[[str, dict[str, str]], str] | None = None
    _last_call_at: float = 0.0

    def __post_init__(self) -> None:
        self.cache_dir = Path(self.cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        if self.rate_limit_per_s is None:
            self.rate_limit_per_s = 8.0 if self.api_key else 2.0

    # --- cache + throttle ---------------------------------------------------

    def _cache_path(self, key: str) -> Path:
        return self.cache_dir / f"{hashlib.sha256(key.encode()).hexdigest()}.json"

    def _cache_read(self, key: str) -> str | None:
        path = self._cache_path(key)
        if not path.exists():
            return None
        try:
            return json.loads(path.read_text(encoding="utf-8")).get("body")
        except (OSError, json.JSONDecodeError):
            return None

    def _cache_write(self, key: str, body: str) -> None:
        self._cache_path(key).write_text(json.dumps({"body": body}), encoding="utf-8")

    def _cache_drop(self, key: str) -> None:
        self._cache_path(key).unlink(missing_ok=True)

    def _throttle(self) -> None:
        gap = 1.0 / (self.rate_limit_per_s or 2.0)
        wait = gap - (time.monotonic() - self._last_call_at)
        if wait > 0:
            time.sleep(wait)
        self._last_call_at = time.monotonic()

    def _get(self, endpoint: str, params: dict[str, str]) -> str:
        if self.api_key:
            params = {**params, "api_key": self.api_key}
        url = f"{ENTREZ_BASE}/{endpoint}"
        if self.fetcher is not None:
            return self.fetcher(url, params)
        self._throttle()
        full = f"{url}?{urllib.parse.urlencode(params)}"
        req = urllib.request.Request(full, headers={"User-Agent": USER_AGENT})
        with urllib.request.urlopen(req, timeout=self.timeout_s) as resp:
            return resp.read().decode("utf-8")

    # --- E-utilities --------------------------------------------------------

    def esearch(self, term: str) -> tuple[list[str], bool, str | None]:
        """Return (id_list, stale, error) for a nuccore title search."""
        key = f"esearch::{term}"
        params = {"db": "nuccore", "term": term, "retmode": "json", "retmax": "5"}
        try:
            body = self._get("esearch.fcgi", params)
            ids = list(json.loads(body).get("esearchresult", {}).get("idlist", []))
            self._cache_write(key, body)  # cache only a parseable response
            return ids, False, None
        except Exception as exc:  # noqa: BLE001 - graceful degrade
            cached = self._cache_read(key)
            if cached is not None:
                # A corrupt cache entry must not abort the pass: parse defensively and,
                # if it is garbage, drop it so a later run refetches instead of re-failing.
                # (This once crashed a full reload; efetch already guarded its cached parse.)
                try:
                    ids = list(json.loads(cached).get("esearchresult", {}).get("idlist", []))
                    return ids, True, f"{exc} (cached)"
                except (json.JSONDecodeError, ValueError, AttributeError):
                    self._cache_drop(key)
            return [], False, str(exc)

    def efetch(self, ident: str) -> tuple[GBRecord | None, bool, str | None]:
        """Fetch one GenBank record by uid or accession; (record, stale, error)."""
        key = f"efetch::{ident}"
        params = {"db": "nuccore", "id": ident, "rettype": "gb", "retmode": "xml"}
        try:
            body = self._get("efetch.fcgi", params)
            self._cache_write(key, body)
            return _parse_gb_xml(body), False, None
        except Exception as exc:  # noqa: BLE001
            cached = self._cache_read(key)
            if cached is not None:
                try:
                    return _parse_gb_xml(cached), True, f"{exc} (cached)"
                except ET.ParseError:
                    pass
            return None, False, str(exc)

    # --- high-level resolution ---------------------------------------------

    def confirm_accession(self, accession: str, query_sequence: str) -> ResolveResult:
        """Confirm a known accession via efetch and grade the match."""
        rec, stale, err = self.efetch(accession)
        return self._grade(rec, accession, query_sequence, stale, err)

    def resolve(self, name: str, query_sequence: str) -> ResolveResult:
        """Resolve an allele *name* to a GenBank record by title search."""
        ids, stale, err = self.esearch(f"{name}[Title]")
        best: ResolveResult | None = None
        for uid in ids:
            rec, s2, e2 = self.efetch(uid)
            graded = self._grade(rec, name, query_sequence, stale or s2, err or e2)
            if graded.confidence == "exact_substring":
                graded.candidates = ids
                return graded
            if best is None or graded.confidence == "name_only":
                best = graded
        if best is None:
            return ResolveResult(None, None, None, None, "none", name, stale, err, ids)
        best.candidates = ids
        return best

    def _grade(
        self,
        rec: GBRecord | None,
        resolved_from: str,
        query_sequence: str,
        stale: bool,
        err: str | None,
    ) -> ResolveResult:
        if rec is None or rec.accession is None:
            return ResolveResult(None, None, None, None, "none", resolved_from, stale, err)
        qn = query_sequence.upper().replace("-", "").replace(".", "")
        rn = (rec.sequence or "").upper()
        if rn and qn and qn in rn:
            conf: MatchConfidence = "exact_substring"
        elif rn and _kmer_overlap(qn, rn) >= _KMER_THRESHOLD:
            conf = "name_only"
        else:
            conf = "name_only" if rn else "none"
        return ResolveResult(
            accession=rec.accession,
            url=nuccore_url(rec.accession),
            gene_id=rec.gene_id,
            definition=rec.definition,
            confidence=conf,
            resolved_from_name=resolved_from,
            stale=stale,
            error=err,
            linked_pubmed_ids=list(rec.pubmed_ids),
        )


def _text(parent: ET.Element, tag: str) -> str | None:
    el = parent.find(tag)
    return el.text if el is not None and el.text else None


def _parse_gb_xml(xml_text: str) -> GBRecord:
    """Parse the first GBSeq from an efetch GB-XML body into a GBRecord."""
    root = ET.fromstring(xml_text)
    gbseq = root.find(".//GBSeq")
    if gbseq is None:
        return GBRecord(None, None, None, None, None)
    accession = _text(gbseq, "GBSeq_accession-version") or _text(
        gbseq, "GBSeq_primary-accession"
    )
    sequence = _text(gbseq, "GBSeq_sequence")
    definition = _text(gbseq, "GBSeq_definition")
    organism = _text(gbseq, "GBSeq_organism")
    gene_id: str | None = None
    for qual in gbseq.findall(".//GBQualifier"):
        if _text(qual, "GBQualifier_name") == "db_xref":
            val = _text(qual, "GBQualifier_value") or ""
            if val.startswith("GeneID:"):
                gene_id = val.split(":", 1)[1]
                break
    # PMIDs the record cites, from its <GBReference_pubmed> reference elements.
    # Deduped, order-preserving — these are the high-precision literature anchor.
    pubmed_ids: list[str] = []
    for ref in gbseq.findall(".//GBReference_pubmed"):
        pmid = (ref.text or "").strip()
        if pmid and pmid not in pubmed_ids:
            pubmed_ids.append(pmid)
    return GBRecord(accession, sequence, definition, organism, gene_id, tuple(pubmed_ids))
