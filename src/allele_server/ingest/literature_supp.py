"""Full-text & supplementary mining — auto-collect allele evidence from PMC OA.

**Status: v1 / experimental.** This is milestone M3 of the literature-mining
plan (``docs/literature-mining-plan.md``): for a candidate article (PMID) we
fetch its **PMC Open-Access supplementary files** *and its article body
full-text XML*, parse them, extract candidate allele **names** and **nucleotide
sequences**, and match those back to the database by name (via
:mod:`allele_server.naming`) and by sequence (exact ``hash_sequence`` + the
existing similarity path for near-misses). Each hit is stored as a ``literature``
PMID annotation so it surfaces through ``queries._references`` alongside the v1
Europe PMC name probe.

**Sequence over name.** A sequence hash-match is the real, cross-species-safe
evidence ("this paper contains this exact allele"); a bare name mention is only a
lead, because the same name string can denote a different allele in another
species. The strongest signal per paper wins (see ``_MATCH_RANK``), and the UI
tiers references accordingly (verified curator > sequence-in-paper > name).

**Honesty (UNVERIFIED):** supplementary tables are exactly where novel-allele
evidence hides, but they are also the noisiest source — a scraped nucleotide run
may be a primer/probe/OCR artefact and an allele name in a cell may be a clone
id or an unrelated table label. Nothing collected here is manually checked: every
hit is written under ``domain == "literature"`` and therefore rendered as
**unverified** (``queries._references`` treats ``domain != "literature"`` as the
only verified provenance). Treat these entries as *leads*, not ground truth. A
companion ``suppmatch:`` annotation records the match type + source file so a
human can adjudicate.

**Scope / coverage.** Only **PMC-OA** (open access, redistributable for
text-mining) is fetched — paywalled publisher HTML/PDF is never touched. Parsing
covers xlsx (``openpyxl``), csv/tsv (stdlib ``csv``) and FASTA/txt (a small
parser); **PDF and docx are deferred** (see ``TODO(pdf-docx)``) — pdfplumber /
python-docx / GROBID are intentionally not added in this phase.

Same rails as the other enrichment passes: disk cache first (idempotent, offline
re-runs), best-effort networking (parse/network failures return empty, never
raise), and every DB write idempotent on the annotation natural key.
"""

from __future__ import annotations

import hashlib
import io
import json
import logging
import os
import re
import threading
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
import zipfile
from collections.abc import Callable, Iterable, Iterator
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path

from sqlalchemy import or_, select
from sqlalchemy.orm import Session

from .. import queries
from ..curie import build as curie_build
from ..db.models import Name, Sequence, SequenceAnnotation
from ..naming import _NAME_RE, AlleleNameError, canonicalise, locus_of, segment_of
from ..sequence import SequenceError, hash_sequence
from . import loader
from .openalex import OpenAlexClient

_LOG = logging.getLogger(__name__)

_EPMC_BASE = "https://www.ebi.ac.uk/europepmc/webservices/rest"
_USER_AGENT = "allele-server/0.1 (contact: aperes7@gmail.com)"

# A candidate sequence must be a contiguous ACGT[N] run of at least this many
# bases — long enough that primer/adapter fragments and short table codes are
# dropped before they ever reach the (expensive) match step.
_MIN_SEQ_LEN = 60
_NT_RUN = re.compile(rf"[ACGTN]{{{_MIN_SEQ_LEN},}}")
# Split a free-text cell / line into whitespace- and punctuation-delimited tokens,
# keeping the characters an allele grammar can contain (letters, digits, / * _ - .).
_TOKEN_SPLIT = re.compile(r"[^A-Za-z0-9/*_.\-]+")

# For a sequence near-miss we demand a *high* identity: literature evidence should
# read as "this paper reports (nearly) this allele", not the generous 0.75 floor
# the neighbour cache uses for exploratory "find similar".
SUPP_SIMILAR_FLOOR = 0.95

# The exact-hash lookup is an indexed point read and runs for every candidate. The
# *near-miss* fallback is a trigram + rapidfuzz scan of the whole table (~1s each), and
# a supplementary FASTA can hold thousands of runs: unbounded, it made the crawl
# database-bound and effectively non-terminating. So the fallback is gated on a
# plausible allele length and rationed per article, and skips are counted, never silent.
_SIMILAR_MIN_LEN = 150
_SIMILAR_MAX_LEN = 800
_SIMILAR_PER_ARTICLE = 20

# Match precedence when a single sequence is hit several ways in one paper: an
# exact sequence beats a near sequence beats a bare name.
_MATCH_RANK = {"sequence_exact": 3, "sequence_similar": 2, "name": 1}

# Supplementary bundles are occasionally raw-data dumps (one observed at 323 MB) that
# contain no allele table and stall the crawl. Skip anything past the cap: an allele
# table is kilobytes to a few megabytes. The deadline is the backstop for a server
# that trickles bytes forever without ever tripping the socket-inactivity timeout.
_MAX_SUPP_BYTES = 64 * 1024 * 1024
_MAX_FULLTEXT_BYTES = 32 * 1024 * 1024
_FETCH_DEADLINE = 180.0

# Total disk the bundle cache may occupy. The cache exists so an interrupted crawl
# resumes cheaply, which is worth a lot for a handful of articles and nothing at
# all if it fills the disk: an article-centric sweep hands this client ~12,600
# PMIDs, and a single observed supplementary bundle is 49 MB. Past the budget the
# client keeps fetching and mining exactly as before and simply stops *retaining*
# what it fetched, so correctness never depends on how much disk is left.
_DEFAULT_CACHE_BUDGET = 20 * 1024 * 1024 * 1024


# --- fetch client ---------------------------------------------------------

@dataclass(frozen=True)
class SuppFile:
    """One supplementary file: its (base) name and raw bytes."""

    name: str
    data: bytes


class PMCSuppClient:
    """Fetch PMC Open-Access supplementary files for a PMID, with a disk cache.

    ``supplementary_files(pmid)`` resolves the PMID to a PMCID (Europe PMC
    search), downloads the article's supplementary-files zip from the Europe PMC
    OA web service, unzips it in memory and returns one :class:`SuppFile` per
    member. The raw zip is cached under ``cache_dir/<pmcid>.zip`` and the
    PMID→PMCID / title lookups under ``cache_dir/<sha256>.json``; a cache hit
    never touches the network. Every network/parse failure degrades to an empty
    result — a missing supplement is not an error (best-effort, like the v1
    Europe PMC name probe).

    ``fetcher`` is an injectable ``url -> bytes`` seam for tests (no network);
    the default uses stdlib ``urllib`` with polite rate limiting.
    """

    def __init__(
        self,
        cache_dir: str | Path,
        *,
        fetcher: Callable[[str], bytes] | None = None,
        min_interval: float = 0.34,
        timeout: int = 60,
        max_cache_bytes: int | None = _DEFAULT_CACHE_BUDGET,
        oa_client: OpenAlexClient | None = None,
    ) -> None:
        self.cache = Path(cache_dir)
        self.cache.mkdir(parents=True, exist_ok=True)
        self.fetcher = fetcher
        # Second opinion on articles Europe PMC has not indexed: it reports no
        # PMCID for records whose PMC copy is live (PMID 42251052 -> PMC13396459).
        self.oa_client = oa_client
        # Counted per *article*, not per call: resolve_pmcid runs up to three times
        # for one PMID (the loop, the bundle fetch, the full text).
        self.pmcid_via_openalex = 0
        self._oa_counted: set[str] = set()
        self.min_interval = min_interval
        self.timeout = timeout
        self.max_cache_bytes = max_cache_bytes
        self._last = 0.0
        self._cache_bytes = self._measure_cache() if max_cache_bytes else 0
        self._budget_reported = False
        # Prefetch runs several fetches at once (see `prefetch`), so the pacer and
        # the cache accounting are shared mutable state. One lock for each: they are
        # held for different durations and must not serialise one another.
        self._pace_lock = threading.Lock()
        self._cache_lock = threading.Lock()

    def _measure_cache(self) -> int:
        return sum(p.stat().st_size for p in self.cache.glob("*") if p.is_file())

    def _wait_turn(self) -> None:
        """Block until this thread may **start** a request, then claim the slot.

        The politeness interval used to be measured from the *end* of one request to
        the start of the next, which serialises by construction: with a 30 s download
        the client made two requests a minute however many threads were waiting.
        Pacing request starts instead keeps the same load on Europe PMC -- at most
        one new request per ``min_interval`` -- while letting the waiting overlap,
        which is the whole cost of this crawl.

        The slot is stamped while the lock is held, so two threads cannot both decide
        the interval has elapsed and leave together.
        """
        with self._pace_lock:
            gap = self.min_interval - (time.monotonic() - self._last)
            if gap > 0:
                time.sleep(gap)
            self._last = time.monotonic()

    def warm(self, pmid: str) -> None:
        """Pull this article's payloads into the disk cache. No parsing, no return.

        Prefetch only: the mining loop calls ``supplementary_files`` /
        ``full_text_fragments`` as before and finds everything local. Splitting it
        this way is what makes concurrency safe -- the parallel part touches only the
        network and the cache, never the database or the match state, so there is no
        session to share and no ordering to preserve.
        """
        pmcid = self.resolve_pmcid(pmid)
        if not pmcid:
            return
        self._zip_bytes(pmcid)
        cp = self.cache / f"{pmcid}.fulltext.xml"
        if not cp.exists():
            raw = self._get(
                f"{_EPMC_BASE}/{pmcid}/fullTextXML", max_bytes=_MAX_FULLTEXT_BYTES
            )
            if raw is not None:
                self._retain(cp, raw)

    def prefetch(self, pmids: Iterable[str], workers: int = 6) -> int:
        """Warm the cache for ``pmids`` concurrently. Returns how many were attempted.

        Best-effort throughout: a worker that raises leaves that article uncached and
        the mining loop fetches it itself, exactly as it would have done anyway.
        """
        def _one(pmid: str) -> None:
            try:
                self.warm(pmid)
            except Exception:  # noqa: BLE001 - prefetch is an optimisation
                pass

        todo = [p for p in pmids if p]
        if todo:
            with ThreadPoolExecutor(max_workers=max(1, workers)) as pool:
                list(pool.map(_one, todo))
        return len(todo)

    # -- resume ledger -------------------------------------------------------
    #
    # A corpus sweep hands this client ~12,600 articles and mining them takes
    # hours, so a run that dies partway must not start over. The bundle cache
    # makes the *fetch* cheap on a second pass but not the parse and the
    # per-candidate database lookups, and an article with no matches leaves no
    # trace in the database to recognise it by. So the ledger is a plain list of
    # PMIDs already mined, appended as each one finishes committing.
    #
    # It records "mined against the alleles we held at the time", which is why
    # nothing consults it unless the caller asks: after an ingest adds alleles,
    # the same articles are worth re-reading, and deleting the file is how you
    # say so.

    def mined(self) -> set[str]:
        try:
            return set((self.cache / "mined.txt").read_text().split())
        except OSError:
            return set()

    def mark_mined(self, pmid: str) -> None:
        with open(self.cache / "mined.txt", "a") as fh:
            fh.write(f"{pmid}\n")

    def _retain(self, path: Path, data: bytes) -> None:
        """Write ``data`` to the cache unless that would break the disk budget.

        Skipping a write costs a re-download on a later run, which is a bounded,
        recoverable cost. Filling the disk during an overnight crawl is neither.
        Small writes (the PMID->PMCID JSONs, and the empty markers that record a
        deliberate skip) always go through: they are what makes a resumed run
        cheap, and they are bytes rather than megabytes.
        """
        if self.max_cache_bytes and data:
            with self._cache_lock:  # concurrent prefetch workers share the running total
                if self._cache_bytes + len(data) > self.max_cache_bytes:
                    if not self._budget_reported:
                        _LOG.warning(
                            "bundle cache at its %.1f GB budget; continuing to fetch "
                            "and mine, but no longer retaining bundles",
                            self.max_cache_bytes / 1024 ** 3,
                        )
                        self._budget_reported = True
                    return
                self._cache_bytes += len(data)
        path.write_bytes(data)

    # -- low-level fetch + cache --------------------------------------------

    def _json_cache(self, key: str) -> Path:
        return self.cache / f"{hashlib.sha256(key.encode()).hexdigest()}.json"

    def _get(
        self, url: str, *, max_bytes: int | None = None, deadline: float = _FETCH_DEADLINE
    ) -> bytes | None:
        """Fetch ``url`` to bytes. ``None`` on failure, ``b""`` when deliberately skipped.

        ``timeout`` only bounds socket *inactivity*: a server that trickles bytes can
        hold the connection open indefinitely (this hung a real crawl for nine hours).
        So we also enforce a total ``deadline`` and a ``max_bytes`` cap, checking the
        advertised Content-Length first and then the running total while streaming.
        A skip returns ``b""`` (a permanent "nothing here", safe to cache) as opposed to
        ``None`` (a transient failure, worth retrying on a later run).
        """
        if self.fetcher is not None:
            try:
                return self.fetcher(url)
            except Exception:  # noqa: BLE001 - injected seam, stay best-effort
                return None
        self._wait_turn()
        req = urllib.request.Request(url, headers={"User-Agent": _USER_AGENT})
        try:
            with urllib.request.urlopen(req, timeout=self.timeout) as resp:
                declared = resp.headers.get("Content-Length")
                if max_bytes and declared and declared.isdigit() and int(declared) > max_bytes:
                    _LOG.info("skipping %s: %s bytes over the %s cap", url, declared, max_bytes)
                    return b""
                started = time.monotonic()
                chunks: list[bytes] = []
                total = 0
                while True:
                    if time.monotonic() - started > deadline:
                        _LOG.info("skipping %s: exceeded %.0fs deadline", url, deadline)
                        return b""
                    chunk = resp.read(1 << 20)
                    if not chunk:
                        break
                    total += len(chunk)
                    if max_bytes and total > max_bytes:
                        _LOG.info("skipping %s: streamed past the %s cap", url, max_bytes)
                        return b""
                    chunks.append(chunk)
                return b"".join(chunks)
        except Exception:  # noqa: BLE001 - network failure -> no supplement
            return None

    def _search_result(self, pmid: str) -> dict | None:
        """First Europe PMC result for a PMID (cached), or None."""
        key = f"epmc-search::{pmid}"
        cp = self._json_cache(key)
        if cp.exists():
            try:
                return json.loads(cp.read_text()) or None
            except json.JSONDecodeError:
                pass
        query = urllib.parse.quote(f"ext_id:{pmid} AND SRC:MED")
        url = f"{_EPMC_BASE}/search?query={query}&format=json&resultType=lite"
        body = self._get(url)
        result: dict | None = None
        if body is not None:
            try:
                hits = json.loads(body.decode()).get("resultList", {}).get("result", [])
                result = hits[0] if hits else None
            except (json.JSONDecodeError, UnicodeDecodeError, IndexError):
                result = None
        cp.write_text(json.dumps(result or {}))
        return result

    def resolve_pmcid(self, pmid: str) -> str | None:
        """Resolve a PMID to its PMCID (e.g. ``PMC1234567``), or None.

        Europe PMC first. When it has no PMCID and an ``oa_client`` was supplied,
        ask OpenAlex, which carries the article's PMC *location* — the PMCID is
        read out of the URL OpenAlex returned, never constructed here. That is the
        only way an article Europe PMC has not indexed becomes fetchable at all.
        """
        res = self._search_result(pmid)
        pmcid = (res or {}).get("pmcid") or None
        if not pmcid and self.oa_client is not None:
            pmcid = self.oa_client.pmcid(pmid)
            if pmcid and pmid not in self._oa_counted:
                self._oa_counted.add(pmid)
                self.pmcid_via_openalex += 1
        return str(pmcid) if pmcid else None

    def article_title(self, pmid: str) -> str | None:
        """The article title for a PMID from the cached search, or None."""
        res = self._search_result(pmid)
        title = (res or {}).get("title")
        return (title or "").strip() or None if title else None

    def _zip_bytes(self, pmcid: str) -> bytes | None:
        """The supplementary-files zip for a PMCID (cached), or None.

        An oversized/stalled bundle comes back as ``b""`` and is cached as such, so the
        skip is permanent and the next run does not re-attempt the same 300 MB download.
        A transient failure (``None``) is not cached and will be retried."""
        cp = self.cache / f"{pmcid}.zip"
        if cp.exists():
            return cp.read_bytes() or None
        url = f"{_EPMC_BASE}/{pmcid}/supplementaryFiles"
        raw = self._get(url, max_bytes=_MAX_SUPP_BYTES)
        if raw is None:
            return None
        self._retain(cp, raw)  # b"" records a deliberate, permanent skip
        return raw or None

    def supplementary_files(self, pmid: str) -> list[SuppFile]:
        """PMC-OA supplementary files for a PMID; ``[]`` if none / not OA."""
        pmcid = self.resolve_pmcid(pmid)
        if not pmcid:
            return []
        raw = self._zip_bytes(pmcid)
        if not raw:
            return []
        try:
            zf = zipfile.ZipFile(io.BytesIO(raw))
        except (zipfile.BadZipFile, OSError):
            return []
        out: list[SuppFile] = []
        with zf:
            for name in zf.namelist():
                if name.endswith("/"):
                    continue
                try:
                    out.append(SuppFile(os.path.basename(name) or name, zf.read(name)))
                except (KeyError, zipfile.BadZipFile, OSError):
                    continue
        return out

    def full_text_fragments(self, pmid: str) -> list[str]:
        """Text blocks of the article **body** (PMC-OA full-text XML), or ``[]``.

        This is the "search inside the text" layer: the paper's own methods,
        results and inline tables — where authors list the alleles they observed
        and, sometimes, print the sequences. Only the OA subset has minable full
        text; a non-OA article (no ``fullTextXML``) caches an empty miss so it is
        never refetched. Returns one string per block element (see ``_jats_fragments``)."""
        pmcid = self.resolve_pmcid(pmid)
        if not pmcid:
            return []
        cp = self.cache / f"{pmcid}.fulltext.xml"
        if cp.exists():
            raw = cp.read_bytes()
        else:
            raw = self._get(f"{_EPMC_BASE}/{pmcid}/fullTextXML", max_bytes=_MAX_FULLTEXT_BYTES)
            if raw is None:
                return []  # transient failure: do not cache, retry on a later run
            self._retain(cp, raw)  # b"" caches a real miss (not OA / no full text)
        return _jats_fragments(raw) if raw else []


# JATS block elements whose text we mine. Leaf-ish blocks (not ``sec``/``body``)
# so a sequence or name is not double-counted through its containers.
_BLOCK_TAGS = frozenset(
    {"p", "td", "th", "title", "caption", "label", "code", "list-item",
     "def", "term", "monospace", "preformat"}
)


def _jats_fragments(raw: bytes) -> list[str]:
    """Split PMC full-text JATS XML into per-block text strings (best-effort).

    Each block's text nodes are joined with spaces (so allele **names** stay
    whitespace-delimited for ``_candidate_names``); the sequence pass despaces
    each block itself, so a nucleotide run split across inline markup rejoins."""
    try:
        root = ET.fromstring(raw)
    except ET.ParseError:
        return []
    frags: list[str] = []
    for el in root.iter():
        tag = el.tag.rsplit("}", 1)[-1].lower()  # drop any XML namespace
        if tag not in _BLOCK_TAGS:
            continue
        txt = " ".join(t.strip() for t in el.itertext() if t and t.strip())
        if txt:
            frags.append(txt)
    return frags


# --- parsers (best-effort: return [] on failure) --------------------------

def _parse_xlsx(data: bytes) -> list[str]:
    """Every non-empty cell of every worksheet, as a string (openpyxl)."""
    try:
        from openpyxl import load_workbook
    except ImportError:  # pragma: no cover - guarded by the `ingest` extra
        _LOG.warning("openpyxl not installed; skipping xlsx (install the ingest extra)")
        return []
    try:
        wb = load_workbook(io.BytesIO(data), read_only=True, data_only=True)
    except Exception:  # noqa: BLE001 - corrupt / unsupported workbook
        return []
    out: list[str] = []
    try:
        for ws in wb.worksheets:
            for row in ws.iter_rows(values_only=True):
                out.extend(str(c) for c in row if c is not None and str(c).strip())
    finally:
        wb.close()
    return out


def _parse_delimited(data: bytes, delimiter: str) -> list[str]:
    """Every non-empty field of a csv/tsv, as a string (stdlib csv)."""
    import csv

    try:
        text = data.decode("utf-8", "replace")
    except Exception:  # noqa: BLE001
        return []
    out: list[str] = []
    try:
        for row in csv.reader(io.StringIO(text), delimiter=delimiter):
            out.extend(c for c in row if c and c.strip())
    except csv.Error:
        return []
    return out


def _parse_fasta(data: bytes) -> list[str]:
    """FASTA as fragments: each header string + each record's joined sequence."""
    text = data.decode("utf-8", "replace")
    out: list[str] = []
    seq: list[str] = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if seq:
                out.append("".join(seq))
                seq = []
            out.append(line[1:])
        else:
            seq.append("".join(line.split()))
    if seq:
        out.append("".join(seq))
    return out


def _parse_text(data: bytes) -> list[str]:
    """Plain text as one fragment per non-empty line."""
    text = data.decode("utf-8", "replace")
    return [ln.strip() for ln in text.splitlines() if ln.strip()]


# Extension -> parser. PDF / doc(x) / legacy .xls are deferred; see TODO(pdf-docx).
_PARSERS: dict[str, Callable[[bytes], list[str]]] = {
    ".xlsx": _parse_xlsx,
    ".xlsm": _parse_xlsx,
    ".csv": lambda d: _parse_delimited(d, ","),
    ".tsv": lambda d: _parse_delimited(d, "\t"),
    ".tab": lambda d: _parse_delimited(d, "\t"),
    ".fasta": _parse_fasta,
    ".fa": _parse_fasta,
    ".fna": _parse_fasta,
    ".faa": _parse_fasta,
    ".txt": _parse_text,
    ".seq": _parse_text,
}
# TODO(pdf-docx): PDF (pdfplumber/GROBID) and docx (python-docx) parsing is
# deferred to a later phase — they need heavy deps not added here. Logged + skipped.
_DEFERRED = {".pdf", ".doc", ".docx", ".xls", ".ppt", ".pptx"}


def _fragments(file: SuppFile) -> list[str]:
    """Text fragments from a supplementary file, dispatched on its extension."""
    ext = os.path.splitext(file.name.lower())[1]
    parser = _PARSERS.get(ext)
    if parser is None:
        if ext in _DEFERRED:
            _LOG.info("deferred file type %s (%s) — not mined in v1", ext, file.name)
        return []
    try:
        return parser(file.data)
    except Exception:  # noqa: BLE001 - any parser blowup is soft-ignored
        _LOG.debug("parse failed for %s", file.name, exc_info=True)
        return []


# --- candidate extraction -------------------------------------------------

def _candidate_names(fragment: str) -> Iterator[str]:
    """Tokens in a fragment that are IG/TR **V** allele names (carry ``*NN``)."""
    for tok in _TOKEN_SPLIT.split(fragment):
        if "*" not in tok or _NAME_RE.match(tok) is None:
            continue
        if segment_of(tok) != "V" or locus_of(tok) is None:
            continue
        yield tok


def _candidate_sequences(fragment: str) -> Iterator[str]:
    """Contiguous ACGT[N] runs of >= ``_MIN_SEQ_LEN`` bases in a fragment."""
    for m in _NT_RUN.finditer(fragment.upper()):
        yield m.group(0)


# --- match back to the DB -------------------------------------------------

def _match_name(
    session: Session, token: str, species_tags: set[str] | None
) -> list[Sequence]:
    """Sequences carrying ``token`` (raw or canonicalised) as a current name."""
    cands = {token}
    try:
        cands.add(canonicalise(token).canonical())
    except AlleleNameError:
        pass
    stmt = (
        select(Sequence)
        .join(Name, Name.sequence_id == Sequence.id)
        .where(
            Name.name.in_(cands),
            Sequence.category == "coding",
            Sequence.archived.is_(False),
        )
        .distinct()
    )
    if species_tags:
        stmt = stmt.where(Sequence.species.in_(species_tags))
    return list(session.scalars(stmt))


class _SeqBudget:
    """Per-article rationing for the expensive near-miss scan (see the constants)."""

    __slots__ = ("seen", "similar_left", "skipped")

    def __init__(self, similar_left: int = _SIMILAR_PER_ARTICLE) -> None:
        self.seen: set[str] = set()
        self.similar_left = similar_left
        self.skipped = 0


def _match_sequence(
    session: Session, cand: str, species_tags: set[str] | None,
    budget: _SeqBudget | None = None,
) -> list[tuple[Sequence, str, float]]:
    """Match a candidate nt run: exact hash (incl. rc), else high-identity similar.

    Reuses the identity contract (``hash_sequence``) and the existing similarity
    path (``queries.search_by_sequence``, the same trigram + rapidfuzz scan the live
    "find similar" uses) so near-misses are ranked exactly as elsewhere.

    The exact lookup always runs. The similarity fallback only runs for a candidate of
    plausible allele length and while the article's budget lasts; otherwise it is
    skipped and counted. Sequence evidence is overwhelmingly exact anyway: a paper
    reproducing a germline reproduces it verbatim.
    """
    try:
        digest = hash_sequence(cand)
    except SequenceError:
        return []
    stmt = select(Sequence).where(
        or_(Sequence.sequence_hash == digest, Sequence.rc_hash == digest),
        Sequence.category == "coding",
        Sequence.archived.is_(False),
    )
    if species_tags:
        stmt = stmt.where(Sequence.species.in_(species_tags))
    exact = list(session.scalars(stmt))
    if exact:
        return [(s, "sequence_exact", 1.0) for s in exact]
    # No exact hit -> the existing similarity path, kept at a high floor and rationed.
    if not (_SIMILAR_MIN_LEN <= len(cand) <= _SIMILAR_MAX_LEN):
        return []
    if budget is not None:
        if budget.similar_left <= 0:
            budget.skipped += 1
            return []
        budget.similar_left -= 1
    try:
        res = queries.search_by_sequence(session, cand, cap=3)
    except ValueError:
        return []
    out: list[tuple[Sequence, str, float]] = []
    for hit in res.hits:
        if hit.score < SUPP_SIMILAR_FLOOR:
            continue
        seq = session.scalar(select(Sequence).where(Sequence.uid == hit.uid))
        if seq is None or (species_tags and seq.species not in species_tags):
            continue
        out.append((seq, "sequence_similar", hit.score))
    return out


# --- article selection ----------------------------------------------------

def _articles(
    session: Session,
    pmids: Iterable[str] | None,
    species_tags: set[str] | None,
    limit: int | None,
) -> list[str]:
    """The PMIDs to mine: the caller's list, else PMIDs already discovered by the
    name probe (``pmid:`` literature annotations), optionally scoped by species."""
    if pmids is not None:
        out = list(dict.fromkeys(str(p).strip() for p in pmids if str(p).strip()))
        return out[:limit] if limit else out
    stmt = (
        select(SequenceAnnotation.key)
        .join(Sequence, Sequence.id == SequenceAnnotation.sequence_id)
        .where(SequenceAnnotation.key.like("pmid:%"))
        .distinct()
    )
    if species_tags:
        stmt = stmt.where(Sequence.species.in_(species_tags))
    seen: list[str] = []
    for (key,) in session.execute(stmt):
        pmid = key.split(":", 1)[1].strip()
        if pmid and pmid not in seen:
            seen.append(pmid)
    return seen[:limit] if limit else seen


def _summary(match_type: str, score: float, filename: str) -> str:
    """Short, human-legible provenance string for the companion annotation."""
    if match_type == "sequence_similar":
        head = f"sequence_similar {score:.0%}"
    else:
        head = match_type
    return f"{head} · {filename}"[:250]


def enrich_literature_supp(
    session: Session,
    client: PMCSuppClient,
    *,
    pmids: Iterable[str] | None = None,
    species_keys: list[str] | None = None,
    limit: int | None = None,
    full_text: bool = True,
    resume: bool = False,
    anchor_names: bool = True,
    prefetch_workers: int = 6,
    prefetch_chunk: int = 64,
) -> dict[str, int]:
    """Mine PMC-OA supplementary files **and article body full text** for allele
    evidence (UNVERIFIED).

    For each article, fetch its supplementary files and (when ``full_text``) the
    body full-text XML, extract candidate allele names + nucleotide sequences,
    match them to the DB (name / exact hash / high-identity similar), and store
    each hit as a ``literature`` ``pmid:`` annotation (surfaced via
    ``queries._references``) plus a companion ``suppmatch:`` annotation recording
    the match type and source (a supplement filename or ``fulltext:PMCID``). When
    a sequence is hit several ways in one paper, the strongest match type wins —
    so a **sequence** match (cross-species-safe, the real evidence) always
    outranks a bare **name** mention (a lead; a same-named allele in another
    species collides on the string but not on the hash).

    ``pmids`` selects the articles; if omitted, PMIDs already discovered by the
    name probe are re-mined (optionally scoped to ``species_keys``). Idempotent
    on the annotation natural key — re-running changes zero rows. With
    ``anchor_names`` (default) a name-only match is stored only in a (species, locus)
    the article also matched by sequence (see ``_anchor``); name links written by
    older runs stay in the database and are filtered at query time. Returns
    ``{"articles", "files", "fulltext_blocks", "matches", "refs", "name",
    "sequence_exact", "sequence_similar", "no_pmcid", "no_content",
    "pmcid_via_openalex", ...}``.
    """
    source = loader.get_or_create_source(
        session, name="literature_supp", version=None, curie_prefix="PMID",
        note="auto-collected PMC-OA supplementary-file matches (unverified)",
    )
    species_tags: set[str] | None = None
    if species_keys:
        from ..species import normalize_species

        species_tags = {normalize_species(s) for s in species_keys}

    # ``no_pmcid`` / ``no_content`` are the honesty counters for the OA-resolution
    # step: an article can be open access (OpenAlex says so, with a URL) and still
    # be unreachable through the PMC-OA services this pass is allowed to fetch --
    # OA at a repository, or a PMC copy Europe PMC has not indexed yet. Those
    # articles are mined for nothing, and that has to show up as a number rather
    # than as an article that quietly contributed no matches.
    stats = {
        "articles": 0, "files": 0, "fulltext_blocks": 0, "matches": 0, "refs": 0,
        "name": 0, "sequence_exact": 0, "sequence_similar": 0, "similar_skipped": 0,
        "already_mined": 0, "name_unanchored": 0,
        "no_pmcid": 0, "no_content": 0, "pmcid_via_openalex": 0,
    }
    # ``resume`` is opt-in because the ledger means "mined against the alleles we
    # held then": after an ingest adds alleles the same articles are worth
    # re-reading, and only the caller knows which run this is.
    done = client.mined() if resume and hasattr(client, "mined") else set()
    selected = _articles(session, pmids, species_tags, limit)
    todo = [p for p in selected if p not in done]
    stats["already_mined"] = len(selected) - len(todo)
    for pmid in _prefetched(client, todo, prefetch_workers, prefetch_chunk):
        stats["articles"] += 1
        pmcid = getattr(client, "resolve_pmcid", lambda _p: None)(pmid)
        if not pmcid:
            # Nothing to fetch: no PMC copy Europe PMC or OpenAlex can point at.
            stats["no_pmcid"] += 1
        budget = _SeqBudget()
        files = client.supplementary_files(pmid)
        # best match per sequence for this article: seq.id -> (seq, type, score, source)
        best: dict[int, tuple[Sequence, str, float, str]] = {}
        for file in files:
            stats["files"] += 1
            for frag in _fragments(file):
                _mine_fragment(session, frag, best, species_tags, file.name, budget=budget)
        ft_getter = getattr(client, "full_text_fragments", None)
        blocks = 0
        if full_text and ft_getter is not None:
            src_label = f"fulltext:{pmcid or pmid}"
            for frag in ft_getter(pmid):
                blocks += 1
                stats["fulltext_blocks"] += 1
                _mine_fragment(session, frag, best, species_tags, src_label,
                               despace=True, budget=budget)
        if pmcid and not files and not blocks:
            # A PMCID that yielded neither a supplementary bundle nor body text:
            # Europe PMC's OA services do not serve this article (404), whatever
            # its licence says. Counted, not swallowed.
            stats["no_content"] += 1
        if budget.skipped:
            stats["similar_skipped"] += budget.skipped
            _LOG.info("pmid %s: skipped %d near-miss scans past the per-article budget",
                      pmid, budget.skipped)
        if anchor_names:
            stats["name_unanchored"] += _anchor(best)
        if not best:
            # Recorded too, and this is the case that matters: an article with no
            # matches leaves nothing in the database, so without the ledger it is
            # re-parsed on every resumed run -- and most articles match nothing.
            _mark(client, pmid)
            continue
        title = _title(client, pmid)
        for seq, mt, score, fname in best.values():
            loader.add_annotation(
                session, seq, source, domain="literature",
                key=f"pmid:{pmid}"[:64],
                value=title or curie_build("PMID", pmid),
            )
            loader.add_annotation(
                session, seq, source, domain="literature",
                key=f"suppmatch:{pmid}:{seq.uid}"[:64],
                value=_summary(mt, score, fname),
            )
            stats["refs"] += 1
            stats["matches"] += 1
            stats[mt] += 1
        session.commit()
        _mark(client, pmid)  # after the commit: the ledger never runs ahead of the DB
    stats["pmcid_via_openalex"] = getattr(client, "pmcid_via_openalex", 0)
    return stats


def _prefetched(
    client: PMCSuppClient, pmids: list[str], workers: int, chunk: int
) -> Iterator[str]:
    """Yield ``pmids`` in order, warming each chunk's downloads concurrently first.

    The mining loop is strictly serial and stays that way: it holds one database
    session, and its per-article match state (``_SeqBudget``, the best-match dict)
    is only coherent one article at a time. What parallelises cleanly is the part
    before it -- fetching a zip and a full-text XML over the network, which is where
    the wall clock actually goes and which touches nothing but the disk cache.

    Chunked rather than prefetching everything up front, for two reasons: mining
    starts after the first chunk instead of after 12,600 downloads, and at most
    ``chunk`` bundles are pulled ahead of what has been consumed, so the disk
    budget is approached gradually rather than in one leap.

    A prefetch failure is invisible here: the article is simply not cached, and the
    mining loop fetches it itself exactly as before.
    """
    warm_many = getattr(client, "prefetch", None)
    if not pmids or workers <= 1 or warm_many is None:
        # No prefetch capability (a stub client in tests, or a caller that asked for
        # none): mine straight through, fetching lazily. Same duck-typing as
        # ``full_text_fragments`` and ``mark_mined`` elsewhere in this pass.
        yield from pmids
        return
    for start in range(0, len(pmids), chunk):
        window = pmids[start:start + chunk]
        t0 = time.monotonic()
        warm_many(window, workers=workers)
        _LOG.info(
            "prefetched %d article(s) [%d/%d] in %.0fs",
            len(window), start + len(window), len(pmids), time.monotonic() - t0,
        )
        yield from window


def _anchor(best: dict[int, tuple[Sequence, str, float, str]]) -> int:
    """Drop name matches the article's own sequence evidence does not support.

    An allele *name* is a string many species' nomenclatures share, so a paper that
    writes ``IGHV1-69*01`` matches that name in every species carrying it. Measured
    on a twelve-article trial of human IGHV papers: 1,806 name matches, of which dog
    (189), cat (170) and sheep (113) each outranked human (112). Worse, a paper that
    merely *lists* names (a gene panel in a spatial-transcriptomics supplement) matched
    1,627 alleles and carried no sequence at all.

    A sequence match cannot collide that way -- it is the hash, or near it -- so it says
    which (species, locus) the article is actually about (``queries.anchor_pairs``), and
    a name match is kept only there. An article with no sequence evidence anchors
    nothing, so all its name matches go: we do not invent a species it never indicated.

    Returns the number dropped. Mutates ``best`` in place.
    """
    anchors = queries.anchor_pairs(
        (seq.species, seq.locus) for seq, mt, _s, _f in best.values() if mt != "name"
    )
    drop = [
        k for k, (seq, mt, _s, _f) in best.items()
        if mt == "name" and (seq.species, seq.locus) not in anchors
    ]
    for k in drop:
        del best[k]
    return len(drop)


def _mark(client: PMCSuppClient, pmid: str) -> None:
    """Note that ``pmid`` has been mined, if the client keeps a ledger."""
    marker = getattr(client, "mark_mined", None)
    if marker is not None:
        marker(pmid)


_WS = re.compile(r"\s+")


def _mine_fragment(
    session: Session,
    frag: str,
    best: dict[int, tuple[Sequence, str, float, str]],
    species_tags: set[str] | None,
    source_label: str,
    *,
    despace: bool = False,
    budget: _SeqBudget | None = None,
) -> None:
    """Fold DB matches for one text fragment into ``best``.

    Names are matched on the fragment as-is (whitespace-delimited tokens).
    Sequences are matched on the fragment optionally despaced (``despace=True``
    for full-text blocks, where a nt run can be split across inline markup);
    supplementary fragments keep their original spacing so adjacent cell runs
    are not merged."""
    for tok in _candidate_names(frag):
        for seq in _match_name(session, tok, species_tags):
            _consider(best, seq, "name", 1.0, source_label)
    seq_src = _WS.sub("", frag) if despace else frag
    for cand in _candidate_sequences(seq_src):
        if budget is not None:
            if cand in budget.seen:   # the same run recurs across cells and blocks
                continue
            budget.seen.add(cand)
        for seq, mt, score in _match_sequence(session, cand, species_tags, budget):
            _consider(best, seq, mt, score, source_label)


def _consider(
    best: dict[int, tuple[Sequence, str, float, str]],
    seq: Sequence,
    match_type: str,
    score: float,
    filename: str,
) -> None:
    """Keep the strongest match per sequence (see ``_MATCH_RANK``)."""
    cur = best.get(seq.id)
    if cur is None or _MATCH_RANK[match_type] > _MATCH_RANK[cur[1]]:
        best[seq.id] = (seq, match_type, score, filename)


def _title(client: PMCSuppClient, pmid: str) -> str | None:
    """Best-effort article title from the client (if it exposes one)."""
    getter = getattr(client, "article_title", None)
    if getter is None:
        return None
    try:
        return getter(pmid)
    except Exception:  # noqa: BLE001 - title is cosmetic; never fail the pass
        return None
