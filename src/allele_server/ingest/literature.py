"""Literature layer — auto-collect PubMed leads by matching the allele name.

For each coding allele we query **Europe PMC** for its exact name (e.g.
``"IGHV1-2*01"``) and record the matching articles as ``literature`` PMID
annotations. This is a *discovery* layer, refreshed on a schedule: it turns the
name every source already gives an allele into a list of papers that mention it.

**Honesty (v1):** matches are by name string only and are **not manually
verified** — a paper mentioning ``IGHV1-2*01`` may not be *about* that allele,
and common gene names collide. The exact-name query is AND'd with an
immunoglobulin/TCR context clause (``_build_query``) so obvious out-of-domain
collisions are dropped, but the result is still a *lead*: every collected
reference is stored/shown as **unverified**; curator-provided references (from
the internal tables) stay verified. Treat unverified entries as leads, not
ground truth.

Same rails as the NCBI pass: disk cache first (idempotent, offline re-runs),
polite rate limiting, graceful on network failure. The DB write is idempotent on
the annotation natural key, so re-running changes zero rows.
"""

from __future__ import annotations

import hashlib
import json
import re
import time
import urllib.parse
import urllib.request
from dataclasses import dataclass
from pathlib import Path

from sqlalchemy import func, select
from sqlalchemy.orm import Session

from ..curie import build as curie_build
from ..db.models import Name, Sequence
from . import loader
from .openalex import OpenAlexClient

_EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
_USER_AGENT = "allele-server/0.1 (contact: aperes7@gmail.com)"

# Immunoglobulin / TCR context clause, AND'd with the exact allele name so a
# bare-name collision (a clone id, an unrelated gene symbol, a random table row)
# in a paper that is *not* about antibody/TCR genetics is dropped before it ever
# becomes a lead. Precision filter for the v1 name probe (still unverified).
_IG_CONTEXT = (
    '(immunoglobulin OR "T cell receptor" OR "T-cell receptor" OR germline '
    'OR antibody OR repertoire OR IGH OR IGK OR IGL OR TRA OR TRB OR TRG OR TRD)'
)


def _build_query(phrase: str, *, scope: str = "title_abstract") -> str:
    """Europe PMC query: exact allele name AND an Ig/TCR context AND a PubMed id.

    ``scope`` controls where the allele name must appear:

    * ``"title_abstract"`` (default, for the displayed reference layer) — the name
      must be in the **title or abstract**, so the paper is much more likely to be
      *about* the allele. Drops the long tail of papers that only mention the name
      in full text or a supplementary table (the main noise source).
    * ``"all"`` — the name anywhere in the record (full text incl. supplements).
      Broader/noisier; this is the right net for **supplementary-file mining**
      (``literature_supp``), where a name in a supp table is exactly the target.

    The context clause is the precision filter (see ``_IG_CONTEXT``); ``SRC:MED``
    restricts to records that carry a PMID so every hit resolves to a reference.
    """
    if scope == "title_abstract":
        core = f'(TITLE:"{phrase}" OR ABSTRACT:"{phrase}")'
    else:
        core = f'"{phrase}"'
    return f'{core} AND {_IG_CONTEXT} AND SRC:MED'


@dataclass(frozen=True)
class Article:
    pmid: str
    title: str | None


class EuropePMCClient:
    """Minimal Europe PMC search client with a JSON disk cache.

    ``search(phrase)`` returns the articles whose text matches the exact phrase,
    capped at ``page_size``. Cached under ``cache_dir/<sha256(phrase)>.json``; a
    cache hit never touches the network. Network/parse failures return ``[]``
    (the pass is best-effort — a missing lead is not an error)."""

    def __init__(self, cache_dir: str | Path, *, page_size: int = 25,
                 min_interval: float = 0.34, timeout: int = 30,
                 scope: str = "title_abstract") -> None:
        self.cache = Path(cache_dir)
        self.cache.mkdir(parents=True, exist_ok=True)
        self.page_size = page_size
        self.min_interval = min_interval
        self.timeout = timeout
        self.scope = scope
        self._last = 0.0

    def _cache_path(self, query: str) -> Path:
        # Key on the full built query (not the bare name) so a change in scope or
        # context self-invalidates the cache — no stale broad hits on a strict re-run.
        return self.cache / f"{hashlib.sha256(query.encode()).hexdigest()}.json"

    def search(self, phrase: str) -> list[Article]:
        query = _build_query(phrase, scope=self.scope)
        cp = self._cache_path(query)
        if cp.exists():
            try:
                raw = json.loads(cp.read_text())
                return [Article(a["pmid"], a.get("title")) for a in raw]
            except (json.JSONDecodeError, KeyError):
                pass
        articles = self._fetch(query)
        cp.write_text(json.dumps([{"pmid": a.pmid, "title": a.title} for a in articles]))
        return articles

    # -- paginated search (article-centric discovery) ------------------------

    def search_query(
        self, query: str, *, max_results: int | None = None, page_size: int = 1000
    ) -> list[Article]:
        """Every article matching a **raw** Europe PMC query, walking all pages.

        The per-allele :meth:`search` builds its own query and reads one page,
        which is the right shape for "papers about this one allele". Discovery is
        the other direction -- a handful of broad queries, each with thousands of
        hits -- so this takes the query as given and pages through it with
        ``cursorMark``. Each page is cached under its own key, so an interrupted
        sweep resumes from where it stopped instead of re-walking the corpus.
        """
        out: list[Article] = []
        seen: set[str] = set()
        cursor = "*"
        while True:
            page, nxt = self._search_page(query, cursor, page_size)
            for art in page:
                if art.pmid not in seen:
                    seen.add(art.pmid)
                    out.append(art)
            if max_results and len(out) >= max_results:
                return out[:max_results]
            # Europe PMC signals exhaustion by returning the cursor unchanged.
            if not page or not nxt or nxt == cursor:
                return out
            cursor = nxt

    def _search_page(
        self, query: str, cursor: str, page_size: int
    ) -> tuple[list[Article], str | None]:
        """One page of a raw query: ``(articles, next_cursor)``, cached on disk."""
        cp = self._cache_path(f"page::{page_size}::{cursor}::{query}")
        if cp.exists():
            try:
                raw = json.loads(cp.read_text())
                return (
                    [Article(a["pmid"], a.get("title")) for a in raw["articles"]],
                    raw.get("next"),
                )
            except (json.JSONDecodeError, KeyError, TypeError):
                pass
        params = urllib.parse.urlencode({
            "query": query, "format": "json", "pageSize": str(page_size),
            "resultType": "lite", "cursorMark": cursor,
        })
        gap = self.min_interval - (time.monotonic() - self._last)
        if gap > 0:
            time.sleep(gap)
        req = urllib.request.Request(
            f"{_EPMC}?{params}", headers={"User-Agent": _USER_AGENT}
        )
        try:
            with urllib.request.urlopen(req, timeout=self.timeout) as resp:
                body = json.loads(resp.read().decode())
        except Exception:  # noqa: BLE001 - best-effort, like the rest of this layer
            return [], None  # not cached: a transient failure should be retried
        finally:
            self._last = time.monotonic()
        articles = [
            Article(str(r["pmid"]).strip(), (r.get("title") or "").strip() or None)
            for r in body.get("resultList", {}).get("result", [])
            if str(r.get("pmid") or "").strip()
        ]
        nxt = body.get("nextCursorMark")
        cp.write_text(json.dumps({
            "articles": [{"pmid": a.pmid, "title": a.title} for a in articles],
            "next": nxt,
        }))
        return articles, nxt

    def _fetch(self, query: str) -> list[Article]:
        # ``query`` is the fully built Europe PMC query (see ``_build_query``):
        # exact allele name in the chosen scope, Ig/TCR context, PMID-bearing records.
        params = urllib.parse.urlencode({
            "query": query, "format": "json", "pageSize": str(self.page_size),
            "resultType": "lite",
        })
        gap = self.min_interval - (time.monotonic() - self._last)
        if gap > 0:
            time.sleep(gap)
        req = urllib.request.Request(f"{_EPMC}?{params}", headers={"User-Agent": _USER_AGENT})
        try:
            with urllib.request.urlopen(req, timeout=self.timeout) as resp:
                body = json.loads(resp.read().decode())
        except Exception:  # best-effort: network/parse failure -> no leads
            return []
        finally:
            self._last = time.monotonic()
        out: list[Article] = []
        for r in body.get("resultList", {}).get("result", []):
            pmid = str(r.get("pmid") or "").strip()
            if pmid:
                out.append(Article(pmid, (r.get("title") or "").strip() or None))
        return out


# --- article-centric discovery --------------------------------------------

# Gene-name prefixes worth probing the corpus with. Derived from the names we
# actually hold, then filtered to this grammar -- the raw prefixes also include
# lab and clone identifiers (MUSI, LJI, LIB, VH, BALB, RP), which as a search
# phrase would drag in whole unrelated literatures.
_GENE_PREFIX_RE = re.compile(r"^(?:IG[HKLI][VDJ]|IGH[GAMED]|IG[KLI]C|TR[ABGD][VDJ])$")

# The mineable corpus is the open-access subset: supplementary bundles and body
# full text come from PMC-OA, so an article outside it yields nothing to match
# against however promising its abstract. Restricting discovery to OA is not a
# quality filter, it is the difference between a work-list and a wish-list.
#
# But Europe PMC's flag is not the same thing as the article's licence, and it is
# wrong often enough to matter: PMID 42251052 is flagged closed with no PMCID
# while a PMC copy (PMC13396459) serves HTTP 200. So the clause selects the
# *cheap* half of the corpus and an optional ``oa_client`` (OpenAlex) is asked
# about the rest — see :mod:`allele_server.ingest.openalex`.
_OA_CLAUSE = "OPEN_ACCESS:Y"


def gene_prefixes(
    session: Session, *, min_sequences: int = 20, species_keys: list[str] | None = None
) -> list[str]:
    """Gene-name prefixes present in the database, commonest first.

    These are the phrases the corpus sweep searches on -- ``IGHV``, ``TRBV``,
    ``IGHG`` -- rather than one phrase per allele. A paper that names any allele
    of a family almost always names the family too, so a couple of dozen prefix
    queries reach the same articles that ~17,800 per-allele queries would.
    """
    prefix = func.upper(func.substring(Name.name, "^[A-Za-z]{2,4}"))
    stmt = (
        select(prefix, func.count(func.distinct(Sequence.id)))
        .join(Sequence, Sequence.id == Name.sequence_id)
        .where(Sequence.category == "coding", Sequence.archived.is_(False))
        .group_by(prefix)
        .having(func.count(func.distinct(Sequence.id)) >= min_sequences)
        .order_by(func.count(func.distinct(Sequence.id)).desc())
    )
    if species_keys:
        from ..species import normalize_species

        stmt = stmt.where(
            Sequence.species.in_([normalize_species(s) for s in species_keys])
        )
    return [p for p, _ in session.execute(stmt) if p and _GENE_PREFIX_RE.match(p)]


def discover_articles(
    session: Session,
    client: EuropePMCClient,
    *,
    prefixes: list[str] | None = None,
    species_keys: list[str] | None = None,
    open_access_only: bool = True,
    max_articles: int | None = None,
    oa_client: OpenAlexClient | None = None,
) -> tuple[list[str], dict[str, int], dict[str, int]]:
    """Find the articles worth mining, in a few dozen queries instead of ~17,800.

    The per-allele probe (:func:`enrich_literature`) asks "which papers mention
    this allele?" once per allele. That is one network round trip per allele
    before a single document has been read, and it re-asks the whole question
    every time a paper is revisited. This asks the question the other way round:
    which papers are about immunoglobulin or TCR genes at all? Each answer is
    then fetched **once** by the supplementary/full-text miner and matched
    against every allele we hold in a single local pass.

    Pass ``oa_client`` (an :class:`~allele_server.ingest.openalex.OpenAlexClient`)
    to second-guess Europe PMC's open-access flag. Each prefix is then queried
    twice — once with ``OPEN_ACCESS:Y``, once without — and the articles only the
    broad query returned are put to OpenAlex; those with a legal OA copy are kept.
    Both queries are page-cached, so the extra sweep costs one pass and nothing on
    a re-run.

    Returns ``(pmids, per_prefix_counts, oa_stats)``. PMIDs are deduplicated
    across prefixes in first-seen order, so the most heavily represented families
    lead. ``oa_stats`` reports what the OA-resolution step did, including how many
    recovered articles have **no** PMC copy — those are kept (a legal OA copy does
    exist) but the PMC-OA miner cannot fetch them, and will report them as
    ``no_pmcid`` rather than dropping them quietly.
    """
    prefixes = prefixes or gene_prefixes(session, species_keys=species_keys)
    seen: dict[str, None] = {}
    counts: dict[str, int] = {}
    oa_stats = {
        "closed_candidates": 0, "openalex_unknown": 0, "confirmed_closed": 0,
        "recovered": 0, "recovered_with_pmc": 0, "recovered_without_pmc": 0,
    }
    resolve = oa_client is not None and open_access_only
    for prefix in prefixes:
        base = f'"{prefix}" AND {_IG_CONTEXT} AND SRC:MED'
        query = f"{base} AND {_OA_CLAUSE}" if open_access_only else base
        found = client.search_query(query)
        counts[prefix] = len(found)
        for art in found:
            seen.setdefault(art.pmid)
        if resolve:
            closed = [a.pmid for a in client.search_query(base) if a.pmid not in seen]
            _recover(oa_client, closed, seen, oa_stats)
        if max_articles and len(seen) >= max_articles:
            break
    pmids = list(seen)
    return (pmids[:max_articles] if max_articles else pmids), counts, oa_stats


def _recover(
    oa_client: OpenAlexClient, closed: list[str], seen: dict[str, None],
    oa_stats: dict[str, int],
) -> None:
    """Add back the Europe-PMC-closed articles OpenAlex confirms are open access."""
    works = oa_client.works(closed)
    for pmid in closed:
        oa_stats["closed_candidates"] += 1
        work = works.get(pmid)
        if work is None:
            oa_stats["openalex_unknown"] += 1
            continue
        if not work.is_oa:
            oa_stats["confirmed_closed"] += 1
            continue
        oa_stats["recovered"] += 1
        oa_stats["recovered_with_pmc" if work.pmcid else "recovered_without_pmc"] += 1
        seen.setdefault(pmid)


def _representative_name(session: Session, seq_id: int) -> str | None:
    """The allele name to search on: prefer the IMGT-nomenclature name."""
    pref = session.scalar(
        select(Name.name)
        .where(Name.sequence_id == seq_id, Name.naming_scheme == "imgt")
        .limit(1)
    )
    return pref or session.scalar(
        select(Name.name).where(Name.sequence_id == seq_id).limit(1)
    )


def enrich_literature(
    session: Session,
    client: EuropePMCClient,
    *,
    species_keys: list[str] | None = None,
    segments: tuple[str, ...] = ("V",),
    limit: int | None = None,
    max_refs: int = 8,
) -> dict[str, int]:
    """Attach unverified PubMed leads to coding alleles by exact-name search.

    Only names that carry an allele designation (``*NN``) are searched — a bare
    gene name is too ambiguous. Returns ``{"searched", "with_refs", "refs"}``.
    """
    source = loader.get_or_create_source(
        session, name="literature", version=None, curie_prefix="PMID",
        note="auto-collected PubMed leads (Europe PMC name match; unverified)",
    )
    stmt = select(Sequence).where(
        Sequence.category == "coding", Sequence.segment.in_(segments),
        Sequence.archived.is_(False),
    )
    if species_keys:
        from ..species import normalize_species

        stmt = stmt.where(Sequence.species.in_([normalize_species(s) for s in species_keys]))
    if limit:
        stmt = stmt.limit(limit)

    stats = {"searched": 0, "with_refs": 0, "refs": 0}
    for seq in session.scalars(stmt):
        name = _representative_name(session, seq.id)
        if not name or "*" not in name:  # need a specific allele, not a bare gene
            continue
        stats["searched"] += 1
        articles = client.search(name)[:max_refs]
        if not articles:
            continue
        stats["with_refs"] += 1
        for art in articles:
            loader.add_annotation(
                session, seq, source, domain="literature",
                key=f"pmid:{art.pmid}"[:64], value=art.title or curie_build("PMID", art.pmid),
            )
            stats["refs"] += 1
        session.commit()
    return stats
