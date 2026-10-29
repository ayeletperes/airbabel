"""REST API routes. Thin projection over :mod:`allele_server.queries`."""

from __future__ import annotations

import threading
from collections.abc import Iterator
from contextlib import contextmanager

from fastapi import APIRouter, Depends, HTTPException, Query, Request
from fastapi.responses import StreamingResponse
from sqlalchemy.orm import Session

from .. import __version__, export, queries
from ..db.session import session_factory
from ..settings import get_settings
from .mcp_gate import RateLimiter
from .schemas import (
    AaSearchRequest,
    Alias,
    AlignmentView,
    AlleleDescription,
    DiffSummary,
    OrthologHint,
    PmidResult,
    ResolveResult,
    RssAllele,
    RssSimilar,
    SearchResults,
    SequenceRecord,
    SequenceSearchRequest,
)

router = APIRouter()


def get_session() -> Iterator[Session]:
    sess = session_factory()()
    try:
        yield sess
    finally:
        sess.close()


def _cap() -> int:
    return get_settings().search_result_cap


# --- de novo partial-search concurrency limiter --------------------------------
# Partial (containment) search is uncacheable and does a full trigram scan +
# local-identity rescore per request. On a small box a burst of them can saturate
# CPU, so we bound how many run at once. Routes are sync (run in the threadpool),
# so a threading semaphore is the right primitive; excess requests wait briefly
# for a slot, then get 503. The semaphore is a lazily-built module singleton sized
# from settings (get_settings() returns a fresh instance, so it can't own it).
_PARTIAL_SEM: threading.BoundedSemaphore | None = None
_PARTIAL_SEM_LOCK = threading.Lock()


def _partial_semaphore() -> threading.BoundedSemaphore:
    global _PARTIAL_SEM
    if _PARTIAL_SEM is None:
        with _PARTIAL_SEM_LOCK:
            if _PARTIAL_SEM is None:
                n = max(1, get_settings().partial_search_concurrency)
                _PARTIAL_SEM = threading.BoundedSemaphore(n)
    return _PARTIAL_SEM


@contextmanager
def maybe_partial_slot(partial: bool) -> Iterator[None]:
    """Hold a concurrency slot for a de novo partial search; no-op otherwise.

    A 503 is raised if the box is saturated and no slot frees within
    ``partial_search_queue_timeout`` seconds. Shared by the REST routes and the
    web UI so both draw on the same CPU budget.
    """
    if not partial:
        yield
        return
    sem = _partial_semaphore()
    timeout = max(0.0, get_settings().partial_search_queue_timeout)
    if not sem.acquire(timeout=timeout):
        raise HTTPException(
            status_code=503,
            detail="Partial search is at capacity; please retry in a moment.",
        )
    try:
        yield
    finally:
        sem.release()


@router.get("/healthz")
def healthz() -> dict[str, str]:
    return {"status": "ok"}


@router.get("/version")
def version(session: Session = Depends(get_session)) -> dict:
    """Code version *and* the stamped dataset release this instance serves.

    Two separate things, kept separate: ``version`` is the software, ``dataset`` is
    the data plus the distance computation it is scored with. An instance that has
    never been stamped reports ``dataset: null`` rather than borrowing the software
    version, because claiming a data version you do not have is the failure this
    endpoint exists to prevent.
    """
    return {
        "name": "allele-server",
        "version": __version__,
        "dataset": queries.dataset_release(session),
        # What this process would compute *now*. Compare against
        # dataset.similarity_digest: a mismatch means the code has moved since the
        # stamp, so the served scores are not the ones the release describes.
        "similarity": queries.similarity_profile(),
        "similarity_digest": queries.similarity_digest(),
    }


@router.get("/api/resolve")
def resolve(
    query: str = Query(..., min_length=1),
    species: str | None = Query(
        None,
        description="Scope resolution to one species, or several comma-separated "
        "(human,rhesus). A name or sequence reused across species is otherwise "
        "resolved to a single, deterministic species; use /api/search/name?exact=1 "
        "to list every species that carries a name.",
    ),
    view: str | None = Query(None, pattern="^airr$"),
    session: Session = Depends(get_session),
) -> ResolveResult | AlleleDescription:
    result = queries.resolve(session, query, species=species)
    if result is None:
        raise HTTPException(status_code=404, detail=f"no match for {query!r}")
    if view == "airr":
        return result.record.airr
    return result


@router.get("/api/sequences/{uid}")
def get_sequence(
    uid: str,
    view: str | None = Query(None, pattern="^airr$"),
    session: Session = Depends(get_session),
) -> SequenceRecord | AlleleDescription:
    result = queries.resolve(session, uid)
    if result is None or result.matched_by not in ("uid", "protein_uid"):
        raise HTTPException(status_code=404, detail=f"no sequence with uid {uid!r}")
    if view == "airr":
        return result.record.airr
    return result.record


@router.get("/api/sequences/{uid}/aliases")
def aliases(uid: str, session: Session = Depends(get_session)) -> list[Alias]:
    result = queries.list_aliases(session, uid)
    if not result:
        raise HTTPException(status_code=404, detail=f"no sequence with uid {uid!r}")
    return result


@router.get("/api/sequences/{uid}/similar")
def similar(
    uid: str,
    min_identity: float = Query(0.85, ge=0.0, le=1.0),
    same_species: bool = Query(True),
    to_species: str | None = Query(
        None,
        description="Restrict neighbours to one species (e.g. 'human') or several "
                    "comma-separated ('human,mouse'); the cross-species "
                    "closest-ortholog lookup. Overrides same_species. A selection "
                    "mixing the query's own species with others is served at the "
                    "0.85 cross-species floor, since below that the precomputed "
                    "cache holds same-species edges the cross-species side lacks.",
    ),
    kind: str = Query("nt", pattern="^(nt|aa)$"),
    functional: bool | None = Query(
        None,
        description="Filter neighbours by functionality: true = functional only, "
        "false = non-functional only, unset = all (incl. unknown).",
    ),
    source: str | None = Query(
        None,
        description="Keep only neighbours contributed by this source (IMGT, OGRDB, "
        "HUSA, MUSA, KIARVA, ...); case-insensitive. Unknown source = no matches.",
    ),
    session: Session = Depends(get_session),
) -> SearchResults:
    fn = queries.find_similar_aa if kind == "aa" else queries.find_similar
    result = fn(
        session, uid, min_identity=min_identity, cap=_cap(),
        same_species=same_species, target_species=to_species,
        functional=functional, source=source,
    )
    if result is None:
        raise HTTPException(status_code=404, detail=f"no sequence with uid {uid!r}")
    return result


@router.get("/api/sequences/{uid}/diff")
def diff(
    uid: str,
    min_identity: float = Query(0.90, ge=0.0, le=1.0),
    session: Session = Depends(get_session),
) -> DiffSummary:
    """Exact nt substitutions to the nearest allele(s), with the protein-level
    synonymous / non-synonymous consequence (INT-1). V substitutions carry an
    IMGT FR/CDR region label when the IMGT-gapped reference is available."""
    result = queries.diff_to_neighbours(
        session, uid, cap=_cap(), min_identity=min_identity
    )
    if result is None:
        raise HTTPException(status_code=404, detail=f"no sequence with uid {uid!r}")
    return result


@router.get("/api/sequences/{uid}/alignment")
def alignment(
    uid: str,
    against: str | None = Query(
        None, description="UID to align against; defaults to the gene's *01, else "
        "the nearest neighbour.",
    ),
    session: Session = Depends(get_session),
) -> AlignmentView:
    """IMGT-gapped AA alignment of a V allele against a reference (SIM-3). 404 for
    a non-V record or when the IMGT-gapped reference is not available."""
    result = queries.imgt_alignment(session, uid, against=against)
    if result is None:
        raise HTTPException(
            status_code=404,
            detail=f"no IMGT alignment for {uid!r} (non-V, or gapped reference "
                   "unavailable)",
        )
    return result


@router.get("/api/sequences/{uid}/orthologs")
def orthologs(
    uid: str, session: Session = Depends(get_session)
) -> list[OrthologHint]:
    """Cross-species nearest allele per species, reciprocal best hits flagged as a
    soft sequence-based orthology hint (SIM-2) — NOT a curated orthology call."""
    if queries.resolve(session, uid) is None:
        raise HTTPException(status_code=404, detail=f"no sequence with uid {uid!r}")
    return queries.reciprocal_orthologs(session, uid)


@router.get("/api/rss/{uid}")
def rss(uid: str, session: Session = Depends(get_session)) -> list[RssAllele]:
    """Coding alleles linked to this RSS, each with the source(s) of that link."""
    found = queries.rss_detail(session, uid)
    if found is None:
        raise HTTPException(status_code=404, detail=f"no RSS with uid {uid!r}")
    return found[1]


@router.get("/api/rss/{uid}/similar")
def rss_similar(
    uid: str,
    min_identity: float = Query(0.85, ge=0.0, le=1.0),
    same_species: bool = Query(True),
    session: Session = Depends(get_session),
) -> RssSimilar:
    """RSS sequences similar to this one (identity-ranked)."""
    result = queries.find_similar_rss(
        session, uid, min_identity=min_identity, cap=_cap(), same_species=same_species
    )
    if result is None:
        raise HTTPException(status_code=404, detail=f"no RSS with uid {uid!r}")
    return result


@router.get("/api/search/name")
def search_name(
    q: str = Query(..., min_length=1),
    species: str | None = Query(
        None, description="One species, or several comma-separated (human,rhesus).",
    ),
    exact: bool = Query(
        False,
        description="Return only exact-name matches, across species unless "
        "species is given. Use to list every species that reuses a name.",
    ),
    partial: bool = Query(
        False,
        description="Browse mode: every allele whose name contains the query, so "
        "a fragment (IGHV1-6, 1-69) returns all variations to refine among. "
        "Default (neither flag) is precise: an exact allele, or a bare gene name "
        "resolved to exactly its own alleles (IGHV1-6 does not pull IGHV1-69).",
    ),
    subgroup: bool = Query(
        False,
        description="Treat the query as an IMGT subgroup (IGHV1) and return every "
        "allele in it, anchored so IGHV1 does not pull IGHV10-1. Returns nothing "
        "for a query that is not a bare subgroup. The default path already falls "
        "back to this when a bare-subgroup query matches nothing else.",
    ),
    functional: bool | None = Query(
        None,
        description="Filter by functionality: true = functional only, false = "
        "non-functional only, unset = all (incl. unknown).",
    ),
    source: str | None = Query(
        None,
        description="Keep only records contributed by this source (IMGT, OGRDB, "
        "HUSA, MUSA, KIARVA, ...); case-insensitive. Unknown source = no matches.",
    ),
    session: Session = Depends(get_session),
) -> SearchResults:
    return queries.search_by_name(
        session, q, cap=_cap(), species=species, exact=exact, partial=partial,
        functional=functional, source=source, subgroup=subgroup,
    )


@router.get("/api/pmid/{pmid}")
def alleles_for_pmid(
    pmid: str,
    species: str | None = Query(
        None, description="One species, or several comma-separated (human,rhesus).",
    ),
    include_name: bool = Query(
        False,
        description="Also list name-only leads (unverified; the paper merely contains "
        "the name string). Default off: they are not evidence. Even when on, only "
        "leads whose species and locus agree with the paper's own sequence/curator "
        "evidence are returned.",
    ),
    session: Session = Depends(get_session),
) -> PmidResult:
    """Every allele a PubMed record attests — the reverse of an allele's reference list.

    Accepts a bare PMID, `PMID:29163486`, or a PubMed URL. Each hit carries the
    evidence tier that attests it (`verified` / `sequence` / `name`), and the per-tier
    totals are reported separately: a mined name-only lead is a weaker claim than a
    curator-verified cross-reference, so they are never summed into one number.

    By default only `verified` and `sequence` hits are listed. A name-only link means
    the paper contains the string, not the allele; pass `include_name=true` to see the
    species/locus-consistent ones, labelled unverified. `n_name_unanchored` counts the
    name links dropped for disagreeing with the paper's own evidence.

    Archived alleles are **included** and flagged, because a paper citing an allele is
    the reason archived records stay reachable at all.
    """
    result = queries.alleles_for_pmid(
        session, pmid, species=species, cap=_cap(), include_name=include_name,
    )
    if result is None:
        raise HTTPException(
            status_code=422,
            detail=f"not a PMID: {pmid!r}. Expected digits, PMID:<digits>, or a "
                   f"PubMed URL.",
        )
    return result


@router.post("/api/search/sequence")
def search_sequence(
    body: SequenceSearchRequest,
    species: str | None = Query(
        None, description="One species, or several comma-separated (human,rhesus).",
    ),
    category: str = Query("coding", pattern="^(coding|leader|rss)$"),
    locus: str | None = Query(
        None,
        description="Restrict the search to one chain (IGH/IGK/IGL/TRA/TRB/TRG/TRD). "
        "Use when the chain is known, e.g. from an existing V call.",
    ),
    partial: bool = Query(
        False,
        description="Local (containment) identity: a truncated or divergent query "
        "resolves against its best-matching window, scoring a contained fragment "
        "near 1.0. Uncacheable and concurrency-limited.",
    ),
    functional: bool | None = Query(
        None,
        description="Filter by functionality: true = functional only, false = "
        "non-functional only, unset = all (incl. unknown).",
    ),
    source: str | None = Query(
        None,
        description="Keep only records contributed by this source (IMGT, OGRDB, "
        "HUSA, MUSA, KIARVA, ...); case-insensitive. Unknown source = no matches.",
    ),
    session: Session = Depends(get_session),
) -> SearchResults:
    try:
        with maybe_partial_slot(partial):
            return queries.search_by_sequence(
                session, body.sequence, cap=_cap(), species=species,
                category=category, partial=partial, locus=locus,
                functional=functional, source=source,
            )
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc


@router.post("/api/search/aa")
def search_aa(
    body: AaSearchRequest,
    species: str | None = Query(
        None, description="One species, or several comma-separated (human,rhesus).",
    ),
    locus: str | None = Query(
        None,
        description="Restrict the search to one chain (IGH/IGK/IGL/TRA/TRB/TRG/TRD). "
        "Use when the chain is known, e.g. from an existing V call.",
    ),
    partial: bool = Query(
        False,
        description="Local (containment) identity: a truncated or divergent query "
        "resolves against its best-matching window. Uncacheable and "
        "concurrency-limited.",
    ),
    functional: bool | None = Query(
        None,
        description="Filter by functionality: true = functional only, false = "
        "non-functional only, unset = all (incl. unknown).",
    ),
    source: str | None = Query(
        None,
        description="Keep only records contributed by this source (IMGT, OGRDB, "
        "HUSA, MUSA, KIARVA, ...); case-insensitive. Unknown source = no matches.",
    ),
    session: Session = Depends(get_session),
) -> SearchResults:
    try:
        with maybe_partial_slot(partial):
            return queries.search_by_aa(
                session, body.sequence_aa, cap=_cap(), species=species,
                partial=partial, locus=locus,
                functional=functional, source=source,
            )
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc


# --- bulk export ---------------------------------------------------------------
# The full dataset is downloadable through the API (IMGT is CC BY 4.0; the other
# sources are from published manuscripts). No per-record cap here; abuse is bounded
# by rate-limiting these endpoints at the proxy, not by hiding the data.

# The bulk-export limiter. Module-level so the window is shared across requests;
# built once from settings, like _PARTIAL_SEM above and for the same reason.
_EXPORT_LIMITER: RateLimiter | None = None
_EXPORT_LIMITER_LOCK = threading.Lock()


def export_rate_limit(request: Request) -> None:
    """Refuse a caller who is pulling the bulk export faster than a human would.

    Applied to the two endpoints that scan the whole table (FASTA and JSON), not
    to /api/export itself -- that one is a 3 KB manifest and is what a client reads
    to discover the others.

    429 with Retry-After, so a well-behaved client waits instead of hammering.
    """
    global _EXPORT_LIMITER
    if _EXPORT_LIMITER is None:
        with _EXPORT_LIMITER_LOCK:
            if _EXPORT_LIMITER is None:
                st = get_settings()
                _EXPORT_LIMITER = RateLimiter(
                    st.export_rate_limit, st.export_rate_window
                )
    limiter = _EXPORT_LIMITER
    # GLOBAL, not per caller, and that is forced by the network rather than chosen.
    # The load balancer in front of this box does not pass the client's address: Caddy
    # sees the balancer's own node, and that rotates over its five nodes, so a
    # "per-caller" bucket silently became five buckets per caller and the cap bit at
    # five times the intended volume. Measured: with the limit at 1, public requests
    # alternated 200/429 while the same requests made on the box gave 200/429/429/429.
    #
    # A global bucket is also the honest instrument for the actual goal, which is not
    # fairness between callers but keeping CPU credits and egress off the bill. One
    # busy client can briefly hold the budget; that is the trade, and the refusal says
    # to cache rather than re-fetch.
    # TODO(export-per-caller): if the balancer is ever set to preserve the client IP
    # (or to send X-Forwarded-For), switch the key back to caller_of(...).
    if limiter.over("__global__"):
        raise HTTPException(
            status_code=429,
            detail=(
                f"Bulk export is limited to {limiter.limit} requests per "
                f"{int(limiter.window)}s across all users. Each one is a full scan "
                f"of the dataset; please cache the file rather than re-fetching it."
            ),
            headers={"Retry-After": str(int(limiter.window))},
        )


@router.get("/api/export")
def export_index(session: Session = Depends(get_session)) -> dict:
    """Describe the bulk-export endpoints and the attribution obligation (CC BY)."""
    counts = export.export_counts(session)
    return {
        "attribution": export.ATTRIBUTION,
        "license": export.LICENCE_NOTE,
        # Cite the release, not just the URL: an export is only reproducible if the
        # reader can say which build it came from and how it was scored.
        "dataset_release": queries.dataset_release(session),
        # Reported, not hidden: the export carries only what a redistributable source
        # attests, and a reader comparing this against the site's own coverage figures
        # needs to see the difference stated rather than infer it from a short file.
        "coding_alleles": counts["exported"],
        "coding_alleles_in_database": counts["matched"],
        "withheld_no_redistribution_licence": counts["withheld"],
        "withheld_sources": sorted(
            m.label for m in export._sources.published()
            if not (m.licence and m.licence.redistributable)
        ),
        "formats": {
            "fasta": "/api/export/alleles.fasta?category=coding&segment=V&species=&kind=nt",
            "json": "/api/export/alleles.json?category=coding&segment=V&species=",
        },
    }


@router.get("/api/export/alleles.fasta",
            dependencies=[Depends(export_rate_limit)])
def export_fasta(
    category: str = Query("coding", pattern="^(coding|leader|rss)$"),
    species: str | None = Query(None),
    segment: str | None = Query(None, pattern="^[VDJ]$"),
    kind: str = Query("nt", pattern="^(nt|aa)$"),
    session: Session = Depends(get_session),
) -> StreamingResponse:
    """Download all matching alleles as FASTA (self-attributing headers)."""
    records = export.export_alleles(
        session, category=category, species=species, segment=segment
    )
    # Provenance travels with the file. A FASTA that leaves here and turns up in
    # someone's pipeline six months later has to be able to say which build it is,
    # and a header is the only place that survives the download.
    rel = queries.dataset_release(session)
    return StreamingResponse(
        export.to_fasta(records, kind=kind), media_type="text/x-fasta",
        headers={
            "Content-Disposition": "attachment; filename=airrbabel_alleles.fasta",
            "X-Data-Attribution": f"{export.LICENCE_NOTE} See /api/export.",
            "X-Dataset-Release": (rel or {}).get("release", "unstamped"),
            "X-Dataset-Digest": (rel or {}).get("content_digest", ""),
        },
    )


@router.get("/api/export/alleles.json",
            dependencies=[Depends(export_rate_limit)])
def export_json(
    category: str = Query("coding", pattern="^(coding|leader|rss)$"),
    species: str | None = Query(None),
    segment: str | None = Query(None, pattern="^[VDJ]$"),
    session: Session = Depends(get_session),
) -> dict:
    """Download all matching alleles as JSON, each with its names + sources."""
    records = export.export_alleles(
        session, category=category, species=species, segment=segment
    )
    return {
        "attribution": export.ATTRIBUTION,
        "dataset_release": queries.dataset_release(session),
        "count": len(records),
        "alleles": records,
    }
