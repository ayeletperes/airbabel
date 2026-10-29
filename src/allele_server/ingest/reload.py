"""Consolidated full reload: truncate everything and load all sources uniformly.

Order matters, and it is a correctness constraint rather than a preference:

1. ``_genedb``: the **all-species** IMGT/GENE-DB pass (IG + TR V/D/J/C, IMGT's
   published V protein, then leaders/RSS). It is the only IMGT input.
2. ``ingest_cache``: OGRDB. After IMGT, so a V allele both carry gets
   IMGT's published protein before OGRDB's translation fallback could fill it.
3. The internal curation tables (HUSA human, MUSA rhesus) and KI/KIARVA, which add
   novel alleles and the RSS network on top of that backbone.
4. ``_enrich``: historical IMGT names, then ``mark_archived``, then the IMGT-name and
   NCBI crawls. ``mark_archived`` flags sequences attested by **no current source**,
   so it must run after *every* current source above. Run it earlier (as it was, before
   the GENE-DB pass existed here) and alleles that pass 2 would have supplied get
   marked archived by mistake.
5. AA backfill, then the neighbour cache, which depends on everything above.

All upserts are idempotent and keyed on the canonical sequence hash, so overlapping
sources corroborate the same sequence rather than duplicating it.
"""

from __future__ import annotations

from pathlib import Path

from sqlalchemy import text
from sqlalchemy.orm import Session

from .curation_table import ingest_husa, ingest_musa
from .ki_xlsx import ingest_ki_xlsx
from .ogrdb import ingest_cache

# Reverse-dependency order. `source_candidates` is deliberately ABSENT: approvals and
# rejections are curatorial state, not data derived from a source, and a reload must
# never destroy them. It carries no FK to sources.id precisely so the CASCADE below
# cannot reach it. See alembic/versions/0005_source_candidates.py and
# tests/test_reload.py::test_reload_preserves_source_candidates.
#
# `dataset_releases` is absent for the same reason: it records what was served and
# when. This reload appends the next row (see stamp_release below); erasing the
# earlier ones would destroy the only account of what an older deployment held.
# It carries no FK to sequences, so the CASCADE cannot reach it either.
_TABLES = "sequence_annotations, sequence_source, sequence_links, names, sequences, sources"


def reload_all(
    session: Session,
    cache_dir: str | Path,
    *,
    musa_path: str | None = None,
    husa_path: str | None = None,
    ki_path: str | None = None,
    species_keys=("human", "mouse", "rhesus"),
    genedb: bool = True,
    enrich: bool = True,
    neighbors: bool = False,
) -> dict[str, object]:
    """Truncate all tables and reload every source. Returns per-source counts.

    ``genedb=True`` (default) runs the all-species IMGT/GENE-DB pass, the only IMGT
    input. Turn it off only for a fast, deliberately partial (IMGT-less) rebuild.

    ``enrich=True`` (default) re-runs the disk-cached enrichment passes. They write
    into tables this reload truncates, so skipping them loses the historical IMGT
    names, the ``archived`` flags, the IMGT clone-names and the NCBI xrefs.

    ``neighbors=True`` rebuilds the precomputed similarity cache after everything
    else (so its denormalised labels reflect the names enrichment adds). The cache is
    **also** rebuilt automatically when it was non-empty before the reload: the
    TRUNCATE cascades into ``sequence_neighbors`` (FK on ``sequences.uid``), so
    without this an existing cache would be silently emptied and every "find similar"
    would fall back to the slow live path until someone re-ran the build.

    See the module docstring: the step order is a correctness constraint.
    """
    # A curated table that is configured but not there is a wrong path, not an
    # opt-out: refuse before the TRUNCATE (the passes commit as they go), so a typo
    # can never leave a half-loaded database. An unset path is the opt-out.
    missing = [p for p in (husa_path, musa_path, ki_path) if p and not Path(p).is_file()]
    if missing:
        raise FileNotFoundError(f"configured curation table(s) not found: {missing}")
    # Check before truncating: afterwards the cascade has already emptied it.
    had_cache = bool(session.scalar(text("SELECT count(*) FROM sequence_neighbors")))
    session.execute(text(f"TRUNCATE {_TABLES} RESTART IDENTITY CASCADE"))
    session.flush()
    counts_invalidated = _invalidate_supp_ledger()

    counts: dict[str, object] = {"supp_ledger_invalidated": counts_invalidated}
    if genedb:
        counts.update(_genedb(session))
    counts.update(ingest_cache(session, cache_dir, species_keys))
    if husa_path:
        counts["husa"] = ingest_husa(session, husa_path)
    if musa_path:
        counts["musa"] = ingest_musa(session, musa_path)
    if ki_path:
        counts["ki"] = ingest_ki_xlsx(session, ki_path)
    if enrich:
        # Must come after every current source above: `mark_archived` inside decides
        # what exists *only* in historical releases.
        counts.update(_enrich(session, species_keys))
    # Fill AA for any V rows whose ingest path didn't translate (before the
    # neighbour build, so the aa cache covers them).
    from .loader import backfill_v_aa

    counts["v_aa_backfill"] = len(backfill_v_aa(session))
    if neighbors or had_cache:
        from .neighbors import build_cross_species, build_neighbors

        counts["neighbors"] = build_neighbors(session)
        counts["neighbors_cross"] = build_cross_species(session)
    # Last, and reading the database rather than the counts above: a release is a
    # statement about what is *stored*, including whether the neighbour cache was
    # rebuilt in this run. See ingest/release_stamp.py.
    session.flush()
    from .release_stamp import stamp_release

    counts["release"] = stamp_release(session)["release"]
    return counts


def _invalidate_supp_ledger() -> int:
    """Drop the supplement miner's resume ledger. Returns the entry count removed.

    ``PMCSuppClient`` records finished PMIDs in ``mined.txt`` so an interrupted
    multi-day sweep resumes instead of re-reading ~12,600 articles. That ledger
    asserts "these articles have been mined **and their matches are in the
    database**" -- and the TRUNCATE above has just destroyed the second half of that
    claim. Left in place, the next resumed sweep would skip every article it had
    already finished and never rewrite a single annotation, so the literature layer
    would come back empty and look complete.

    Only the ledger goes. The downloaded bundles and full-text XML stay cached, so
    re-mining costs CPU and database time, not another pass over the network -- which
    is the expensive half.

    Best-effort: a missing or unreadable ledger is the normal case (no sweep has run
    on this host), not an error.
    """
    from ..settings import get_settings

    ledger = Path(get_settings().literature_supp_cache_dir) / "mined.txt"
    try:
        n = sum(1 for line in ledger.read_text().splitlines() if line.strip())
        ledger.unlink()
    except OSError:
        return 0
    return n


def _genedb(session: Session) -> dict[str, object]:
    """The all-species IMGT/GENE-DB pass: IG + TR V/J coding, then leaders + V-RS.

    ``fetch_reference`` is idempotent (it reuses the bulk FASTA already on disk), so
    the only network cost here is ``ingest_leaders_rss``, which queries IMGT GENElect
    per species. Both write through the same hash-keyed upserts as every other source,
    so an allele OGRDB also carries is corroborated, not duplicated.
    """
    from ..settings import get_settings
    from ..species import known_tags
    from .imgt_constant import ingest_constant
    from .imgt_genedb import (
        fetch_reference,
        fetch_reference_aa,
        fetch_reference_gapped,
        ingest_leaders_rss,
        ingest_reference,
        ingest_reference_aa,
    )

    st = get_settings()
    out: dict[str, object] = {}
    fasta = fetch_reference(st.imgt_genedb_fasta)
    fetch_reference_gapped(st.imgt_genedb_fasta_gapped)  # diff & align's reference
    out["imgt_genedb_v"] = ingest_reference(session, fasta, segments=("V", "D", "J"))
    # Same file, no second download: the constant domains are already in it.
    out["imgt_genedb_c"] = ingest_constant(session, fasta)
    session.flush()
    # IMGT's published V proteins, so aa_origin distinguishes a sourced protein from
    # one we translated. Must run after the coding ingest: it joins on name.
    out["imgt_genedb_aa"] = ingest_reference_aa(
        session, fetch_reference_aa(st.imgt_genedb_fasta_aa)
    )
    session.flush()  # leaders/RSS match V coding by name, so it must see the rows above
    out["imgt_leaders_rss"] = ingest_leaders_rss(
        session, sorted(known_tags()), cache_dir=st.imgt_genelect_cache_dir
    )
    session.flush()
    return out


def _enrich(session: Session, species_keys) -> dict[str, object]:
    """Run the disk-cached, idempotent enrichment passes (IMGT names + NCBI).

    Call only after **all** current sources are loaded: ``mark_archived`` marks the
    sequences no current source attests, so anything not yet ingested would be
    mislabelled archived and dropped from coverage and default search.
    """
    from ..settings import get_settings
    from ._ncbi_client import NCBIClient
    from .imgt_crawl import crawl_imgt_names
    from .imgt_history import ingest_imgt_history, mark_archived
    from .ncbi import enrich_ncbi

    st = get_settings()
    out: dict[str, object] = {}
    # Former/retired official IMGT names from historical release snapshots — reload
    # truncates the names table, so these must be re-attached here or they're lost.
    out["imgt_history"] = ingest_imgt_history(session, st.imgt_history_dir)
    # Flag alleles that exist only in historical releases (gone from every current
    # source) so they don't inflate current coverage / pollute search results.
    out["archived"] = mark_archived(session)
    out["imgt_names"] = crawl_imgt_names(
        session, cache_dir=st.imgt_crawl_cache_dir, species_keys=tuple(species_keys)
    )
    client = NCBIClient(cache_dir=st.ncbi_cache_dir, api_key=st.ncbi_api_key)
    out["ncbi"] = enrich_ncbi(session, client, species_keys=list(species_keys))
    return out
