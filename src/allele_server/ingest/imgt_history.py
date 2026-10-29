"""Recover *former* official IMGT allele names via sequence-keyed reconciliation
against historical IMGT/GENE-DB reference releases.

IMGT renames alleles over time — notably across the mouse IGH renomenclature —
and only serves the current release. Loading an *older* release through the
normal FASTA path attaches that release's official name as an ``imgt`` alias to
the **same sequence** (identity is the sequence hash), so a search for a retired
name resolves to today's allele. Each release is recorded as its own ``IMGT``
source *version*, so former-vs-current is recoverable from the source version.

Provenance of the snapshots: historical IMGT/GENE-DB releases are taken from
Jamie Heather's ``genedb-releases`` archive
(https://github.com/JamieHeather/genedb-releases), which recovered pre-2022
releases via the Internet Archive Wayback Machine and snapshots GENE-DB weekly
since 2023. Credit and thanks to Jamie for assembling them. The data itself is
IMGT's: the public ``source`` record attributes it to **IMGT** (with the release
version), so the reference provenance shown to users and consumers is IMGT.
"""

from __future__ import annotations

from pathlib import Path

from sqlalchemy import delete, insert, select, update
from sqlalchemy.orm import Session

from ..curie import build as curie_build
from ..naming import locus_of, scheme_label, segment_of
from ..sequence import hash_sequence
from ..species import SpeciesError, normalize_species
from . import loader
from .fasta import parse_record, read_fasta
# The IG + TR coding loci, shared with the current-release pass so the two cannot
# drift: a locus the main ingest loads but the history pass skips is exactly how
# TR ended up with no former names at all.
from .imgt_genedb import _CODING_LOCI, _base_species


#: The IMGT labels that publish a segment as ONE sequence — the same set
#: ``imgt_genedb.ingest_reference`` admits. Everything else in the reference FASTA is
#: a part (a constant-region domain, a leader exon, a transmembrane region), and a
#: part cannot hash-match the whole we store for that allele.
_HISTORY_LABELS = ("V-REGION", "D-REGION", "J-REGION")
#: Segments reconciled from the snapshots. C is absent because IMGT publishes IGH
#: constants per domain while we store them assembled; see the note in
#: ``ingest_imgt_release``.
_HISTORY_SEGMENTS = ("V", "D", "J")


def _species_key(imgt_species: str) -> str | None:
    """IMGT's species field (pipe column 3) → our canonical tag, or ``None``.

    Delegates to the **same pair the current-release pass uses** (``_base_species``
    then ``normalize_species``) rather than carrying its own table. This module used
    to hold a three-entry prefix map — human, mouse, rhesus — written when that was
    the whole database. The scope grew to 30-odd species and the map did not follow,
    so the history layer stayed silently three-species: no former names, no archived
    records and no version chain for dog, cow, sheep, zebrafish and the rest, even
    though the snapshots carry them. Same class of omission as the IG-only locus
    filter this file also had.

    ``_base_species`` takes the first two words after dropping IMGT's strain suffix,
    which is what makes the awkward fields work: ``Gorilla gorilla gorilla_Kamilah``
    → GORILLA, ``Mus musculus castaneus`` → MOUSE, ``Ovis aries_Rambouillet`` →
    SHEEP. Subspecies therefore fold into the parent species — which is not a new
    decision here, it is agreement with how the current data was already loaded.
    """
    try:
        return normalize_species(_base_species(imgt_species))
    except SpeciesError:
        return None


# One IMGT/GENE-DB release per year (earliest available each year), spanning the
# renomenclature history — a broad spread for recovering former/retired official
# names by sequence hash. Directory names from JamieHeather/genedb-releases.
_RELEASES = (
    "2013-01-20_WBM_201311-0", "2015-03-19_WBM_201511-0",
    "2016-12-19_WBM_201821-0", "2018-09-27_WBM_201900-0",
    "2019-03-01_WBM_201908-0", "2020-01-31_WBM_202001-0",
    "2021-12-27_JH-GeneDL_202151-7", "2022-01-17_JH-GeneDL_20222-7",
    "2023-02-20_GENEDB_20237-7", "2024-01-07_GENEDB_20241-7",
    "2025-01-05_GENEDB_20251-7", "2026-01-04_GENEDB_20261-7",
)
_RAW_BASE = (
    "https://raw.githubusercontent.com/JamieHeather/genedb-releases/main/releases/"
    "{dirname}/{fname}"
)
# The reference-FASTA filename prefix changed over time (early releases drop the
# "IMGT"); try both. Ungapped, F+ORF+all-pseudogenes nucleotide set.
_FASTA_NAMES = (
    "IMGTGENEDB-ReferenceSequences.fasta-nt-WithoutGaps-F%2BORF%2BallP",
    "GENEDB-ReferenceSequences.fasta-nt-WithoutGaps-F%2BORF%2BallP",
)


def fetch_imgt_history(snapshot_dir: str | Path, *, releases=_RELEASES) -> list[Path]:
    """Download the curated IMGT/GENE-DB release snapshots into ``snapshot_dir`` as
    ``<release-date>.fasta``. Idempotent — skips any already present. Network."""
    import urllib.error
    import urllib.request

    out = Path(snapshot_dir)
    out.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    for dirname in releases:
        dest = out / f"{dirname[:10]}.fasta"  # date prefix is the release label
        if not dest.exists() or dest.stat().st_size == 0:
            err: Exception | None = None
            for fname in _FASTA_NAMES:
                try:
                    urllib.request.urlretrieve(
                        _RAW_BASE.format(dirname=dirname, fname=fname), dest
                    )
                    err = None
                    break
                except urllib.error.HTTPError as exc:
                    err = exc
            if err is not None:
                raise err
        paths.append(dest)
    return paths


def ingest_imgt_history(session: Session, snapshot_dir: str | Path) -> dict[str, int]:
    """Ingest every ``*.fasta`` snapshot in ``snapshot_dir`` (release label = file
    stem). Idempotent, so safe to re-run on every reload. A missing/empty directory
    is a no-op (returns zero counts) — former names just aren't recovered."""
    sd = Path(snapshot_dir)
    agg = {"releases": 0, "names_added": 0, "new_sequences": 0}
    if not sd.is_dir():
        return agg
    for fasta in sorted(sd.glob("*.fasta")):
        c = ingest_imgt_release(session, fasta, release=fasta.stem)
        agg["releases"] += 1
        agg["names_added"] += c["names_added"]
        agg["new_sequences"] += c["new_sequences"]
    return agg


def mark_archived(session: Session) -> dict[str, int]:
    """Flag coding sequences whose names come **only** from historical IMGT release
    snapshots — gone from every current source — as ``archived``, and record their
    first/last-seen release dates as ``history`` annotations.

    Order-independent and idempotent: it re-derives the flag and dates from scratch
    each run (reset all, then set), so it's correct regardless of which release
    first created a sequence and survives reloads. These archived alleles are then
    excluded from current coverage stats and default search/similarity results, but
    stay reachable directly and surface their similar *current* alleles.
    """
    from ..db.models import Name, Sequence, SequenceAnnotation, Source

    hist = dict(session.execute(
        select(Source.id, Source.version).where(Source.note.like("historical IMGT%"))
    ).all())
    session.execute(update(Sequence).values(archived=False))  # reset, then re-derive
    if not hist:
        session.commit()
        return {"archived": 0}
    hist_ids = list(hist)
    has_hist = select(Name.sequence_id).where(Name.source_id.in_(hist_ids))
    has_current = select(Name.sequence_id).where(Name.source_id.notin_(hist_ids))
    archived_ids = [i for (i,) in session.execute(
        select(Sequence.id).where(
            Sequence.category == "coding",
            Sequence.id.in_(has_hist),
            Sequence.id.notin_(has_current),
        )
    )]
    if not archived_ids:
        session.commit()
        return {"archived": 0}
    session.execute(
        update(Sequence).where(Sequence.id.in_(archived_ids)).values(archived=True)
    )
    # first/last seen release date, from the historical sources naming each sequence
    seen: dict[int, list[str]] = {}
    for sid, src_id in session.execute(
        select(Name.sequence_id, Name.source_id).where(
            Name.sequence_id.in_(archived_ids), Name.source_id.in_(hist_ids)
        )
    ).all():
        if hist.get(src_id):
            seen.setdefault(sid, []).append(hist[src_id])
    session.execute(delete(SequenceAnnotation).where(  # idempotent: rewrite history
        SequenceAnnotation.domain == "history",
        SequenceAnnotation.sequence_id.in_(archived_ids),
    ))
    src0 = hist_ids[0]
    rows: list[dict] = []
    for sid, vers in seen.items():
        vers = sorted(vers)
        rows.append({"sequence_id": sid, "source_id": src0, "domain": "history",
                     "key": "first_seen", "value_text": vers[0]})
        rows.append({"sequence_id": sid, "source_id": src0, "domain": "history",
                     "key": "last_seen", "value_text": vers[-1]})
    if rows:
        session.execute(insert(SequenceAnnotation), rows)
    session.commit()
    return {"archived": len(archived_ids)}


def ingest_imgt_release(
    session: Session, fasta_path: str | Path, *, release: str,
) -> dict[str, int]:
    """Ingest one historical IMGT/GENE-DB FASTA snapshot (multi-species), attaching
    its (possibly former) official allele names to the matching sequences.

    ``release`` is the IMGT release label (date or number) recorded as the source
    version — so a name that no longer matches the current release is a *former*
    official name, recoverable from its source version. Each record's species is
    read from the FASTA header (pipe column 3). Idempotent: re-running attaches the
    same names with zero new rows.

    Returns counts: ``names_added`` (attached to existing or new sequences),
    ``new_sequences`` (present in this release but not in the current data — i.e.
    retired alleles), and skip tallies.
    """
    source = loader.get_or_create_source(
        session, name="IMGT", version=str(release), curie_prefix="IMGT/GENE-DB",
        note=f"historical IMGT/GENE-DB reference (release {release})",
    )
    names_added = new_sequences = skipped_non_v = skipped_species = 0
    for header, seq in read_fasta(fasta_path):
        parts = header.split("|")
        sp_key = _species_key(parts[2]) if len(parts) > 2 else None
        if sp_key is None:
            skipped_species += 1
            continue
        rec = parse_record(header, seq)
        if not rec.allele or not rec.sequence:
            continue
        # Same admission rule as the current-release pass (imgt_genedb.ingest_reference):
        # the whole-unit labels, IG **and** TR. Matching it is the point — history
        # reconciles to current data by sequence hash, so a record admitted by one pass
        # and refused by the other can never meet its own counterpart.
        #
        # TR was excluded here while every other pass already covered it, which made
        # the history layer silently IG-only: a renamed or re-cut TRBV allele got no
        # former name, no archived record and no "Earlier records" entry, even though
        # the snapshots carry it and the whole downstream chain (mark_archived,
        # prev_versions, the allele page) is segment- and locus-agnostic. Reported by
        # Jamie Heather (TRBV15*02 extended after publication; TRBV10-2*03 added then
        # withdrawn once *02 was completed, so *03 is properly a former name of *02).
        #
        # C is deliberately absent, and it is the *assembly* that excludes it, not a
        # policy: IMGT publishes an IGH constant per domain (CH1 / H / CH2 / CH3-CHS /
        # EX1… — ~40 labels across the snapshots), while we store the assembled
        # secreted region as one sequence. A historical CH1 record therefore cannot
        # hash-match its own current allele; ingesting it would mint a spurious
        # archived "C allele" that is really one domain. Only C-REGION (IGK/IGL/IGI,
        # published whole) could reconcile, and taking just those would put a new
        # asymmetry inside C. See ingest/imgt_constant.py.
        # TODO(history-c): revisit if the per-domain composition is ever hashed too.
        label = parts[4].strip() if len(parts) > 4 else ""
        segment = segment_of(rec.allele)
        if (label not in _HISTORY_LABELS
                or segment not in _HISTORY_SEGMENTS
                or locus_of(rec.allele) not in _CODING_LOCI):
            skipped_non_v += 1
            continue
        tag = sp_key  # already canonical: _species_key normalises
        digest = hash_sequence(rec.sequence)
        seq_row = loader.get_coding(session, tag, digest)
        if seq_row is None:  # present in this release but not current → retired
            seq_row = loader.upsert_coding(
                session, species=tag, segment=segment,
                locus=locus_of(rec.allele), nt=rec.sequence,
            )
            new_sequences += 1
        loader.add_name(session, seq_row, source, rec.allele, scheme_label(rec.allele))
        loader.add_attribution(
            session, seq_row, source,
            original_id=curie_build("IMGT/GENE-DB", rec.allele),
            original_sequence=rec.sequence,
        )
        # AA is V-only by design (decision #6): D, J and C carry no protein.
        if segment == "V":
            loader.set_v_aa(session, seq_row)
        names_added += 1
    return {
        "names_added": names_added, "new_sequences": new_sequences,
        "skipped_non_v": skipped_non_v, "skipped_species": skipped_species,
    }
