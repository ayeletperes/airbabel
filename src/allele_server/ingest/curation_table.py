"""Ingest the internal curation tables (MUSA rhesus macaque / HUSA human).

Both are the same shape: one row per allele with the coding ``seq``, cross-
reference name columns, and the **RSS** components for every segment. This pass
populates V, D and J coding alleles + their names + the RSS sequences (as
first-class ``rss`` rows linked via ``sequence_links``), so the RSS-variation
network works.

**All four RSS are ingested** (``_RSS_PARTS``): the V-RS, both D-RS, and the J-RS.
They are one concept -- heptamer, spacer, nonamer -- differing only in spacer
length (the 12/23 rule, which follows the locus rather than the segment) and in
which way round the parts sit relative to the coding sequence. A D allele has two,
one on each flank. The consequence is that RSS of different kinds come out the same
length -- a D-RS and an IGKV V-RS are both 28 nt -- so they must never share a
similarity pool; see ``queries.comparison_class``.

The leader block stays V-shaped on purpose: it keys on the ``l_part*`` columns,
which are empty on a D or J row (a D and a J have no leader), so those rows
contribute a coding sequence, names and RSS and nothing else.

Leaders (``l_part1``/``l_part2``) are intentionally **not** ingested: the full
leader (exon1+intron+exon2) needs a separate curation (see the plan).
"""

from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import NamedTuple

from sqlalchemy.orm import Session

from ..curie import CurieError
from ..curie import build as curie_build
from ..naming import locus_of, scheme_label
from ..sequence import SequenceError
from ..species import normalize_species
from . import loader

# Per-table name columns: (column, default scheme). An internal set contributes
# ONLY its own allele / IgLabel name — never an IMGT-nomenclature name (decision:
# HUSA/MUSA must not appear as a *source* of an IMGT name; if a sequence also has
# an official IMGT name, that comes from the IMGT source itself). So we deliberately
# drop each table's IMGT-name column (MUSA ``imgt`` / HUSA ``husa``).
#
# The column is ``allele``, NOT MUSA's ``iglabel`` / HUSA's ``iglabel_allele``.
# Those two are **gene**-level IgLabels carrying a ``*00`` placeholder — every
# value in both ends in ``*00``, and there are fewer distinct values than there
# are distinct sequences (HUSA 1,349 labels for 1,353 sequences), so they cannot
# name an allele. Using them published ~5,100 gene labels as allele names.
# ``allele`` is allele-level and exactly one per distinct sequence.
_NAME_COLUMN = "allele"
MUSA_NAMES = ((_NAME_COLUMN, "iglabel"),)
HUSA_NAMES = ((_NAME_COLUMN, "iglabel"),)

# Segments these tables carry, taken from ``gene_type`` (IGHV / IGHD / IGKJ ...).
# The pass was V-only, which silently dropped ~900 curated D and J alleles that
# were sitting in the same files -- 406 in HUSA, 513 in MUSA. The segment must be
# read per row and never hardcoded: a widened filter with a hardcoded segment
# would store a J allele labelled V, which is worse than not loading it at all.
_SEGMENTS = frozenset({"V", "D", "J"})


def _segment_of(gene_type: str) -> str | None:
    """``IGHD`` -> ``"D"``; None for anything this pass does not ingest.

    C is deliberately absent: the tables carry no constant-region rows, and IMGT's
    C comes assembled from per-domain records (``ingest/imgt_constant.py``).
    """
    seg = gene_type[3:4].upper() if len(gene_type) >= 4 else ""
    return seg if seg in _SEGMENTS else None

# Underlying rhesus-macaque germline databases the MUSA set aggregates. Each
# becomes a first-class Source; an allele's name in that column is recorded as a
# per-source alias + a provenance row, so a record shows exactly which databases
# attest it. (column, source name, naming scheme, CURIE prefix.)
MUSA_SOURCES = (
    ("kimdb", "KIMDB", "kimdb", "KIMDB"),
    ("rhgldb+", "RhGLDB+", "rhgldb", "RHGLDB"),
    ("vrc", "VRC", "vrc", "VRC"),
    ("guo", "Guo", "guo", "GUO"),
)


class _RssPart(NamedTuple):
    segment: str
    subtype: str
    relation: str
    components: tuple[tuple[str, str], ...]  # (component name, table column)


# The recombination signal sequences these tables carry, one entry per RS. Both
# tables hold every column for all three segments; only the V-RS was ever read.
#
# An RS is a heptamer, a spacer and a nonamer, and the same three parts in every
# case -- what differs is the spacer length and which way round the parts sit.
# Components are listed in genomic 5'->3' order **as they lie beside the coding
# sequence**, so the stored string reads the way the locus does:
#
#     [ 5'D-RS ] D [ 3'D-RS ]              nonamer-spacer-heptamer, then heptamer-spacer-nonamer
#     V [ V-RS ]                           heptamer-spacer-nonamer, 3' of V
#     [ J-RS ] J                           nonamer-spacer-heptamer, 5' of J
#
# Spacer length follows the 12/23 rule and so varies by locus, not by segment: in
# the tables IGHV/IGLV V-RS and IGHJ/IGKJ J-RS are 23, while IGKV V-RS, IGLJ J-RS
# and both D-RS are 12. Nothing here needs to encode that -- the spacer is taken as
# published -- but it is why a D-RS and an IGKV V-RS come out the same length, and
# therefore why the two must not share a similarity pool
# (``queries.comparison_class``).
_RSS_PARTS: tuple[_RssPart, ...] = (
    _RssPart("V", "v_rs", "has_rss_vrs", (
        ("heptamer", "v_heptamer"), ("spacer", "spacer_3"), ("nonamer", "v_nonamer"),
    )),
    _RssPart("D", "d_rs_5", "has_rss_drs5", (
        ("nonamer", "d_5_nonamer"), ("spacer", "spacer_5"), ("heptamer", "d_5_heptamer"),
    )),
    _RssPart("D", "d_rs_3", "has_rss_drs3", (
        ("heptamer", "d_3_heptamer"), ("spacer", "spacer_3"), ("nonamer", "d_3_nonamer"),
    )),
    _RssPart("J", "j_rs", "has_rss_jrs", (
        ("nonamer", "j_nonamer"), ("spacer", "spacer_5"), ("heptamer", "j_heptamer"),
    )),
)

# Plausible component lengths. Not the canonical 7 / 23-or-12 / 9: these are curated
# tables and a heptamer of 6 or a spacer of 22 is a real observation, not a typo.
# The bounds exist to catch **corruption**, which is present -- a whole 31,361 nt
# sequence sits in one table's `v_heptamer`, and three such rows were ingested as
# "RSS" of 31,382 and 5,405 nt before this guard. Rejections are counted, never
# silent.
_RSS_BOUNDS = {"heptamer": (5, 9), "nonamer": (7, 11), "spacer": (9, 26)}


def _plausible(components) -> bool:
    """True when every component's length is within ``_RSS_BOUNDS``."""
    return all(
        _RSS_BOUNDS[name][0] <= len(value) <= _RSS_BOUNDS[name][1]
        for name, value in components
    )


def _functional(value: str) -> bool | None:
    v = (value or "").strip().lower()
    if not v:
        return None
    return v in ("f", "functional")


# Column widths (migration 0007), mirrored here so an over-long value is refused
# deliberately and counted, instead of aborting a whole table with a DataError
# partway through. These used to read 64 and 128, which was invisible while the name
# column was a 13-character gene label and silently stripped 446 MUSA rows of their
# name the moment it became `allele`: a novel allele there is named by its seed plus
# the mutations separating it from that seed and reaches 784 characters. The head of
# that string is not a substitute -- 649 heads cover 2,153 distinct sequences -- so
# the identifier is stored whole. See 0007 for the measurement.
_MAX_NAME = 1024
_MAX_ORIGINAL_ID = 1024


def _add_names(session, seq, source, raw: str, default_scheme: str) -> int:
    """Attach each comma-separated name. Returns the count refused as over-long."""
    refused = 0
    for token in (raw or "").split(","):
        name = token.split("|")[0].strip()  # drop a |Mus_musculus_<strain> qualifier
        if not name:
            continue
        if len(name) > _MAX_NAME:
            refused += 1
            continue
        loader.add_name(session, seq, source, name, scheme_label(name, default_scheme))
    return refused


def _attribute(session, seq, source, prefix, ident, original_sequence) -> bool:
    """Add a provenance row (source attests seq) keyed by a CURIE. Returns added."""
    ident = (ident or "").strip()
    if not ident:
        return False
    if len(ident) + len(prefix) + 1 > _MAX_ORIGINAL_ID:
        return False
    try:
        loader.add_attribution(
            session, seq, source,
            original_id=curie_build(prefix, ident), original_sequence=original_sequence,
        )
    except CurieError:
        return False
    return True


def ingest_curation_table(
    session: Session,
    path: str | Path,
    *,
    source_name: str,
    species_key: str,
    delimiter: str,
    name_columns,
    genbank_column: str | None = None,
    source_columns=(),
    pubid_column: str | None = None,
    attribution_id_column: str | None = None,
    limit: int | None = None,
) -> dict[str, int]:
    """Ingest a MUSA/HUSA-shaped table. Returns counts.

    ``source_columns`` = iterable of ``(column, source_name, scheme, curie_prefix)``
    for underlying databases this set aggregates (e.g. KIMDB/RhGLDB+ for rhesus):
    each becomes a first-class Source, and a populated cell yields a per-source
    alias + provenance row. ``attribution_id_column`` adds a provenance row for the
    aggregating source itself (e.g. MUSA) on every row; ``pubid_column`` records
    PMIDs as xref annotations.
    """
    tag = normalize_species(species_key)
    source = loader.get_or_create_source(
        session, name=source_name, version=None, curie_prefix=source_name,
        note=f"internal curated allele table ({source_name})",
    )
    # Register the underlying databases once.
    sub_sources = {
        col: (
            loader.get_or_create_source(
                session, name=sname, version=None, curie_prefix=prefix,
                note=f"rhesus macaque germline source ({sname})",
            ),
            scheme, prefix,
        )
        for col, sname, scheme, prefix in source_columns
    }
    stats = {
        "coding": 0, "rss_links": 0, "leader_links": 0, "skipped": 0,
        "attributions": 0,
        # Names the columns cannot hold, and rows where the aggregating source
        # ended up with no provenance row. Reported rather than dropped in
        # silence: a nameless allele is invisible to every search we offer.
        "refused_long_name": 0, "unattributed": 0, "rss_rejected": 0,
    }
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter=delimiter):
            if limit and stats["coding"] >= limit:
                break
            gene_type = (row.get("gene_type") or "").strip()
            segment = _segment_of(gene_type)
            if segment is None:
                continue
            allele = (row.get("allele") or "").strip()
            seq_nt = (row.get("seq") or "").strip()
            try:
                coding = loader.upsert_coding(
                    session, species=tag, segment=segment,
                    locus=locus_of(allele) or gene_type[:3],
                    nt=seq_nt,
                    functional=_functional(row.get("functional", "")),
                )
            except (SequenceError, KeyError):
                stats["skipped"] += 1
                continue
            stats["coding"] += 1

            for col, default_scheme in name_columns:
                stats["refused_long_name"] += _add_names(
                    session, coding, source, row.get(col, ""), default_scheme
                )
            # Aggregating-source provenance on every row (e.g. MUSA = curated set).
            if attribution_id_column and source.curie_prefix:
                if _attribute(session, coding, source, source.curie_prefix,
                              row.get(attribution_id_column) or allele, seq_nt):
                    stats["attributions"] += 1
                else:
                    stats["unattributed"] += 1
            # Underlying databases: per-source alias + provenance.
            for col, (sub_source, scheme, prefix) in sub_sources.items():
                for token in (row.get(col) or "").split(","):
                    name = token.split("|")[0].strip()
                    if not name:
                        continue
                    if len(name) > _MAX_NAME:
                        stats["refused_long_name"] += 1
                        continue
                    loader.add_name(session, coding, sub_source, name, scheme)
                    if _attribute(session, coding, sub_source, prefix, name, seq_nt):
                        stats["attributions"] += 1
            if pubid_column:
                for pid in re.split(r"[,\s]+", (row.get(pubid_column) or "").strip()):
                    if pid.isdigit():
                        loader.add_annotation(
                            session, coding, source, domain="xref",
                            key=f"pmid:{pid}"[:64], value=curie_build("PMID", pid),
                        )
            if genbank_column:
                for acc in re.split(r"[,\s]+", (row.get(genbank_column) or "").strip()):
                    if acc:
                        loader.add_annotation(
                            session, coding, source, domain="xref",
                            key=f"genbank:{acc}"[:64], value=curie_build("GENBANK", acc),
                        )
            notes = (row.get("rss_notes") or "").strip()
            if notes:
                loader.add_annotation(
                    session, coding, source, domain="rss", key="notes", value=notes
                )
            # Recombination signal sequences. One entity per RS, stored complete.
            for part in _RSS_PARTS:
                if part.segment != segment:
                    continue
                comps = [(name, (row.get(col) or "").strip())
                         for name, col in part.components]
                if not all(v for _n, v in comps):
                    continue
                if not _plausible(comps):
                    stats["rss_rejected"] += 1
                    continue
                try:
                    rss = loader.upsert_noncoding(
                        session, species=tag, category="rss", subtype=part.subtype,
                        nt="".join(v for _n, v in comps),
                    )
                except SequenceError:
                    stats["rss_rejected"] += 1
                    continue
                loader.add_link(session, coding, rss, part.relation, source)
                # The component lengths, in order. A heptamer is not always 7 nor a
                # nonamer always 9 in these tables, so the boundaries are not
                # recoverable from the total length -- and the order differs between
                # RS kinds. Same shape as the leader's exon1_len and the constant
                # region's domain list: store the composition, keep the sequence whole.
                loader.add_annotation(
                    session, rss, source, domain="rss", key="composition",
                    value=",".join(f"{n}:{len(v)}" for n, v in comps),
                )
                stats["rss_links"] += 1
            # Spliced leader = exon1 (l_part1) + exon2 (l_part2): the signal-peptide
            # CDS as it appears in the mRNA. Stored complete (not split), mirroring
            # the RSS decision; the exon1 length is kept so the splice boundary is
            # recoverable.
            # TODO(leader-intron): the genomic leader (l_part1 + intron + l_part2) is
            # not in these tables. When a genomic source is wired in (e.g. IMGT
            # L-PART genomic), store it as a second representation
            # (subtype "l_genomic") linked with the same has_leader relation.
            l1 = (row.get("l_part1") or "").strip()
            l2 = (row.get("l_part2") or "").strip()
            if l1 and l2:
                try:
                    leader = loader.upsert_noncoding(
                        session, species=tag, category="leader", subtype="l_spliced",
                        nt=l1 + l2,
                    )
                    loader.add_link(session, coding, leader, "has_leader", source)
                    loader.add_annotation(
                        session, leader, source, domain="leader",
                        key="exon1_len", value=len(l1),
                    )
                    stats["leader_links"] += 1
                except SequenceError:
                    pass
    return stats


def ingest_musa(session: Session, path: str | Path, **kw) -> dict[str, int]:
    # MUSA = Macaca (rhesus macaque) curated set, NOT mouse. It aggregates several
    # underlying rhesus databases (KIMDB / RhGLDB+ / VRC / Guo); each is credited.
    return ingest_curation_table(
        session, path, source_name="MUSA", species_key="rhesus", delimiter=",",
        name_columns=MUSA_NAMES, genbank_column="genbank",
        source_columns=MUSA_SOURCES, pubid_column="pubid",
        attribution_id_column=_NAME_COLUMN, **kw,
    )


def ingest_husa(session: Session, path: str | Path, **kw) -> dict[str, int]:
    # attribution_id_column is what writes the sequence_source row. Without it HUSA
    # contributed sequences and names but was credited for nothing, which left 1,112
    # human alleles -- the HUSA-only novel ones nothing else publishes -- with no
    # source at all, against the "attribution on every record" rule.
    return ingest_curation_table(
        session, path, source_name="HUSA", species_key="human", delimiter="\t",
        name_columns=HUSA_NAMES, genbank_column=None,
        attribution_id_column=_NAME_COLUMN, **kw,
    )
