"""Direct IMGT/GENE-DB ingestion across ALL species, straight from imgt.org.

Two sources, both from IMGT/GENE-DB:
  * **Bulk V/D/J coding** — the ``ReferenceSequences`` FASTA (one download; every
    species; V-/D-/J-REGION only — it carries no leader or RSS).
  * **Leader + RSS** — recovered per species/locus via **GENElect** label queries,
    matched to the coding allele by name: ``L-PART1`` + ``L-PART2`` spliced into
    the signal-peptide leader (V only — a D and a J have no leader), and every
    recombination signal IMGT publishes (``V-RS``, ``5'D-RS``, ``3'D-RS``,
    ``J-RS`` — see ``_RS_LABELS``).

A static FASTA + a templated query URL. IMGT
data is CC BY 4.0 (attribution); redistributable with attribution.
"""
from __future__ import annotations

import hashlib
import json
import logging
import urllib.request
from pathlib import Path

from sqlalchemy import select
from sqlalchemy.orm import Session

from ..curie import build as curie_build
from ..db.models import Name, Sequence
from ..naming import dual_locus_of, locus_of, scheme_label, segment_of
from ..species import SpeciesError, normalize_species, scientific_name
from . import loader

_BULK_URL = (
    "https://www.imgt.org/download/GENE-DB/"
    "IMGTGENEDB-ReferenceSequences.fasta-nt-WithoutGaps-F+ORF+allP"
)
# The same records projected as protein. Narrower than the nt file by design: IMGT
# publishes AA only for F/ORF/in-frame-P, so an out-of-frame pseudogene has nt here
# and no AA anywhere — which is why a V row can legitimately end up 'translated'.
_BULK_AA_URL = (
    "https://www.imgt.org/download/GENE-DB/"
    "IMGTGENEDB-ReferenceSequences.fasta-AA-WithoutGaps-F+ORF+inframeP"
)
# IMGT's one-line release tag for the current GENE-DB build (e.g. ``202631-7``).
_RELEASE_URL = "https://www.imgt.org/download/GENE-DB/RELEASE"
# The only IMGT-gapped nucleotide set IMGT publishes (no allP variant exists).
_BULK_GAPPED_URL = (
    "https://www.imgt.org/download/GENE-DB/"
    "IMGTGENEDB-ReferenceSequences.fasta-nt-WithGaps-F+ORF+inframeP"
)
_GENELECT = (
    "https://www.imgt.org/genedb/GENElect?query=8.1+{gene}"
    "&species={species}&IMGTlabel={label}"
)
# The recombination signal sequences GENElect publishes, one entry per RS:
# (segment, IMGTlabel, subtype, link relation). Every label here was confirmed
# against the live service (Homo sapiens, 2026-08-13): V-RS 470 IGHV records,
# 5'D-RS and 3'D-RS 38 IGHD records each at 28 nt, J-RS 20 IGHJ records.
#
# Only V-RS was ever requested, so the D and J signals were absent for every
# species -- the same reader-narrower-than-source gap that hid the curated D and J
# alleles. IMGT publishes each RS whole rather than as heptamer/spacer/nonamer, so
# unlike the curation-table path there are no component lengths to record; the
# extent is IMGT's and is stored as published.
_RS_LABELS = (
    ("V", "V-RS", "v_rs", "has_rss_vrs"),
    ("D", "5'D-RS", "d_rs_5", "has_rss_drs5"),
    ("D", "3'D-RS", "d_rs_3", "has_rss_drs3"),
    ("J", "J-RS", "j_rs", "has_rss_jrs"),
)

_IG_LOCI = ("IGH", "IGK", "IGL")
_TR_LOCI = ("TRA", "TRB", "TRG", "TRD")
# Coding, leader and RSS recovery all cover IG + TR. Leader/RSS used to stop at IG
# "for now", which was a gap rather than a decision: GENElect serves the TR signals
# too, and TRBD is the clearest demonstration of the 12/23 rule in the whole set --
# TRBD1*01 carries a 28 nt 5'D-RS (12 nt spacer) and a 39 nt 3'D-RS (23 nt), because
# it rearranges to a V on one side and a J on the other.
_CODING_LOCI = _IG_LOCI + _TR_LOCI

_LOG = logging.getLogger(__name__)

#: Key prefix used in :func:`ingest_reference` counts for records whose IMGT species
#: is not in ``allele_server.species._SPECIES``. Species tags are uppercase short
#: tags (``HUMAN``), so a lowercase prefixed key can never collide with a real one.
UNMAPPED_PREFIX = "unmapped:"


def _segment_locus(allele: str, label: str) -> tuple[str | None, str | None]:
    """Segment and locus of an IMGT record, from its name or, failing that, the record.

    ``naming`` does not parse every name IMGT publishes -- orphons
    (``IGHD1/OR15-1a*01``), lettered genes (``IGKV12-e*01``), ``IGLV1S6*01`` -- and
    returns None for them. The record still says what it is: the label gives the
    segment (``V-REGION`` -> ``V``) and the name's first three letters the locus.
    Without this fallback those alleles were dropped without a count.
    """
    return segment_of(allele) or label[:1] or None, locus_of(allele) or allele[:3] or None


def _base_species(imgt_field: str) -> str:
    """'Mus musculus_BALB/c' -> 'Mus musculus' (drop the IMGT strain suffix)."""
    return " ".join(imgt_field.split("_")[0].split()[:2])


def _species_source(
    session: Session, tag: str, cache: dict, release: str | None = None,
) -> object:
    """A **per-species** IMGT/GENE-DB source. Generic allele names (e.g.
    ``IGHV1-1*01``) are reused across many species, and ``names`` is unique on
    ``(name, scheme, source_id)`` — so one shared source would let only the first
    species claim each name. A source per species keeps the name uniqueness
    per-species while still grouping as "IMGT" in coverage (grouped by source name).
    """
    if tag not in cache:
        cache[tag] = loader.get_or_create_source(
            session, name="IMGT", version=f"GENE-DB:{tag}", curie_prefix="IMGT/GENE-DB",
            note=f"IMGT/GENE-DB reference (all species); release {release or 'unknown'}",
        )
    return cache[tag]


def _read_fasta_text(text: str) -> list[tuple[str, str]]:
    """Parse FASTA records from a text blob -> [(header_without_'>', seq), ...]."""
    out: list[tuple[str, str]] = []
    h, s = None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if h is not None:
                out.append((h, "".join(s)))
            h, s = line[1:].strip(), []
        elif h is not None:
            s.append(line.strip())
    if h is not None:
        out.append((h, "".join(s)))
    return out


def _release_path(fasta_path: str | Path) -> Path:
    return Path(f"{fasta_path}.release")


def read_release(fasta_path: str | Path) -> str | None:
    """The IMGT release recorded beside the bulk FASTA when it was downloaded, or None."""
    p = _release_path(fasta_path)
    return (p.read_text().strip() or None) if p.is_file() else None


def fetch_reference(dest: str | Path) -> Path:
    """Download the bulk ungapped ReferenceSequences FASTA (idempotent).

    The release tag is read in the same breath and written beside the file
    (``<dest>.release``): IMGT serves only the current build, so the tag is
    knowable at download time and never afterwards.
    """
    dest = Path(dest)
    if not dest.exists() or dest.stat().st_size == 0:
        with urllib.request.urlopen(_RELEASE_URL, timeout=60) as resp:  # noqa: S310
            release = resp.read().decode().strip()
        urllib.request.urlretrieve(_BULK_URL, dest)
        _release_path(dest).write_text(release + "\n")
    return dest


def fetch_reference_gapped(dest: str | Path) -> Path:
    """Download the bulk IMGT-gapped nucleotide FASTA (idempotent).

    Not ingested: it is the V numbering reference the query server reads for diff &
    align (``queries._imgt_gapped_refs``), so it is fetched with the rest of GENE-DB.
    """
    dest = Path(dest)
    if not dest.exists() or dest.stat().st_size == 0:
        urllib.request.urlretrieve(_BULK_GAPPED_URL, dest)
    return dest


def fetch_reference_aa(dest: str | Path) -> Path:
    """Download the bulk ungapped **amino-acid** ReferenceSequences FASTA (idempotent)."""
    dest = Path(dest)
    if not dest.exists() or dest.stat().st_size == 0:
        urllib.request.urlretrieve(_BULK_AA_URL, dest)
    return dest


def ingest_reference_aa(
    session: Session, fasta_path: str | Path, *, loci=_CODING_LOCI,
) -> dict[str, int]:
    """Attach IMGT's **published** V-REGION protein to the V alleles we hold.

    Without this pass every AA in the database is our own translation, and
    ``aa_origin`` is NULL — meaning "we do not know which path produced this", which
    is the honest value but not a useful one. ``loader.set_v_aa`` already knows how to
    reconcile a sourced protein with the nucleotides (and what to do when IMGT
    publishes a reference protein its own nucleotides do not encode: keep ours, record
    theirs as an ``aa/imgt_reference_aa`` annotation). It was simply never given the
    sourced side on the all-species path — only the removed per-species airrdb
    route supplied it.

    Matched by ``(species, allele name)`` rather than by hash: the AA file is a
    different projection of the same records, so the name is the join, and the
    reconciliation against our stored nucleotides happens inside ``set_v_aa``.

    Idempotent. Returns counts keyed on the resulting ``aa_origin`` value, plus
    ``"unmatched"`` for AA records naming an allele we do not hold.
    """
    from .fasta import read_fasta, require_fasta

    fasta_path = require_fasta(fasta_path)
    counts: dict[str, int] = {}
    src_cache: dict[str, object] = {}
    # (species_tag, name) -> Sequence, for every V allele we hold under an IMGT name.
    index: dict[tuple[str, str], Sequence] = {}
    for row, name in session.execute(
        select(Sequence, Name.name)
        .join(Name, Name.sequence_id == Sequence.id)
        .where(Sequence.category == "coding", Sequence.segment == "V")
    ).all():
        index.setdefault((row.species, name), row)

    for header, aa in read_fasta(fasta_path):
        p = header.split("|")
        if len(p) < 5 or not aa or p[4].strip() != "V-REGION":
            continue
        allele, sp_field = p[1], p[2]
        if _segment_locus(allele, "V-REGION")[1] not in loci:
            continue
        try:
            tag = normalize_species(_base_species(sp_field))
        except SpeciesError:
            continue
        seq = index.get((tag, allele))
        if seq is None:
            counts["unmatched"] = counts.get("unmatched", 0) + 1
            continue
        loader.set_v_aa(
            session, seq, aa_source=aa,
            source=_species_source(session, tag, src_cache),
        )
        key = seq.aa_origin or "unset"
        counts[key] = counts.get(key, 0) + 1
    session.commit()
    return counts


def ingest_reference(
    session: Session, fasta_path: str | Path, *,
    segments=("V", "D", "J"), loci=_CODING_LOCI,
) -> dict[str, int]:
    """Ingest IG + TR V/D/J coding sequences for every species from the bulk FASTA.

    Header = ``accession|allele|species|functionality|label|...``. Keeps only
    V-/D-/J-REGION of the requested loci (IG + TR by default); species mapped
    through the registry (strain suffix stripped).

    Returns per-species counts keyed on the canonical tag (``HUMAN``). Records whose
    IMGT species is **not** in the registry cannot be ingested — the tag is part of the
    UID, and registration is deliberately a code change (``species._SPECIES``), never
    something ingest may invent. They are therefore not loaded, but they are **not**
    silently dropped either: each unmapped scientific name is counted under
    ``"unmapped:<Scientific name>"`` in the returned dict and logged as a warning. A
    new species must not abort the load (that would block every other species), but
    "IMGT added a species" has to be visible — it is exactly the event change tracking
    needs to see.
    """
    src_cache: dict[str, object] = {}
    counts: dict[str, int] = {}
    unmapped: dict[str, int] = {}
    from .fasta import parse_record, read_fasta

    release = read_release(fasta_path)
    if release is None:
        _LOG.warning(
            "IMGT/GENE-DB: no release recorded beside %s; the source will say "
            "'release unknown'. Re-download it (delete the file) to record one.",
            fasta_path,
        )
    for header, seq in read_fasta(fasta_path):
        p = header.split("|")
        if len(p) < 5 or not seq:
            continue
        allele, sp_field, label = p[1], p[2], p[4].strip()
        if label not in ("V-REGION", "D-REGION", "J-REGION"):
            continue
        seg, locus = _segment_locus(allele, label)
        if seg not in segments or locus not in loci:
            continue
        try:
            tag = normalize_species(_base_species(sp_field))
        except SpeciesError:
            name = _base_species(sp_field)
            unmapped[name] = unmapped.get(name, 0) + 1
            continue
        rec = parse_record(header, seq)
        coding = loader.upsert_coding(
            session, species=tag, segment=seg, locus=locus, nt=seq,
            functional=rec.functional, secondary_locus=dual_locus_of(allele),
        )
        # One sequence can carry several IMGT names with different calls (IGHV1-116*01
        # P, IGHV1-119*01 F). Functional if IMGT calls it F under any of them, so the
        # result does not depend on which record the file lists first.
        if rec.functional and coding.functional is False:
            coding.functional = True
        src = _species_source(session, tag, src_cache, release)
        loader.add_name(session, coding, src, allele, scheme_label(allele))
        # The GenBank accession is the primary public record (it drives the GenBank
        # xref on the record page); the IMGT name is the fallback when there is none.
        loader.add_attribution(
            session, coding, src,
            original_id=(curie_build("GENBANK", rec.accession) if rec.accession
                         else curie_build("IMGT/GENE-DB", allele)),
            original_sequence=seq,
        )
        # The raw IMGT tokens, kept because `functional` collapses ORF and P together.
        if rec.functionality:
            loader.add_annotation(session, coding, src, domain="functionality",
                                  key="imgt", value=rec.functionality)
        if rec.label:
            loader.add_annotation(session, coding, src, domain="imgt",
                                  key="label", value=rec.label)
        counts[tag] = counts.get(tag, 0) + 1
    session.commit()
    for name, n in sorted(unmapped.items(), key=lambda kv: (-kv[1], kv[0])):
        _LOG.warning(
            "IMGT/GENE-DB: %d coding record(s) skipped for unregistered species %r — "
            "add it to allele_server.species._SPECIES to ingest them.", n, name,
        )
        counts[UNMAPPED_PREFIX + name] = n
    return counts


def _genelect_cache_path(cache_dir: str | Path, gene: str, species_sci: str, label: str) -> Path:
    key = f"{gene}|{species_sci}|{label}"
    return Path(cache_dir) / f"{hashlib.sha256(key.encode()).hexdigest()}.json"


def genelect(
    gene: str, species_sci: str, label: str, *, timeout: int = 90,
    cache_dir: str | Path | None = None,
) -> dict[str, str]:
    """Query GENElect for one IMGT label; return ``{allele: sequence}``.

    Needs ``requests`` (``verify=False``: IMGT's TLS cert). Returns ``{}`` on no data
    (GENElect cannot distinguish 'no such label' from 'no data for species').

    With ``cache_dir`` the response is memoised on disk, keyed on gene/species/label.
    This is the only pass that hit imgt.org live on every reload (~350 queries), which
    made a rebuild hostage to IMGT being up. An empty ``{}`` **is** cached: it is a real
    answer for most species/locus pairs. A network failure is **not** cached and is left
    to propagate, so a rebuild fails loudly rather than silently losing leaders/RSS.
    """
    cp = _genelect_cache_path(cache_dir, gene, species_sci, label) if cache_dir else None
    if cp is not None and cp.exists():
        try:
            return json.loads(cp.read_text())
        except json.JSONDecodeError:
            pass  # corrupt entry: refetch below and overwrite

    import urllib3
    from requests import get

    urllib3.disable_warnings()
    url = _GENELECT.format(
        gene=gene, species=species_sci.replace(" ", "+"), label=label
    )
    txt = get(url, verify=False, timeout=timeout).text
    out: dict[str, str] = {}
    if "\n>" in txt:
        block = txt[txt.index("\n>") + 1:]
        block = block.split("</pre>")[0].split("</b>")[0]  # FASTA sits inside <pre>
        for h, s in _read_fasta_text(block):
            parts = h.split("|")
            if len(parts) > 1 and s:
                out[parts[1]] = s
    if cp is not None:
        cp.parent.mkdir(parents=True, exist_ok=True)
        cp.write_text(json.dumps(out))
    return out


def ingest_leaders_rss(
    session: Session, species_tags, *, loci=_CODING_LOCI, release: str = "current",
    progress=None, cache_dir: str | Path | None = None,
) -> dict[str, int]:
    """Recover + ingest the leader (L-PART1+L-PART2 spliced) and **every** RSS
    GENElect publishes, matched to their coding allele by name.

    Leader stored as a ``leader`` sequence (subtype ``l_spliced``) linked
    ``has_leader`` with an ``exon1_len`` annotation, and it stays V-only: a D and a
    J have no leader. RSS are stored one per signal (``_RS_LABELS``) as an ``rss``
    sequence linked by its own relation — the V-RS, both D-RS and the J-RS, the same
    shape as the curation tables.

    A label is only requested for a (locus, segment) we actually hold alleles for,
    so no query is spent asking IMGT for the D genes of IGK.
    """
    src_cache: dict[str, object] = {}
    out = {"leaders": 0, "rss": 0, "queries": 0}
    for tag in species_tags:
        sci = scientific_name(tag)
        source = _species_source(session, tag, src_cache)
        # allele-name -> coding sequence row, per segment, for this species
        by_segment: dict[str, dict[str, Sequence]] = {}
        for nm, seq in session.execute(
            select(Name.name, Sequence)
            .join(Sequence, Sequence.id == Name.sequence_id)
            .where(Sequence.species == tag, Sequence.category == "coding",
                   Sequence.segment.in_({s for s, _l, _st, _r in _RS_LABELS}))
        ):
            by_segment.setdefault(seq.segment, {}).setdefault(nm, seq)
        if not by_segment:
            continue
        for locus in loci:
            vgene = locus + "V"
            v_alleles = by_segment.get("V", {})
            if any(a.startswith(vgene) for a in v_alleles):
                l1 = genelect(vgene, sci, "L-PART1", cache_dir=cache_dir)
                l2 = genelect(vgene, sci, "L-PART2", cache_dir=cache_dir)
                out["queries"] += 2
                for allele, coding in v_alleles.items():
                    if not allele.startswith(vgene):
                        continue
                    a, b = l1.get(allele, "").strip(), l2.get(allele, "").strip()
                    if not (a and b):
                        continue
                    leader = loader.upsert_noncoding(
                        session, species=tag, category="leader",
                        subtype="l_spliced", nt=a + b,
                    )
                    loader.add_link(session, coding, leader, "has_leader", source)
                    loader.add_annotation(
                        session, leader, source, domain="leader",
                        key="exon1_len", value=len(a),
                    )
                    out["leaders"] += 1
            for segment, label, subtype, relation in _RS_LABELS:
                gene = locus + segment
                alleles = by_segment.get(segment, {})
                # IGK and IGL have no D gene; asking for one wastes a request and
                # teaches nothing, since GENElect cannot distinguish "no such label"
                # from "no data for this species".
                if not any(a.startswith(gene) for a in alleles):
                    continue
                found = genelect(gene, sci, label, cache_dir=cache_dir)
                out["queries"] += 1
                for allele, coding in alleles.items():
                    if not allele.startswith(gene):
                        continue
                    rs = found.get(allele, "").strip()
                    if not rs:
                        continue
                    rss = loader.upsert_noncoding(
                        session, species=tag, category="rss", subtype=subtype, nt=rs,
                    )
                    loader.add_link(session, coding, rss, relation, source)
                    out["rss"] += 1
        session.commit()
        if progress:
            progress(tag, out)
    return out
