"""Constant-region (C) alleles from IMGT/GENE-DB.

IMGT does not publish a heavy-chain constant region as one record. It publishes the
**domains** — ``CH1``, the hinge (``H``, or ``H1``..``H5`` where the hinge is
multi-exon), ``CH2``, ``CH3``, sometimes ``CH4``..``CH10``, and the secreted tail
``CHS`` (often merged into its predecessor as ``CH3-CHS``). Light chains are the easy
case: ``IGKC``/``IGLC`` ship as a real, whole ``C-REGION`` record.

So heavy-chain C is a **composition**, and the question is whose composition it is.
It is IMGT's: ``IGHG1*04`` (JN582178) has its four domains at 11..304, 305..349,
350..679 and 680..1000 — contiguous — and IMGT publishes for that one allele an
assembled ``C-REGION`` at 11..1000, 990 nt. Concatenating the domains in domain order
reproduces that record exactly. That is the only assembled IGH C-REGION in the whole
reference set, and it is enough: it pins the rule to IMGT rather than to us. Six human
isotypes then land on their published lengths independently — IGHG1 990 (330 aa),
IGHG3 1131 (the four-exon hinge), IGHG4 981, IGHM 1356 (452 aa), IGHE 1284, IGHA1
1059, IGHD 1152 — which is corroboration the single anchor could not give.

**The secreted form, not the membrane form.** Every heavy-chain constant gene also
carries membrane exons (``M``/``M1``/``M2``, plus the transmembrane and cytoplasmic
labels), spliced in place of the secreted tail on the membrane-bound receptor. They
are a different transcript, not a longer version of this one, so they are excluded
here rather than appended. Storing the membrane form as a second representation
(``subtype`` ``c_membrane``) is the natural extension — see ``TODO(c-membrane)`` —
and follows the leader precedent exactly: one sequence per representation, with the
composition recorded so the boundaries stay recoverable.

**No new download.** The domains are already in the bulk ReferenceSequences FASTA
that ``imgt_genedb`` fetches. This module re-reads that file; it never goes to the
network.

**TR constants are deferred, not forgotten.** TRA/TRB/TRD ship as ``EX1``..``EX4``,
which would compose the same way, but ``EX4UTR`` (TRA/TRD) mixes coding sequence with
untranslated tail, and TRG carries ``EX2A``/``EX2B``/``EX2C``/``EX2R``/``EX2T`` for
the duplicated TRGC exons — a real choice about which copy belongs in the region, and
one IMGT does not make for us. Guessing there would be adjudication, which is not
this service's job. See ``TODO(tr-constant)``.
"""

from __future__ import annotations

import logging
from pathlib import Path

from sqlalchemy.orm import Session

from ..curie import build as curie_build
from ..naming import scheme_label
from ..species import SpeciesError, normalize_species
from . import loader

_LOG = logging.getLogger(__name__)

#: Domain order within a secreted constant region. Merged labels sit at the rank of
#: their first component (``CH3-CHS`` where ``CH3`` would be), so a gene that ships
#: ``CH3``+``CHS`` separately and one that ships them merged compose identically.
#: The hinge variants ``H1``..``H5`` are the multi-exon hinges (human IGHG3 has four)
#: and sort between ``CH1`` and ``CH2``. The ``D``-suffixed labels (``CH2D2``) are
#: IMGT's duplicated domains in a few non-human species.
DOMAIN_ORDER: dict[str, int] = {
    "CH1": 10,
    "H": 20, "H1": 21, "H2": 22, "H3": 23, "H4": 24, "H5": 25,
    "H-CH2": 20,
    "CH2": 30, "CH2D": 31, "CH2D2": 32, "CH2D3": 33,
    "CH3": 40, "CH3D": 41, "CH3D2": 42, "CH3D3": 43, "CH3-CHS": 40,
    "CH4": 50, "CH4D": 51, "CH4D2": 52, "CH4D3": 53, "CH4-CHS": 50,
    "CH5": 60, "CH6": 70, "CH7": 80, "CH7-CHS": 80,
    "CH8": 90, "CH9": 100, "CH10": 110, "CHX": 115,
    "CHS": 120,
}

#: Labels that terminate a secreted constant region. Their absence does not make a
#: record unusable, but it does make it partial, and that has to be visible.
_TERMINAL = frozenset({"CHS", "CH3-CHS", "CH4-CHS", "CH7-CHS"})

#: The domain that starts a constant region. Without it there is no constant region
#: to assemble, only a fragment of one.
_FIRST = "CH1"

#: Whole-region records that need no assembly at all (light chains).
_WHOLE = "C-REGION"

#: Loci whose whole ``C-REGION`` records we take directly. IGI is the IgI locus IMGT
#: registers for a few non-mammalian species.
_WHOLE_LOCI = ("IGK", "IGL", "IGI")

#: Loci we assemble from domains.
_DOMAIN_LOCI = ("IGH",)

SUBTYPE = "c_secreted"


def _locus(allele: str) -> str:
    """``IGHG1*01`` -> ``IGH``. The locus is the first three characters of every
    IMGT gene name, which is why this does not go through ``naming.locus_of``: it
    must work before we have decided the name is one we can parse."""
    return allele[:3]


def assemble(domains: dict[str, str]) -> tuple[str, list[tuple[str, int]]] | None:
    """Concatenate ``{label: sequence}`` into the secreted constant region.

    Returns ``(sequence, [(label, length), ...])`` in domain order, or ``None`` when
    ``CH1`` is absent — a record starting at CH2 or CH3 is a fragment, and letting it
    through would put a 300 nt piece beside a 1,019 nt whole region under one segment,
    each looking equally complete.
    """
    usable = {k: v for k, v in domains.items() if k in DOMAIN_ORDER}
    if _FIRST not in usable:
        return None
    ordered = sorted(usable.items(), key=lambda kv: (DOMAIN_ORDER[kv[0]], kv[0]))
    return "".join(v for _, v in ordered), [(k, len(v)) for k, v in ordered]


def is_partial(labels) -> bool:
    """True when the region has no secreted tail, so it stops mid-gene."""
    return not (set(labels) & _TERMINAL)


def ingest_constant(
    session: Session, fasta_path: str | Path, *,
    loci=_DOMAIN_LOCI + _WHOLE_LOCI, release: str = "current",
) -> dict[str, int]:
    """Ingest constant-region alleles from the bulk reference FASTA.

    Two paths, one output shape: IGH is assembled from its domains, IGK/IGL/IGI are
    taken whole from their ``C-REGION`` records. Both land as ``segment='C'``,
    ``subtype='c_secreted'``, with the composition recorded as annotations so nothing
    about how the sequence was built is lost.

    Returns per-species counts keyed on the canonical tag, plus two diagnostic keys:
    ``"partial"`` (assembled, but with no secreted tail) and ``"skipped:no-CH1"``
    (domain records that do not start a region). Unregistered species are counted
    under ``UNMAPPED_PREFIX`` exactly as in :func:`imgt_genedb.ingest_reference` —
    a species IMGT has added must be visible, not silently absent.
    """
    from .fasta import read_fasta, require_fasta
    from .imgt_genedb import UNMAPPED_PREFIX, _base_species, _species_source

    fasta_path = require_fasta(fasta_path)

    # (species_field, allele) -> {label: seq}; plus the header fields we need once.
    domains: dict[tuple[str, str], dict[str, str]] = {}
    meta: dict[tuple[str, str], tuple[str, str]] = {}   # -> (functionality, accession)
    whole: list[tuple[str, str, str, str, str]] = []    # sp, allele, seq, func, acc

    for header, seq in read_fasta(fasta_path):
        p = header.split("|")
        if len(p) < 5 or not seq:
            continue
        acc, allele, sp_field, func, label = p[0], p[1], p[2], p[3], p[4].strip()
        locus = _locus(allele)
        if locus not in loci:
            continue
        if locus in _WHOLE_LOCI and label == _WHOLE:
            whole.append((sp_field, allele, seq, func, acc))
        elif locus in _DOMAIN_LOCI and label in DOMAIN_ORDER:
            key = (sp_field, allele)
            domains.setdefault(key, {})[label] = seq
            meta.setdefault(key, (func, acc))

    counts: dict[str, int] = {}
    unmapped: dict[str, int] = {}
    src_cache: dict[str, object] = {}
    partial = skipped = 0

    def _load(sp_field, allele, seq, func, acc, composition):
        nonlocal partial
        try:
            tag = normalize_species(_base_species(sp_field))
        except SpeciesError:
            name = _base_species(sp_field)
            unmapped[name] = unmapped.get(name, 0) + 1
            return
        row = loader.upsert_coding(
            session, species=tag, segment="C", locus=_locus(allele), nt=seq,
            functional=(True if func.strip().startswith("F") else None),
            subtype=SUBTYPE,
        )
        src = _species_source(session, tag, src_cache)
        loader.add_name(session, row, src, allele, scheme_label(allele))
        loader.add_attribution(
            session, row, src,
            original_id=curie_build("IMGT/GENE-DB", allele), original_sequence=seq,
        )
        if composition is not None:
            # The composition is the audit trail for a sequence we built rather than
            # read. Without it, an assembled region is indistinguishable from one
            # IMGT published whole, and the domain boundaries are unrecoverable.
            loader.add_annotation(
                session, row, src, domain="constant", key="domains",
                value="+".join(f"{lab}:{n}" for lab, n in composition),
            )
            loader.add_annotation(
                session, row, src, domain="constant", key="assembled_from",
                value=str(len(composition)),
            )
            if is_partial(lab for lab, _ in composition):
                loader.add_annotation(
                    session, row, src, domain="constant", key="partial",
                    value=True,
                )
                partial += 1
        counts[tag] = counts.get(tag, 0) + 1

    for (sp_field, allele), doms in domains.items():
        built = assemble(doms)
        if built is None:
            skipped += 1
            continue
        seq, composition = built
        func, acc = meta[(sp_field, allele)]
        _load(sp_field, allele, seq, func, acc, composition)

    for sp_field, allele, seq, func, acc in whole:
        _load(sp_field, allele, seq, func, acc, None)

    session.commit()
    counts["partial"] = partial
    counts["skipped:no-CH1"] = skipped
    if skipped:
        _LOG.info(
            "IMGT/GENE-DB C: %d allele(s) had constant-region domains but no CH1 "
            "(fragmentary records, mostly pseudogenes) and were not assembled.",
            skipped,
        )
    for name, n in sorted(unmapped.items(), key=lambda kv: (-kv[1], kv[0])):
        _LOG.warning(
            "IMGT/GENE-DB C: %d record(s) skipped for unregistered species %r — "
            "add it to allele_server.species._SPECIES to ingest them.", n, name,
        )
        counts[UNMAPPED_PREFIX + name] = n
    return counts


# TODO(c-membrane): the membrane-bound form (…CH3 + M1 + M2, replacing the secreted
# tail) is a separate transcript of the same gene. Add it as subtype 'c_membrane'
# linked to the secreted form, rather than lengthening this one.
# TODO(tr-constant): TRA/TRB/TRD EX1..EX4 compose the same way, but EX4UTR mixes CDS
# with UTR and TRG's EX2A/B/C/R/T duplicates need a rule IMGT does not supply.
