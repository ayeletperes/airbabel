"""Allele name parsing and scheme detection across naming schemes.

The same parsing as the HUSA curation pipeline (``husa_core``), so aliases bucket
the same way in both. Names are labels, never keys; always carry the scheme.
Parse into ``(segment, family, gene, allele, suffix)``, zero-pad the allele to
two digits, preserve the original string verbatim.

Schemes: ``IUIS`` (canonical), ``VDJBASE`` (``_<variant>`` suffix), ``IGDISCOVER``
(``_S<n>``), ``OGRDB`` (``*i<n>`` inferred marker), ``IGLABEL``/``HUSA`` (4-char
ASC-style gene label), ``ASC`` (bare 4-char code), ``NCBI`` (GenBank accession).
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from enum import StrEnum

__all__ = [
    "AlleleName",
    "AlleleNameError",
    "Scheme",
    "canonicalise",
    "detect_scheme",
    "parse_name",
    "scheme_for_vdjbase_name",
]


class Scheme(StrEnum):
    IUIS = "iuis"
    VDJBASE = "vdjbase"
    IGDISCOVER = "igdiscover"
    OGRDB = "ogrdb"
    HUSA = "husa"
    IGLABEL = "iglabel"
    ASC = "asc"
    NCBI = "ncbi"


class AlleleNameError(ValueError):
    """Raised when a name cannot be parsed under the requested scheme."""


_NAME_RE: re.Pattern[str] = re.compile(
    # IG (IGH/IGK/IGL) and TR (TRA/TRB/TRG/TRD) loci, each V/D/J/C.
    r"^(?P<segment>(?:IG[HKL]|TR[ABGD])[VDJC])"
    # Family is normally numeric, but the TRG J-P genes are not: IMGT names them
    # TRGJP, TRGJP1 and TRGJP2 -- a letter where a digit is otherwise required. With
    # a digits-only family they failed to parse at all, so segment_of and locus_of
    # both returned None and EVERY ingest dropped them silently, IMGT's included: the
    # database held zero TRGJP alleles from any source. They are real functional
    # joining genes (TRGJP is preferentially used in Vgamma9Vdelta2 T cells), so this
    # was data loss, not filtering. Found because the KI-TCR ingest reports its skips.
    r"(?P<family>P?[0-9]+(?:/[0-9]+)?[A-Z]?|P)"
    r"(?:-(?P<gene>[A-Z0-9]+(?:-[A-Z0-9]+)*))?"
    # IMGT dual TRA/TRD designation for the shared V genes rearranged in both loci,
    # e.g. TRAV14/DV4, TRAV15-1/DV6-1, TRAV15D-1/DV6D-1. Primary locus stays TRA
    # (from <segment>); the second locus (TRD) is recovered by dual_locus_of().
    r"(?P<dual>/DV[0-9]+[A-Z]?(?:-[0-9]+)?)?"
    r"\*(?P<allele>i?[0-9]+)"
    r"(?:_(?P<suffix>.+))?$"
)

# Constant genes that ``_NAME_RE`` cannot reach, because IMGT does not name them
# with a "C" segment letter plus a family number. Two shapes:
#
#   * IGH is named by **isotype**: IGHM, IGHD, IGHE, IGHG1..IGHG7, IGHA1, IGHGP,
#     and across species IGHY (avian IgY), IGHZ/IGHT (teleost), IGHO, IGHG2A/2B/2C,
#     IGHG5-1, IGHDD1P, IGHEP1 ...
#   * IGKC / TRAC / TRDC carry the C but no family number, which ``_NAME_RE``
#     requires. (IGLC2, TRBC1 do have one and parse there already.)
#
# **Tried only after ``_NAME_RE`` fails**, and that ordering is what disambiguates
# IgD-the-isotype from D-the-diversity-segment: ``IGHD3-10*01`` matches ``_NAME_RE``
# (family 3, gene 10) and stays a D segment, while ``IGHD*01`` has no family, fails
# there, and lands here as the constant delta gene. Verified against the bulk
# reference FASTA: zero D-REGION records are named ``IGHD*NN``, so nothing is
# reclassified by this rule.
_CONSTANT_RE: re.Pattern[str] = re.compile(
    r"^(?P<locus>IG[HKL]|TR[ABGD])"
    r"(?P<isotype>C|[ADEGMOTYZ][A-Z0-9]*(?:-[0-9]+)?)"
    r"\*(?P<allele>i?[0-9]+)"
    r"(?:_(?P<suffix>.+))?$"
)

_HUSA_LABEL_RE: re.Pattern[str] = re.compile(r"^[A-Z0-9]{4}$")
_NCBI_ACCESSION_RE: re.Pattern[str] = re.compile(r"^[A-Z]{1,3}_?[0-9]+(?:\.[0-9]+)?$")
_IGDISCOVER_SUFFIX_RE: re.Pattern[str] = re.compile(r"_S\d+", re.IGNORECASE)
_OGRDB_INFERRED_RE: re.Pattern[str] = re.compile(r"\*i\d+", re.IGNORECASE)


@dataclass(frozen=True, slots=True)
class AlleleName:
    """Parsed, canonicalised allele name. ``raw`` preserves the exact input."""

    segment: str
    family: str
    gene: str | None
    allele: str
    suffix: str | None
    scheme: Scheme
    raw: str
    dual: str | None = None  # IMGT TRA/TRD dual designation, e.g. "/DV4"

    def canonical(self) -> str:
        """Rebuild the canonical string form from the parsed components."""
        if self.scheme in (Scheme.ASC, Scheme.NCBI):
            return self.allele
        base = f"{self.segment}{self.family}"
        if self.gene is not None:
            base = f"{base}-{self.gene}"
        if self.dual is not None:
            base = f"{base}{self.dual}"  # already includes the leading '/'
        base = f"{base}*{self.allele}"
        if self.suffix is not None:
            base = f"{base}_{self.suffix}"
        return base


def _pad_allele(digits: str) -> str:
    prefix = ""
    rest = digits
    if rest.startswith("i"):
        prefix = "i"
        rest = rest[1:]
    stripped = rest.lstrip("0") or "0"
    return prefix + stripped.zfill(2)


def _match(raw: str) -> re.Match[str]:
    if not isinstance(raw, str):
        raise AlleleNameError(f"name must be str, got {type(raw).__name__}")
    stripped = raw.strip()
    if not stripped:
        raise AlleleNameError("name is empty")
    m = _NAME_RE.match(stripped)
    if m is None:
        raise AlleleNameError(f"name does not match any known allele grammar: {raw!r}")
    return m


def parse_name(raw: str, scheme: Scheme) -> AlleleName:
    """Parse a raw allele name under the requested scheme."""
    if scheme is Scheme.ASC:
        stripped = raw.strip()
        if not _HUSA_LABEL_RE.match(stripped):
            raise AlleleNameError(
                f"ASC name must be a 4-char alphanumeric cluster code: {raw!r}"
            )
        return AlleleName("", "", None, stripped, None, scheme, raw)
    if scheme is Scheme.NCBI:
        stripped = raw.strip()
        if not stripped:
            raise AlleleNameError("NCBI accession is empty")
        if not _NCBI_ACCESSION_RE.match(stripped):
            raise AlleleNameError(f"not a GenBank accession: {raw!r}")
        return AlleleName("", "", None, stripped, None, scheme, raw)
    m = _match(raw)
    suffix = m["suffix"]
    gene = m["gene"]
    allele_raw = m["allele"]
    if scheme is Scheme.IUIS and suffix is not None:
        raise AlleleNameError(f"IUIS names cannot carry a '_' variant suffix: {raw!r}")
    iglabel_like = scheme in (Scheme.HUSA, Scheme.IGLABEL)
    if (scheme is Scheme.IUIS or iglabel_like) and allele_raw.startswith("i"):
        raise AlleleNameError(
            f"{scheme.value.upper()} names cannot carry the VDJbase "
            f"inferred-allele 'i' marker: {raw!r}"
        )
    if iglabel_like:
        if gene is None:
            raise AlleleNameError(
                f"{scheme.value.upper()} name requires a gene/label after "
                f"the family: {raw!r}"
            )
        if suffix is not None:
            raise AlleleNameError(
                f"{scheme.value.upper()} names cannot carry a '_' variant "
                f"suffix: {raw!r}"
            )
    return AlleleName(
        segment=m["segment"],
        family=m["family"],
        gene=gene,
        allele=_pad_allele(allele_raw),
        suffix=suffix,
        scheme=scheme,
        raw=raw,
        dual=m["dual"],
    )


def detect_scheme(raw: str) -> Scheme | None:
    """Guess the naming scheme for a raw allele name from grammar alone."""
    if not isinstance(raw, str):
        return None
    stripped = raw.strip()
    if not stripped:
        return None
    m = _NAME_RE.match(stripped)
    if m is None:
        return None
    gene = m["gene"]
    suffix = m["suffix"]
    if gene is not None and _HUSA_LABEL_RE.match(gene) and suffix is None:
        return Scheme.IGLABEL
    if _IGDISCOVER_SUFFIX_RE.search(stripped):
        return Scheme.IGDISCOVER
    if _OGRDB_INFERRED_RE.search(stripped):
        return Scheme.OGRDB
    if suffix is not None:
        return Scheme.VDJBASE
    return Scheme.IUIS


def scheme_for_vdjbase_name(name: str) -> str:
    """Return the canonical scheme bucket for a VDJbase-style name (grammar, not source)."""
    if _IGDISCOVER_SUFFIX_RE.search(name):
        return "igdiscover"
    if _OGRDB_INFERRED_RE.search(name):
        return "ogrdb"
    if "_" in name:
        return "vdjbase"
    return "iuis"


def scheme_label(name: str, default: str = "imgt") -> str:
    """The scheme bucket to store/display for an allele name. Canonical
    IUIS-nomenclature names are labelled ``imgt`` (IMGT is the de-facto authority
    for that nomenclature); special grammars (``_suffix`` -> vdjbase, ``*i`` ->
    ogrdb, ``_S`` -> igdiscover) keep their detected scheme. Single source of truth
    so every ingest path labels the same name the same way (no iuis/imgt split)."""
    sch = detect_scheme(name)
    if sch is None or sch is Scheme.IUIS:
        return default
    return sch.value


def canonicalise(raw: str, scheme: Scheme | None = None) -> AlleleName:
    """Parse a raw name, detecting the scheme if one is not supplied."""
    resolved = scheme if scheme is not None else detect_scheme(raw)
    if resolved is None:
        raise AlleleNameError(f"cannot detect scheme for name: {raw!r}")
    return parse_name(raw, resolved)


def segment_of(raw: str) -> str | None:
    """Return the single-letter V/D/J/C segment for a name, or None.

    Convenience for ingestion: ``IGHV1-2*02`` -> ``V``, ``IGHG1*01`` -> ``C``.
    Returns ``None`` for ASC/NCBI-style names that carry no IG grammar.
    """
    if not isinstance(raw, str):
        return None
    stripped = raw.strip()
    m = _NAME_RE.match(stripped)
    if m is not None:
        return m["segment"][-1]
    return "C" if _CONSTANT_RE.match(stripped) else None


def locus_of(raw: str) -> str | None:
    """Return the locus (``IGH``/``IGK``/``IGL``) for a name, or None.

    ``IGHV1-2*02`` -> ``IGH``, ``IGHG1*01`` -> ``IGH``. The AIRR ``locus`` field.
    Returns ``None`` for names without IG grammar.
    """
    if not isinstance(raw, str):
        return None
    stripped = raw.strip()
    m = _NAME_RE.match(stripped)
    if m is not None:
        return m["segment"][:3]
    c = _CONSTANT_RE.match(stripped)
    return c["locus"] if c else None


_SUBGROUP_RE: re.Pattern[str] = re.compile(
    # The head of _NAME_RE and nothing else: segment + family, no gene, no allele.
    r"^(?:IG[HKL]|TR[ABGD])[VDJC][0-9]+(?:/[0-9]+)?[A-Z]?$"
)


def as_subgroup(raw: str) -> str | None:
    """Return the canonical subgroup if ``raw`` *is* a bare subgroup, else None.

    IMGT calls this the subgroup and :class:`AlleleName` calls it ``family``:
    ``IGHV1`` is the subgroup of ``IGHV1-69*01``. A bare subgroup is not a name —
    :func:`parse_name` rejects it, because ``_NAME_RE`` requires ``*<allele>`` — so
    it needs its own recogniser rather than a loosened one.

    Strict by design: ``IGHV1-69`` is a gene, not a subgroup, and returns None. That
    keeps subgroup search from quietly reinterpreting a precise gene query.
    """
    if not isinstance(raw, str):
        return None
    stripped = raw.strip().upper()
    return stripped if _SUBGROUP_RE.match(stripped) else None


def subgroup_of(raw: str) -> str | None:
    """Return the subgroup a full allele name belongs to, or None.

    ``IGHV1-69*01`` -> ``IGHV1``. The inverse direction from :func:`as_subgroup`:
    that one asks "is this a subgroup?", this one asks "which subgroup is this in?".
    """
    m = _NAME_RE.match(raw.strip()) if isinstance(raw, str) else None
    if m is None:
        return None
    return f"{m['segment']}{m['family']}"


def dual_locus_of(raw: str) -> str | None:
    """Return the SECOND locus of an IMGT dual TRA/TRD V gene, or None.

    The shared V genes (IMGT ``TRAV14/DV4*01`` and similar) rearrange in both the
    TRA and TRD loci. :func:`locus_of` returns the primary locus (``TRA`` from the
    ``TRAV`` segment); this returns the second (``TRD``) so the allele can be made
    discoverable under both. ``None`` for ordinary single-locus names.
    """
    m = _NAME_RE.match(raw.strip()) if isinstance(raw, str) else None
    if m is None or not m["dual"]:
        return None
    return "TRD"
