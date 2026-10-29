"""Bulk data export, gated on each source's declared licence.

IMGT/GENE-DB is **CC BY 4.0** as of the 2026 terms of use
(https://www.imgt.org/about/termsofuse.php), and OGRDB, MUSA, Guo and NCBI are
likewise redistributable with attribution. Export is therefore ENABLED and served
through the API; abuse is controlled by rate-limiting the endpoints, not by hiding
the data behind gated download buttons.

**But not every source is redistributable, and the export must not pretend otherwise.**
Several sets publish no licence at all (HUSA, KIARVA, KIMDB, RhGLDB+, VRC). We may
ingest and credit those — that is what ``docs/source-policy.md`` permits — but
republishing them in a public bulk download is a different act, and one their terms do
not cover. KI-TCR states this explicitly: *"For other uses, please contact the
corresponding authors."*

So :func:`export_alleles` keeps only what a **redistributable** source attests, and
reports the count it withheld rather than silently shrinking. A sequence survives if
*any* redistributable source names or attests it — an allele IMGT also publishes is
IMGT's to give, whatever else corroborates it — and within a surviving record the
names and source credits are likewise restricted, so the export never credits a source
whose contribution it withheld.
"""

from __future__ import annotations

from collections.abc import Iterator

from sqlalchemy import func, select
from sqlalchemy.orm import Session

from . import sources as _sources
from .db.models import Name, Sequence, SequenceSource, Source
from .settings import get_settings
from .species import normalize_species


def _published_unshelved() -> list:
    """Published manifests minus the shelved ones. User-facing attribution and licence
    strings must not name a source that is shelved from the site (e.g. HUSA while its
    manuscript is unpublished), or they leak its existence."""
    shelved = {s.lower() for s in get_settings().shelved_sources}
    return [m for m in _sources.published() if m.name.lower() not in shelved]

# Export is enabled; what it may contain is decided per source, from each manifest's
# declared licence, by `redistributable_sources()` below.
EXPORT_ENABLED = True


def attribution_text() -> str:
    """The attribution line, derived from the source manifests.

    Built from every published source so it cannot drift from what we actually hold.
    The hand-written constant this replaces had already gone stale: it named HUSA and
    KIARVA but not MUSA, VRC, Guo or NCBI, and nothing would have caught that.

    Sources whose licence is undeclared are still *credited* here — attribution is owed
    regardless of redistribution terms — but their data is what
    ``sources.redistributable`` gates.
    """
    published = _published_unshelved()
    if not published:  # pragma: no cover - only reachable with an empty registry
        return "Data aggregated by AIRBabel. Attribute the originating sources on reuse."
    parts = []
    for manifest in published:
        licence = manifest.licence.name
        parts.append(f"{manifest.label} ({licence})" if licence else manifest.label)
    return (
        f"Data aggregated by AIRBabel from {', '.join(parts)}. "
        f"Attribute the originating sources when reusing this data; see /sources for "
        f"each one's citation."
    )


def licence_note() -> str:
    """One line stating the terms the exported data is held under.

    Also derived, for the same reason: the two places that used to state this in the API
    said different things and both were hand-written.

    It reports what the export *contains*. Undeclared-licence sources are now withheld
    by :func:`export_alleles`, so naming them here as "included on the basis of their
    publication" would be false; they are named as withheld instead.
    """
    declared = sorted({
        m.licence.name for m in _published_unshelved() if m.licence.name
    })
    undeclared = sorted(
        m.label for m in _published_unshelved() if m.licence.status == "none_declared"
    )
    parts = []
    if declared:
        parts.append(f"Includes data licensed under {', '.join(declared)}.")
    if undeclared:
        parts.append(
            f"The following publish no explicit licence, so their sequences are "
            f"**withheld from this export** and are available only through the site "
            f"and the per-record API, with credit: {', '.join(undeclared)}."
        )
    parts.append(
        "Attribution to the originating sources is required; see /sources for each "
        "source's citation."
    )
    return " ".join(parts)


def redistributable_sources() -> set[str]:
    """Source *names* whose manifest declares the data redistributable.

    Names, not slugs or labels: this is compared against ``sources.name`` in the
    database (``RhGLDB+``, ``ncbi_nucleotide``), which is what the ingest wrote.
    ``resolve`` handles the spelling differences between the two.

    A source with no manifest is treated as **not** redistributable. That is the safe
    direction and it cannot happen in practice — ``loader.get_or_create_source``
    refuses a source without an approved manifest — but the export is the wrong place
    to assume an invariant holds.
    """
    out: set[str] = set()
    for manifest in _sources.registry().values():
        licence = manifest.licence
        if not (licence and licence.redistributable):
            continue
        out.add(manifest.name)
        out.update(manifest.aliases or ())
        out.add(manifest.label)
    return out


# Module constants for the importers that read them directly. Sources are static per
# process, so computing these once at import is correct and keeps the old call shape.
ATTRIBUTION = attribution_text()
LICENCE_NOTE = licence_note()


def _conds(category: str, species: str | None, segment: str | None, archived: bool):
    conds = [Sequence.category == category]
    if not archived:
        conds.append(Sequence.archived.is_(False))
    if species:
        conds.append(Sequence.species == normalize_species(species))
    if segment:
        conds.append(Sequence.segment == segment)
    return conds


def export_counts(
    session: Session, *, category: str = "coding", species: str | None = None,
    segment: str | None = None, archived: bool = False,
) -> dict[str, int]:
    """``{matched, exported, withheld}`` for the same filters as :func:`export_alleles`.

    A filter that drops records must say how many it dropped; otherwise a shrunken
    export is indistinguishable from a shrunken database. Served on ``/api/export`` so
    the number is visible without diffing two downloads.
    """
    conds = _conds(category, species, segment, archived)
    matched = session.scalar(
        select(func.count()).select_from(Sequence).where(*conds)
    ) or 0
    exported = session.scalar(
        select(func.count()).select_from(Sequence)
        .where(*conds, Sequence.id.in_(_redistributable_ids()))
    ) or 0
    return {"matched": matched, "exported": exported, "withheld": matched - exported}


def _redistributable_ids():
    """Subquery: sequence ids attested by at least one redistributable source.

    Both routes count, because they carry different information and an allele may have
    only one of them: a ``Name`` row means that source *calls* the sequence something,
    a ``SequenceSource`` row means it *contributed* the sequence. Either is a licence
    to redistribute if the source grants one.
    """
    ok = redistributable_sources()
    by_name = (
        select(Name.sequence_id)
        .join(Source, Source.id == Name.source_id)
        .where(Source.name.in_(ok))
    )
    by_attrib = (
        select(SequenceSource.sequence_id)
        .join(Source, Source.id == SequenceSource.source_id)
        .where(Source.name.in_(ok))
    )
    return by_name.union(by_attrib)


def export_alleles(
    session: Session, *, category: str = "coding", species: str | None = None,
    segment: str | None = None, archived: bool = False,
) -> list[dict]:
    """Return every **redistributable** matching allele as a self-attributing record.

    Uses three bulk queries (sequences, names, attributions) rather than the
    per-record graph, so a full-dataset dump stays fast. Each record carries its
    names (with scheme + source) and the set of contributing sources.

    Sequences no redistributable source attests are **omitted** — see the module
    docstring. Within a record, names and source credits are restricted the same way,
    so the export never credits a source whose contribution it withheld. Use
    :func:`export_counts` to report how many were held back.
    """
    conds = _conds(category, species, segment, archived)
    conds = [*conds, Sequence.id.in_(_redistributable_ids())]
    ok_sources = redistributable_sources()
    seqs = list(session.scalars(select(Sequence).where(*conds).order_by(Sequence.uid)))

    # Sources for attribution come from BOTH a sequence's names (Name -> Source) and
    # its provenance rows (SequenceSource): many alleles are attested only by a name
    # from their source, with no separate attribution row. Merge both so `sources` is
    # complete — an empty attribution would defeat the CC BY requirement.
    names_by: dict[int, list[dict]] = {}
    src_by: dict[int, set] = {}
    seen: dict[int, set] = {}  # dedup (name, scheme, source) — one name per source-version
    for sid, name, scheme, src in session.execute(
        select(Name.sequence_id, Name.name, Name.naming_scheme, Source.name)
        .join(Sequence, Sequence.id == Name.sequence_id)
        .join(Source, Source.id == Name.source_id)
        .where(*conds, Source.name.in_(ok_sources))
    ):
        src_by.setdefault(sid, set()).add(src)
        key = (name, scheme, src)
        if key not in seen.setdefault(sid, set()):
            seen[sid].add(key)
            names_by.setdefault(sid, []).append({"name": name, "scheme": scheme, "source": src})
    for sid, sname in session.execute(
        select(SequenceSource.sequence_id, Source.name)
        .join(Sequence, Sequence.id == SequenceSource.sequence_id)
        .join(Source, Source.id == SequenceSource.source_id)
        .where(*conds, Source.name.in_(ok_sources)).distinct()
    ):
        src_by.setdefault(sid, set()).add(sname)

    out = []
    for s in seqs:
        out.append({
            "uid": s.uid, "species": s.species, "locus": s.locus,
            "secondary_locus": s.secondary_locus, "segment": s.segment,
            "category": s.category, "functional": s.functional, "length": s.length,
            "names": names_by.get(s.id, []),
            "sources": sorted(src_by.get(s.id, set())),
            "sequence": s.sequence, "sequence_aa": s.sequence_aa,
        })
    return out


def to_fasta(records: list[dict], *, kind: str = "nt") -> Iterator[str]:
    """Stream records as FASTA. ``kind`` = ``nt`` (coding sequence) or ``aa``.

    Header: ``>UID name1,name2|species|locus|segment|F=functional|src=…`` so the
    export is self-describing and attributed at the record level.
    """
    for r in records:
        seq = r["sequence"] if kind == "nt" else r["sequence_aa"]
        if not seq:
            continue
        names = ",".join(dict.fromkeys(n["name"] for n in r["names"])) or "-"  # distinct
        yield (
            f">{r['uid']} {names}|{r['species']}|{r['locus'] or ''}|{r['segment'] or ''}"
            f"|F={r['functional']}|src={','.join(r['sources'])}\n{seq}\n"
        )
