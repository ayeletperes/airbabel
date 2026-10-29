"""NCBI cross-reference enrichment pass.

For each coding sequence: if it already carries a GenBank accession (from the
IMGT header), confirm it via ``efetch`` and grade the match; otherwise resolve
by allele name via ``esearch``. Writes an ``ncbi`` alias (the accession), an
``ncbi_nucleotide`` source attribution (only when the allele nt is confirmed
contained in the GenBank record), and ``xref`` annotations (URL, NCBI Gene,
match confidence). When the sequence is confirmed, the PMIDs the GenBank record
cites are also written as ``xref`` ``pmid:<id>`` annotations — a high-precision
literature anchor that ``queries._references`` surfaces as **verified**.
Idempotent and cache-backed.
"""

from __future__ import annotations

import logging

from sqlalchemy import select
from sqlalchemy.orm import Session

from ..curie import build as curie_build
from ..curie import split as curie_split
from ..db.models import Name, Sequence, SequenceSource
from . import loader
from ._ncbi_client import NCBIClient

_LOG = logging.getLogger(__name__)


def _representative_name(session: Session, seq_id: int) -> str | None:
    pref = session.scalar(
        select(Name.name)
        .where(Name.sequence_id == seq_id, Name.naming_scheme.in_(("iuis", "imgt")))
        .limit(1)
    )
    if pref:
        return pref
    return session.scalar(select(Name.name).where(Name.sequence_id == seq_id).limit(1))


def _existing_accession(session: Session, seq_id: int) -> str | None:
    curie = session.scalar(
        select(SequenceSource.original_id)
        .where(
            SequenceSource.sequence_id == seq_id,
            SequenceSource.original_id.like("GENBANK:%"),
        )
        .limit(1)
    )
    return curie_split(curie)[1] if curie else None


def enrich_ncbi(
    session: Session,
    client: NCBIClient,
    *,
    species_keys: list[str] | None = None,
    segments: tuple[str, ...] = ("V",),
    limit: int | None = None,
) -> dict[str, int]:
    """Walk coding sequences and attach NCBI cross-references.

    Returns counts: ``confirmed`` (nt contained in GenBank record), ``name_only``
    (name matched, sequence not confirmed), ``unresolved``, ``stale``, and
    ``pubmed_refs`` (verified PMIDs cited by confirmed GenBank records).
    """
    source = loader.get_or_create_source(
        session, name="ncbi_nucleotide", version=None, curie_prefix="GENBANK",
        note="NCBI nuccore via E-utilities",
    )
    stmt = select(Sequence).where(
        Sequence.category == "coding", Sequence.segment.in_(segments)
    )
    if species_keys:
        from ..species import normalize_species

        tags = [normalize_species(s) for s in species_keys]
        stmt = stmt.where(Sequence.species.in_(tags))
    if limit:
        stmt = stmt.limit(limit)

    stats = {
        "confirmed": 0, "name_only": 0, "unresolved": 0, "stale": 0,
        "pubmed_refs": 0, "errors": 0,
    }
    # NCBI enrichment is best-effort cross-reference data. It runs inside a full
    # reload, so a single malformed response or record must not abort the rebuild
    # (one corrupt cache entry once did). Isolate each sequence: log and carry on.
    for seq in session.scalars(stmt):
        try:
            acc = _existing_accession(session, seq.id)
            if acc:
                res = client.confirm_accession(acc, seq.sequence)
            else:
                name = _representative_name(session, seq.id)
                if not name:
                    stats["unresolved"] += 1
                    continue
                res = client.resolve(name, seq.sequence)
            if res.stale:
                stats["stale"] += 1
            if res.accession is None:
                stats["unresolved"] += 1
                continue

            loader.add_name(session, seq, source, res.accession, "ncbi")
            if res.url:
                loader.add_annotation(
                    session, seq, source, domain="xref", key="ncbi_url", value=res.url
                )
            loader.add_annotation(
                session, seq, source, domain="xref", key="match_confidence",
                value=res.confidence,
            )
            if res.gene_id:
                loader.add_annotation(
                    session, seq, source, domain="xref", key="ncbi_gene",
                    value=curie_build("NCBIGENE", res.gene_id),
                )
            if res.confidence == "exact_substring":
                loader.add_attribution(
                    session, seq, source,
                    original_id=curie_build("GENBANK", res.accession),
                    original_sequence=seq.sequence,
                )
                # High-precision literature anchor: papers the confirmed GenBank
                # record cites. Stored under domain ``xref`` so ``_references`` marks
                # them verified. Idempotent on ``(seq, source, "xref", pmid:<id>)``.
                for pmid in res.linked_pubmed_ids:
                    loader.add_annotation(
                        session, seq, source, domain="xref",
                        key=f"pmid:{pmid}"[:64], value=curie_build("PMID", pmid),
                    )
                    stats["pubmed_refs"] += 1
                stats["confirmed"] += 1
            else:
                stats["name_only"] += 1
        except Exception as exc:  # noqa: BLE001 - best-effort; never abort a reload
            stats["errors"] += 1
            _LOG.warning("NCBI enrich skipped %s: %s", seq.uid, exc)
    return stats
