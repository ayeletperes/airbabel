"""Idempotent upsert helpers shared by every ingestion path.

Every function is safe to call repeatedly: re-running ingestion produces zero
row changes. UID minting flows through :mod:`allele_server.uid` with a DB-backed
collision check.
"""

from __future__ import annotations

import logging

from sqlalchemy import select
from sqlalchemy.orm import Session

from ..db.models import (
    AA_ORIGIN_SOURCED,
    AA_ORIGIN_TRANSLATED,
    AA_ORIGIN_TRANSLATED_OVER_CONFLICT,
    Name,
    Sequence,
    SequenceAnnotation,
    SequenceLink,
    SequenceSource,
    Source,
)
from ..sequence import (
    SequenceError,
    hash_aa,
    hash_sequence,
    is_ambiguous,
    normalize,
    reverse_complement_hash,
    translate_v,
)
from ..settings import get_settings
from ..sources import SourceError, SourceManifest, resolve_strict
from ..uid import Category, mint_protein_uid, mint_uid

_log = logging.getLogger(__name__)


def _apply_manifest(src: Source, manifest: SourceManifest) -> None:
    """Copy the manifest's metadata onto a source row.

    Called on create *and* on every subsequent fetch, so editing a manifest propagates
    on the next ingest. This is the one documented exception to "re-running an ingest
    changes zero rows" (``docs/source-policy.md``): the alternative is a separate sync
    step someone has to remember, which is the drift the manifests exist to remove.
    """
    src.kind = manifest.role
    src.manifest_id = manifest.slug
    src.license = manifest.licence.name
    src.licence_url = manifest.licence.url
    src.redistributable = manifest.licence.redistributable
    src.citation = manifest.citation.text or None
    src.citation_status = manifest.citation.status
    src.doi = manifest.citation.doi
    src.url = manifest.homepage_url
    if manifest.approval is not None:
        src.approved_by = manifest.approval.approved_by
        src.approved_at = manifest.approval.approved_at


def get_or_create_source(
    session: Session,
    *,
    name: str,
    version: str | None = None,
    curie_prefix: str | None = None,
    germline_set_ref: str | None = None,
    note: str | None = None,
) -> Source:
    """Fetch or create a source row (natural key ``(name, version)``).

    **The gate.** Every ingest path routes through here, so this is where a source with
    no approved manifest is refused. The manifest is truth for kind, licence, citation,
    DOI, homepage and redistributability; those are projected onto the row on every
    call. See ``docs/source-policy.md``.

    ``kind`` / ``redistributable`` / ``license`` used to be caller-supplied and are
    gone: every call site passed ``redistributable=True`` and only one ever passed a
    licence, which is how the export gate came to mean nothing.

    Raises:
        SourceError: if no manifest resolves ``name`` (or one of its aliases), if that
            manifest carries no ``[approval]`` block, or if a caller-supplied
            ``curie_prefix`` contradicts the manifest's. Suppressed to a warning by
            ``settings.source_gate``, which exists for a broken deploy, not for
            development.
    """
    gate = get_settings().source_gate
    manifest: SourceManifest | None = None
    try:
        manifest = resolve_strict(name)
    except SourceError as exc:
        if gate == "strict":
            raise
        if gate == "warn":
            _log.warning("source gate (warn): %s", exc)
        elif gate != "off":
            raise SourceError(
                f"invalid source_gate {gate!r}; expected 'strict', 'warn' or 'off'"
            ) from exc

    if manifest is not None and curie_prefix and manifest.curie_prefix:
        if curie_prefix != manifest.curie_prefix:
            raise SourceError(
                f"source {name!r}: caller passed curie_prefix {curie_prefix!r} but its "
                f"manifest declares {manifest.curie_prefix!r}. One of the two is wrong; "
                f"the manifest is the one to fix."
            )

    stmt = select(Source).where(Source.name == name)
    stmt = stmt.where(Source.version == version) if version is not None else stmt.where(
        Source.version.is_(None)
    )
    src = session.scalar(stmt)
    if src is None:
        src = Source(
            name=name,
            version=version,
            # Overwritten by _apply_manifest below; only reachable as-is when the gate
            # is off, and then the honest default is "we know nothing, so no export".
            kind="external",
            redistributable=False,
            curie_prefix=curie_prefix,
            germline_set_ref=germline_set_ref,
            note=note,
        )
        session.add(src)
    elif germline_set_ref and not src.germline_set_ref:
        src.germline_set_ref = germline_set_ref
    if manifest is not None:
        _apply_manifest(src, manifest)
    session.flush()
    return src


def _uid_exists(session: Session):
    def f(uid: str) -> str | None:
        return session.scalar(select(Sequence.sequence_hash).where(Sequence.uid == uid))

    return f


def get_coding(session: Session, species: str, sequence_hash: str) -> Sequence | None:
    """Return the existing coding sequence for ``(species, hash)`` or None."""
    return session.scalar(
        select(Sequence).where(
            Sequence.species == species,
            Sequence.category == "coding",
            Sequence.sequence_hash == sequence_hash,
        )
    )


def upsert_coding(
    session: Session,
    *,
    species: str,
    segment: str,
    locus: str | None,
    nt: str,
    functional: bool | None = None,
    secondary_locus: str | None = None,
    subtype: str | None = None,
) -> Sequence:
    """Create (or return) a coding sequence row, minting its UID.

    Idempotent on ``(species, 'coding', sequence_hash)``. ``nt`` is normalised
    (gaps stripped) before storage; the verbatim original belongs in
    ``sequence_source``. ``secondary_locus`` is the second locus for a shared
    TRA/TRD V gene (see :func:`allele_server.naming.dual_locus_of`); it is
    backfilled on an existing row so a reload can add it. ``subtype`` names the
    representation where a segment has more than one (constant regions:
    ``c_secreted``); backfilled the same way.
    """
    digest = hash_sequence(nt)
    existing = get_coding(session, species, digest)
    if existing is not None:
        if functional is not None and existing.functional is None:
            existing.functional = functional
        if secondary_locus is not None and existing.secondary_locus is None:
            existing.secondary_locus = secondary_locus
        if subtype is not None and existing.subtype is None:
            existing.subtype = subtype
        return existing
    norm = normalize(nt)
    uid = mint_uid(species, Category.CODING, segment, digest, _uid_exists(session))
    seq = Sequence(
        uid=uid,
        species=species,
        locus=locus,
        secondary_locus=secondary_locus,
        category="coding",
        segment=segment,
        subtype=subtype,
        sequence=norm,
        sequence_hash=digest,
        rc_hash=reverse_complement_hash(nt),
        functional=functional,
        length=len(norm),
        has_ambiguity=is_ambiguous(nt),
    )
    session.add(seq)
    session.flush()
    return seq


_CATEGORY = {"rss": Category.RSS, "leader": Category.LEADER}


def upsert_noncoding(
    session: Session, *, species: str, category: str, subtype: str, nt: str
) -> Sequence:
    """Create (or return) an RSS or leader sequence row, minting its UID.

    Idempotent on ``(species, category, sequence_hash)``. ``category`` is
    ``rss`` or ``leader``; ``subtype`` is e.g. ``heptamer`` / ``l_part1``.
    """
    digest = hash_sequence(nt)
    existing = session.scalar(
        select(Sequence).where(
            Sequence.species == species,
            Sequence.category == category,
            Sequence.sequence_hash == digest,
        )
    )
    if existing is not None:
        return existing
    norm = normalize(nt)
    uid = mint_uid(species, _CATEGORY[category], None, digest, _uid_exists(session))
    seq = Sequence(
        uid=uid,
        species=species,
        category=category,
        subtype=subtype,
        sequence=norm,
        sequence_hash=digest,
        rc_hash=reverse_complement_hash(nt),
        length=len(norm),
        has_ambiguity=is_ambiguous(nt),
    )
    session.add(seq)
    session.flush()
    return seq


def add_link(
    session: Session, from_seq: Sequence, to_seq: Sequence, relation: str, source: Source
) -> SequenceLink:
    """Link two sequences (e.g. coding -> its RSS). Idempotent on the natural key."""
    existing = session.scalar(
        select(SequenceLink).where(
            SequenceLink.from_sequence_id == from_seq.id,
            SequenceLink.to_sequence_id == to_seq.id,
            SequenceLink.relation == relation,
            SequenceLink.source_id == source.id,
        )
    )
    if existing is not None:
        return existing
    link = SequenceLink(
        from_sequence_id=from_seq.id, to_sequence_id=to_seq.id,
        relation=relation, source_id=source.id,
    )
    session.add(link)
    session.flush()
    return link


def _translate_best_frame(nt: str) -> str:
    """Frame-1 translation (the V-REGION reading frame). If frame 1 carries a stop
    but a shifted forward frame is clean, use that frame — ORF alleles are often
    read in another frame. No clean frame (a true pseudogene) keeps frame 1, whose
    stops are stripped downstream by ``normalize_aa``."""
    f1 = translate_v(nt)
    if "*" not in f1:
        return f1
    for off in (1, 2):
        try:
            shifted = translate_v(nt[off:])
        except SequenceError:
            continue
        if "*" not in shifted:
            return shifted
    return f1


# Strip gaps + whitespace but KEEP '*' — so a non-functional allele's stop codons
# stay visible in the stored AA. Hashing still strips '*' (hash_aa), so protein
# identity is unaffected by whether stops are shown.
_AA_KEEP_STOPS = dict.fromkeys(ord(c) for c in ".-\t\n\r\x0b\x0c \xa0﻿")


def _clean_keep_stops(seq_aa: str) -> str:
    cleaned = seq_aa.upper().translate(_AA_KEEP_STOPS)
    if not cleaned:
        raise SequenceError("AA sequence is empty after normalisation")
    return cleaned


def _aa_structural_conflict(a: str, b: str) -> bool:
    """True if two AAs disagree *beyond* nucleotide ambiguity — a wrong-protein /
    frameshift difference, not a lone ``X`` (ambiguous codon) vs a resolved residue.
    Compared on the stop-stripped forms; length divergence or several non-``X``
    mismatches count as structural.
    """
    a, b = a.replace("*", ""), b.replace("*", "")
    if abs(len(a) - len(b)) > 2:
        return True
    n = min(len(a), len(b))
    if n == 0:
        return False
    diffs = sum(1 for i in range(n) if a[i] != b[i] and a[i] != "X" and b[i] != "X")
    return diffs >= 3 or diffs > 0.02 * n


def set_v_aa(
    session: Session, seq: Sequence, aa_source: str | None = None,
    source: Source | None = None,
) -> None:
    """Attach the V-REGION AA + aa_hash + protein UID to a V coding row.

    The stored AA must stay **consistent with the stored nucleotides** — it is what
    the nt encodes. A directly-sourced AA (e.g. IMGT's V-AA) is trusted only when it
    agrees with the sequence's own translation. IMGT sometimes publishes a *reference*
    protein its own nucleotides don't encode (silently healing a deposit indel — e.g.
    ``IGHV5-10-1*02``, a 1-nt deletion); storing that beside our nt would make the row
    internally inconsistent, so on a **structural conflict** we keep the translation
    and record the sourced reference AA as an ``aa/imgt_reference_aa`` annotation
    (when a ``source`` is given). No-op for non-V/non-coding or nothing translatable.

    Which of the three paths produced the stored AA is recorded in ``aa_origin``
    (:data:`AA_ORIGINS`), so a reader can tell a published protein from one we
    computed. It is a column, not an annotation, because it describes the single
    stored AA rather than any source's claim about it.
    """
    if seq.category != "coding" or seq.segment != "V":
        return
    # The sequence's own translation (best frame; '*' kept so pseudogene stops show).
    try:
        translated = _clean_keep_stops(_translate_best_frame(seq.sequence))
    except SequenceError:
        translated = None

    if aa_source:
        sourced = _clean_keep_stops(aa_source)
        if translated is not None and _aa_structural_conflict(sourced, translated):
            aa = translated
            origin = AA_ORIGIN_TRANSLATED_OVER_CONFLICT
            if source is not None:
                add_annotation(
                    session, seq, source, domain="aa",
                    key="imgt_reference_aa", value=sourced,
                )
        else:
            aa = sourced
            origin = AA_ORIGIN_SOURCED
    elif seq.sequence_aa:
        return  # already set by a prior pass
    elif translated is not None:
        aa = translated
        origin = AA_ORIGIN_TRANSLATED
    else:
        return

    digest = hash_aa(aa)  # hash_aa strips '*' → protein grouping unaffected by stops
    if seq.sequence_aa == aa and seq.aa_hash == digest and seq.aa_origin == origin:
        return
    seq.sequence_aa = aa
    seq.aa_hash = digest
    seq.aa_origin = origin
    seq.protein_uid = mint_protein_uid(seq.species, seq.segment, digest)
    session.flush()


def backfill_v_aa(session: Session) -> list[str]:
    """Translate + attach AA to every V coding row that still lacks it (rows
    created by paths that didn't call :func:`set_v_aa`). Idempotent — only fills
    NULLs. Returns the uids that were filled (for an incremental aa-cache update).

    Rows this fills get ``aa_origin = 'translated'`` from :func:`set_v_aa`, because
    that is what happened. Rows that already carry an AA are left alone, origin
    included: a NULL origin means the row predates the column and we do **not** know
    which path produced it. Stamping those ``translated`` would be a guess, and on a
    database loaded before the column existed it would relabel every genuinely
    published IMGT protein as ours. Only a reload re-derives them honestly."""
    filled: list[str] = []
    for seq in session.scalars(
        select(Sequence).where(
            Sequence.category == "coding", Sequence.segment == "V",
            Sequence.sequence_aa.is_(None),
        )
    ):
        set_v_aa(session, seq)
        if seq.sequence_aa is not None:
            filled.append(seq.uid)
    session.flush()
    return filled


def add_name(
    session: Session, seq: Sequence, source: Source, name: str, scheme: str
) -> Name | None:
    """Attach an alias. Idempotent on ``(sequence, name, scheme, source)``.

    A name is a *label*: the same allele name can apply to several distinct
    sequences (IMGT strain variants reuse e.g. ``IGHV1-2*01``), so the key includes
    the sequence — each variant keeps its own name. Re-attaching the identical name
    to the SAME sequence is a no-op.
    """
    existing = session.scalar(
        select(Name).where(
            Name.name == name,
            Name.naming_scheme == scheme,
            Name.source_id == source.id,
            Name.sequence_id == seq.id,
        )
    )
    if existing is not None:
        return existing
    row = Name(
        sequence_id=seq.id, source_id=source.id, name=name, naming_scheme=scheme
    )
    session.add(row)
    session.flush()
    return row


def backfill_names_from_attributions(session: Session) -> int:
    """Give a name to any coding allele that has a source attribution but no Name.

    These are strain-variant sequences that reused an allele name and lost it to
    the old ``(name, scheme, source)`` unique key (now keyed per sequence too). The
    name lives in the attribution's CURIE ``original_id`` (e.g.
    ``IMGT/GENE-DB:IGHV1-2*01``) — re-attach it. Idempotent; returns names added.
    """
    from ..naming import scheme_label

    named = select(Name.sequence_id).distinct()
    unnamed = session.execute(
        select(Sequence).where(Sequence.category == "coding", Sequence.id.notin_(named))
    ).scalars().all()
    added = 0
    for seq in unnamed:
        attr = session.execute(
            select(SequenceSource.original_id, Source)
            .join(Source, Source.id == SequenceSource.source_id)
            .where(SequenceSource.sequence_id == seq.id).limit(1)
        ).first()
        if not attr:
            continue
        original_id, source = attr
        name = original_id.split(":", 1)[-1].strip()  # CURIE value after the prefix
        if not name:
            continue
        add_name(session, seq, source, name, scheme_label(name))
        added += 1
    session.commit()
    return added


def add_attribution(
    session: Session,
    seq: Sequence,
    source: Source,
    *,
    original_id: str,
    original_sequence: str,
) -> SequenceSource:
    """Record that ``source`` attests ``seq``. Idempotent on the natural key."""
    existing = session.scalar(
        select(SequenceSource).where(
            SequenceSource.sequence_id == seq.id,
            SequenceSource.source_id == source.id,
            SequenceSource.original_id == original_id,
        )
    )
    if existing is not None:
        return existing
    row = SequenceSource(
        sequence_id=seq.id,
        source_id=source.id,
        original_id=original_id,
        original_sequence=original_sequence,
    )
    session.add(row)
    session.flush()
    return row


def add_annotation(
    session: Session,
    seq: Sequence,
    source: Source,
    *,
    domain: str,
    key: str,
    value: str | int | float | bool,
) -> SequenceAnnotation:
    """Attach one key/value annotation. Idempotent on ``(seq, source, domain, key)``."""
    existing = session.scalar(
        select(SequenceAnnotation).where(
            SequenceAnnotation.sequence_id == seq.id,
            SequenceAnnotation.source_id == source.id,
            SequenceAnnotation.domain == domain,
            SequenceAnnotation.key == key,
        )
    )
    if existing is not None:
        return existing
    row = SequenceAnnotation(
        sequence_id=seq.id, source_id=source.id, domain=domain, key=key
    )
    # bool must be checked before int (bool is a subclass of int).
    if isinstance(value, bool):
        row.value_bool = value
    elif isinstance(value, int):
        row.value_int = value
    elif isinstance(value, float):
        row.value_float = value
    else:
        row.value_text = str(value)
    session.add(row)
    session.flush()
    return row
