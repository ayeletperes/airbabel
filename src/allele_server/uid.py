"""Deterministic per-sequence UID ("iglabel") minting.

A unique sequence's identity is ``(species, category, sequence_hash)``. Its UID
is derived from that, never randomly allocated, so the same sequence yields the
same UID in any database and the family of databases federates without a
central registry::

    UID = f"{SPECIES_TAG}-{CATEGORY_CODE}{base32(sha256(seq))[:N]}"
        HUMAN-V-A3F2C10B   coding, N=8 (40 bits, collision-safe past ~100k seqs)
        HUMAN-L-Q9X4Z      leader, N=5
        HUMAN-R-K7M2Q      rss,    N=5

The trailing token is always a prefix of ``base32(sequence_hash)``, so any UID
is verifiable straight from the sequence. Coding is the most numerous category,
hence the longest token. ``CATEGORY_CODE`` is the segment letter (V/D/J) for
coding, ``L`` for leader, ``R`` for rss.

Collision backstop: at the default lengths a clash needs an astronomically
unlikely hash prefix collision within one ``(species, category)`` namespace. If
it ever happens, the *newcomer* extends its token (taking more of its own hash)
until free; the incumbent keeps its shorter UID. The full ``sequence_hash`` is
always stored, so it remains the authoritative cross-database key regardless.
"""

from __future__ import annotations

from base64 import b32encode
from enum import StrEnum

from .sequence import hash_sequence

__all__ = [
    "Category",
    "UidError",
    "base32_of_hash",
    "category_code",
    "default_length",
    "mint_protein_uid",
    "mint_uid",
]

# base32 chars taken from aa_hash for the (non-unique) protein UID token.
_PROTEIN_LENGTH = 6

# base32 character budget per category (5 bits each).
_DEFAULT_LENGTH: dict["Category", int] = {}
_MAX_LENGTH = 52  # 256-bit sha256 -> 52 base32 chars of real entropy


class UidError(ValueError):
    """Raised when a UID cannot be constructed."""


class Category(StrEnum):
    CODING = "coding"
    LEADER = "leader"
    RSS = "rss"


_DEFAULT_LENGTH.update({Category.CODING: 8, Category.LEADER: 5, Category.RSS: 5})

_SEGMENT_CODES = frozenset({"V", "D", "J", "C"})


def default_length(category: Category) -> int:
    """Default base32 token length for a category."""
    return _DEFAULT_LENGTH[category]


def category_code(category: Category, segment: str | None) -> str:
    """Single-letter code embedded in the UID.

    Coding uses the segment letter (V/D/J/C); leader is ``L``, rss is ``R``.

    Raises:
        UidError: If a coding category has a missing/invalid segment.
    """
    if category is Category.CODING:
        seg = (segment or "").upper()
        if seg not in _SEGMENT_CODES:
            raise UidError(
                f"coding UID needs segment in {sorted(_SEGMENT_CODES)}, got {segment!r}"
            )
        return seg
    if category is Category.LEADER:
        return "L"
    if category is Category.RSS:
        return "R"
    raise UidError(f"unknown category {category!r}")  # pragma: no cover


def base32_of_hash(sequence_hash: str) -> str:
    """Return the uppercase, unpadded base32 of a sha256 hex digest."""
    try:
        raw = bytes.fromhex(sequence_hash)
    except ValueError as exc:
        raise UidError(f"sequence_hash is not hex: {sequence_hash!r}") from exc
    if len(raw) != 32:
        raise UidError(
            f"sequence_hash must be a 32-byte sha256 (64 hex chars), got {len(raw)} bytes"
        )
    return b32encode(raw).decode("ascii").rstrip("=")


def mint_uid(
    species_tag: str,
    category: Category,
    segment: str | None,
    sequence_hash: str,
    exists_hash=None,
) -> str:
    """Build the UID for a sequence.

    Args:
        species_tag: Canonical tag from :mod:`allele_server.species` (``HUMAN``).
        category: One of :class:`Category`.
        segment: ``V``/``D``/``J`` for coding, ignored otherwise.
        sequence_hash: ``sha256`` hex digest of the normalised sequence
            (from :func:`allele_server.sequence.hash_sequence`).
        exists_hash: Optional callable ``uid -> sequence_hash | None`` used to
            detect collisions. Returns the hash already bound to ``uid``, or
            ``None`` if free. When omitted, the default-length UID is returned
            without a collision check.

    Returns:
        The UID string.

    Raises:
        UidError: On bad inputs, or if the hash space is exhausted (impossible
            in practice).
    """
    code = category_code(category, segment)
    token = base32_of_hash(sequence_hash)
    n = default_length(category)
    while n <= _MAX_LENGTH:
        uid = f"{species_tag}-{code}{token[:n]}"
        if exists_hash is None:
            return uid
        bound = exists_hash(uid)
        if bound is None or bound == sequence_hash:
            return uid
        # Different sequence already holds this UID: newcomer takes more chars.
        n += 1
    raise UidError(
        f"could not allocate a unique UID for {sequence_hash} under "
        f"{species_tag}-{code} (hash space exhausted)"
    )


def mint_protein_uid(species_tag: str, segment: str, aa_hash: str) -> str:
    """Build the V-indicating protein UID for a coding allele's AA.

    ``f"{SPECIES_TAG}-{SEGMENT}P-{base32(aa_hash)[:6]}"`` (e.g. ``HUMAN-VP-Q9X4Z3``).
    Derived purely from ``aa_hash``, so every nucleotide allele encoding the same
    protein shares it. Deliberately **not** uniqueness-checked: collisions are the
    point (grouping), and the full ``aa_hash`` stays the authoritative key.

    Raises:
        UidError: If ``segment`` is not V/D/J or ``aa_hash`` is not a valid digest.
    """
    seg = (segment or "").upper()
    if seg not in _SEGMENT_CODES:
        raise UidError(
            f"protein UID needs segment in {sorted(_SEGMENT_CODES)}, got {segment!r}"
        )
    token = base32_of_hash(aa_hash)[:_PROTEIN_LENGTH]
    return f"{species_tag}-{seg}P-{token}"


def mint_uid_for_sequence(
    species_tag: str,
    category: Category,
    segment: str | None,
    sequence: str,
    exists_hash=None,
) -> tuple[str, str]:
    """Convenience: hash a raw sequence and mint its UID.

    Returns ``(uid, sequence_hash)``.
    """
    digest = hash_sequence(sequence)
    return mint_uid(species_tag, category, segment, digest, exists_hash), digest
