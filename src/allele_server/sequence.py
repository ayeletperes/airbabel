"""Sequence identity primitives: normalisation, hashing, ambiguity detection.

The same identity contract as the HUSA curation pipeline (``husa_core``), so
databases built on it federate. Locked rules:

    * Canonical hash = ``sha256(normalize(seq))`` as a hex digest.
    * Normalisation uppercases, strips gaps (``.`` ``-``) and stray whitespace,
      validates against the IUPAC DNA alphabet, rejects empty and RNA input.
    * Ambiguity codes (N, R, Y, S, W, K, M, B, D, H, V) are preserved, never
      collapsed. ``is_ambiguous`` flags them so the caller can warn users.
    * Reverse complement is a *different* sequence. Its hash is exposed as a
      side attribute via ``reverse_complement_hash``.
"""

from __future__ import annotations

from hashlib import sha256

__all__ = [
    "SequenceError",
    "aa_hash",
    "hash_aa",
    "hash_sequence",
    "is_ambiguous",
    "normalize",
    "normalize_aa",
    "reverse_complement_hash",
    "translate_v",
]

_UNAMBIGUOUS: frozenset[str] = frozenset("ACGT")
_AMBIGUOUS: frozenset[str] = frozenset("NRYSWKMBDHV")
_IUPAC: frozenset[str] = _UNAMBIGUOUS | _AMBIGUOUS

_STRIP_CHARS: str = ".-\t\n\r\x0b\x0c \xa0﻿"
_STRIP_TABLE: dict[int, None] = dict.fromkeys(ord(c) for c in _STRIP_CHARS)

_COMPLEMENT_TABLE: dict[int, int] = {
    ord(a): ord(b)
    for a, b in zip("ACGTNRYSWKMBDHV", "TGCANYRSWMKVHDB", strict=True)
}


class SequenceError(ValueError):
    """Raised when a sequence cannot be normalised into a canonical form."""


def normalize(seq: str) -> str:
    """Return the canonical form of a DNA sequence.

    Canonical form is uppercase, gapless, and IUPAC-validated. Trailing
    ``N`` padding is preserved because stripping it would change identity.

    Raises:
        SequenceError: If ``seq`` is not a string, is empty after
            normalisation, contains ``U`` (RNA), or contains any character
            outside the IUPAC DNA alphabet.
    """
    if not isinstance(seq, str):
        raise SequenceError(f"sequence must be str, got {type(seq).__name__}")
    cleaned = seq.upper().translate(_STRIP_TABLE)
    if not cleaned:
        raise SequenceError("sequence is empty after normalisation")
    if "U" in cleaned:
        raise SequenceError("RNA sequence (contains 'U') rejected; supply DNA")
    invalid = sorted({c for c in cleaned if c not in _IUPAC})
    if invalid:
        raise SequenceError(f"sequence contains non-IUPAC characters: {invalid!r}")
    return cleaned


def hash_sequence(seq: str) -> str:
    """Return the canonical sha256 hex digest of ``normalize(seq)``.

    Any two inputs that normalise to the same string share an identity. This
    is the value stored in ``sequences.sequence_hash`` and the basis of the
    deterministic UID.
    """
    return sha256(normalize(seq).encode("ascii")).hexdigest()


def is_ambiguous(seq: str) -> bool:
    """Return whether the normalised sequence carries any IUPAC ambiguity."""
    return not _AMBIGUOUS.isdisjoint(normalize(seq))


# Standard genetic code, codon -> single-letter AA. ``*`` is stop. Codons
# containing any IUPAC ambiguity code (or trailing N padding) collapse to
# ``X``: we don't try to resolve them, but we don't drop them either, so
# the AA hash stays stable across re-runs of the same record.
_CODON_TABLE: dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def translate_v(seq: str) -> str:
    """Translate an in-frame ungapped V-REGION DNA sequence to amino acids.

    The DB stores V-REGION ungapped from codon 1, so frame is always 0. A
    partial trailing codon is dropped rather than rejected.

    Raises:
        SequenceError: If the normalised sequence is shorter than one full
            codon, or :func:`normalize` rejects the input.
    """
    dna = normalize(seq)
    if len(dna) < 3:
        raise SequenceError(
            f"length {len(dna)} is shorter than one codon; cannot translate"
        )
    full_codon_len = (len(dna) // 3) * 3
    out = []
    for i in range(0, full_codon_len, 3):
        codon = dna[i:i + 3]
        out.append(_CODON_TABLE.get(codon, "X"))
    return "".join(out)


def aa_hash(seq: str) -> str:
    """Return the sha256 hex digest of the translated V-REGION.

    Side-attribute analogue of :func:`hash_sequence`: the join key for finding
    sequences that encode the same protein despite differing at the nucleotide
    level. NEVER mixed with ``sequence_hash``. Routed through :func:`hash_aa`
    so it matches a directly-sourced AA string (IMGT V-AA) of the same protein.
    """
    return hash_aa(translate_v(seq))


_AA_STRIP_TABLE: dict[int, None] = dict.fromkeys(ord(c) for c in ".-*\t\n\r\x0b\x0c \xa0﻿")


def normalize_aa(seq_aa: str) -> str:
    """Return the canonical form of an amino-acid sequence.

    Uppercases and strips gaps (``.`` ``-``), stop markers (``*``), and
    whitespace. No alphabet validation beyond non-emptiness — translations may
    legitimately carry ``X`` for ambiguous codons.

    Raises:
        SequenceError: If empty after normalisation.
    """
    if not isinstance(seq_aa, str):
        raise SequenceError(f"AA sequence must be str, got {type(seq_aa).__name__}")
    cleaned = seq_aa.upper().translate(_AA_STRIP_TABLE)
    if not cleaned:
        raise SequenceError("AA sequence is empty after normalisation")
    return cleaned


def hash_aa(seq_aa: str) -> str:
    """Return the sha256 hex digest of a *given* amino-acid string.

    Distinct from :func:`aa_hash` (which translates nt first). Use this for an
    AA sequence sourced directly (e.g. IMGT's V-AA) so it hashes into the same
    ``aa_hash`` as ``aa_hash(nt)`` when they agree.
    """
    return sha256(normalize_aa(seq_aa).encode("ascii")).hexdigest()


def reverse_complement_hash(seq: str) -> str:
    """Return the sha256 hex digest of the reverse complement.

    Stored as a side attribute so resolvers can answer "did you mean the
    reverse complement of X?" without treating the two as the same entry.
    """
    rc = normalize(seq).translate(_COMPLEMENT_TABLE)[::-1]
    return sha256(rc.encode("ascii")).hexdigest()
