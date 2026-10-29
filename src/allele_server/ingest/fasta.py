"""Minimal FASTA reading + IMGT header parsing for ingestion.

Handles both IMGT pipe-delimited headers
(``>X02897|IGHV1-2*02|Homo_sapiens|F|V-REGION|...``) and plain
``>IGHV1-2*02`` headers (OGRDB / internal sets).
"""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path

# IMGT FASTA header positional fields (the subset we use).
_IMGT_FIELDS = (
    "accession", "allele", "species", "functionality", "label",
    "position", "length_nt",
)


def require_fasta(path: str | Path) -> Path:
    """Return ``path``, raising if it is not a readable, non-empty file.

    :func:`read_fasta` yields nothing for a missing file, which suits callers where
    the input is genuinely optional. For an ingest it is a trap: the pass returns
    zero counts and reads as "nothing to do" rather than "the file is not there".
    That is the silent fallback this codebase forbids, so a required input goes
    through here first.
    """
    p = Path(path)
    if not p.is_file():
        raise FileNotFoundError(f"required FASTA not found: {p}")
    if p.stat().st_size == 0:
        raise ValueError(f"required FASTA is empty: {p}")
    return p


def read_fasta(path: str | Path) -> Iterator[tuple[str, str]]:
    """Yield ``(header, sequence)`` from a FASTA file (header without ``>``).

    Sequence whitespace is collapsed. Empty/non-existent files yield nothing.
    """
    p = Path(path)
    if not p.is_file():
        return
    header: str | None = None
    seq: list[str] = []
    for raw in p.read_text().splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(seq)
            header = line[1:]
            seq = []
        else:
            seq.append("".join(line.split()))
    if header is not None:
        yield header, "".join(seq)


@dataclass(frozen=True, slots=True)
class FastaRecord:
    """One parsed allele record from a germline FASTA."""

    allele: str
    sequence: str
    accession: str | None
    species: str | None
    functionality: str | None
    label: str | None

    @property
    def functional(self) -> bool | None:
        """AIRR ``functional`` boolean from the IMGT token, or None if absent.

        Strict: only ``F`` is functional. ``ORF`` and ``P`` (pseudogene) map to
        ``False``. The raw token is preserved separately as an annotation.
        """
        if self.functionality is None or self.functionality == "":
            return None
        return self.functionality.strip().strip("[]()") == "F"


def parse_record(header: str, sequence: str) -> FastaRecord:
    """Parse one FASTA header+sequence into a :class:`FastaRecord`.

    Pipe-delimited headers are read as IMGT; a header with no ``|`` is treated
    as a bare allele name.
    """
    if "|" in header:
        parts = [p.strip() for p in header.split("|")]
        fields = {
            name: (parts[i] if i < len(parts) else "")
            for i, name in enumerate(_IMGT_FIELDS)
        }
        return FastaRecord(
            allele=fields["allele"] or parts[0],
            sequence=sequence,
            accession=fields["accession"] or None,
            species=(fields["species"] or None),
            functionality=fields["functionality"] or None,
            label=fields["label"] or None,
        )
    return FastaRecord(
        allele=header.split()[0],
        sequence=sequence,
        accession=None,
        species=None,
        functionality=None,
        label=None,
    )
