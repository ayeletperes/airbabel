"""KI-TCR: population-based human TR germline alleles from Karolinska Institutet.

The TR counterpart of KIARVA, from the same laboratory (Karlsson Hedestam). Ten
per-segment FASTA files, one per (locus, segment), served at fixed URLs under
https://gkhlab.gitlab.io/tcr/sequences/.

**Ingested and credited, not republished.** The site publishes no licence; what it
grants is narrower than one and says so: *"Please cite the following paper when using
these data for academic purposes ... For other uses, please contact the corresponding
authors."* The manifest records that as ``licence.redistributable = false``, and
``export.export_alleles`` filters on that flag, so these sequences reach the site and
the per-record API but never the bulk download. See ``src/allele_server/sources/ki_tcr.toml``.

The site also publishes four **archaic hominin** sets (Altai, Chagyrskaya, Denisova,
Vindija). Those are deliberately not ingested: they are archaic genomes, not modern
human germline, and filing them under species HUMAN would assert something false about
every one of them.
TODO(ki-tcr-archaic): they would need their own species tags to be representable at all.
"""

from __future__ import annotations

from pathlib import Path

from sqlalchemy.orm import Session

from ..curie import build as curie_build
from ..naming import locus_of, scheme_label, segment_of
from ..sequence import SequenceError
from ..species import normalize_species
from . import loader
from .fasta import read_fasta

_BASE = "https://gkhlab.gitlab.io/tcr/sequences/{stem}.fasta"

#: One file per (locus, segment) the site publishes. Named explicitly rather than
#: discovered, so a silently added file cannot enter the database unreviewed and a
#: silently removed one fails loudly instead of shrinking the ingest.
_FILES = (
    "TRAV", "TRAJ",
    "TRBV", "TRBD", "TRBJ",
    "TRDV", "TRDD", "TRDJ",
    "TRGV", "TRGJ",
)

_DOI = "10.1016/j.immuni.2023.01.026"


def fetch_ki_tcr(cache_dir: str | Path) -> list[Path]:
    """Download the ten FASTA files into ``cache_dir``. Idempotent; network.

    A file that comes back empty or without a FASTA header is an error, not an empty
    set: the URLs are fixed, so a blank response means the site moved or an error page
    was served, and silently ingesting nothing is how a source disappears unnoticed.
    """
    import urllib.request

    out = Path(cache_dir)
    out.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    for stem in _FILES:
        dest = out / f"{stem}.fasta"
        if not dest.exists() or dest.stat().st_size == 0:
            urllib.request.urlretrieve(_BASE.format(stem=stem), dest)
        text = dest.read_text(errors="replace").lstrip()
        if not text.startswith(">"):
            raise ValueError(
                f"KI-TCR {stem}: no FASTA content at {_BASE.format(stem=stem)} "
                f"(got {text[:80]!r}). Refusing to ingest an empty set."
            )
        paths.append(dest)
    return paths


def ingest_ki_tcr(session: Session, cache_dir: str | Path) -> dict[str, int]:
    """Ingest the cached KI-TCR FASTAs. Idempotent; returns per-outcome counts.

    Every sequence is upserted on its hash, so the ~70% that IMGT already publishes
    are corroborated with a second source and a KI-TCR name rather than duplicated.
    """
    tag = normalize_species("human")
    source = loader.get_or_create_source(
        session, name="KI-TCR", version=_DOI, curie_prefix="KI-TCR",
        note="KI-TCR — Karolinska Institutet human T-cell receptor database",
    )
    stats = {"coding": 0, "skipped_segment": 0, "skipped_sequence": 0, "files": 0}
    for path in sorted(Path(cache_dir).glob("*.fasta")):
        if path.stem not in _FILES:
            continue
        stats["files"] += 1
        for header, seq in read_fasta(path):
            name = header.split()[0].strip() if header else ""
            if not name or "*" not in name or not seq:
                stats["skipped_sequence"] += 1
                continue
            segment = segment_of(name)
            locus = locus_of(name)
            # Guard both, not just one: segment_of must not default an unrecognised
            # name to "V", and a TR file must not smuggle in a non-TR locus.
            if segment not in ("V", "D", "J") or locus not in (
                    "TRA", "TRB", "TRG", "TRD"):
                stats["skipped_segment"] += 1
                continue
            try:
                coding = loader.upsert_coding(
                    session, species=tag, segment=segment, locus=locus, nt=seq,
                )
            except SequenceError:
                stats["skipped_sequence"] += 1
                continue
            loader.add_name(session, coding, source, name, scheme_label(name))
            loader.add_attribution(
                session, coding, source,
                original_id=curie_build(source.curie_prefix or "KI-TCR", name),
                original_sequence=seq,
            )
            stats["coding"] += 1
    return stats
