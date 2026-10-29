"""The review report: what the automation hands a human before they decide.

One markdown file per candidate, structured so the decision fits on one screen — the
nine checklist items from ``docs/source-policy.md`` with the machine's answer pre-filled,
the proposed citation, and the commands to act on it.

Regenerable from ``source_candidates.evidence``, so the file is a cache and the
committed audit record is the manifest's ``[approval]`` block.

**Where the machine could not answer, the dossier says so.** An absent number must never
read as a zero: "no sequence file was reachable" and "nothing overlaps" are opposite
conclusions, and a reviewer skimming a table cannot tell them apart unless we say which.
"""

from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

from sqlalchemy.orm import Session

from ..db.models import SourceCandidate
from ..settings import get_settings
from ..species import SpeciesError, normalize_species

__all__ = ["dossier_path", "render", "write_dossier"]

_UNKNOWN = "_not determined — the automation could not answer this; check by hand_"


def _slug(value: str) -> str:
    out = "".join(c if c.isalnum() or c in "-_." else "_" for c in value)
    return out.strip("_") or "candidate"


def dossier_path(candidate: SourceCandidate) -> Path:
    root = Path(get_settings().source_dossier_dir)
    return root / candidate.channel / f"{_slug(candidate.external_id)}.md"


def _species_line(candidate: SourceCandidate) -> str:
    if not candidate.species_hint:
        return _UNKNOWN
    parts = []
    for raw in (s.strip() for s in candidate.species_hint.split(",") if s.strip()):
        try:
            parts.append(f"{normalize_species(raw)} (recognised)")
        except SpeciesError:
            parts.append(f"**{raw} — not in species._SPECIES; add it first**")
    return " · ".join(parts)


def _citation_block(candidate: SourceCandidate) -> str:
    if candidate.doi:
        return (
            f"DOI [{candidate.doi}](https://doi.org/{candidate.doi})\n\n"
            f"**Verify the full reference before approving.** The manifest needs the "
            f"citation text, not just the identifier."
        )
    return (
        "No DOI was found on this channel. Look for the paper or the deposit's own "
        "'how to cite' statement. If there genuinely is none, the manifest records "
        "`citation.status = \"none_available\"` **with a written justification** — see "
        "docs/source-policy.md."
    )


def render(candidate: SourceCandidate, *, generated_at: datetime | None = None) -> str:
    """The dossier markdown for a candidate. Pure; ``write_dossier`` persists it."""
    stamp = (generated_at or datetime.now(UTC)).strftime("%Y-%m-%d")
    evidence = json.loads(candidate.evidence) if candidate.evidence else {}
    title = candidate.title or candidate.external_id
    score = f"{candidate.score:.2f}" if candidate.score is not None else "n/a"

    lines = [
        f"# {title}",
        "",
        "| | |",
        "|---|---|",
        f"| candidate | `{candidate.ref}` |",
        f"| status | `{candidate.status}` (score {score}) |",
        f"| first seen | {candidate.first_seen_at:%Y-%m-%d} via {candidate.channel} |",
        f"| url | {candidate.url or '—'} |",
        f"| dossier generated | {stamp} |",
        "",
        "## Decide",
        "",
        "Work these in order; items 1-3 can end the review on their own.",
        "",
        f"- [ ] **1. Provenance** — {title}. "
        "Is this a primary germline set, or a repackaging of one we already hold? "
        "(Mirrors of IMGT and OGRDB are the common false positive.)",
        "- [ ] **2. Licence** — "
        + (f"declared as **{candidate.licence}**" if candidate.licence
           else "**none found by the crawler.** Check the deposit page yourself; if "
                "there really is none, that is `none_declared`, not a guess."),
        "- [ ] **3. Redistributable** — follows from (2). An undeclared licence means "
        "`redistributable = false`: we may ingest and credit it, not republish it.",
        f"- [ ] **4. Species / loci** — {_species_line(candidate)}"
        + (f" · loci {candidate.loci_hint}" if candidate.loci_hint else ""),
        "- [ ] **5. Overlap with what we hold** — "
        "_no sequence file was fetched for this candidate, so overlap is unmeasured._ "
        "Fetch the deposit and compare on `sequence.hash_sequence` before deciding: a "
        "source that is entirely already-held is a corroboration, not a new source.",
        "- [ ] **6. Sequence quality** — do the sequences pass `sequence.normalize` "
        "(IUPAC, no RNA)? Is the length distribution sane for the segment? Any "
        "duplicates inside the deposit?",
        "- [ ] **7. Naming scheme** — what does `naming.detect_scheme` say? A source "
        "claiming IMGT nomenclature it does not own gets bucketed to its own scheme.",
        "- [ ] **8. Citation** — see below. This is the item that most often blocks an "
        "approval, and it should.",
        "- [ ] **9. Ingest cost** — which module loads it; does it need a new species, "
        "naming scheme or parser?",
        "",
        "## Proposed citation",
        "",
        _citation_block(candidate),
        "",
        "## Evidence trail",
        "",
    ]

    fetches = evidence.get("fetches") or []
    if fetches:
        lines += ["| fetched | url | cache |", "|---|---|---|"]
        for f in fetches:
            lines.append(
                f"| {f.get('at', '—')} | {f.get('url', '—')} | `{f.get('cache', '—')}` |"
            )
    else:
        lines.append(
            "_No fetch record._ This candidate was recorded by hand, so there is no "
            "crawl trail to show."
        )
    if evidence.get("raw"):
        lines += [
            "",
            "<details><summary>Raw channel record</summary>",
            "",
            "```json",
            json.dumps(evidence["raw"], indent=2, sort_keys=True),
            "```",
            "",
            "</details>",
        ]

    slug = _slug(candidate.external_id).replace("-", "_").replace(".", "_").casefold()
    lines += [
        "",
        "## Act",
        "",
        "```bash",
        f"allele-sources draft   {candidate.ref} --manifest {slug}",
        f"$EDITOR src/allele_server/sources/{slug}.toml   # licence + citation",
        f"allele-sources approve {candidate.ref} --manifest {slug} "
        f"--by 'Your Name' --note '...'",
        f"allele-sources reject  {candidate.ref} --by 'Your Name' --note 'why'",
        f"allele-sources defer   {candidate.ref} --by 'Your Name' --until YYYY-MM-DD "
        f"--note 'why'",
        "```",
        "",
    ]
    return "\n".join(lines)


def write_dossier(
    session: Session, candidate: SourceCandidate, *, generated_at: datetime | None = None
) -> Path:
    """Write the dossier and record its path on the candidate.

    Promotes a ``discovered``/``triaged`` candidate to ``dossier`` — the first state a
    human is asked to look at. A decided candidate keeps its status: regenerating the
    report must not reopen a closed decision.
    """
    # Promote before rendering, so the status printed in the report is the status the
    # database holds. A dossier that says "discovered" about a row that says "dossier"
    # sends the reviewer looking for a bug that is not there.
    if candidate.status in ("discovered", "triaged"):
        candidate.status = "dossier"
    path = dossier_path(candidate)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(render(candidate, generated_at=generated_at), encoding="utf-8")
    candidate.dossier_path = str(path)
    session.flush()
    return path
