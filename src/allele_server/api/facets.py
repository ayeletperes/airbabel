"""Result-set facets for the UI.

Pure, transport-free aggregation over an already-computed :class:`SearchResults`
(no DB, no search logic) so the results table can show a per-species distribution
and drive in-page filtering. Lives apart from ``queries.py`` deliberately: it only
*reads* search output and never touches how similarity is computed or ranked.
"""

from __future__ import annotations

from collections import Counter

from ..species import common_name
from .schemas import SearchResults


def _type_of(hit) -> str | None:
    """Display 'type' for a hit: the V/D/J segment for coding, else the category."""
    if hit.category and hit.category != "coding":
        return hit.category
    return hit.sequence_type or None


def result_facets(results: SearchResults | None) -> dict:
    """Summarize ``results.hits`` for the results-table summary panel.

    Returns counts by species / locus / type, the per-species list of identity
    scores (for the strip plot), and the identity axis bounds. Safe on empty or
    ``None`` input (returns ``total == 0`` with empty collections).
    """
    hits = results.hits if results is not None else []
    if not hits:
        return {
            "total": 0,
            "by_species": [],
            "by_locus": [],
            "by_type": [],
            "points_by_species": {},
            "identity_min": 0.8,
            "identity_max": 1.0,
        }

    sp_counts: Counter[str] = Counter()
    locus_counts: Counter[str] = Counter()
    type_counts: Counter[str] = Counter()
    # Per species, one point per hit (uid + label + score) so the strip plot can
    # show a hover tooltip and jump to the matching table row on click.
    points_by_species: dict[str, list[dict]] = {}

    for h in hits:
        sp_counts[h.species] += 1
        points_by_species.setdefault(h.species, []).append(
            {"uid": h.uid, "label": h.label or h.uid, "score": h.score}
        )
        if h.locus:
            locus_counts[h.locus] += 1
        t = _type_of(h)
        if t:
            type_counts[t] += 1

    def _ranked(counter: Counter[str], *, with_label: bool = False) -> list[dict]:
        # Sorted by count desc, then key asc for a stable order.
        out = []
        for key, count in sorted(counter.items(), key=lambda kv: (-kv[1], kv[0])):
            entry = {"key": key, "count": count}
            if with_label:
                entry["label"] = common_name(key)
            out.append(entry)
        return out

    scores = [h.score for h in hits]
    return {
        "total": len(hits),
        "by_species": _ranked(sp_counts, with_label=True),
        "by_locus": _ranked(locus_counts),
        "by_type": _ranked(type_counts),
        "points_by_species": points_by_species,
        "identity_min": min(scores),
        "identity_max": max(scores),
    }
