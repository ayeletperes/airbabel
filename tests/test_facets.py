"""Unit tests for result_facets (pure aggregation over SearchResults; no DB)."""

from __future__ import annotations

from allele_server.api.facets import result_facets
from allele_server.api.schemas import SearchHit, SearchResults
from allele_server.species import common_name


def _hit(species, locus, score, label, seqtype="V"):
    return SearchHit(
        uid=f"U-{label}", label=label, species=species, category="coding",
        locus=locus, sequence_type=seqtype, score=score, matched_by="similar",
    )


def test_result_facets_basic():
    res = SearchResults(query="x", count=4, capped=False, hits=[
        _hit("HUMAN", "IGH", 0.99, "a"),
        _hit("HUMAN", "IGH", 0.95, "b"),
        _hit("RHESUS", "IGH", 0.90, "c"),
        _hit("MOUSE", "IGK", 0.82, "d"),
    ])
    f = result_facets(res)
    assert f["total"] == 4
    # species counts sum to total; most-common first; labels via common_name
    assert sum(s["count"] for s in f["by_species"]) == 4
    assert f["by_species"][0] == {"key": "HUMAN", "count": 2, "label": common_name("HUMAN")}
    # locus + per-species points (uid/label/score) + axis bounds
    assert {loc["key"] for loc in f["by_locus"]} == {"IGH", "IGK"}
    human_pts = f["points_by_species"]["HUMAN"]
    assert len(human_pts) == 2
    assert human_pts[0] == {"uid": "U-a", "label": "a", "score": 0.99}
    assert f["identity_min"] == 0.82 and f["identity_max"] == 0.99


def test_result_facets_empty():
    assert result_facets(None)["total"] == 0
    empty = result_facets(SearchResults(query="x", count=0, capped=False, hits=[]))
    assert empty["total"] == 0 and empty["by_species"] == []
