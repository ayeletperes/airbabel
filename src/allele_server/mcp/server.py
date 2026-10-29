"""MCP tools over the allele reference DB.

Every tool is a thin wrapper: it opens a read-only session, calls one existing
:mod:`allele_server.queries` function, and returns the resulting Pydantic model
as plain JSON. No new resolution/similarity logic lives here — the query core is
already transport-agnostic (locked decision #10), so REST, the web UI and this
server all share exactly the same code path.

All tools are read-only. Resolution is one indexed lookup; the similarity tools
read the precomputed neighbour cache and can fall back to a live scan. When mounted
on the site, ``api.mcp_gate.McpGate`` rate-limits callers, with a tighter budget on
the similarity tools (``HEAVY_TOOLS`` there must match the tool names here).
"""

from __future__ import annotations

from typing import Any

from fastapi import HTTPException
from mcp.server.mcpserver import MCPServer
from sqlalchemy import select

from .. import queries
from ..api.routes import maybe_partial_slot
from ..db.models import Sequence
from ..db.session import session_scope
from ..settings import get_settings

server = MCPServer("allele-server")


def _dump(model: Any) -> Any:
    """Pydantic model (or list of them) -> JSON-ready plain data."""
    if isinstance(model, list):
        return [m.model_dump(mode="json") for m in model]
    return model.model_dump(mode="json")


def _resolve(query: str, species: str | None) -> dict:
    with session_scope() as session:
        result = queries.resolve(session, query, species=species)
        if result is None:
            raise ValueError(f"no match for {query!r}")
        return _dump(result)


def _similar(
    uid: str, min_identity: float, kind: str, same_species: bool, to_species: str | None
) -> dict:
    """One similarity call, holding a slot on this box's expensive-query budget.

    ``find_similar`` is served from the precomputed neighbour cache, but a miss
    falls through to the live ranker — the same full-scan-and-rescore cost that
    ``maybe_partial_slot`` exists to bound for partial search. There is one CPU on
    the serve box, so MCP callers draw on the same semaphore as the REST routes and
    the web UI rather than getting a private, unbounded lane around it. The flag is
    named for that caller; passing True here means "this is one of the expensive
    ones", which every similarity call is.
    """
    fn = queries.find_similar_aa if kind == "aa" else queries.find_similar
    try:
        with maybe_partial_slot(True), session_scope() as session:
            result = fn(
                session, uid, min_identity=min_identity,
                cap=get_settings().search_result_cap,
                same_species=same_species, target_species=to_species,
            )
            if result is None:
                raise ValueError(f"no sequence with uid {uid!r}")
            return _dump(result)
    except HTTPException as exc:  # 503: no slot freed within the queue timeout
        raise ValueError(f"server busy, retry shortly ({exc.detail})") from exc


@server.tool()
def resolve_allele(query: str, species: str | None = None) -> dict:
    """Resolve an allele name, UID or protein UID to its one canonical record.

    Names are labels, not keys: the same name can be reused across species, so
    pass ``species`` (e.g. "human", or "human,rhesus") to pin an ambiguous query.
    Returns the matched record plus how it matched (name / uid / protein_uid).
    """
    return _resolve(query, species)


@server.tool()
def resolve_sequence(sequence: str, species: str | None = None) -> dict:
    """Resolve a raw nucleotide or amino-acid sequence to its canonical record.

    Identity is by sequence hash (and reverse-complement / AA hash), so this is
    the authoritative "is this sequence already known, and under what name" call.
    ``species`` scopes the match (identity is per-species). Same underlying
    resolver as ``resolve_allele`` — kind is auto-detected from the input.
    """
    return _resolve(sequence, species)


@server.tool()
def find_similar(
    uid: str,
    min_identity: float = 0.85,
    kind: str = "nt",
    same_species: bool = True,
    to_species: str | None = None,
) -> dict:
    """Nearest coding alleles to ``uid`` by identity (self excluded).

    ``kind`` is "nt" (default) or "aa". Set ``same_species=False`` or name
    ``to_species`` for cross-species neighbours. Served from the precomputed
    neighbour cache with a live fallback; results are capped by the server's
    configured search cap.
    """
    return _similar(uid, min_identity, kind, same_species, to_species)


@server.tool()
def get_cross_species_neighbors(
    uid: str,
    min_identity: float = 0.85,
    kind: str = "nt",
    to_species: str | None = None,
) -> dict:
    """Nearest alleles to ``uid`` in *other* species — the closest-ortholog lookup.

    Optionally restrict to one or more target species via ``to_species``
    ("gorilla", "human,mouse"). This is a sequence-based signal, not a curated
    orthology call.

    The query's own species is excluded from ``hits``; ``dropped_same_species``
    reports how many were removed, so the scan is never silently narrowed.
    """
    result = _similar(uid, min_identity, kind, same_species=False, to_species=to_species)
    with session_scope() as session:
        own = session.scalar(select(Sequence.species).where(Sequence.uid == uid))
    kept = [h for h in result["hits"] if h["species"] != own]
    # The underlying search caps BEFORE this filter, so a query whose same-species
    # neighbours fill the cap can return few (or no) cross-species hits even though
    # more exist. `capped` + `dropped_same_species` together make that visible.
    # ponytail: post-filter; push the exclusion into queries.find_similar only if
    #           capped-out cross-species queries turn out to be common in practice.
    result["dropped_same_species"] = len(result["hits"]) - len(kept)
    result["hits"] = kept
    result["count"] = len(kept)
    return result


@server.tool()
def get_aliases(uid: str) -> list[dict]:
    """Every known alias (name + scheme + source) for the sequence with this UID.

    The aliases path behind ``GET /api/sequences/{uid}/aliases``. Raises if the
    UID is unknown.
    """
    with session_scope() as session:
        result = queries.list_aliases(session, uid)
        if not result:
            raise ValueError(f"no sequence with uid {uid!r}")
        return _dump(result)


def main() -> None:
    """Run the MCP server over stdio (the default MCP transport)."""
    server.run()


if __name__ == "__main__":
    main()
