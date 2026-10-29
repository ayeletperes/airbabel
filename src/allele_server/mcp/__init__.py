"""MCP server: a thin, read-only wrapper over :mod:`allele_server.queries`.

Locked decision #10 made the query core transport-agnostic, so this package is a
transport shim only — it adds no resolution or similarity logic of its own. See
:mod:`allele_server.mcp.server`.
"""

from .server import main, server

__all__ = ["main", "server"]
