"""Token check + per-caller rate limit in front of the ``/mcp`` endpoint.

``allele_server.mcp.server`` was written for "a trusted local/desktop MCP client,
not a public endpoint", and said so, with a note to add a token check and a
per-caller limit if that ever changed. It changed: the endpoint is live on
airbabel.org and the URL is being handed to collaborators.

The public endpoint is **open by design**: the data is public and read-only, and a
token nobody can get without emailing the maintainer would keep AIRBabel out of the
agents it is meant to be reachable from. What needs guarding is CPU, so the controls
are rate limits, with the token kept only for a private deployment:

* **Token.** ``Authorization: Bearer <token>``, compared with
  :func:`hmac.compare_digest` so a wrong guess costs the same time as a right one.
  Unset (the default, and airbabel.org's setting) means open; setting
  ``MCP_TOKEN`` turns the check on, so the deployment decides rather than the code.
* **Similarity budget.** The similarity tools (:data:`HEAVY_TOOLS`) are the ones
  that can fall through to a live scan. They get a tighter per-caller limit plus a
  global one, so many callers together cannot starve the website either. Every
  other tool is an indexed lookup and only counts against the general limit.
* **Rate limit.** A fixed window per caller. Not because an attacker is the worry
  -- the data is public and read-only -- but because this box is a two-core
  burstable instance sharing CPU credits with the website, and one agent looping
  over similarity queries can throttle airbabel.org for everyone. The same reason
  the Caddyfile already blocks crawlers from the FASTA export.

Refusals are JSON-RPC errors, not bare HTTP codes: an MCP client shows the
``message`` to its user, so "set MCP_TOKEN" reaches a person instead of a 401
that surfaces as "the server is broken".
"""

from __future__ import annotations

import hmac
import json
import time
from collections import defaultdict, deque
from threading import Lock

#: JSON-RPC codes. -32000 is the generic server error the spec reserves for
#: implementation-defined conditions, which is what both of these are.
_ERR_AUTH = -32001
_ERR_RATE = -32002

#: MCP tools that may run a live similarity scan (a neighbour-cache miss falls
#: through to the ranker). Names must match ``allele_server.mcp.server``.
HEAVY_TOOLS = frozenset({"find_similar", "get_cross_species_neighbors"})
#: A JSON-RPC request carries a name or a sequence; anything bigger is not one.
_MAX_BODY = 1 << 20


class RateLimiter:
    """A fixed-window request count per caller. Shared by /mcp and the bulk export.

    Extracted rather than written twice: the two endpoints have different budgets
    but the same question, and two copies of a counter drift.

    Not distributed, and deliberately so — one process serves this site, so a dict
    and a lock are the whole mechanism. A second worker would need something
    shared; ``ponytail: in-process only, move to Redis if this is ever > 1 worker``.
    """

    def __init__(self, limit: int, window: float):
        self.limit = max(0, limit)        # 0 disables
        self.window = max(1.0, window)
        self._hits: dict[str, deque[float]] = defaultdict(deque)
        self._lock = Lock()

    def over(self, caller: str) -> bool:
        """True when this caller has spent its budget; records the hit otherwise."""
        if not self.limit:
            return False
        now = time.monotonic()
        cutoff = now - self.window
        with self._lock:
            hits = self._hits[caller]
            while hits and hits[0] < cutoff:
                hits.popleft()
            if len(hits) >= self.limit:
                return True
            hits.append(now)
            # Callers that have gone quiet are dropped, so a long-running process
            # cannot grow this without bound.
            if len(self._hits) > 1024:
                for key in [k for k, v in self._hits.items() if not v]:
                    del self._hits[key]
        return False


def caller_of(scope_or_headers, client=None) -> str:
    """Who to count against: the forwarded client IP, else the socket peer.

    Behind Caddy and the load balancer every socket peer is the proxy, so counting
    those makes one bucket for the whole internet — and the first busy caller then
    locks everyone out. The left-most X-Forwarded-For entry is the original client.
    It is spoofable, but this guards CPU rather than secrets, and the alternative
    fails closed for everybody.
    """
    if isinstance(scope_or_headers, dict) and "headers" in scope_or_headers:
        headers = {k.lower(): v for k, v in scope_or_headers.get("headers", [])}
        fwd = headers.get(b"x-forwarded-for", b"").decode(errors="replace")
        client = scope_or_headers.get("client")
    else:  # a Starlette request's headers (str -> str)
        fwd = scope_or_headers.get("x-forwarded-for", "") or ""
    if fwd.strip():
        return fwd.split(",")[0].strip()
    if client:
        return client[0] if isinstance(client, (tuple, list)) else str(client)
    return "unknown"


class McpGate:
    """ASGI middleware: authenticate, then rate-limit, then delegate.

    Wraps the MCP sub-app rather than living in FastAPI middleware, so it cannot
    be bypassed by a route added later and does not run on any other path.
    """

    def __init__(self, app, *, token: str = "", limit: int = 0, window: float = 60.0,
                 heavy_limit: int = 0, heavy_global_limit: int = 0):
        self.app = app
        self.token = token or ""
        self.limiter = RateLimiter(limit, window)
        self.heavy = RateLimiter(heavy_limit, window)
        self.heavy_global = RateLimiter(heavy_global_limit, window)

    def _authorised(self, scope) -> bool:
        if not self.token:
            return True
        headers = {k.lower(): v for k, v in scope.get("headers", [])}
        raw = headers.get(b"authorization", b"").decode(errors="replace")
        prefix = "bearer "
        if not raw.lower().startswith(prefix):
            return False
        return hmac.compare_digest(raw[len(prefix):].strip(), self.token)

    # --- ASGI ----------------------------------------------------------------
    async def __call__(self, scope, receive, send):
        if scope["type"] != "http":
            await self.app(scope, receive, send)
            return
        if not self._authorised(scope):
            await _refuse(
                send, 401, _ERR_AUTH,
                "AIRBabel MCP requires a token. Set it in your MCP client as the "
                "header 'Authorization: Bearer <token>'. Ask the maintainer for one.",
            )
            return
        caller = caller_of(scope)
        if self.limiter.over(caller):
            await _refuse(
                send, 429, _ERR_RATE,
                f"Rate limit reached ({self.limiter.limit} requests per "
                f"{int(self.limiter.window)}s). AIRBabel runs on a small shared "
                f"instance; "
                f"please slow down and retry shortly.",
            )
            return
        calls: list[str] = []
        if scope.get("method") == "POST":
            body, receive = await _buffer(receive)
            if body is None:
                await _refuse(send, 413, _ERR_RATE, "Request body too large.")
                return
            calls = _mcp_calls(body)
            send = _tag_response(send, calls)
        heavy = sum(c.removeprefix("tools/call:") in HEAVY_TOOLS
                    for c in calls if c.startswith("tools/call:"))
        if self.heavy.limit or self.heavy_global.limit:
            for _ in range(heavy):
                # Per caller first, so one busy agent spends its own budget and not
                # the shared one.
                if self.heavy.over(caller) or self.heavy_global.over("*"):
                    await _refuse(
                        send, 429, _ERR_RATE,
                        f"Similarity-search limit reached ({self.heavy.limit} per "
                        f"caller, {self.heavy_global.limit} in total, per "
                        f"{int(self.heavy.window)}s). Name and sequence lookups are "
                        f"not affected; retry the similarity call shortly.",
                    )
                    return
        await self.app(scope, receive, send)


async def _buffer(receive):
    """Read the whole request body, and a ``receive`` that replays it downstream.

    ``(None, receive)`` when the body exceeds :data:`_MAX_BODY`.
    """
    chunks, size = [], 0
    while True:
        msg = await receive()
        if msg["type"] != "http.request":  # client went away; let the app see it
            break
        chunks.append(msg.get("body", b""))
        size += len(chunks[-1])
        if size > _MAX_BODY:
            return None, receive
        if not msg.get("more_body"):
            break
    body = b"".join(chunks)
    replayed = False

    async def replay():
        nonlocal replayed
        if not replayed:
            replayed = True
            return {"type": "http.request", "body": body, "more_body": False}
        return await receive()

    return body, replay


def _mcp_calls(body: bytes) -> list[str]:
    """The JSON-RPC methods in this body, with the tool name for ``tools/call``:
    e.g. ``["initialize"]`` or ``["tools/call:find_similar"]``. Arguments are never
    included, so nothing a caller searched for ends up in a log.

    A body that does not parse is not ours to judge: the MCP app rejects it.
    """
    try:
        msg = json.loads(body)
    except ValueError:
        return []
    calls = []
    for m in msg if isinstance(msg, list) else [msg]:
        if not isinstance(m, dict) or not isinstance(m.get("method"), str):
            continue
        params = m.get("params") if isinstance(m.get("params"), dict) else {}
        name = params.get("name") if m["method"] == "tools/call" else None
        calls.append(f"tools/call:{name}" if isinstance(name, str) else m["method"])
    return calls


#: Response header naming the MCP calls in the request. Caddy's JSON access log
#: records response headers, so this is how ``scripts/ai_traffic.py`` sees which
#: tools are used without the app keeping a log of its own.
CALL_HEADER = b"x-mcp-call"


def _tag_response(send, calls: list[str]):
    """``send`` that adds :data:`CALL_HEADER` to the response start."""
    if not calls:
        return send
    value = ",".join(calls)[:200].encode()

    async def tagged(msg):
        if msg["type"] == "http.response.start":
            msg = {**msg, "headers": [*msg.get("headers", []), (CALL_HEADER, value)]}
        await send(msg)

    return tagged


async def _refuse(send, status: int, code: int, message: str) -> None:
    """One JSON-RPC error response, so the client shows the reason to its user."""
    body = json.dumps(
        {"jsonrpc": "2.0", "id": None, "error": {"code": code, "message": message}}
    ).encode()
    await send({
        "type": "http.response.start",
        "status": status,
        "headers": [(b"content-type", b"application/json"),
                    (b"content-length", str(len(body)).encode())],
    })
    await send({"type": "http.response.body", "body": body})
