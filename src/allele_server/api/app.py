"""FastAPI application: REST API + a minimal Jinja2 web UI."""

from __future__ import annotations

import time
from pathlib import Path
from urllib.parse import quote, urlencode

from fastapi import Depends, FastAPI, HTTPException, Query, Request
from fastapi.responses import FileResponse, HTMLResponse, RedirectResponse, Response
from fastapi.templating import Jinja2Templates
from sqlalchemy.orm import Session

from .. import __version__, queries
from ..examples_meta import SEARCH_EXAMPLES
from ..settings import get_settings
from .facets import result_facets
from .mcp_gate import McpGate
from .routes import get_session, maybe_partial_slot, router


_TEMPLATES = Jinja2Templates(directory=str(Path(__file__).parent / "web" / "templates"))
# Shown inline on /tutorial and served as a download from /tutorial/airbabel_example.py.
# One file, read at request time, so the page and the download can never drift apart --
# and so the script stays a real, runnable, lintable module rather than a string in a
# template.
_EXAMPLE_SCRIPT = Path(__file__).parent / "web" / "static" / "airbabel_example.py"
# Available in every template (e.g. the header nav) without per-route context.
_TEMPLATES.env.globals["feedback_url"] = get_settings().feedback_url
_TEMPLATES.env.globals["version"] = __version__
_TEMPLATES.env.globals["public_url"] = get_settings().public_url.rstrip("/")
_STATIC = Path(__file__).parent / "web" / "static"
_RELEASE: list = [float("-inf"), None]


def _data_release() -> dict | None:
    """The stamped dataset release, for the footer of every page.

    Looked up at most once a minute: one indexed row, but on every page view.
    ``None`` on a database never stamped; the footer then says so rather than
    showing a date it does not have.
    """
    now = time.monotonic()
    if now - _RELEASE[0] > 60:
        from ..db.session import session_scope

        with session_scope() as s:
            _RELEASE[:] = [now, queries.dataset_release(s)]
    return _RELEASE[1]


_TEMPLATES.env.globals["data_release"] = _data_release
# "Try an example" chips. The URL is built here rather than in the template so the
# encoding is done once, at import, by the same urlencode the /examples page uses.
_TEMPLATES.env.globals["search_examples"] = [
    {**ex, "url": ex.get("url") or ("/?" + urlencode(ex["params"]))}
    for ex in SEARCH_EXAMPLES
]


def _species_sel(species: list[str] | None) -> list[str]:
    """The species tags a multi-select should show as chosen.

    Accepts what either form of the parameter delivers — repeated
    ``?species=human&species=rhesus`` or the older comma form
    ``?species=human,rhesus`` — so an existing link keeps selecting the right
    options rather than silently reverting the picker to "all species".
    Upper-cased because the options carry canonical tags; genuine normalisation
    (aliases, common names) happens in ``queries._species_tags``, which is the one
    that must raise on an unknown species.
    """
    return [
        part.strip().upper()
        for raw in (species or [])
        for part in raw.split(",")
        if part.strip()
    ]


def _here(request: Request) -> str:
    """This page as a same-site relative URL, minus any ``back`` it was given.

    Used as the ``back`` target of links leaving the page. Dropping the incoming
    ``back`` is what keeps the chain finite: embedding it made every hop a longer,
    never-seen URL, and a crawler ignoring robots.txt walked that tree at ~70k
    expensive /similar requests a day.
    """
    query = urlencode([(k, v) for k, v in request.query_params.multi_items() if k != "back"])
    return request.url.path + (f"?{query}" if query else "")


def _safe_back(back: str | None) -> str | None:
    """``back`` only when it is a same-site path; it is rendered as a link's href,
    so anything else (``javascript:``, another host, ``//host``) is dropped."""
    if back and back.startswith("/") and not back.startswith("//"):
        return back
    return None


def _mcp_mount():
    """The MCP streamable-HTTP sub-app + the lifespan that runs its session
    manager, or ``(None, None)`` when the optional ``mcp`` extra is not installed
    (the site then serves normally, just without ``/mcp``). Mounting it here gives
    a remote MCP endpoint on the same host and DB as the API, so other people only
    need the URL — no clone, no local DB, no keys."""
    try:
        from allele_server.mcp.server import server as mcp_server
    except ImportError:
        return None, None
    from contextlib import asynccontextmanager

    from mcp.server.transport_security import TransportSecuritySettings

    # path="/" so the mount point (/mcp) is the whole endpoint, not /mcp/mcp.
    #
    # DNS-rebinding protection defaults to a 127.0.0.1-only Host allowlist, which
    # 421s every request once Caddy fronts this as airbabel.org. Name the hostnames
    # instead of turning the guard off: the endpoint is public and read-only, but it
    # is also the CPU-expensive one, and with the guard off any page a visitor loads
    # can drive POSTs at it from their browser. CORS would stop that page reading the
    # response; it does not stop the request running, so the cost lands anyway.
    #
    # allowed_origins stays empty deliberately: a missing Origin passes (real MCP
    # clients send none), while anything browser-originated is refused.
    mcp_app = mcp_server.streamable_http_app(  # + session manager
        streamable_http_path="/",
        transport_security=TransportSecuritySettings(
            allowed_hosts=get_settings().mcp_allowed_hosts,
        ),
    )

    @asynccontextmanager
    async def lifespan(_app: FastAPI):
        async with mcp_server.session_manager.run():
            yield

    return mcp_app, lifespan


class _HeadAsGet:
    """Answer HEAD on the website pages by running GET and dropping the body.

    FastAPI routes answer GET only, so HEAD got a 405, and SEO tools and link
    checkers that probe with HEAD reported the site as broken. /api/ and /mcp are
    left alone: HEAD on /api/export/* would run a full-table scan for nothing, and
    GET on /mcp opens an event stream.
    """

    def __init__(self, app):
        self.app = app

    async def __call__(self, scope, receive, send):
        if (scope["type"] != "http" or scope["method"] != "HEAD"
                or scope["path"].startswith(("/api/", "/mcp"))):
            await self.app(scope, receive, send)
            return

        async def headers_only(msg):
            if msg["type"] == "http.response.body":
                if msg.get("more_body"):
                    return
                msg = {"type": "http.response.body", "body": b"", "more_body": False}
            await send(msg)

        await self.app({**scope, "method": "GET"}, receive, headers_only)


def create_app() -> FastAPI:
    mcp_app, mcp_lifespan = _mcp_mount()
    app = FastAPI(
        title="AIRBabel",
        version=__version__,
        description=(
            "Allele / leader / RSS reference with deterministic UIDs, alias "
            "resolution, AA + similarity search, AIRR-aligned responses."
        ),
        lifespan=mcp_lifespan,
    )
    app.include_router(router)
    app.add_middleware(_HeadAsGet)

    # dataset_stats is ~4 s of counting and backs the landing page; the data only
    # changes on a reload, so hold results per process, keyed on the release stamp
    # (a one-row lookup) so a stamped reload is picked up on the next request.
    # ponytail: a reload that is never stamped stays cached until restart.
    _cache: dict[str, tuple[object, object]] = {}

    def _cached(key, session, fn):
        rel = queries.dataset_release(session)
        stamp = rel and rel["built_at"]
        hit = _cache.get(key)
        if hit is None or hit[0] != stamp:
            hit = _cache[key] = (stamp, fn())
        return hit[1]

    @app.get("/sitemap.xml")
    def sitemap(session: Session = Depends(get_session)):
        """Every current record page plus the static pages, for search engines."""
        from xml.sax.saxutils import escape

        base = get_settings().public_url.rstrip("/")
        prefix = {"coding": "/seq/", "rss": "/rss/", "leader": "/leader/"}

        def build() -> str:
            paths = ["/", "/overview", "/sources", "/examples", "/tutorial", "/help"]
            paths += [prefix[cat] + uid for uid, cat in queries.sitemap_uids(session)
                      if cat in prefix]
            # 50,000 URLs is the protocol's per-file limit; split into an index past it.
            if len(paths) > 50_000:
                raise RuntimeError(f"sitemap has {len(paths)} URLs; split it")
            urls = "".join(f"<url><loc>{escape(base + p)}</loc></url>" for p in paths)
            return ('<?xml version="1.0" encoding="UTF-8"?>'
                    '<urlset xmlns="http://www.sitemaps.org/schemas/sitemap/0.9">'
                    f"{urls}</urlset>")

        return Response(_cached("sitemap", session, build), media_type="application/xml")

    @app.get("/llms.txt")
    def llms_txt():
        """Plain-text site guide for LLM crawlers and agents (llmstxt.org)."""
        return FileResponse(_STATIC / "llms.txt", media_type="text/plain; charset=utf-8")

    indexnow_key = get_settings().indexnow_key
    if indexnow_key:
        @app.get(f"/{indexnow_key}.txt", include_in_schema=False)
        def indexnow_key_file():
            """IndexNow ownership proof; see allele_server.indexnow."""
            return Response(indexnow_key, media_type="text/plain")

    @app.get("/favicon.svg")
    def favicon():
        return FileResponse(_STATIC / "favicon.svg", media_type="image/svg+xml")

    @app.get("/favicon.ico", include_in_schema=False)
    def favicon_ico():
        return RedirectResponse(url="/favicon.svg", status_code=301)
    if mcp_app is not None:
        # Wrap, don't mount bare: the gate authenticates and rate-limits before the
        # MCP app sees a byte, and wrapping means no route added later can bypass it.
        settings = get_settings()
        app.mount("/mcp", McpGate(  # remote MCP over Streamable HTTP
            mcp_app,
            token=settings.mcp_token,
            limit=settings.mcp_rate_limit,
            window=settings.mcp_rate_window,
            heavy_limit=settings.mcp_similar_rate_limit,
            heavy_global_limit=settings.mcp_similar_global_limit,
        ))

    @app.get("/help", response_class=HTMLResponse)
    def help_page(request: Request):
        """Human-readable guide to the website and REST API."""
        return _TEMPLATES.TemplateResponse(request, "docs.html", {})

    @app.get("/tutorial", response_class=HTMLResponse)
    def tutorial_page(request: Request):
        """Worked questions, plus the runnable Python client shown inline."""
        return _TEMPLATES.TemplateResponse(
            request, "tutorial.html", {"script": _EXAMPLE_SCRIPT.read_text()},
        )

    @app.get("/tutorial/airbabel_example.py")
    def tutorial_script():
        """The same file the tutorial page shows, as a download."""
        return FileResponse(
            _EXAMPLE_SCRIPT, media_type="text/x-python",
            filename=_EXAMPLE_SCRIPT.name,
        )

    @app.get("/", response_class=HTMLResponse)
    def home(
        request: Request,
        q: str | None = None,
        mode: str = "auto",
        # A list, so the picker can offer several species. Repeated params
        # (?species=human&species=rhesus) and the older comma form
        # (?species=human,rhesus) both arrive here and both end up comma-joined,
        # which is what queries._species_tags already parses.
        species: list[str] | None = Query(None),
        min: float = 0.85,
        cat: str = "coding",
        kind: str = "nt",
        partial: int = 0,
        functional: str | None = None,
        source: str | None = None,
        names: int = 0,   # PMID mode: also list name-only (unverified) leads
        session: Session = Depends(get_session),
    ):
        results = None
        error = None
        pmid_result = None   # PMID -> the alleles that paper attests
        choices = None       # cross-species / multi-record disambiguation buttons
        choices_msg = None
        # Advanced filters. Functionality is a tri-state select ("" = any, "1" =
        # functional, "0" = non-functional); source is a plain name (or "" = any).
        func_val = {"1": True, "0": False}.get(functional)
        src = source or None
        # Whether the results carry a meaningful sequence-identity score. A name
        # lookup does not (its match is on the label, not the sequence), so the
        # Identity column is shown only for sequence / similarity results.
        show_identity = False
        # Species scope. Selecting nothing already meant "all species", but getting
        # *back* to all meant Ctrl-clicking every chosen species off one at a time.
        # "*" is an explicit "All species" option, and because a plain click in a
        # multi-select replaces the whole selection, it clears the picker in one
        # click. It wins over any tag picked alongside it: asking for everything is
        # not a narrowing.
        sp_sel = _species_sel(species)
        sp = None if "*" in sp_sel else (",".join(sp_sel) or None)
        cat = cat if cat in ("coding", "leader", "rss") else "coding"
        kind = "aa" if kind == "aa" else "nt"
        if q and q.strip():
            query = q.strip()
            try:
                if mode == "pmid":
                    pmid_result = queries.alleles_for_pmid(
                        session, query, species=sp, include_name=bool(names),
                    )
                    if pmid_result is None:
                        error = (
                            "Enter a PubMed ID (digits), PMID:29163486, or a PubMed URL."
                        )
                elif mode == "auto":
                    # A bare 4-9 digit string cannot be an allele name (they carry
                    # letters) or a UID (SPECIES-prefixed), so digits are unambiguous
                    # and worth answering as a PMID. Only when it actually matches:
                    # otherwise fall through, rather than replacing a name search
                    # with an empty citation page.
                    maybe = queries.alleles_for_pmid(session, query, species=sp)
                    if maybe is not None and maybe.count:
                        pmid_result = maybe
                    else:
                        # Exact allele name first — surfaces cross-species overlap
                        # instead of silently picking one species.
                        matches = queries.resolve_exact_name(session, query, species=sp)
                        if len(matches) == 1:
                            return RedirectResponse(url=f"/seq/{matches[0].uid}", status_code=303)
                        if len(matches) > 1:
                            nsp = len({m.species for m in matches})
                            choices = [
                                {"species": m.species, "uid": m.uid, "href": f"/seq/{m.uid}"}
                                for m in matches
                            ]
                            choices_msg = (
                                f"“{query}” matches {len(matches)} records"
                                + (f" across {nsp} species" if nsp > 1 else "")
                                + " — the sequences differ; pick one:"
                            )
                        else:
                            hit = queries.resolve(session, query)
                            if hit is not None:
                                return RedirectResponse(
                                    url=f"/seq/{hit.record.native.uid}", status_code=303
                                )
                            results = queries.search_by_name(
                                session, query, species=sp,
                                functional=func_val, source=src,
                            )
                elif mode == "similar":
                    # Show, inline on the front page, the alleles whose SEQUENCE is
                    # most similar to what was typed, ranked by % identity — just
                    # like the other modes, no redirect. A name is resolved to its
                    # sequence first; a genuine cross-species name needs a pick;
                    # anything else (a pasted sequence) is compared directly.
                    matches = queries.resolve_exact_name(session, query, species=sp)
                    target_uid = None
                    if matches:
                        if len({m.species for m in matches}) > 1:
                            choices = [
                                {"species": m.species, "uid": m.uid,
                                 "href": f"/?q={quote(query)}&mode=similar"
                                         f"&species={m.species.lower()}&min={min}"
                                         f"&kind={kind}"}
                                for m in matches
                            ]
                            choices_msg = (
                                f"“{query}” exists in {len(choices)} species — pick one "
                                "to see its similar alleles, or open its allele page:"
                            )
                        else:
                            target_uid = max(matches, key=lambda m: m.length).uid
                    else:
                        hit = queries.resolve(session, query)
                        if hit is not None:
                            target_uid = hit.record.native.uid
                    if target_uid is not None:
                        sim_fn = (
                            queries.find_similar_aa if kind == "aa"
                            else queries.find_similar
                        )
                        # Respect the species selector: "any species" (no filter)
                        # → search across species; a chosen species → restrict to
                        # *that* species. It has to go in as target_species:
                        # same_species=True means "the query's own species", which
                        # for a rhesus allele with human selected returned rhesus
                        # neighbours — the opposite of what was asked. Several may
                        # be selected (comma-separated); target_species takes a list.
                        results = sim_fn(
                            session, target_uid, min_identity=min,
                            same_species=bool(sp), target_species=sp,
                            functional=func_val, source=src,
                        )
                        show_identity = True
                    elif choices is None:
                        # not a known allele name/UID — compare the text as a sequence
                        try:
                            with maybe_partial_slot(bool(partial)):
                                if kind == "aa":
                                    results = queries.search_by_aa(
                                        session, query, species=sp,
                                        partial=bool(partial),
                                        functional=func_val, source=src,
                                    )
                                else:
                                    results = queries.search_by_sequence(
                                        session, query, species=sp, category=cat,
                                        partial=bool(partial),
                                        functional=func_val, source=src,
                                    )
                            show_identity = True
                        except ValueError:
                            # strictly sequence-similarity: never fall back to a
                            # fuzzy name search.
                            error = (
                                "Enter an allele name, Server UID, or a nucleotide "
                                "sequence to find similar sequences."
                            )
                            results = None
                elif mode == "name":
                    results = queries.search_by_name(
                        session, query, species=sp, partial=bool(partial),
                        functional=func_val, source=src,
                    )
                elif mode == "sequence":
                    with maybe_partial_slot(bool(partial)):
                        results = queries.search_by_sequence(
                            session, query, species=sp, category=cat,
                            partial=bool(partial),
                            functional=func_val, source=src,
                        )
                    show_identity = True
                elif mode == "aa":
                    with maybe_partial_slot(bool(partial)):
                        results = queries.search_by_aa(
                            session, query, species=sp, partial=bool(partial),
                            functional=func_val, source=src,
                        )
                    show_identity = True
            except ValueError as exc:
                error = str(exc)
        # Apply the min-identity threshold to the sequence-similarity modes
        # (find_similar already honours it; this also covers a pasted sequence).
        if results is not None and mode in ("sequence", "aa", "similar") and 0 < min <= 1:
            kept = [h for h in results.hits if h.score >= min]
            results.hits = kept
            results.count = len(kept)
        multi_species = bool(
            results and len({h.species for h in results.hits}) > 1
        )
        # Landing / no-query branch: build the dataset overview summary. The
        # per-species length + concordance charts live on /overview.
        landing = not (q and q.strip())
        stats = (_cached("stats", session, lambda: queries.dataset_stats(session))
                 if landing else None)
        release = queries.dataset_release(session) if landing else None
        return _TEMPLATES.TemplateResponse(
            request, "search.html",
            {"q": q or "", "mode": mode, "species": sp_sel,
             "min": min, "cat": cat,
             "kind": kind, "partial": partial, "show_identity": show_identity,
             "functional": functional or "", "source": source or "",
             "results": results, "error": error, "multi_species": multi_species,
             "pmid_result": pmid_result, "names": names,
             "choices": choices, "choices_msg": choices_msg,
             "back_url": _here(request), "facets": result_facets(results),
             "stats": stats, "release": release,
             "species_options": queries.species_options(session),
             "source_options": queries.list_sources(session)},
        )

    @app.get("/overview", response_class=HTMLResponse)
    def overview(
        request: Request,
        species: str | None = None,
        session: Session = Depends(get_session),
    ):
        sp = species or None
        stats = _cached("stats", session, lambda: queries.dataset_stats(session))
        return _TEMPLATES.TemplateResponse(
            request, "overview.html",
            {"species": species or "",
             "stats": stats,
             "length_dists": {
                 s: d for s in stats["segments"]
                 if (d := queries.length_distribution(session, species=sp, segment=s))},
             "source_conc": queries.source_concordance(session, species=sp)},
        )

    @app.get("/examples", response_class=HTMLResponse)
    def examples_page(request: Request, session: Session = Depends(get_session)):
        """Real records that show why names alone can't identify an allele."""
        return _TEMPLATES.TemplateResponse(
            request, "examples.html", {"examples": queries.example_showcase(session)},
        )

    @app.get("/sources", response_class=HTMLResponse)
    def sources_page(request: Request, session: Session = Depends(get_session)):
        """Every data source, its licence and citation — the paper trail for the dataset."""
        from .. import export

        return _TEMPLATES.TemplateResponse(
            request, "sources.html",
            {"sources": queries.source_catalog(session),
             "release": queries.dataset_release(session),
             "attribution": export.ATTRIBUTION},
        )

    @app.get("/seq/{uid}", response_class=HTMLResponse)
    def detail(request: Request, uid: str, session: Session = Depends(get_session)):
        hit = queries.resolve(session, uid)
        if hit is None or hit.matched_by not in ("uid", "protein_uid"):
            raise HTTPException(status_code=404, detail=f"no sequence with uid {uid!r}")
        if hit.record.native.category == "rss":
            return RedirectResponse(url=f"/rss/{hit.record.native.uid}", status_code=303)
        if hit.record.native.category == "leader":
            return RedirectResponse(url=f"/leader/{hit.record.native.uid}", status_code=303)
        # The record page is a lean identity record. Everything comparative (diff,
        # alignment, cross-species, protein siblings) lives on the similarity query,
        # reached via the "Find similar" entry point — a comparison only exists in a
        # query, never on a lone record.
        return _TEMPLATES.TemplateResponse(
            request, "detail.html", {"record": hit.record},
        )

    @app.get("/rss/{uid}", response_class=HTMLResponse)
    def rss_page(request: Request, uid: str, session: Session = Depends(get_session)):
        found = queries.rss_detail(session, uid)
        if found is None:
            raise HTTPException(status_code=404, detail=f"no RSS with uid {uid!r}")
        rss, alleles, sources = found
        kind = queries.RSS_KINDS.get(rss.subtype or "")
        return _TEMPLATES.TemplateResponse(
            request, "rss.html",
            {
                "rss": rss, "alleles": alleles, "sources": sources,
                # The label and the part order both come from the subtype: a J-RS
                # and a 5'D-RS carry their nonamer first, so slicing the head as a
                # heptamer would mislabel half the signals we now hold.
                "kind_label": kind[0] if kind else (rss.subtype or "RSS"),
                "parts": queries.rss_parts(rss.subtype, rss.sequence),
            },
        )

    @app.get("/leader/{uid}", response_class=HTMLResponse)
    def leader_page(request: Request, uid: str, session: Session = Depends(get_session)):
        found = queries.leader_detail(session, uid)
        if found is None:
            raise HTTPException(status_code=404, detail=f"no leader with uid {uid!r}")
        leader, alleles, sources, exon1_len = found
        return _TEMPLATES.TemplateResponse(
            request, "leader.html",
            {"leader": leader, "alleles": alleles, "sources": sources,
             "exon1_len": exon1_len},
        )

    @app.get("/rss/{uid}/similar", response_class=HTMLResponse)
    def rss_similar_page(
        request: Request, uid: str, min: float = 0.85, cross_species: int = 0,
        session: Session = Depends(get_session),
    ):
        sim = queries.find_similar_rss(
            session, uid, min_identity=min, same_species=not cross_species
        )
        if sim is None:
            raise HTTPException(status_code=404, detail=f"no RSS with uid {uid!r}")
        profile = queries.conservation_profile(
            session, [sim.uid] + [h.uid for h in sim.hits[:60]], kind="nt"
        )
        return _TEMPLATES.TemplateResponse(
            request, "rss_similar.html",
            {"sim": sim, "min": min, "cross_species": bool(cross_species),
             "profile": profile},
        )

    @app.get("/seq/{uid}/similar", response_class=HTMLResponse)
    def similar(
        request: Request,
        uid: str,
        min: float = 0.85,
        cross_species: int = 0,
        to_species: list[str] | None = Query(None),
        kind: str = "nt",
        back: str | None = None,
        session: Session = Depends(get_session),
    ):
        if back is not None and _safe_back(back) is None:
            # An unusable back (the old absolute, ever-nesting form, or anything
            # off-site) goes to the same page without it, before any similarity
            # work: a crawler replaying its backlog of old URLs then costs a
            # redirect and lands on a page it has already seen.
            return RedirectResponse(url=_here(request), status_code=301)
        kind = "aa" if kind == "aa" else "nt"
        fn = queries.find_similar_aa if kind == "aa" else queries.find_similar
        # Species scope, in one selector: named species restrict the comparison to
        # those species' closest orthologs — several may be named, which is the point
        # of the selection ("compare me against human and mouse, nothing else"). "*"
        # (or the legacy cross_species flag) spans all species; otherwise same-species.
        sel = _species_sel(to_species)
        tags = [t for t in sel if t != "*"]
        target = ",".join(tags) or None
        # "*" wins over a tag selection: asking for everything is not a narrowing.
        span_all = "*" in sel or bool(cross_species and not target)
        results = fn(
            session, uid, min_identity=min,
            same_species=not span_all, target_species=target,
        )
        if results is None:
            raise HTTPException(status_code=404, detail=f"no sequence with uid {uid!r}")
        rec = queries.resolve(session, uid)
        label = rec.record.airr.label if rec else None
        category = rec.record.native.category if rec else "coding"
        seed_length = rec.record.native.length if rec else None
        multi = len({h.species for h in results.hits}) > 1
        profile = queries.conservation_profile(
            session, [uid] + [h.uid for h in results.hits[:60]], kind=kind
        )
        # Current selection for the species picker: a list, so several tags stay lit.
        # Empty = same species; ["*"] = all.
        species_sel = ["*"] if span_all else tags
        # Cross-species query: flag hits that are reciprocal best hits (a soft,
        # sequence-based orthology hint, never a curated call). Only computed when
        # the query actually crosses species, so same-species search pays nothing.
        reciprocal_uids: set[str] = set()
        if multi and category == "coding":
            reciprocal_uids = {
                o.uid for o in queries.reciprocal_orthologs(session, uid) if o.reciprocal
            }
        return _TEMPLATES.TemplateResponse(
            request, "similar.html",
            {"uid": uid, "label": label, "kind": kind, "min": min,
             "category": category, "profile": profile,
             "cross_species": bool(span_all), "species_sel": species_sel,
             "species_options": queries.species_options(session), "back": back,
             "results": results, "multi_species": multi,
             "seed_uid": uid, "seed_length": seed_length,
             "reciprocal_uids": reciprocal_uids,
             # Links from here return to the original search, not to this page:
             # see _here for why the chain must not grow.
             "back_url": back or _here(request), "facets": result_facets(results)},
        )

    @app.get("/seq/{uid}/compare/{other}", response_class=HTMLResponse)
    def compare_fragment(
        request: Request, uid: str, other: str,
        session: Session = Depends(get_session),
    ):
        """HTML fragment: the nt diff + IMGT alignment of one similarity hit against
        the query. Fetched on demand when a hit row is expanded, so the results page
        stays lean (nothing computed until asked)."""
        pair = queries.compare_pair(session, uid, other)
        if pair is None:
            raise HTTPException(status_code=404, detail="unknown sequence")
        diff, alignment = pair
        return _TEMPLATES.TemplateResponse(
            request, "_compare.html", {"diff": diff, "alignment": alignment},
        )

    return app


app = create_app()


def main() -> None:
    import uvicorn

    # Bind to localhost only; public exposure is a later, gated decision.
    uvicorn.run(app, host="127.0.0.1", port=8000)


if __name__ == "__main__":
    main()
