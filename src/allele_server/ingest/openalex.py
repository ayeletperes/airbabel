"""OpenAlex OA resolution — a second opinion on "is this article mineable?".

Europe PMC's ``OPEN_ACCESS:Y`` flag is the only gate the corpus sweep used to
apply, and it is demonstrably wrong on real records. Measured, both verified by
hand:

* PMID 42251052 (*Nat Commun* 2026, ``10.1038/s41467-026-73730-9``) — Europe PMC
  reports ``isOpenAccess:N``, ``inPMC:N``, no PMCID. The article is open access
  (CC BY-NC-ND at nature.com) **and** has a PMC copy, ``PMC13396459``, which
  serves HTTP 200. Europe PMC simply has not indexed it.
* PMID 36796364 (Corcoran *et al.*, *Immunity* 2023) — genuinely absent from PMC,
  but OpenAlex reports ``is_oa=True, oa_status=green`` with a free repository
  copy.

So the sweep asks OpenAlex about the articles Europe PMC calls closed, and keeps
the ones that have a legal OA copy. OpenAlex also carries the *locations* of that
copy, which is where a missing PMCID comes back from: the PMC landing page URL
the API returns contains it. **No URL here is constructed by us** — a PMCID is
only ever read out of a URL OpenAlex actually returned.

**What this actually recovers, measured.** Random 200 of the 1,935 IGHV articles
Europe PMC flags closed (2026-09-03): 146 are open access per OpenAlex, 18 of
those have a PMC copy at all, and **0** of those 18 are served by Europe PMC's OA
web services — the copies exist but sit outside the redistributable OA subset
(``fullTextXML`` 404s, ``supplementaryFiles`` answers HTTP 200 with a 165-byte
``errorBean``). The flag is wrong; almost none of what it wrongly excluded is
reachable through the fetch path this layer is licensed to use. So the pass is
opt-in, and every recovered-but-unfetchable article is *counted*
(``recovered_without_pmc`` here, ``no_pmcid`` / ``no_content`` in the miner)
rather than left looking like an article that was read and matched nothing.
See ``docs/literature-mining-plan.md``, addendum 2026-09-03.

Same rails as the other enrichment clients: one small JSON per PMID on disk (a
cache hit never touches the network, so a re-run re-queries nothing), polite
pacing, and best-effort networking — a failed batch is not cached, so it is
retried, while "OpenAlex has no record of this PMID" is cached as an empty
object because that answer is stable.

No API key is required; OpenAlex asks only for a mailto in the User-Agent.
"""

from __future__ import annotations

import json
import re
import time
import urllib.parse
import urllib.request
from collections.abc import Callable, Iterable
from dataclasses import dataclass
from pathlib import Path

_OPENALEX = "https://api.openalex.org/works"
_MAILTO = "aperes7@gmail.com"
_USER_AGENT = f"allele-server/0.1 (mailto:{_MAILTO})"

# Only the fields we act on. Asking for the whole record costs OpenAlex (and us)
# a few hundred KB per batch for data nothing reads.
_SELECT = "ids,open_access,locations"

# OpenAlex ORs up to 50 filter values per request, which is also its page size.
_BATCH = 50

# A PMCID as it appears inside a URL OpenAlex returned, e.g.
# ``https://pmc.ncbi.nlm.nih.gov/articles/PMC13396459/``.
_PMCID_RE = re.compile(r"\bPMC\d{4,}\b")


@dataclass(frozen=True)
class OAWork:
    """What OpenAlex says about one article's open-access status."""

    pmid: str
    is_oa: bool
    oa_status: str | None
    pmcid: str | None
    oa_url: str | None
    oa_source: str | None


class OpenAlexClient:
    """Look up articles by PMID, one cached JSON per PMID, fetched in batches.

    ``works(pmids)`` returns ``{pmid: OAWork}`` for every PMID OpenAlex knows;
    PMIDs it has no record of are absent (and that absence is cached).
    ``fetcher`` is an injectable ``url -> bytes`` seam for tests, matching
    :class:`~allele_server.ingest.literature_supp.PMCSuppClient`.
    """

    def __init__(
        self,
        cache_dir: str | Path,
        *,
        fetcher: Callable[[str], bytes] | None = None,
        min_interval: float = 0.11,  # OpenAlex's polite pool allows 10 req/s
        timeout: int = 30,
    ) -> None:
        self.cache = Path(cache_dir)
        self.cache.mkdir(parents=True, exist_ok=True)
        self.fetcher = fetcher
        self.min_interval = min_interval
        self.timeout = timeout
        self._last = 0.0

    # -- cache ---------------------------------------------------------------

    def _path(self, pmid: str) -> Path:
        return self.cache / f"{pmid}.json"

    def _cached(self, pmid: str) -> dict | None:
        """The stored record for ``pmid``: a dict, ``{}`` for a known miss, or None."""
        cp = self._path(pmid)
        if not cp.exists():
            return None
        try:
            return json.loads(cp.read_text())
        except json.JSONDecodeError:
            return None

    # -- fetch ---------------------------------------------------------------

    def _wait_turn(self) -> None:
        gap = self.min_interval - (time.monotonic() - self._last)
        if gap > 0:
            time.sleep(gap)
        self._last = time.monotonic()

    def _get(self, url: str) -> bytes | None:
        if self.fetcher is not None:
            try:
                return self.fetcher(url)
            except Exception:  # noqa: BLE001 - injected seam, stay best-effort
                return None
        self._wait_turn()
        req = urllib.request.Request(url, headers={"User-Agent": _USER_AGENT})
        try:
            with urllib.request.urlopen(req, timeout=self.timeout) as resp:
                return resp.read()
        except Exception:  # noqa: BLE001 - network failure -> no verdict, retry later
            return None

    def _fetch_batch(self, pmids: list[str]) -> None:
        """Fetch up to ``_BATCH`` PMIDs and write one cache entry each.

        A PMID the response does not mention gets ``{}`` — OpenAlex has no record
        of it, which is a stable answer worth caching. A failed request writes
        nothing, so the next run retries it."""
        params = urllib.parse.urlencode({
            "filter": "pmid:" + "|".join(pmids),
            "select": _SELECT,
            "per-page": str(_BATCH),
            "mailto": _MAILTO,
        })
        body = self._get(f"{_OPENALEX}?{params}")
        if body is None:
            return
        try:
            results = json.loads(body.decode()).get("results") or []
        except (json.JSONDecodeError, UnicodeDecodeError, AttributeError):
            return
        got: set[str] = set()
        for work in results:
            pmid = _pmid_of(work)
            if pmid:
                got.add(pmid)
                self._path(pmid).write_text(json.dumps(work))
        for pmid in pmids:
            if pmid not in got:
                self._path(pmid).write_text("{}")

    # -- public --------------------------------------------------------------

    def works(self, pmids: Iterable[str]) -> dict[str, OAWork]:
        """``{pmid: OAWork}`` for the PMIDs OpenAlex knows. Cached; idempotent."""
        wanted = list(dict.fromkeys(p for p in pmids if p))
        missing = [p for p in wanted if self._cached(p) is None]
        for start in range(0, len(missing), _BATCH):
            self._fetch_batch(missing[start:start + _BATCH])
        out: dict[str, OAWork] = {}
        for pmid in wanted:
            raw = self._cached(pmid)
            if raw:
                out[pmid] = _to_work(pmid, raw)
        return out

    def pmcid(self, pmid: str) -> str | None:
        """The PMCID of this article's PMC copy, read out of an OpenAlex URL.

        This is the recovery path for a record Europe PMC has not indexed: the
        PMC copy exists and is reachable, Europe PMC's search just does not know
        its PMCID yet."""
        work = self.works([pmid]).get(pmid)
        return work.pmcid if work else None


def _pmid_of(work: dict) -> str | None:
    """The bare PMID from a work's ``ids.pmid`` URL, or None."""
    url = (work.get("ids") or {}).get("pmid") or ""
    tail = str(url).rstrip("/").rsplit("/", 1)[-1]
    return tail if tail.isdigit() else None


def _to_work(pmid: str, raw: dict) -> OAWork:
    oa = raw.get("open_access") or {}
    pmcid: str | None = None
    oa_url: str | None = None
    oa_source: str | None = None
    for loc in raw.get("locations") or []:
        if not loc.get("is_oa"):
            continue
        urls = [u for u in (loc.get("landing_page_url"), loc.get("pdf_url")) if u]
        if pmcid is None:
            for url in urls:
                hit = _PMCID_RE.search(url)
                if hit:
                    pmcid = hit.group(0)
                    break
        if oa_url is None and urls:
            oa_url = urls[0]
            oa_source = (loc.get("source") or {}).get("display_name")
    return OAWork(
        pmid=pmid,
        is_oa=bool(oa.get("is_oa")),
        oa_status=oa.get("oa_status"),
        pmcid=pmcid,
        oa_url=oa.get("oa_url") or oa_url,
        oa_source=oa_source,
    )
