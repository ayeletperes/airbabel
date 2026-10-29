"""Runtime configuration, loaded from the environment / ``.env``."""

from __future__ import annotations

import os
from pathlib import Path

from pydantic_settings import BaseSettings, SettingsConfigDict


def _cache_home(name: str) -> str:
    """A cache directory outside the working tree, per the XDG base-dir spec.

    The literature caches used to default to ``./<name>_cache``, i.e. inside the
    checkout -- which here lives in a synced Dropbox folder, so every downloaded
    supplementary bundle was queued for upload to the cloud and mirrored onto
    every other machine. That is wrong for a cache by definition (regenerable,
    write-heavy, and in a full-corpus sweep tens of gigabytes), and it competes
    for the same disk as the data. Still overridable through ``.env`` like every
    other path; only the default moved.
    """
    root = os.environ.get("XDG_CACHE_HOME") or str(Path.home() / ".cache")
    return str(Path(root) / "allele_server" / name)


class Settings(BaseSettings):
    model_config = SettingsConfigDict(env_file=".env", extra="ignore")

    # Sources hidden from the served (interactive/API) site. Distinct from a source's
    # `redistributable` flag: KI-TCR is non-redistributable but *published*, so it
    # stays visible; HUSA is *unpublished* internal curation, so it is shelved — its
    # names, its source chip and any sequence it alone attests never reach the site.
    # Removing "HUSA" here un-shelves it (queries gate on this set). Names are matched
    # case-insensitively against `sources.name`.
    shelved_sources: frozenset[str] = frozenset({"HUSA"})

    # SQLAlchemy / psycopg3 connection string.
    database_url: str = "postgresql+psycopg://allele:allele@localhost:5433/allele_server"

    # Separate database for the test suite, so `pytest` never touches served data.
    database_url_test: str = (
        "postgresql+psycopg://allele:allele@localhost:5433/allele_server_test"
    )

    # Where ingestion caches OGRDB downloads.
    ogrdb_cache_dir: str = "./ogrdb_cache"

    # Curation tables + KI workbook for the consolidated reload. Unset by default:
    # these are per-maintainer absolute paths to files that are not in this repo, so a
    # baked-in default is both wrong for every other checkout and a disclosure of the
    # maintainer's filesystem layout in a published tree. Set them in `.env`:
    #
    #   HUSA_TABLE=/path/to/husa_pipeline/data/husa/husa.tsv
    #   MUSA_TABLE=/path/to/MUSA/MUSA_with_rss_and_leader_and_asc_<date>.csv
    #   KI_XLSX=/path/to/kirsanov_2026_immunity_mmc3.xlsx
    #
    # Unset is the opt-out: `reload_all` skips that source, so an unconfigured clone
    # still reloads everything else. A path that IS set but missing is a typo, and it
    # refuses before the TRUNCATE rather than half-loading (see `reload_all`).
    # For HUSA, point at the husa_pipeline output (post
    # 2587-gate-fix + post 48-donor removal), NOT the results_v3/asc_final/current/
    # symlink — that one is stale (still 2026-09-10) and was never repointed after the
    # rebuild moved to husa_pipeline. (master, 2026-09-18)
    husa_table: str = ""
    musa_table: str = ""
    ki_xlsx: str = ""

    # SQLAlchemy connection-pool size. max_overflow is set to 2x this in
    # db.session; raise for higher request concurrency against the cache.
    db_pool_size: int = 20

    # Hard cap on rows any SEARCH endpoint returns. Bulk export is a separate,
    # uncapped endpoint (allele_server.export / /api/export) — rate-limit it upstream.
    search_result_cap: int = 100

    # De novo *partial* (containment) search is uncacheable and CPU-heavy — a full
    # trigram scan + local-identity rescore per request. On a small instance a burst
    # of them can saturate CPU, so the API bounds how many run at once; excess
    # requests wait up to partial_search_queue_timeout seconds for a slot, then 503.
    partial_search_concurrency: int = 2
    partial_search_queue_timeout: float = 2.0

    # Host header allowlist for the /mcp endpoint (MCP's DNS-rebinding guard).
    # The guard stays ON; this names the hostnames it should accept, so the guard
    # survives being fronted by a reverse proxy instead of being switched off.
    # "host:*" matches that host on any port. Set MCP_ALLOWED_HOSTS (JSON list) on a
    # box reachable under another name -- e.g. the lab box behind its ngrok domain.
    mcp_allowed_hosts: list[str] = [
        "127.0.0.1:*", "localhost:*", "airbabel.org", "www.airbabel.org",
    ]

    # --- /mcp access control -------------------------------------------------
    # Bearer token for the MCP endpoint. Empty means open, and that is the public
    # setting: the data is public and read-only, and a token would keep AIRBabel
    # out of the agents it should be reachable from. Set MCP_TOKEN only on a
    # private deployment; CPU is guarded by the limits below.
    mcp_token: str = ""
    # Requests per caller per mcp_rate_window seconds; 0 disables. This guards CPU,
    # not secrets: airbabel.org is a two-core burstable instance and an agent
    # looping over queries can spend the credits the website needs.
    mcp_rate_limit: int = 120
    # The similarity tools (mcp_gate.HEAVY_TOOLS) can run a live scan, so they get
    # a tighter budget: per caller, and in total across every caller. On
    # airbabel.org the balancer hides the client (see export_rate_limit), so
    # "caller" is one of its nodes and the global figure is the one that binds.
    mcp_similar_rate_limit: int = 20
    mcp_similar_global_limit: int = 60
    mcp_rate_window: float = 60.0

    # Bulk export, per caller. Much tighter than /mcp because each request is a
    # different animal: a full-table scan taking 6-8 seconds and returning 5-12 MB,
    # against a two-core burstable instance. Egress is nearly free (AWS gives 100
    # GB/month), so what this protects is CPU credits -- spend those and the whole
    # site throttles, which is the same failure the Caddyfile's crawler rule avoids.
    # This budget is GLOBAL, not per caller: the load balancer does not pass the
    # client address, so a per-caller key silently splits into one bucket per balancer
    # node. Twenty full scans a minute across everyone is already a lot of CPU, and
    # far more than real use -- nobody downloads a reference set twenty times a
    # minute. See the note on export_rate_limit in api/routes.py.
    export_rate_limit: int = 20
    export_rate_window: float = 60.0

    # "Report a bug" link shown in the UI. Empty by default: the site is distributed
    # without pointing at the (private) GitHub issue tracker, which outside readers
    # cannot file to anyway. Set FEEDBACK_URL to show the link again.
    feedback_url: str = ""

    # The site's public origin, for absolute URLs in the sitemap, canonical links
    # and share cards. Not taken from the request: behind Caddy + the load balancer
    # the app sees http://127.0.0.1. Set PUBLIC_URL on another deployment.
    public_url: str = "https://airbabel.org"

    # IndexNow (Bing, Yandex, Seznam, ...): a key the site serves at /<key>.txt to
    # prove it owns the URLs it submits. Public by design, so not a secret; empty
    # turns the route and `python -m allele_server.indexnow` off.
    indexnow_key: str = ""

    # Enrichment passes (opt-in; network + disk cache).
    ncbi_api_key: str | None = None
    ncbi_email: str = "aperes7@gmail.com"
    ncbi_cache_dir: str = "./ncbi_cache"
    imgt_crawl_cache_dir: str = "./imgt_crawl_cache"
    # Historical IMGT/GENE-DB release snapshots (former/retired official names),
    # downloaded via --fetch-history and re-ingested on every --reload --enrich.
    imgt_history_dir: str = "./imgt_history"
    # Cached bulk IMGT/GENE-DB ReferenceSequences FASTA for the all-species ingest.
    imgt_genedb_fasta: str = "./imgt_genedb_ref.fasta"
    # The same records as protein: supplies the *published* V-REGION AA, so aa_origin
    # can say 'sourced' instead of leaving every protein looking like our translation.
    imgt_genedb_fasta_aa: str = "./imgt_genedb_ref_aa.fasta"
    # The same records IMGT-gapped (F+ORF+in-frame P, the only gapped nt file IMGT
    # publishes). Not ingested: the query server reads it as the V numbering reference
    # for diff & align, so it ships with every deployment.
    imgt_genedb_fasta_gapped: str = "./imgt_genedb_ref_gapped.fasta"
    # GENElect leader/RSS responses, memoised so a reload never depends on imgt.org
    imgt_genelect_cache_dir: str = "./imgt_genelect_cache"
    # Literature layer (opt-in; Europe PMC search by allele name -> unverified
    # PubMed leads). Disk-cached; re-runs are idempotent.
    literature_cache_dir: str = _cache_home("literature")
    # Supplementary-file mining (opt-in; PMC-OA supplements -> unverified name +
    # sequence matches). Disk-cached (raw files + PMID->PMCID lookups); idempotent.
    literature_supp_cache_dir: str = _cache_home("literature_supp")

    # Source onboarding (docs/source-policy.md).
    # Extra directory of source manifests, overlaid on the ones shipped in
    # allele_server/sources/. For a deployment that holds a source the repo does not.
    source_manifest_dir: str | None = None
    # Where review dossiers are written. Gitignored and regenerable, like every other
    # cache: the durable audit record is the manifest's [approval] block.
    source_dossier_dir: str = "./source_dossiers"
    # The ingest gate: refuse ("strict"), complain ("warn"), or ignore ("off") a source
    # with no approved manifest. An operator escape hatch for a broken deploy, NOT a
    # development convenience — tests/test_sources_completeness.py always runs strict,
    # so "warn" cannot quietly become the permanent state.
    source_gate: str = "strict"


def get_settings() -> Settings:
    """Return a freshly-loaded Settings instance."""
    return Settings()
