"""Shared fixtures.

Tests run against a **separate** Postgres database (``database_url_test``), created
and migrated on demand, so the suite never touches the served dev data.

Preparing that database can fail in two very different ways, and they are treated
differently (CLAUDE.md: *fail loudly, no silent fallbacks*):

*absent* -- nothing is listening: connection refused, host unreachable, DNS failure.
    CI without Postgres, or a contributor without Docker. The DB-backed tests are
    **skipped**; that is exactly what the skip is for.

*broken* -- Postgres answers, but the test database cannot be brought to head: it is
    stamped at an Alembic revision this checkout does not contain (the branch-switch
    trap), a migration errors, an extension is missing, credentials are refused.
    That is a broken environment, not a missing one, so the run is **aborted** with the
    reason and the remedy. Skipping here used to print ``40 passed, 61 skipped`` -- a
    green-looking run that tested almost nothing, with the cause visible only under
    ``-rs``.
"""

from __future__ import annotations

import os
import socket
from typing import NoReturn

import pytest
from sqlalchemy import create_engine, inspect, text
from sqlalchemy.engine import make_url
from sqlalchemy.exc import DBAPIError

from allele_server.db.session import configure_engine, session_factory
from allele_server.settings import get_settings

_TABLES = (
    # Truncated here but deliberately NOT in ingest/reload.py::_TABLES. Test isolation
    # needs a clean table; a production reload must preserve approvals and rejections.
    # source_candidates carries no FK to sources, so the CASCADE below cannot reach it
    # and it has to be named explicitly. The divergence is intended — see
    # tests/test_reload.py::test_reload_preserves_source_candidates.
    "source_candidates",
    "sequence_annotations",
    "sequence_source",
    "sequence_links",
    "names",
    "sequences",
    "sources",
)


def _port_is_open(host: str, port: int, timeout: float = 3.0) -> bool:
    """Cheap "is anything listening there" probe, independent of the driver."""
    try:
        with socket.create_connection((host, port), timeout=timeout):
            return True
    except OSError:
        return False


def _is_server_absent(test_url: str, exc: BaseException) -> bool:
    """True only when there is genuinely no database to test against.

    Two signals, neither of which is message text:

    * a SQLSTATE means the *server* reported the error, so Postgres is up (``28xxx``
      bad credentials, ``42501`` insufficient privilege, ``3D000`` no such database);
    * otherwise psycopg failed before any conversation -- but that covers both
      "connection refused" and "the server said FATAL", because psycopg3 attaches no
      SQLSTATE to connect-time failures. So probe the socket directly: a port that
      accepts a TCP connection means something *is* there and the failure is a broken
      configuration, not an absent database.
    """
    orig = getattr(exc, "orig", exc)
    if getattr(orig, "sqlstate", None):
        return False
    url = make_url(test_url)
    if not url.host:
        return True  # unix socket: nothing cheap to probe, treat as absent
    return not _port_is_open(url.host, url.port or 5432)


def _abort(
    config: pytest.Config, url, headline: str, detail: str = "", remedy: list[str] | None = None
) -> NoReturn:
    """Stop the whole run with a banner naming the problem and the remedy.

    Aborting rather than failing one fixture is deliberate: a test database that
    cannot reach head invalidates every DB-backed test at once, and a partial run is
    the very thing that used to read as green.
    """
    u = make_url(url)
    banner = [
        "",
        "=" * 79,
        "TEST DATABASE UNUSABLE -- aborting (these tests are NOT being skipped)",
        "=" * 79,
        f"  database : {u.render_as_string(hide_password=True)}",
        f"  problem  : {headline}",
    ]
    banner += [f"             {line}" for line in detail.splitlines() if line.strip()]
    banner += [""]
    banner += remedy or [
        "  Postgres is reachable, so the environment is broken rather than missing.",
        "  The test database is disposable -- recreate it and re-run:",
        "",
        f"      psql -h {u.host} -p {u.port} -U {u.username} -d postgres \\",
        f"          -c 'DROP DATABASE IF EXISTS \"{u.database}\"' \\",
        f"          -c 'CREATE DATABASE \"{u.database}\"'",
        "",
        "  pytest recreates the schema (alembic upgrade head) on the next run.",
    ]
    banner += ["=" * 79, ""]
    reporter = config.pluginmanager.getplugin("terminalreporter")
    if reporter is not None:
        reporter.write_line("\n".join(banner), red=True, bold=True)
    pytest.exit(f"test database unusable: {headline}", returncode=1)


def _ensure_database(test_url: str) -> None:
    """Create the test database if it doesn't exist (connect via 'postgres')."""
    url = make_url(test_url)
    admin = url.set(database="postgres")
    eng = create_engine(admin, isolation_level="AUTOCOMMIT", future=True)
    try:
        with eng.connect() as conn:
            exists = conn.execute(
                text("SELECT 1 FROM pg_database WHERE datname = :n"),
                {"n": url.database},
            ).scalar()
            if not exists:
                conn.execute(text(f'CREATE DATABASE "{url.database}"'))
    finally:
        eng.dispose()


def _stamped_revision(test_url: str) -> str | None:
    """The revision the test database currently claims, or None if never migrated."""
    eng = create_engine(test_url, future=True)
    try:
        if not inspect(eng).has_table("alembic_version"):
            return None
        with eng.connect() as conn:
            return conn.execute(text("SELECT version_num FROM alembic_version")).scalar()
    finally:
        eng.dispose()


def _revision_exists(script, revision: str) -> bool:
    """Whether this checkout's migration tree knows ``revision``."""
    try:
        return script.get_revision(revision) is not None
    except Exception:  # noqa: BLE001 - alembic raises CommandError/ResolutionError
        return False


@pytest.fixture(scope="session")
def engine(pytestconfig: pytest.Config):
    test_url = get_settings().database_url_test
    safe_url = make_url(test_url).render_as_string(hide_password=True)

    # --- absent vs broken -------------------------------------------------------
    # The first connection is the only place a skip can be earned.
    try:
        _ensure_database(test_url)
    except DBAPIError as exc:
        if _is_server_absent(test_url, exc):
            pytest.skip(
                f"No Postgres at {safe_url}: {exc.orig or exc}. "
                "Start it with: docker compose -f docker-compose.dev.yml up -d"
            )
        _abort(
            pytestconfig,
            test_url,
            "a server is listening, but the connection was rejected",
            str(exc.orig or exc),
            remedy=[
                "  Something is answering on that host/port, so this is a configured",
                "  server refusing us -- wrong credentials, pg_hba, or another service on",
                "  the port. Fix DATABASE_URL_TEST (or the server), or point it at a host",
                "  with no Postgres at all if you meant to run without a database.",
            ],
        )
    except Exception as exc:  # noqa: BLE001 - a broken URL/driver is not an absent DB
        _abort(
            pytestconfig,
            test_url,
            "could not create the test database",
            f"{type(exc).__name__}: {exc}",
        )

    # Postgres is up from here on: every remaining failure is loud.
    os.environ["DATABASE_URL"] = test_url  # alembic/env.py reads it through settings
    from alembic.config import Config
    from alembic.script import ScriptDirectory

    from alembic import command

    cfg = Config(str(pytestconfig.rootpath / "alembic.ini"))
    cfg.set_main_option("script_location", str(pytestconfig.rootpath / "alembic"))
    script = ScriptDirectory.from_config(cfg)
    head = script.get_current_head()

    # The branch-switch trap: the test DB was migrated by a branch whose migrations
    # this checkout does not have, so `upgrade head` can't even find where it stands.
    stamped = _stamped_revision(test_url)
    if stamped is not None and not _revision_exists(script, stamped):
        _abort(
            pytestconfig,
            test_url,
            f"stamped at Alembic revision {stamped!r}, which does not exist in this checkout",
            f"This checkout's head is {head!r}. The database was migrated by another "
            f"branch (a newer/other migration head) and cannot be downgraded from here.",
        )

    try:
        command.upgrade(cfg, "head")
    except Exception as exc:  # noqa: BLE001
        _abort(
            pytestconfig,
            test_url,
            "'alembic upgrade head' failed",
            f"{type(exc).__name__}: {exc}",
        )

    now_at = _stamped_revision(test_url)
    if now_at != head:
        _abort(
            pytestconfig,
            test_url,
            f"migrations ran but the database is at {now_at!r}, not head {head!r}",
        )

    try:
        eng = configure_engine(test_url)
        with eng.connect() as conn:
            conn.execute(text("SELECT 1"))
    except Exception as exc:  # noqa: BLE001
        _abort(
            pytestconfig,
            test_url,
            "the migrated test database is not usable",
            f"{type(exc).__name__}: {exc}",
        )
    return eng


@pytest.fixture
def session(engine):
    """A clean session: all tables truncated before the test (on the TEST db)."""
    factory = session_factory()
    sess = factory()
    sess.execute(text(f"TRUNCATE {', '.join(_TABLES)} RESTART IDENTITY CASCADE"))
    sess.commit()
    try:
        yield sess
    finally:
        sess.rollback()
        sess.close()


@pytest.fixture
def manifest_dir(tmp_path, monkeypatch):
    """A temporary overlay directory of source manifests, plus a registry reset.

    Lets a test exercise the ingest gate against a fabricated source without shipping a
    real ``.toml``. The registry is memoised and process-wide, so it is refreshed on the
    way in *and* on the way out — otherwise a fabricated source would leak into every
    later test in the session.
    """
    import allele_server.sources as sources_pkg

    monkeypatch.setenv("SOURCE_MANIFEST_DIR", str(tmp_path))
    sources_pkg.registry(refresh=True)
    try:
        yield tmp_path
    finally:
        monkeypatch.undo()
        sources_pkg.registry(refresh=True)
