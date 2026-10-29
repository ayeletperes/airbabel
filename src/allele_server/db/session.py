"""Engine and session factory, configured from :mod:`allele_server.settings`."""

from __future__ import annotations

from collections.abc import Iterator
from contextlib import contextmanager

from sqlalchemy import create_engine
from sqlalchemy.engine import Engine
from sqlalchemy.orm import Session, sessionmaker

from ..settings import get_settings

_engine: Engine | None = None
_Session: sessionmaker[Session] | None = None


def get_engine() -> Engine:
    """Return the process-wide SQLAlchemy engine (lazily created)."""
    global _engine, _Session
    if _engine is None:
        st = get_settings()
        _engine = create_engine(
            st.database_url, future=True,
            pool_size=st.db_pool_size, max_overflow=2 * st.db_pool_size,
            pool_pre_ping=True,
        )
        _Session = sessionmaker(bind=_engine, future=True, expire_on_commit=False)
    return _engine


def configure_engine(url: str) -> Engine:
    """Rebind the process-wide engine + session factory to ``url``.

    Used by the test suite to point everything (including the app's request
    sessions) at the throwaway test database instead of the served dev DB.
    """
    global _engine, _Session
    if _engine is not None:
        _engine.dispose()
    st = get_settings()
    _engine = create_engine(
        url, future=True,
        pool_size=st.db_pool_size, max_overflow=2 * st.db_pool_size,
        pool_pre_ping=True,
    )
    _Session = sessionmaker(bind=_engine, future=True, expire_on_commit=False)
    return _engine


def session_factory() -> sessionmaker[Session]:
    """Return the configured session factory."""
    if _Session is None:
        get_engine()
    assert _Session is not None
    return _Session


@contextmanager
def session_scope() -> Iterator[Session]:
    """Transactional session context: commit on success, rollback on error."""
    sess = session_factory()()
    try:
        yield sess
        sess.commit()
    except Exception:
        sess.rollback()
        raise
    finally:
        sess.close()
