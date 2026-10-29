"""Alembic environment. URL + metadata come from allele_server."""

from __future__ import annotations

from logging.config import fileConfig

from sqlalchemy import create_engine

from alembic import context
from allele_server.db.models import Base
from allele_server.settings import get_settings

config = context.config
if config.config_file_name is not None:
    # `disable_existing_loggers=False`: alembic.ini configures alembic's own loggers,
    # it must not silence the application's. The default (True) disables every logger
    # already imported, so anything logged afterwards by `allele_server.*` — e.g. the
    # ingest warning about an unregistered IMGT species — would vanish in any process
    # that ran a migration first (the test suite does exactly this).
    fileConfig(config.config_file_name, disable_existing_loggers=False)

target_metadata = Base.metadata
DB_URL = get_settings().database_url


def run_migrations_offline() -> None:
    context.configure(
        url=DB_URL,
        target_metadata=target_metadata,
        literal_binds=True,
        dialect_opts={"paramstyle": "named"},
        compare_type=True,
    )
    with context.begin_transaction():
        context.run_migrations()


def run_migrations_online() -> None:
    engine = create_engine(DB_URL, future=True)
    with engine.connect() as connection:
        context.configure(
            connection=connection,
            target_metadata=target_metadata,
            compare_type=True,
        )
        with context.begin_transaction():
            context.run_migrations()
    engine.dispose()


if context.is_offline_mode():
    run_migrations_offline()
else:
    run_migrations_online()
