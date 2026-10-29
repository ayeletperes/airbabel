"""sources: manifest link, citation/licence/approval columns; kind gains 'derived'

The source manifest (``allele_server.sources``, one reviewed TOML per source) becomes
the single source of truth for a source's licence, citation and redistribution terms.
Those values are *projected* onto the ``sources`` row by
``ingest/loader.get_or_create_source``, so a row is self-describing: the export
attribution line and the citation-completeness test read the database, not the files,
and nothing has to remember to run a sync step. See ``docs/source-policy.md``.

Three changes beyond the new columns:

* ``license`` widens 64 -> 128. Real licence strings do not fit in 64 characters.
* ``kind`` gains ``'derived'``. ``docs/ingestion.md`` already classifies the mined
  literature layers that way while the CHECK constraint forbade it, so the constraint
  was describing a schema we do not have.
* ``citation_status`` is constrained to the two values the policy allows. There is
  deliberately no third "unknown" value: that is how nine ``TODO(citation)`` entries
  survived in the code this replaces.

Additive; no data migration. The values are projected on the next ingest, and
``tests/test_sources_completeness.py`` fails if any row is still missing them.

Revision ID: 0004
Revises: 0003
Create Date: 2026-07-29
"""
from __future__ import annotations

from collections.abc import Sequence

from alembic import op

revision: str = "0004"
down_revision: str | None = "0003"
branch_labels: str | Sequence[str] | None = None
depends_on: str | Sequence[str] | None = None


_UP: tuple[str, ...] = (
    "ALTER TABLE sources ALTER COLUMN license TYPE VARCHAR(128)",
    "ALTER TABLE sources ADD COLUMN manifest_id VARCHAR(64)",
    "ALTER TABLE sources ADD COLUMN citation TEXT",
    "ALTER TABLE sources ADD COLUMN citation_status VARCHAR(16)",
    "ALTER TABLE sources ADD COLUMN doi VARCHAR(128)",
    "ALTER TABLE sources ADD COLUMN url VARCHAR(512)",
    "ALTER TABLE sources ADD COLUMN licence_url VARCHAR(512)",
    "ALTER TABLE sources ADD COLUMN approved_by VARCHAR(128)",
    "ALTER TABLE sources ADD COLUMN approved_at TIMESTAMP WITH TIME ZONE",
    "ALTER TABLE sources DROP CONSTRAINT ck_sources_kind",
    "ALTER TABLE sources ADD CONSTRAINT ck_sources_kind "
    "CHECK (kind IN ('internal', 'external', 'derived'))",
    "ALTER TABLE sources ADD CONSTRAINT ck_sources_citation_status "
    "CHECK (citation_status IS NULL OR citation_status IN ('verified', 'none_available'))",
    "CREATE INDEX ix_sources_manifest_id ON sources (manifest_id)",
)

_DOWN: tuple[str, ...] = (
    "DROP INDEX IF EXISTS ix_sources_manifest_id",
    "ALTER TABLE sources DROP CONSTRAINT IF EXISTS ck_sources_citation_status",
    "ALTER TABLE sources DROP CONSTRAINT ck_sources_kind",
    # Any 'derived' row would violate the narrowed constraint; fold those back first
    # rather than letting the ALTER fail on data this migration itself allowed in.
    "UPDATE sources SET kind = 'external' WHERE kind = 'derived'",
    "ALTER TABLE sources ADD CONSTRAINT ck_sources_kind "
    "CHECK (kind IN ('internal', 'external'))",
    "ALTER TABLE sources DROP COLUMN approved_at",
    "ALTER TABLE sources DROP COLUMN approved_by",
    "ALTER TABLE sources DROP COLUMN licence_url",
    "ALTER TABLE sources DROP COLUMN url",
    "ALTER TABLE sources DROP COLUMN doi",
    "ALTER TABLE sources DROP COLUMN citation_status",
    "ALTER TABLE sources DROP COLUMN citation",
    "ALTER TABLE sources DROP COLUMN manifest_id",
    # Last, and deliberately able to fail: a licence string longer than 64 characters
    # cannot be re-narrowed without losing it, and silent truncation is worse than a
    # failed downgrade.
    "ALTER TABLE sources ALTER COLUMN license TYPE VARCHAR(64)",
)


def upgrade() -> None:
    for stmt in _UP:
        op.execute(stmt)


def downgrade() -> None:
    for stmt in _DOWN:
        op.execute(stmt)
