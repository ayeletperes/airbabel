"""sequence_neighbors: integer keys, natural PK, denormalised to_uid

Re-key the similarity cache from text ``uid`` columns to integer ``sequences.id``
columns and use the natural key ``(from_id, kind, to_id)`` as the primary key (no
surrogate ``id``). Text-uid keys made the two big btree indexes ~1.76 GB and the
surrogate-id pkey another 307 MB; integer keys + dropping the surrogate cut the whole
cache from ~3.1 GB to ~1.0 GB. The ``to_*`` attributes stay denormalised (``to_uid``
added) so the cached read still assembles a hit with no per-row alias join.

The table is **rebuilt from scratch** after this migration (the neighbour build pass),
so this revision simply drops and recreates it — no row-level data migration. Safe
because the cache is derived, not source-of-truth.

Revision ID: 0002
Revises: 0001
Create Date: 2026-07-17
"""
from __future__ import annotations

from collections.abc import Sequence

from alembic import op

revision: str = "0002"
down_revision: str | None = "0001"
branch_labels: str | Sequence[str] | None = None
depends_on: str | Sequence[str] | None = None


_NEW = (
    """
CREATE TABLE sequence_neighbors (
	from_id INTEGER NOT NULL,
	kind VARCHAR(8) NOT NULL,
	to_id INTEGER NOT NULL,
	identity REAL NOT NULL,
	to_uid VARCHAR(32),
	label VARCHAR(256),
	to_species VARCHAR(16),
	to_locus VARCHAR(8),
	to_segment VARCHAR(8),
	to_category VARCHAR(16),
	PRIMARY KEY (from_id, kind, to_id),
	FOREIGN KEY(from_id) REFERENCES sequences (id) ON DELETE CASCADE,
	FOREIGN KEY(to_id) REFERENCES sequences (id) ON DELETE CASCADE
)
    """,
    "CREATE INDEX ix_neighbor_topk ON sequence_neighbors (from_id, kind, identity)",
    "CREATE INDEX ix_neighbor_to_id ON sequence_neighbors (to_id)",
)

# The former (0001) table, for downgrade().
_OLD = (
    """
CREATE TABLE sequence_neighbors (
	id SERIAL NOT NULL,
	from_uid VARCHAR(32) NOT NULL,
	to_uid VARCHAR(32) NOT NULL,
	kind VARCHAR(8) NOT NULL,
	identity FLOAT NOT NULL,
	label VARCHAR(256),
	to_species VARCHAR(16),
	to_locus VARCHAR(8),
	to_segment VARCHAR(8),
	to_category VARCHAR(16),
	PRIMARY KEY (id),
	CONSTRAINT uq_neighbor_edge UNIQUE (from_uid, kind, to_uid),
	FOREIGN KEY(from_uid) REFERENCES sequences (uid) ON DELETE CASCADE,
	FOREIGN KEY(to_uid) REFERENCES sequences (uid) ON DELETE CASCADE
)
    """,
    "CREATE INDEX ix_neighbor_topk ON sequence_neighbors (from_uid, kind, identity)",
)


def upgrade() -> None:
    op.execute("DROP TABLE IF EXISTS sequence_neighbors CASCADE")
    for stmt in _NEW:
        op.execute(stmt)


def downgrade() -> None:
    op.execute("DROP TABLE IF EXISTS sequence_neighbors CASCADE")
    for stmt in _OLD:
        op.execute(stmt)
