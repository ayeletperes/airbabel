"""dataset_releases: one stamped, verifiable build of the served dataset

What a served instance can currently say about itself is "the newest row was written
on <date>". That is a load timestamp, not a release: it does not say *which* data, and
it says nothing at all about the **distance computation**, even though every identity
score, every neighbour list and every threshold on the site is a function of it. Two
instances can hold the same alleles and rank them differently — different floors, a
different scorer for a segment, or a neighbour cache built before a parameter changed —
and nothing on either page would reveal it.

So a release row records both halves and digests each:

* ``content_digest``  — sha256 over the sequence identity set (species | category |
  sequence_hash) **and** the name set (name | scheme | source). Recomputable from the
  database, so two instances can be compared by one string.
* ``similarity_digest`` — sha256 over ``queries.similarity_profile()``: the metric,
  the per-segment equal-length policy, the cache floors and cap. Changing
  ``NEIGHBOR_FLOOR`` changes this digest, which is the point.

**No foreign key to sequences, and not in reload's TRUNCATE list** — deliberately, and
for the same reason as ``source_candidates`` (migration 0005): the release history is
the record of what was served and when. A reload writes the next row; it must never
erase the previous ones.

Revision ID: 0009
Revises: 0008
Create Date: 2026-09-01
"""
from __future__ import annotations

from collections.abc import Sequence

from alembic import op

revision: str = "0009"
down_revision: str | None = "0008"
branch_labels: str | Sequence[str] | None = None
depends_on: str | Sequence[str] | None = None


_UP: tuple[str, ...] = (
    """
    CREATE TABLE dataset_releases (
        id SERIAL NOT NULL,
        release VARCHAR(32) NOT NULL,
        built_at TIMESTAMP WITH TIME ZONE NOT NULL,
        content_digest VARCHAR(64) NOT NULL,
        similarity_digest VARCHAR(64) NOT NULL,
        code_version VARCHAR(32) NOT NULL,
        neighbors_built BOOLEAN DEFAULT false NOT NULL,
        neighbor_edges BIGINT,
        details TEXT,
        PRIMARY KEY (id),
        CONSTRAINT uq_dataset_releases_release UNIQUE (release)
    )
    """,
    "CREATE INDEX ix_dataset_releases_built_at ON dataset_releases (built_at)",
)

_DOWN: tuple[str, ...] = (
    "DROP TABLE IF EXISTS dataset_releases",
)


def upgrade() -> None:
    for stmt in _UP:
        op.execute(stmt)


def downgrade() -> None:
    for stmt in _DOWN:
        op.execute(stmt)
