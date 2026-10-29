"""source_candidates: the discovery -> review -> decision queue

Where a candidate source lives between "a crawler found it" and "a human decided".
The automation may set any status up to ``dossier``; the three decision states are
reachable only with a named person and a timestamp, and that is enforced here rather
than in Python:

    CONSTRAINT ck_source_candidates_decided CHECK (
        status NOT IN ('approved','rejected','deferred')
        OR (decided_by IS NOT NULL AND decided_at IS NOT NULL))

That constraint *is* the human-review requirement. A rule that lives only in the code
meant to obey it is a convention, not a gate.

**There is deliberately no foreign key to sources.id.** ``ingest/reload.py::reload_all``
runs ``TRUNCATE ... sources RESTART IDENTITY CASCADE``; an FK would put this table
inside that cascade, and a routine ``--reload`` would silently destroy every approval
and every rejection — curatorial state that is not derived from any source and cannot
be rebuilt. The link to a manifest is the plain text ``manifest_id``, which survives
precisely because it is not referential. See ``docs/source-policy.md``.

Revision ID: 0005
Revises: 0004
Create Date: 2026-07-29
"""
from __future__ import annotations

from collections.abc import Sequence

from alembic import op

revision: str = "0005"
down_revision: str | None = "0004"
branch_labels: str | Sequence[str] | None = None
depends_on: str | Sequence[str] | None = None


_UP: tuple[str, ...] = (
    """
    CREATE TABLE source_candidates (
        id SERIAL NOT NULL,
        channel VARCHAR(24) NOT NULL,
        external_id VARCHAR(256) NOT NULL,
        dedupe_key VARCHAR(256),
        title TEXT,
        url VARCHAR(1024),
        doi VARCHAR(128),
        licence VARCHAR(128),
        species_hint VARCHAR(256),
        loci_hint VARCHAR(128),
        status VARCHAR(16) DEFAULT 'discovered' NOT NULL,
        score REAL,
        suppressed BOOLEAN DEFAULT false NOT NULL,
        manifest_id VARCHAR(64),
        dossier_path VARCHAR(512),
        decided_by VARCHAR(128),
        decided_at TIMESTAMP WITH TIME ZONE,
        decision_note TEXT,
        deferred_until TIMESTAMP WITH TIME ZONE,
        evidence TEXT,
        first_seen_at TIMESTAMP WITH TIME ZONE NOT NULL,
        last_seen_at TIMESTAMP WITH TIME ZONE NOT NULL,
        PRIMARY KEY (id),
        CONSTRAINT uq_source_candidates_natkey UNIQUE (channel, external_id),
        CONSTRAINT ck_source_candidates_status CHECK (status IN
            ('discovered', 'triaged', 'dossier', 'in_review', 'approved', 'rejected',
             'deferred', 'ingested')),
        CONSTRAINT ck_source_candidates_channel CHECK (channel IN
            ('europepmc', 'ogrdb', 'imgt', 'zenodo', 'figshare', 'github',
             'ncbi_bioproject', 'ncbi_sra', 'manual')),
        CONSTRAINT ck_source_candidates_decided CHECK (
            status NOT IN ('approved', 'rejected', 'deferred')
            OR (decided_by IS NOT NULL AND decided_at IS NOT NULL))
    )
    """,
    "CREATE INDEX ix_source_candidates_status ON source_candidates (status)",
    "CREATE INDEX ix_source_candidates_dedupe_key ON source_candidates (dedupe_key)",
    "CREATE INDEX ix_source_candidates_manifest_id ON source_candidates (manifest_id)",
)

_DOWN: tuple[str, ...] = (
    "DROP TABLE IF EXISTS source_candidates",
)


def upgrade() -> None:
    for stmt in _UP:
        op.execute(stmt)


def downgrade() -> None:
    for stmt in _DOWN:
        op.execute(stmt)
