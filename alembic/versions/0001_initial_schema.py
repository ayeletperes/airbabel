"""initial schema (squashed): sources, sequences, names, links, source, annotations, neighbours

Squashes the former 0001-0005. Those revisions are gone: the old 0001 built its tables
with ``Base.metadata.create_all`` against the *live* models, so a fresh database already
had ``archived`` and the newer unique key, and 0004/0005 then failed on it. ``alembic
upgrade head`` could not create a database at all.

The DDL below is therefore **frozen literal SQL**, generated once from the models and
never re-read from them. That is the point of the squash: a baseline that cannot drift.
A future schema change is a new revision with its own explicit DDL. Do not reintroduce
``create_all`` here, or the next additive migration breaks in exactly the same way.

Revision ID: 0001
Revises:
Create Date: 2026-07-10
"""
from __future__ import annotations

from collections.abc import Sequence

from alembic import op

revision: str = "0001"
down_revision: str | None = None
branch_labels: str | Sequence[str] | None = None
depends_on: str | Sequence[str] | None = None


# GIN trigram indexes for fuzzy name / sequence search. They reference gin_trgm_ops,
# so they come after the tables and after the pg_trgm extension.
_TRGM_INDEXES = (
    ("ix_names_name_trgm", "names", "name"),
    ("ix_sequences_sequence_trgm", "sequences", "sequence"),
    ("ix_sequences_sequence_aa_trgm", "sequences", "sequence_aa"),
)

_DDL: tuple[str, ...] = (
    """
CREATE TABLE sequences (
	id SERIAL NOT NULL,
	uid VARCHAR(32) NOT NULL,
	species VARCHAR(16) NOT NULL,
	locus VARCHAR(8),
	category VARCHAR(16) NOT NULL,
	segment VARCHAR(8),
	subtype VARCHAR(16),
	sequence TEXT NOT NULL,
	sequence_hash VARCHAR(64) NOT NULL,
	rc_hash VARCHAR(64) NOT NULL,
	sequence_aa TEXT,
	aa_hash VARCHAR(64),
	protein_uid VARCHAR(32),
	functional BOOLEAN,
	length INTEGER NOT NULL,
	has_ambiguity BOOLEAN NOT NULL,
	archived BOOLEAN DEFAULT 'false' NOT NULL,
	created_at TIMESTAMP WITH TIME ZONE NOT NULL,
	PRIMARY KEY (id),
	CONSTRAINT uq_sequences_uid UNIQUE (uid),
	CONSTRAINT uq_sequences_species_category_hash UNIQUE (species, category, sequence_hash),
	CONSTRAINT ck_sequences_length_positive CHECK (length > 0),
	CONSTRAINT ck_sequences_category CHECK (category IN ('coding', 'leader', 'rss'))
)
    """,
    """
CREATE INDEX ix_sequences_aa_hash ON sequences (aa_hash)
    """,
    """
CREATE INDEX ix_sequences_archived ON sequences (archived)
    """,
    """
CREATE INDEX ix_sequences_category_segment ON sequences (category, segment)
    """,
    """
CREATE INDEX ix_sequences_protein_uid ON sequences (protein_uid)
    """,
    """
CREATE INDEX ix_sequences_rc_hash ON sequences (rc_hash)
    """,
    """
CREATE INDEX ix_sequences_segment_length ON sequences (segment, length)
    """,
    """
CREATE INDEX ix_sequences_sequence_hash ON sequences (sequence_hash)
    """,
    """
CREATE INDEX ix_sequences_species_locus ON sequences (species, locus)
    """,
    """
CREATE TABLE sources (
	id SERIAL NOT NULL,
	name VARCHAR(64) NOT NULL,
	version VARCHAR(64),
	curie_prefix VARCHAR(32),
	germline_set_ref VARCHAR(128),
	kind VARCHAR(16) NOT NULL,
	license VARCHAR(64),
	redistributable BOOLEAN DEFAULT false NOT NULL,
	note VARCHAR(256),
	loaded_at TIMESTAMP WITH TIME ZONE NOT NULL,
	PRIMARY KEY (id),
	CONSTRAINT uq_sources_name_version UNIQUE NULLS NOT DISTINCT (name, version),
	CONSTRAINT ck_sources_kind CHECK (kind IN ('internal', 'external'))
)
    """,
    """
CREATE TABLE names (
	id SERIAL NOT NULL,
	sequence_id INTEGER NOT NULL,
	source_id INTEGER NOT NULL,
	name VARCHAR(256) NOT NULL,
	naming_scheme VARCHAR(16) NOT NULL,
	confidence FLOAT,
	first_seen_at TIMESTAMP WITH TIME ZONE NOT NULL,
	last_seen_at TIMESTAMP WITH TIME ZONE NOT NULL,
	PRIMARY KEY (id),
	CONSTRAINT uq_names_name_scheme_source_seq UNIQUE (name, naming_scheme, source_id, sequence_id),
	FOREIGN KEY(sequence_id) REFERENCES sequences (id) ON DELETE CASCADE,
	FOREIGN KEY(source_id) REFERENCES sources (id) ON DELETE RESTRICT
)
    """,
    """
CREATE INDEX ix_names_sequence_id ON names (sequence_id)
    """,
    """
CREATE TABLE sequence_annotations (
	id SERIAL NOT NULL,
	sequence_id INTEGER NOT NULL,
	source_id INTEGER NOT NULL,
	domain VARCHAR(32) NOT NULL,
	key VARCHAR(64) NOT NULL,
	value_text TEXT,
	value_float FLOAT,
	value_int INTEGER,
	value_bool BOOLEAN,
	first_seen_at TIMESTAMP WITH TIME ZONE NOT NULL,
	last_seen_at TIMESTAMP WITH TIME ZONE NOT NULL,
	PRIMARY KEY (id),
	CONSTRAINT uq_sequence_annotations_natkey UNIQUE (sequence_id, source_id, domain, key),
	CONSTRAINT ck_sequence_annotations_exactly_one_value CHECK ((CASE WHEN value_text  IS NOT NULL THEN 1 ELSE 0 END) + (CASE WHEN value_float IS NOT NULL THEN 1 ELSE 0 END) + (CASE WHEN value_int   IS NOT NULL THEN 1 ELSE 0 END) + (CASE WHEN value_bool  IS NOT NULL THEN 1 ELSE 0 END) = 1),
	FOREIGN KEY(sequence_id) REFERENCES sequences (id) ON DELETE CASCADE,
	FOREIGN KEY(source_id) REFERENCES sources (id) ON DELETE RESTRICT
)
    """,
    """
CREATE INDEX ix_sequence_annotations_seq_domain ON sequence_annotations (sequence_id, domain)
    """,
    """
CREATE TABLE sequence_links (
	id SERIAL NOT NULL,
	from_sequence_id INTEGER NOT NULL,
	to_sequence_id INTEGER NOT NULL,
	relation VARCHAR(32) NOT NULL,
	source_id INTEGER NOT NULL,
	first_seen_at TIMESTAMP WITH TIME ZONE NOT NULL,
	last_seen_at TIMESTAMP WITH TIME ZONE NOT NULL,
	PRIMARY KEY (id),
	CONSTRAINT uq_sequence_links_natkey UNIQUE (from_sequence_id, to_sequence_id, relation, source_id),
	CONSTRAINT ck_sequence_links_no_self CHECK (from_sequence_id <> to_sequence_id),
	FOREIGN KEY(from_sequence_id) REFERENCES sequences (id) ON DELETE CASCADE,
	FOREIGN KEY(to_sequence_id) REFERENCES sequences (id) ON DELETE CASCADE,
	FOREIGN KEY(source_id) REFERENCES sources (id) ON DELETE RESTRICT
)
    """,
    """
CREATE INDEX ix_sequence_links_from ON sequence_links (from_sequence_id)
    """,
    """
CREATE INDEX ix_sequence_links_to ON sequence_links (to_sequence_id)
    """,
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
    """
CREATE INDEX ix_neighbor_topk ON sequence_neighbors (from_uid, kind, identity)
    """,
    """
CREATE TABLE sequence_source (
	id SERIAL NOT NULL,
	sequence_id INTEGER NOT NULL,
	source_id INTEGER NOT NULL,
	original_id VARCHAR(128) NOT NULL,
	original_sequence TEXT NOT NULL,
	first_seen_at TIMESTAMP WITH TIME ZONE NOT NULL,
	PRIMARY KEY (id),
	CONSTRAINT uq_sequence_source_natkey UNIQUE (sequence_id, source_id, original_id),
	FOREIGN KEY(sequence_id) REFERENCES sequences (id) ON DELETE CASCADE,
	FOREIGN KEY(source_id) REFERENCES sources (id) ON DELETE RESTRICT
)
    """,
    """
CREATE INDEX ix_sequence_source_sequence_id ON sequence_source (sequence_id)
    """,
    """
CREATE INDEX ix_sequence_source_source_id ON sequence_source (source_id)
    """,
)


# Reverse dependency order, for downgrade().
_TABLES = ('sequence_source', 'sequence_neighbors', 'sequence_links', 'sequence_annotations', 'names', 'sources', 'sequences')


def upgrade() -> None:
    op.execute("CREATE EXTENSION IF NOT EXISTS pg_trgm")
    for stmt in _DDL:
        op.execute(stmt)
    for name, table, column in _TRGM_INDEXES:
        op.execute(
            f"CREATE INDEX IF NOT EXISTS {name} ON {table} "
            f"USING gin ({column} gin_trgm_ops)"
        )


def downgrade() -> None:
    for name, _table, _column in _TRGM_INDEXES:
        op.execute(f"DROP INDEX IF EXISTS {name}")
    for table in _TABLES:
        op.execute(f"DROP TABLE IF EXISTS {table} CASCADE")
