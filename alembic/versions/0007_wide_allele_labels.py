"""Widen names.name and sequence_source.original_id so a curated allele label fits

MUSA names a novel allele by its seed plus the mutations that separate it from that
seed -- ``IGLV2-66YQ*02_i0gg_c1g_a2g_g3c_g7a_c8t_10tgaggg15_c20a_...`` -- which runs
to 784 characters. That string is not decoration: measured over the whole table it is
a perfect identifier, 3,777 distinct values against 3,777 distinct sequences, one
each way. Truncating it to the head before the first underscore is what it looks
like it should be, and is wrong: 649 heads collide, covering 2,153 sequences, so
``IGHV1-F5R3*02`` alone would name ten different alleles. A label that does not
identify an allele is the defect we removed when we stopped naming by the gene-level
``*00`` IgLabel; reintroducing it in allele-level clothing would be worse, because it
would look right.

So the identifier stays whole and the columns grow to hold it. 1024 rather than
unbounded ``Text``: both columns sit inside unique btree indexes, where an unbounded
value can exceed the index-entry limit and fail at insert time on whichever row
happens to be pathological. A bound of 1024 clears the widest real value (784, plus a
``MUSA:`` prefix) with headroom, keeps every index entry far below the btree limit,
and leaves an over-long value something the ingest guards can refuse deliberately and
count, rather than discover as a DataError halfway through a table.

Widening only. No value is truncated, so the downgrade can only be attempted, and it
raises if any row would not survive it.

Revision ID: 0007
Revises: 0006
Create Date: 2026-08-13
"""
from __future__ import annotations

from collections.abc import Sequence

from alembic import op

revision: str = "0007"
down_revision: str | None = "0006"
branch_labels: str | Sequence[str] | None = None
depends_on: str | Sequence[str] | None = None


def upgrade() -> None:
    op.execute("ALTER TABLE names ALTER COLUMN name TYPE VARCHAR(1024)")
    op.execute(
        "ALTER TABLE sequence_source ALTER COLUMN original_id TYPE VARCHAR(1024)"
    )


def downgrade() -> None:
    # Narrowing loses data. Refuse rather than silently truncate an identifier.
    op.execute(
        """
        DO $$
        DECLARE n bigint;
        BEGIN
            SELECT count(*) INTO n FROM names WHERE length(name) > 256;
            IF n > 0 THEN
                RAISE EXCEPTION
                    '% names exceed 256 chars; downgrading would truncate them', n;
            END IF;
            SELECT count(*) INTO n
                FROM sequence_source WHERE length(original_id) > 128;
            IF n > 0 THEN
                RAISE EXCEPTION
                    '% original_ids exceed 128 chars; downgrading would truncate them',
                    n;
            END IF;
        END $$;
        """
    )
    op.execute("ALTER TABLE names ALTER COLUMN name TYPE VARCHAR(256)")
    op.execute(
        "ALTER TABLE sequence_source ALTER COLUMN original_id TYPE VARCHAR(128)"
    )
