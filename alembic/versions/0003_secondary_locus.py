"""sequences.secondary_locus: dual TRA/TRD membership for shared TCR V genes

The IMGT V genes shared between the TRA and TRD loci carry a combined name
(``TRAV14/DV4*01``) and rearrange in both loci. Identity is one record per
``(species, category, sequence_hash)``, so the sequence cannot be duplicated across
loci; instead the primary ``locus`` holds TRA (from the ``TRAV`` segment) and this
new nullable ``secondary_locus`` holds TRD, making the allele discoverable under both.
NULL for every ordinary single-locus record.

Additive column + index only; no data migration (the value is set on reload).

Revision ID: 0003
Revises: 0002
Create Date: 2026-07-27
"""
from __future__ import annotations

from collections.abc import Sequence

from alembic import op

revision: str = "0003"
down_revision: str | None = "0002"
branch_labels: str | Sequence[str] | None = None
depends_on: str | Sequence[str] | None = None


def upgrade() -> None:
    op.execute("ALTER TABLE sequences ADD COLUMN secondary_locus VARCHAR(8)")
    op.execute(
        "CREATE INDEX ix_sequences_secondary_locus ON sequences (secondary_locus)"
    )


def downgrade() -> None:
    op.execute("DROP INDEX IF EXISTS ix_sequences_secondary_locus")
    op.execute("ALTER TABLE sequences DROP COLUMN secondary_locus")
