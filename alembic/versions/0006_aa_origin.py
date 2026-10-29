"""sequences.aa_origin: whether the stored AA was published or computed by us

``loader.set_v_aa`` produces the stored AA three different ways — a sourced AA that
agrees with our translation, our translation where no source published one, and our
translation kept over a sourced AA that structurally conflicts with the nucleotides
(the conflicting reference is demoted to an ``aa/imgt_reference_aa`` annotation). Only
the third case left any trace, so a reader could not tell a published protein from one
we computed.

A column rather than a ``sequence_annotations`` row: annotations are keyed
``(sequence, source, domain, key)``, and "we translated this" is not a claim any source
makes — recording it against IMGT would invert its meaning, and the pure-translation
path has no source at all. Provenance here describes the one stored AA, so it is
single-valued, like the AA itself.

Additive nullable column + index, and deliberately no data migration. NULL means "this
row predates the column and we do not know which path produced its AA" — not
"translated". Guessing would relabel every genuinely published IMGT protein as ours on
any database loaded before the column existed. Origins are re-derived honestly by a
reload, which calls ``set_v_aa`` down every path.

Revision ID: 0006
Revises: 0005
Create Date: 2026-07-31
"""
from __future__ import annotations

from collections.abc import Sequence

from alembic import op

revision: str = "0006"
down_revision: str | None = "0005"
branch_labels: str | Sequence[str] | None = None
depends_on: str | Sequence[str] | None = None


def upgrade() -> None:
    op.execute("ALTER TABLE sequences ADD COLUMN aa_origin VARCHAR(24)")
    op.execute("CREATE INDEX ix_sequences_aa_origin ON sequences (aa_origin)")


def downgrade() -> None:
    op.execute("DROP INDEX IF EXISTS ix_sequences_aa_origin")
    op.execute("ALTER TABLE sequences DROP COLUMN aa_origin")
