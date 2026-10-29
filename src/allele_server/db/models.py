"""SQLAlchemy 2.x models for the allele_server schema.

The network of tables:

    sources          provenance + redistribution + CURIE prefix
    sequences        every unique sequence (coding/leader/rss) + its UID
    names            aliases (labels) per source, never keys
    sequence_links   coding allele <-> its leader/RSS (and future variants)
    sequence_source  which source(s) attest each sequence (CURIE original_id)
    sequence_annotations   per-source key/value metadata

Identity rules mirror ``allele_server.sequence`` / ``uid``: a sequence's
identity is ``(species, category, sequence_hash)``; names are labels paired
with a scheme; AA + protein UID live on V coding rows.

The ``pg_trgm`` extension and the GIN trigram indexes on ``sequences.sequence``,
``sequences.sequence_aa`` and ``names.name`` are created in Alembic migration
``0001`` (they need the extension to exist first), not via ``create_all``.
"""

from __future__ import annotations

from datetime import UTC, datetime

from sqlalchemy import (
    REAL,
    BigInteger,
    Boolean,
    CheckConstraint,
    DateTime,
    Float,
    ForeignKey,
    Index,
    Integer,
    String,
    Text,
    UniqueConstraint,
    text,
)
from sqlalchemy.orm import DeclarativeBase, Mapped, mapped_column, relationship


def _utcnow() -> datetime:
    return datetime.now(UTC)


# Provenance of ``sequences.aa_origin``. The stored AA is always consistent with the
# stored nucleotides by design, which is why a reference AA that structurally conflicts
# with them is demoted to an ``aa/imgt_reference_aa`` annotation and we keep our own
# translation — ``translated_over_conflict`` is that case, and it is the one a reader
# most needs flagged.
# published by a source, and it agrees with our translation
AA_ORIGIN_SOURCED = "sourced"
# no source published one, so we translated the nt
AA_ORIGIN_TRANSLATED = "translated"
# a source published one the nt does not encode; ours is stored, theirs annotated
AA_ORIGIN_TRANSLATED_OVER_CONFLICT = "translated_over_conflict"
AA_ORIGINS = (
    AA_ORIGIN_SOURCED, AA_ORIGIN_TRANSLATED, AA_ORIGIN_TRANSLATED_OVER_CONFLICT,
)


class Base(DeclarativeBase):
    """Shared declarative base; one metadata for the whole schema."""


class Source(Base):
    """A data source, as actually loaded.

    The metadata columns (``manifest_id`` through ``approved_at``) are a **projection**
    of the source's manifest in ``allele_server.sources``, refreshed by
    ``ingest/loader.get_or_create_source`` on every ingest. The manifest is the thing to
    edit; these columns exist so a row is self-describing and so export attribution and
    the completeness test can read the database rather than the files. See
    ``docs/source-policy.md``.
    """

    __tablename__ = "sources"

    id: Mapped[int] = mapped_column(Integer, primary_key=True)
    name: Mapped[str] = mapped_column(String(64), nullable=False)
    version: Mapped[str | None] = mapped_column(String(64), nullable=True)
    # AIRR/standard CURIE prefix for this source (OGRDB, IMGT/GENE-DB, GENBANK,
    # HUSA, MUSA). See allele_server.curie.
    curie_prefix: Mapped[str | None] = mapped_column(String(32), nullable=True)
    # OGRDB germline-set CURIE ref, e.g. "OGRDB:Human_IGH:2021.11".
    germline_set_ref: Mapped[str | None] = mapped_column(String(128), nullable=True)
    # internal | external | derived ('derived' = a layer we compute, e.g. the mined
    # literature evidence, as opposed to a germline set someone published).
    kind: Mapped[str] = mapped_column(String(16), nullable=False)

    # --- projected from the manifest ----------------------------------------------
    # Manifest slug; NULL only for a row written before the manifest existed.
    manifest_id: Mapped[str | None] = mapped_column(String(64), nullable=True)
    license: Mapped[str | None] = mapped_column(String(128), nullable=True)
    licence_url: Mapped[str | None] = mapped_column(String(512), nullable=True)
    citation: Mapped[str | None] = mapped_column(Text, nullable=True)
    # 'verified' (citation holds the reference) | 'none_available' (none exists; the
    # manifest carries the justification). No third value: see migration 0004.
    citation_status: Mapped[str | None] = mapped_column(String(16), nullable=True)
    doi: Mapped[str | None] = mapped_column(String(128), nullable=True)
    url: Mapped[str | None] = mapped_column(String(512), nullable=True)
    # Gates bulk export. False whenever the licence is undeclared: an undeclared
    # licence is not permission.
    redistributable: Mapped[bool] = mapped_column(
        Boolean, nullable=False, server_default=text("false")
    )
    approved_by: Mapped[str | None] = mapped_column(String(128), nullable=True)
    approved_at: Mapped[datetime | None] = mapped_column(
        DateTime(timezone=True), nullable=True
    )

    note: Mapped[str | None] = mapped_column(String(256), nullable=True)
    loaded_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), nullable=False, default=_utcnow
    )

    __table_args__ = (
        UniqueConstraint(
            "name", "version", name="uq_sources_name_version",
            postgresql_nulls_not_distinct=True,
        ),
        CheckConstraint(
            "kind IN ('internal', 'external', 'derived')", name="ck_sources_kind"
        ),
        CheckConstraint(
            "citation_status IS NULL OR citation_status IN ('verified', 'none_available')",
            name="ck_sources_citation_status",
        ),
        # Declared here as well as in migration 0004 so autogenerate never proposes
        # dropping it.
        Index("ix_sources_manifest_id", "manifest_id"),
    )


class Sequence(Base):
    __tablename__ = "sequences"

    id: Mapped[int] = mapped_column(Integer, primary_key=True)
    uid: Mapped[str] = mapped_column(String(32), nullable=False)
    species: Mapped[str] = mapped_column(String(16), nullable=False)
    locus: Mapped[str | None] = mapped_column(String(8), nullable=True)  # IGH/IGK/IGL
    # Second locus for the shared TCR V genes that rearrange in both TRA and TRD
    # (IMGT ``TRAV14/DV4``): locus=TRA, secondary_locus=TRD. NULL for everything else.
    secondary_locus: Mapped[str | None] = mapped_column(String(8), nullable=True)
    category: Mapped[str] = mapped_column(String(16), nullable=False)  # coding|leader|rss
    # AIRR sequence_type for coding (V/D/J); NULL for leader/rss.
    segment: Mapped[str | None] = mapped_column(String(8), nullable=True)
    # l_part1/l_part2/heptamer/spacer/nonamer for non-coding; NULL for coding.
    subtype: Mapped[str | None] = mapped_column(String(16), nullable=True)
    sequence: Mapped[str] = mapped_column(Text, nullable=False)  # AIRR coding_sequence (nt)
    sequence_hash: Mapped[str] = mapped_column(String(64), nullable=False)
    rc_hash: Mapped[str] = mapped_column(String(64), nullable=False)
    # AA representation: V coding only.
    sequence_aa: Mapped[str | None] = mapped_column(Text, nullable=True)
    aa_hash: Mapped[str | None] = mapped_column(String(64), nullable=True)
    protein_uid: Mapped[str | None] = mapped_column(String(32), nullable=True)
    # Where ``sequence_aa`` came from — a property of the one stored AA, so a column
    # rather than a (sequence, source, domain, key) annotation: "we translated this"
    # is not a claim any source makes. One of AA_ORIGINS; NULL where there is no AA.
    aa_origin: Mapped[str | None] = mapped_column(String(24), nullable=True)
    functional: Mapped[bool | None] = mapped_column(Boolean, nullable=True)
    length: Mapped[int] = mapped_column(Integer, nullable=False)
    has_ambiguity: Mapped[bool] = mapped_column(Boolean, nullable=False)
    # True for sequences present only in historical IMGT release snapshots and gone
    # from every current source (retired/archived alleles). Excluded from current
    # coverage stats and default search/similarity results; reachable directly.
    archived: Mapped[bool] = mapped_column(
        Boolean, nullable=False, default=False, server_default="false"
    )
    created_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), nullable=False, default=_utcnow
    )

    names: Mapped[list["Name"]] = relationship(
        back_populates="sequence", cascade="all, delete-orphan"
    )
    attributions: Mapped[list["SequenceSource"]] = relationship(
        back_populates="sequence", cascade="all, delete-orphan"
    )
    annotations: Mapped[list["SequenceAnnotation"]] = relationship(
        back_populates="sequence", cascade="all, delete-orphan"
    )

    __table_args__ = (
        UniqueConstraint("uid", name="uq_sequences_uid"),
        UniqueConstraint(
            "species", "category", "sequence_hash",
            name="uq_sequences_species_category_hash",
        ),
        CheckConstraint("length > 0", name="ck_sequences_length_positive"),
        CheckConstraint(
            "category IN ('coding', 'leader', 'rss')", name="ck_sequences_category"
        ),
        Index("ix_sequences_sequence_hash", "sequence_hash"),
        Index("ix_sequences_aa_hash", "aa_hash"),
        Index("ix_sequences_aa_origin", "aa_origin"),
        Index("ix_sequences_protein_uid", "protein_uid"),
        Index("ix_sequences_rc_hash", "rc_hash"),
        Index("ix_sequences_category_segment", "category", "segment"),
        Index("ix_sequences_segment_length", "segment", "length"),
        Index("ix_sequences_species_locus", "species", "locus"),
        Index("ix_sequences_secondary_locus", "secondary_locus"),
        # Declared here so the model matches the schema: it existed only in the old
        # migration 0004, which meant autogenerate would have proposed dropping it.
        Index("ix_sequences_archived", "archived"),
    )


class Name(Base):
    __tablename__ = "names"

    id: Mapped[int] = mapped_column(Integer, primary_key=True)
    sequence_id: Mapped[int] = mapped_column(
        ForeignKey("sequences.id", ondelete="CASCADE"), nullable=False
    )
    source_id: Mapped[int] = mapped_column(
        ForeignKey("sources.id", ondelete="RESTRICT"), nullable=False
    )
    # 1024, not 64: a curated set may name a novel allele by its seed plus the
    # mutations separating it from that seed (MUSA runs to 784 chars), and that
    # whole string is the identifier -- its head is not unique. See 0007.
    name: Mapped[str] = mapped_column(String(1024), nullable=False)
    naming_scheme: Mapped[str] = mapped_column(String(16), nullable=False)
    confidence: Mapped[float | None] = mapped_column(Float, nullable=True)
    first_seen_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), nullable=False, default=_utcnow
    )
    last_seen_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), nullable=False, default=_utcnow
    )

    sequence: Mapped[Sequence] = relationship(back_populates="names")
    source: Mapped[Source] = relationship()

    __table_args__ = (
        # Includes sequence_id: a name is a *label* and the same allele name can
        # legitimately apply to several distinct sequences (e.g. IMGT strain
        # variants reuse IGHV1-2*01). Idempotency is per (sequence, name, scheme,
        # source); only re-attaching the identical name to the SAME sequence is a
        # no-op. (Was (name, scheme, source), which silently dropped the name from
        # every strain variant after the first.)
        UniqueConstraint(
            "name", "naming_scheme", "source_id", "sequence_id",
            name="uq_names_name_scheme_source_seq",
        ),
        Index("ix_names_sequence_id", "sequence_id"),
    )


class SequenceLink(Base):
    __tablename__ = "sequence_links"

    id: Mapped[int] = mapped_column(Integer, primary_key=True)
    from_sequence_id: Mapped[int] = mapped_column(
        ForeignKey("sequences.id", ondelete="CASCADE"), nullable=False
    )
    to_sequence_id: Mapped[int] = mapped_column(
        ForeignKey("sequences.id", ondelete="CASCADE"), nullable=False
    )
    relation: Mapped[str] = mapped_column(String(32), nullable=False)
    source_id: Mapped[int] = mapped_column(
        ForeignKey("sources.id", ondelete="RESTRICT"), nullable=False
    )
    first_seen_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), nullable=False, default=_utcnow
    )
    last_seen_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), nullable=False, default=_utcnow
    )

    __table_args__ = (
        UniqueConstraint(
            "from_sequence_id", "to_sequence_id", "relation", "source_id",
            name="uq_sequence_links_natkey",
        ),
        CheckConstraint(
            "from_sequence_id <> to_sequence_id", name="ck_sequence_links_no_self"
        ),
        Index("ix_sequence_links_from", "from_sequence_id"),
        Index("ix_sequence_links_to", "to_sequence_id"),
    )


class SequenceSource(Base):
    __tablename__ = "sequence_source"

    id: Mapped[int] = mapped_column(Integer, primary_key=True)
    sequence_id: Mapped[int] = mapped_column(
        ForeignKey("sequences.id", ondelete="CASCADE"), nullable=False
    )
    source_id: Mapped[int] = mapped_column(
        ForeignKey("sources.id", ondelete="RESTRICT"), nullable=False
    )
    # CURIE where a standard prefix exists (GENBANK:X02897, OGRDB:..., ...).
    original_id: Mapped[str] = mapped_column(String(1024), nullable=False)
    original_sequence: Mapped[str] = mapped_column(Text, nullable=False)
    first_seen_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), nullable=False, default=_utcnow
    )

    sequence: Mapped[Sequence] = relationship(back_populates="attributions")
    source: Mapped[Source] = relationship()

    __table_args__ = (
        UniqueConstraint(
            "sequence_id", "source_id", "original_id",
            name="uq_sequence_source_natkey",
        ),
        Index("ix_sequence_source_sequence_id", "sequence_id"),
        Index("ix_sequence_source_source_id", "source_id"),
    )


class SequenceAnnotation(Base):
    __tablename__ = "sequence_annotations"

    id: Mapped[int] = mapped_column(Integer, primary_key=True)
    sequence_id: Mapped[int] = mapped_column(
        ForeignKey("sequences.id", ondelete="CASCADE"), nullable=False
    )
    source_id: Mapped[int] = mapped_column(
        ForeignKey("sources.id", ondelete="RESTRICT"), nullable=False
    )
    domain: Mapped[str] = mapped_column(String(32), nullable=False)
    key: Mapped[str] = mapped_column(String(64), nullable=False)
    value_text: Mapped[str | None] = mapped_column(Text, nullable=True)
    value_float: Mapped[float | None] = mapped_column(Float, nullable=True)
    value_int: Mapped[int | None] = mapped_column(Integer, nullable=True)
    value_bool: Mapped[bool | None] = mapped_column(Boolean, nullable=True)
    first_seen_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), nullable=False, default=_utcnow
    )
    last_seen_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), nullable=False, default=_utcnow
    )

    sequence: Mapped[Sequence] = relationship(back_populates="annotations")

    __table_args__ = (
        UniqueConstraint(
            "sequence_id", "source_id", "domain", "key",
            name="uq_sequence_annotations_natkey",
        ),
        CheckConstraint(
            "(CASE WHEN value_text  IS NOT NULL THEN 1 ELSE 0 END) + "
            "(CASE WHEN value_float IS NOT NULL THEN 1 ELSE 0 END) + "
            "(CASE WHEN value_int   IS NOT NULL THEN 1 ELSE 0 END) + "
            "(CASE WHEN value_bool  IS NOT NULL THEN 1 ELSE 0 END) = 1",
            name="ck_sequence_annotations_exactly_one_value",
        ),
        Index("ix_sequence_annotations_seq_domain", "sequence_id", "domain"),
    )


class SequenceNeighbor(Base):
    """Precomputed top-k nearest neighbours per sequence (a similarity cache).

    Directed edges keyed on integer ``sequences.id`` (``from_id`` -> ``to_id``),
    with the natural key ``(from_id, kind, to_id)`` as the primary key — no surrogate
    ``id``. Integer keys make the two big indexes ~4x smaller than the former text-uid
    keys (3.1 GB cache -> ~1 GB). The ``to_*`` attributes stay denormalised so a cached
    read still assembles each ``SearchHit`` with no per-hit alias join.

    Safe to key on ``id`` (renumbered on reload) rather than the content-stable ``uid``
    **because the cache is always fully rebuilt on reload** — the ``TRUNCATE`` cascades
    into this table (FK on ``sequences.id``), so the ids are never stale. Do not
    renumber ``sequences.id`` without rebuilding this table.
    """

    __tablename__ = "sequence_neighbors"

    from_id: Mapped[int] = mapped_column(
        ForeignKey("sequences.id", ondelete="CASCADE"), primary_key=True
    )
    kind: Mapped[str] = mapped_column(String(8), primary_key=True)  # nt | aa
    to_id: Mapped[int] = mapped_column(
        ForeignKey("sequences.id", ondelete="CASCADE"), primary_key=True
    )
    identity: Mapped[float] = mapped_column(REAL, nullable=False)  # 0..1

    # Denormalised attributes of the *to* sequence, copied at build time so the
    # cached read assembles each SearchHit with no Sequence fetch / alias join.
    to_uid: Mapped[str | None] = mapped_column(String(32), nullable=True)
    # Must track names.name (1024): this is a denormalised copy of it, so a narrower
    # column here fails the rebuild on the first long curated label. See 0008.
    label: Mapped[str | None] = mapped_column(String(1024), nullable=True)
    to_species: Mapped[str | None] = mapped_column(String(16), nullable=True)
    to_locus: Mapped[str | None] = mapped_column(String(8), nullable=True)
    to_segment: Mapped[str | None] = mapped_column(String(8), nullable=True)
    to_category: Mapped[str | None] = mapped_column(String(16), nullable=True)

    __table_args__ = (
        # top-k read: WHERE from_id=? AND kind=? ORDER BY identity DESC
        Index("ix_neighbor_topk", "from_id", "kind", "identity"),
        # look up / exclude by target sequence (archived filter joins on to_id)
        Index("ix_neighbor_to_id", "to_id"),
    )


class SourceCandidate(Base):
    """A prospective data source, between "discovered" and "decided".

    Discovery writes rows here and nothing else; the three decision states
    (``approved`` / ``rejected`` / ``deferred``) require a named person and a timestamp,
    enforced by ``ck_source_candidates_decided`` rather than by the code that would
    otherwise be trusted to check. See ``docs/source-policy.md``.

    **No foreign key to sources.id, on purpose.** ``reload_all`` truncates ``sources``
    with CASCADE; an FK would drag every approval and rejection down with it. The link
    is the plain-text ``manifest_id``. See migration 0005.
    """

    __tablename__ = "source_candidates"

    STATUSES = (
        "discovered", "triaged", "dossier", "in_review",
        "approved", "rejected", "deferred", "ingested",
    )
    #: Statuses only a human may set (and only with decided_by + decided_at).
    DECIDED = ("approved", "rejected", "deferred")
    CHANNELS = (
        "europepmc", "ogrdb", "imgt", "zenodo", "figshare", "github",
        "ncbi_bioproject", "ncbi_sra", "manual",
    )

    id: Mapped[int] = mapped_column(Integer, primary_key=True)
    channel: Mapped[str] = mapped_column(String(24), nullable=False)
    # Stable per channel: a PMID, an OGRDB set CURIE, an IMGT release, a DOI, owner/repo.
    external_id: Mapped[str] = mapped_column(String(256), nullable=False)
    # Cross-channel identity: the normalised DOI when there is one, else
    # "<channel>:<external_id>". A rejection on any row suppresses every row sharing it,
    # so a rejected deposit does not return as the paper announcing it.
    dedupe_key: Mapped[str | None] = mapped_column(String(256), nullable=True)

    title: Mapped[str | None] = mapped_column(Text, nullable=True)
    url: Mapped[str | None] = mapped_column(String(1024), nullable=True)
    doi: Mapped[str | None] = mapped_column(String(128), nullable=True)
    licence: Mapped[str | None] = mapped_column(String(128), nullable=True)
    species_hint: Mapped[str | None] = mapped_column(String(256), nullable=True)
    loci_hint: Mapped[str | None] = mapped_column(String(128), nullable=True)

    status: Mapped[str] = mapped_column(
        String(16), nullable=False, server_default=text("'discovered'")
    )
    score: Mapped[float | None] = mapped_column(REAL, nullable=True)
    suppressed: Mapped[bool] = mapped_column(
        Boolean, nullable=False, server_default=text("false")
    )
    manifest_id: Mapped[str | None] = mapped_column(String(64), nullable=True)
    dossier_path: Mapped[str | None] = mapped_column(String(512), nullable=True)

    decided_by: Mapped[str | None] = mapped_column(String(128), nullable=True)
    decided_at: Mapped[datetime | None] = mapped_column(
        DateTime(timezone=True), nullable=True
    )
    decision_note: Mapped[str | None] = mapped_column(Text, nullable=True)
    deferred_until: Mapped[datetime | None] = mapped_column(
        DateTime(timezone=True), nullable=True
    )

    # Raw channel record + the cache paths it came from, as JSON, so the dossier can be
    # regenerated without re-crawling.
    evidence: Mapped[str | None] = mapped_column(Text, nullable=True)
    first_seen_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), nullable=False, default=_utcnow
    )
    last_seen_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), nullable=False, default=_utcnow
    )

    __table_args__ = (
        UniqueConstraint("channel", "external_id", name="uq_source_candidates_natkey"),
        CheckConstraint(
            "status IN ('discovered', 'triaged', 'dossier', 'in_review', 'approved', "
            "'rejected', 'deferred', 'ingested')",
            name="ck_source_candidates_status",
        ),
        CheckConstraint(
            "channel IN ('europepmc', 'ogrdb', 'imgt', 'zenodo', 'figshare', 'github', "
            "'ncbi_bioproject', 'ncbi_sra', 'manual')",
            name="ck_source_candidates_channel",
        ),
        # The human-in-the-loop requirement, as DDL. Nothing can approve anonymously.
        CheckConstraint(
            "status NOT IN ('approved', 'rejected', 'deferred') "
            "OR (decided_by IS NOT NULL AND decided_at IS NOT NULL)",
            name="ck_source_candidates_decided",
        ),
        Index("ix_source_candidates_status", "status"),
        Index("ix_source_candidates_dedupe_key", "dedupe_key"),
        Index("ix_source_candidates_manifest_id", "manifest_id"),
    )

    @property
    def ref(self) -> str:
        """The CLI handle for this candidate: ``<channel>:<external_id>``."""
        return f"{self.channel}:{self.external_id}"


class DatasetRelease(Base):
    """One stamped build of the served dataset: what was loaded, and how it is scored.

    Written by ``ingest.reload.stamp_release`` at the end of a reload (and by
    ``allele-ingest --stamp-release`` for a database built before this table existed).
    Read by ``queries.dataset_release`` for ``/version``, ``/api/export`` and the
    landing page.

    Two digests, because a served answer depends on two things that move independently:

    * ``content_digest`` — the data. sha256 over the sequence identity set and the
      name set, so an instance holding different alleles, or the same alleles under
      different names, reports a different string.
    * ``similarity_digest`` — the **distance computation**. sha256 over
      ``queries.similarity_profile()``. Identity scores, neighbour lists and every
      threshold on the site are functions of those parameters, so a release that did
      not pin them would not be reproducible.

    ``neighbors_built`` is not cosmetic: with the cache absent, "find similar" falls
    through to the live ranker. Same metric and same answers, but a different
    performance and completeness story, so a release states which one it is.

    **No FK to sequences, and not in ``reload._TABLES``** (see migration 0009): the
    release history must outlive the data it describes.
    """

    __tablename__ = "dataset_releases"

    id: Mapped[int] = mapped_column(Integer, primary_key=True)
    #: Human-facing version, ``YYYY.MM.DD`` (``.N`` for the Nth build of a day).
    release: Mapped[str] = mapped_column(String(32), nullable=False)
    built_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), nullable=False, default=_utcnow
    )
    content_digest: Mapped[str] = mapped_column(String(64), nullable=False)
    similarity_digest: Mapped[str] = mapped_column(String(64), nullable=False)
    code_version: Mapped[str] = mapped_column(String(32), nullable=False)
    neighbors_built: Mapped[bool] = mapped_column(
        Boolean, nullable=False, server_default=text("false")
    )
    neighbor_edges: Mapped[int | None] = mapped_column(BigInteger, nullable=True)
    #: JSON: counts by segment/category, sources with their versions, and the full
    #: similarity profile the digest was taken over. Text, matching how the rest of
    #: this schema stores JSON (see SourceCandidate.evidence).
    details: Mapped[str | None] = mapped_column(Text, nullable=True)

    __table_args__ = (
        UniqueConstraint("release", name="uq_dataset_releases_release"),
        Index("ix_dataset_releases_built_at", "built_at"),
    )
