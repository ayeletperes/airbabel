"""Pydantic request/response schemas.

Every record is returned with two views: a ``native`` block (our UIDs, hashes,
source tags) and an ``airr`` block shaped as an AIRR ``AlleleDescription`` for
semantic interoperability.
"""

from __future__ import annotations

from pydantic import BaseModel, Field


class Ontology(BaseModel):
    """AIRR Ontology value (CURIE id + label), e.g. NCBITaxon for species."""

    id: str
    label: str


class SourceTag(BaseModel):
    """One source attesting a sequence, with its CURIE cross-reference."""

    name: str
    version: str | None = None
    curie_prefix: str | None = None
    redistributable: bool
    original_id: str | None = Field(
        default=None, description="CURIE cross-reference, e.g. GENBANK:X02897"
    )
    url: str | None = None


class Alias(BaseModel):
    name: str
    scheme: str
    source: str
    url: str | None = None
    # A former-release IMGT name (recovered from a historical snapshot), and the
    # IMGT release it came from — not the current official designation.
    historical: bool = False
    release: str | None = None


class AliasGroup(BaseModel):
    """One distinct name (a ``(name, scheme)`` pair), with **every source** that
    attributes it collapsed into ``sources``. The same designation often comes from
    many sources (and IMGT spans many release versions); grouping keeps the names
    table to one row per name instead of one per source-version."""

    name: str
    scheme: str
    url: str | None = None
    historical: bool = False
    release: str | None = None
    sources: list[str] = []


class RecordVersion(BaseModel):
    """Another **record of the same allele**: a different sequence carrying the same
    IMGT designation in a different release range.

    IMGT revises the extent of a record without changing its name — a partial
    V-REGION later completed, or a longer record trimmed. Because a sequence's
    identity is its hash, each extent is a separate row, and the older one loses its
    current name. Nothing was renamed and no allele went away; only the span of
    sequence IMGT publishes under that name changed.

    ``relation`` says how this record relates to the one being viewed:
    ``"extends"`` (it contains the viewed sequence), ``"trims"`` (it is contained in
    the viewed sequence), or ``"differs"`` (neither, so the bases changed too).

    Only ever populated from **disjoint** release ranges. Two sequences that held one
    name in the *same* release are a genuine ambiguity, not a version chain, and are
    deliberately excluded.
    """

    uid: str
    name: str
    length: int
    relation: str
    current: bool = False
    first_seen: str | None = None
    last_seen: str | None = None


class Xref(BaseModel):
    """A cross-reference to an external resource (CURIE + resolved URL)."""

    label: str
    curie: str
    url: str | None = None


class Reference(BaseModel):
    """A literature reference (PubMed) attesting an allele.

    ``kind`` tiers the evidence, strongest first:

    * ``"verified"`` — curator-provided (from the internal tables / GenBank xref).
    * ``"sequence"`` — auto-mined, but the allele's **sequence** was found in the
      paper's full text or a supplementary file (hash match; cross-species-safe).
    * ``"name"`` — auto-mined by **name** only (a lead; the same name can denote a
      different allele in another species, so treat as a pointer, not proof).

    ``verified`` (kept for back-compat) is true iff ``kind == "verified"``.
    ``detail`` carries the mining provenance for a sequence/name hit (e.g.
    ``"sequence_exact · fulltext:PMC1234567"``)."""

    pmid: str
    url: str | None = None
    title: str | None = None
    verified: bool = False
    kind: str = "name"
    detail: str | None = None
    source: str | None = None


class RssLink(BaseModel):
    """A linked RSS (or leader) sequence flanking a coding allele."""

    uid: str
    subtype: str | None = None
    sequence: str
    relation: str
    # For a spliced leader (exon1+exon2): length of exon1, so the splice boundary
    # is known. None for RSS / when not recorded.
    exon1_len: int | None = None
    # Sources that assert this RSS/leader links to the allele (e.g. IMGT, MUSA).
    sources: list[str] = []


class RssAllele(BaseModel):
    """A coding allele linked to an RSS, with the source(s) of that link."""

    uid: str
    label: str | None = None
    species: str
    locus: str | None = None
    sources: list[str] = []


class RssHit(BaseModel):
    """One RSS sequence in a similarity result."""

    uid: str
    sequence: str
    species: str
    subtype: str | None = None
    score: float
    allele_count: int


class RssSimilar(BaseModel):
    """RSS sequences similar to a query RSS."""

    uid: str
    sequence: str
    count: int
    capped: bool = False
    hits: list[RssHit] = []


class AlleleDescription(BaseModel):
    """AIRR-aligned view of a coding allele (subset of AIRR AlleleDescription)."""

    label: str | None
    aliases: list[str] = []
    coding_sequence: str
    sequence_aa: str | None = None
    functional: bool | None = None
    locus: str | None = None
    sequence_type: str | None = None
    species: Ontology
    germline_set_ref: str | None = None


class NativeRecord(BaseModel):
    """Our native record shape."""

    uid: str
    species: str
    category: str
    segment: str | None = None
    locus: str | None = None
    # Second locus for shared TRA/TRD V genes (IMGT TRAV14/DV4): locus=TRA,
    # secondary_locus=TRD. NULL for ordinary single-locus records.
    secondary_locus: str | None = None
    subtype: str | None = None
    sequence: str
    sequence_hash: str
    sequence_aa: str | None = None
    aa_hash: str | None = None
    # Provenance of ``sequence_aa``: "sourced" (published by a source and consistent
    # with the nucleotides), "translated" (no source published one, so we translated
    # the nt), or "translated_over_conflict" (a source published one that the stored
    # nucleotides do not encode; ours is stored and theirs kept as an
    # ``aa/imgt_reference_aa`` annotation). NULL where there is no AA.
    aa_origin: str | None = None
    protein_uid: str | None = None
    functional: bool | None = None
    length: int
    has_ambiguity: bool
    # Official allele designations (IUIS/IMGT/OGRDB/VDJbase/...). >1 distinct
    # name = the same sequence is known under multiple allele names. (Raw, one per
    # source; ``allele_name_groups`` is the display view, collapsed by name.)
    allele_names: list[Alias] = []
    # Alternative names + identifiers (IMGT clone/former names, NCBI accession).
    other_names: list[Alias] = []
    # Display views: one row per distinct name, all its sources collapsed together.
    allele_name_groups: list[AliasGroup] = []
    other_name_groups: list[AliasGroup] = []
    # Former official IMGT designations for this exact sequence, recovered from
    # historical IMGT/GENE-DB release snapshots (each with the release it's from).
    former_names: list[Alias] = []
    # True when the allele exists only in historical releases **and** no other record
    # carries its name today — genuinely retired. A record IMGT merely re-cut (see
    # ``newer_version``) is not retired: its allele is still published, at a different
    # extent, so `build_record` leaves this false.
    retired: bool = False
    # Earlier records of this same allele, superseded by this one (oldest first).
    prev_versions: list[RecordVersion] = []
    # The record that supersedes this one: same allele name, later release range.
    newer_version: RecordVersion | None = None
    # Archived = present only in historical IMGT snapshots, gone from every current
    # source. Excluded from current coverage + default search; reachable directly.
    archived: bool = False
    # Set when this record's segment has a restricted similarity, and says how (see
    # queries.SIMILARITY_EQUAL_LENGTH). D is compared only against D alleles of the
    # same length, scored positionally. NULL = the ordinary edit-distance ranking.
    # A property of the segment, surfaced per record so a client can explain a short
    # or empty neighbour list rather than treating it as "nothing is similar".
    similarity_note: str | None = None
    # For a constant region we assembled from separately-published domains: the
    # ordered composition, e.g. "CH1:294+H:45+CH2:330+CH3-CHS:321". NULL when the
    # source published the sequence whole (all light chains, every V/D/J).
    composition: str | None = None
    # Release dates this archived allele was first / last seen in (IMGT snapshots).
    first_seen: str | None = None
    last_seen: str | None = None
    # External cross-references (NCBI Gene, HGNC, accessions) with links.
    xrefs: list[Xref] = []
    # Literature references (PubMed). Curator-provided ones are verified; the
    # literature layer auto-collects unverified leads (may contain wrong entries).
    references: list[Reference] = []
    # Linked RSS (heptamer/spacer/nonamer) sequences flanking this allele.
    rss: list[RssLink] = []
    # Linked leader sequences (spliced exon1+exon2 signal peptide).
    leaders: list[RssLink] = []
    sources: list[SourceTag] = []


class SequenceRecord(BaseModel):
    """A full record: native + AIRR views."""

    native: NativeRecord
    airr: AlleleDescription


class ResolveResult(BaseModel):
    matched_by: str = Field(description="uid | protein_uid | name | sequence_hash | rc_hash | aa_hash")
    record: SequenceRecord


class SearchHit(BaseModel):
    uid: str
    label: str | None = None
    species: str
    category: str | None = None  # coding / leader / rss
    locus: str | None = None
    sequence_type: str | None = None
    # Sequence length in nt. Surfaced so a hit that differs in length from the query
    # is visible — those are exactly the ones the conservation logo excludes.
    length: int | None = None
    score: float
    matched_by: str
    aliases: list[str] = []
    # UID of a later IMGT record for this same allele name, at a different extent.
    # Set so a superseded record never sits in a result list looking current.
    superseded_by: str | None = None


class ArchivedMatch(BaseModel):
    """A name that matched nothing current, but does match a record we still hold.

    Archived records are excluded from search on purpose — the current allele is what
    a name query wants. Returning plain silence is the wrong way to do it: the record
    is here, the query was not a mistake, and a paper citing that name deserves an
    answer. So a zero-result name search reports what it deliberately skipped.
    """

    uid: str
    name: str
    first_seen: str | None = None
    last_seen: str | None = None
    # Set when the name is not gone at all — a later record carries it at a
    # different extent, and that record is what the searcher actually wants.
    superseded_by: str | None = None
    superseded_by_name: str | None = None
    superseded_by_length: int | None = None
    superseded_relation: str | None = None


class SearchResults(BaseModel):
    query: str
    count: int
    capped: bool
    hits: list[SearchHit] = []
    # Only populated when ``count`` is 0: archived records the query does match.
    archived_matches: list[ArchivedMatch] = []


class PmidHit(BaseModel):
    """One allele attested by a paper, carrying **how** the paper attests it."""

    uid: str
    label: str | None = None
    species: str
    locus: str | None = None
    sequence_type: str | None = None
    # Evidence tier, same vocabulary and ranking as ``Reference.kind``:
    # "verified" (curator/GenBank xref), "sequence" (the paper's text or supplement
    # contains this exact sequence), "name" (the paper mentions the name — a lead,
    # and the same string denotes other sequences in other species).
    kind: str
    detail: str | None = None
    source: str | None = None
    # An archived allele is excluded from ordinary search but must appear here: a
    # paper citing it is precisely the reason the record is kept reachable.
    archived: bool = False
    superseded_by: str | None = None


class PmidResult(BaseModel):
    """Every allele one PubMed record attests, tiered by evidence.

    The reverse of the allele page's reference list. Tiers are kept separate and never
    collapsed into a flat count: a mined name-only lead is not the same claim as a
    curator-verified cross-reference, and reporting one total would present the weaker
    as the stronger.
    """

    pmid: str
    url: str
    title: str | None = None
    count: int
    capped: bool = False
    # Per-tier totals over the whole match set, not just the returned page.
    n_verified: int = 0
    n_sequence: int = 0
    # Name-only leads whose (species, locus) agrees with the paper's own evidence.
    # They are returned in ``hits`` only when the caller asked (``include_name``).
    n_name: int = 0
    # Name-only links dropped because their species/locus disagrees with the paper
    # (or the paper has no sequence/curator evidence to say what it is about).
    n_name_unanchored: int = 0
    include_name: bool = False
    hits: list[PmidHit] = []


# --- diff with consequence (INT-1) / IMGT alignment (SIM-3) ------------------

class NtSub(BaseModel):
    """One nucleotide substitution between a record and a neighbour, with its
    protein-level consequence. Positions are 1-based in the **query** sequence."""

    pos: int
    from_nt: str
    to_nt: str
    # IMGT unique-numbering codon + FR/CDR region label — V only, and only when the
    # IMGT-gapped reference was available to number against; None otherwise.
    imgt_codon: int | None = None
    region: str | None = None
    codon_from: str | None = None
    codon_to: str | None = None
    aa_from: str | None = None
    aa_to: str | None = None
    synonymous: bool | None = None


class NeighbourDiff(BaseModel):
    """The diff from a record to one nearby allele."""

    uid: str
    label: str | None = None
    identity: float
    same_length: bool
    length_delta: int  # len(neighbour) - len(query)
    n_subs: int
    n_nonsyn: int  # codons whose amino acid changed
    n_syn: int     # substitutions that left the amino acid unchanged
    subs: list[NtSub] = []
    note: str | None = None


class DiffSummary(BaseModel):
    """A record's exact differences to its nearest allele(s), with protein
    consequence (INT-1). ``region_labeled`` says whether FR/CDR labels were
    assigned (V, with the IMGT-gapped reference available)."""

    uid: str
    label: str | None = None
    segment: str | None = None
    region_labeled: bool = False
    region_source: str | None = None
    neighbours: list[NeighbourDiff] = []


class AlignmentRow(BaseModel):
    """One aligned allele in the IMGT-gapped view (SIM-3): the gapped nucleotide
    (primary) and amino-acid rows, each already padded to the view length."""

    uid: str
    label: str | None = None
    is_reference: bool = False
    gapped_nt: str
    gapped_aa: str


class AlignmentView(BaseModel):
    """IMGT-gapped alignment of the query against one hit, mismatches resolvable by
    region. Nucleotide is the primary view, amino acid a secondary toggle. V only."""

    uid: str
    label: str | None = None
    against_uid: str
    against_label: str | None = None
    nt_length: int
    aa_length: int
    # region label per position, aligned to the gapped nt / aa rows respectively
    nt_regions: list[str] = []
    aa_regions: list[str] = []
    rows: list[AlignmentRow] = []


class OrthologHint(BaseModel):
    """A cross-species nearest allele, flagged when the hit is reciprocal
    (each is the other's nearest in its species). A sequence-based hint only —
    NOT a curated orthology assignment."""

    species: str
    uid: str
    label: str | None = None
    identity: float
    reciprocal: bool = False
    back_uid: str | None = None
    back_label: str | None = None


class SequenceSearchRequest(BaseModel):
    sequence: str


class AaSearchRequest(BaseModel):
    sequence_aa: str
