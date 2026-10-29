# Data model

Seven tables. In the maintainer repository, read this before `docs/ingestion.md`,
which covers how they are filled (that document is not part of the public subset).

The single idea the schema encodes: **a sequence is the identity; a name is a label
attached to it by a source.** Everything else follows from that.

```
                         ┌──────────┐
                         │ sources  │  who says so
                         └────┬─────┘
        ┌─────────────┬───────┼────────────┬──────────────┐
        │             │       │            │              │
   ┌────▼────┐  ┌─────▼─────┐ │  ┌─────────▼────────┐ ┌───▼──────────────┐
   │  names  │  │ sequence_ │ │  │    sequence_     │ │   sequence_      │
   │         │  │  source   │ │  │   annotations    │ │     links        │
   └────┬────┘  └─────┬─────┘ │  └─────────┬────────┘ └───┬──────────┬───┘
        │             │       │            │              │          │
        └─────────────┴───────┴────────────┴──────────────┘          │
                              │                                      │
                        ┌─────▼──────┐                               │
                        │ sequences  │◄──────────────────────────────┘
                        └─────┬──────┘   (coding ──has_leader──► leader)
                              │          (coding ──has_rss_*───► rss)
                     ┌────────▼──────────┐
                     │ sequence_neighbors│  derived cache (FK on uid)
                     └───────────────────┘
```

---

## `sequences` — the identity table

Every nucleotide sequence we know, of any kind.

| column | notes |
|---|---|
| `uid` | **unique.** Deterministic: `SPECIES-CATEGORYCODE+base32(sha256)`. Minted from the sequence, so it is stable across reloads. This is the key the site and API lead with. |
| `species`, `category`, `sequence_hash` | **the natural key** (`uq_sequences_species_category_hash`) |
| `sequence`, `length`, `has_ambiguity` | normalised sequence (uppercase, ungapped, IUPAC-validated) |
| `rc_hash` | hash of the reverse complement. A different sequence, tracked as a side attribute. |
| `sequence_aa`, `aa_hash` | V coding rows only |
| `protein_uid` | `SPECIES-VP-<base32[:6]>`, derived from `aa_hash`. **Not unique** — it groups the alleles whose nucleotides differ but whose V-REGION amino acids are identical. |
| `locus`, `segment`, `subtype` | `IGH`/`IGK`/`IGL`/`TRA`/`TRB`/`TRG`/`TRD`; `V`/`D`/`J`; see below |
| `functional` | as reported by the source; may be NULL |
| `archived` | present only in historical IMGT releases, in no current source. Excluded from coverage counts and default search, still reachable directly. Set by `mark_archived`. |

**Identity rule:** `(species, category, sequence_hash)`. `sequence_hash = sha256(normalize(seq))`
(`allele_server.sequence`). Two sources that report the same DNA corroborate one row;
they never duplicate it. This is why one row can carry several allele names.

**Controlled values**

| `category` | `subtype` | meaning |
|---|---|---|
| `coding` | *(null)* | a V/D/J allele |
| `coding` | `c_secreted` | a constant region, secreted form (IG only) |
| `leader` | `l_spliced` | signal peptide, `exon1 + exon2` spliced (the genomic form, with the intron, is not yet ingested) |
| `rss` | `v_rs` | the recombination signal 3' of a V |
| `rss` | `d_rs_5` | the signal 5' of a D |
| `rss` | `d_rs_3` | the signal 3' of a D |
| `rss` | `j_rs` | the signal 5' of a J |

Each RSS is one entity: heptamer + spacer + nonamer, stored complete, with the
parts in the order they lie beside the coding sequence — so a `j_rs` and a
`d_rs_5` read nonamer-first and a `v_rs` and a `d_rs_3` heptamer-first. Where
the components were ingested separately (the curation tables) the exact
composition is kept as an `rss/composition` annotation, because a heptamer is
not always 7 nor a nonamer always 9 and the boundaries are otherwise
unrecoverable. IMGT publishes each signal whole, so those rows carry no
composition.

**The subtype is load-bearing, not a label.** Spacer length follows the 12/23
rule and therefore the locus, so signals of different kinds come out the same
length — an IGKV `v_rs` and either D signal are all 28 nt. Similarity blocks and
the live ranker key on `coalesce(segment, subtype)` (`queries.comparison_class`)
precisely so those never share a pool: every RSS carries `segment IS NULL`, so
keying on segment alone put all of them together.

## `sources` — who says so

Natural key `(name, version)`. `kind` is `external` or `internal`.

`version` carries real meaning and varies by source: an IMGT release (`202624-7`), a
per-species marker (`GENE-DB:HUMAN`), or a DOI (KIARVA). **IMGT appears as many source
rows**, one per species and per historical release — grouped by `name` for display.
That is deliberate: `names` is unique on `(name, scheme, source_id, sequence_id)`, and
generic allele names like `IGHV1-1*01` recur across species, so a single shared IMGT
source row would let only the first species claim each name.

`redistributable` and `license` record the terms. IMGT is CC BY 4.0.

## `names` — labels, never keys

`(sequence_id, source_id, name, naming_scheme)` is the natural key. One sequence may
carry many names, from many sources, in many schemes (`imgt`, `iuis`, `ogrdb`,
`iglabel`, `vdjbase`, `kimdb`, …).

Two invariants worth knowing:

* **The internal sets contribute only their own IgLabel.** HUSA and MUSA must never
  appear as the *source* of an IMGT designation. Their ingest reads only the iglabel
  column (`curation_table.HUSA_NAMES` / `MUSA_NAMES`).
* **A name may be historical.** Names attached to a source whose `note` starts with
  `historical IMGT` are former designations, surfaced separately from current ones.

## `sequence_source` — provenance and accessions

What each source recorded for this sequence: `original_id` (a CURIE, e.g.
`GENBANK:X62106` or `IMGT/GENE-DB:IGHV1-2*02`) plus `original_sequence` as submitted.
Natural key `(sequence_id, source_id, original_id)`.

The GenBank accession here is the sequence's **primary public deposit**; it belongs to
the original submitters, not to any curator database.

## `sequence_annotations` — a typed key/value sidecar

Natural key `(sequence_id, source_id, domain, key)`. Values land in exactly one of
`value_text` / `value_float` / `value_int` / `value_bool`.

| `domain` | example keys | meaning |
|---|---|---|
| `functionality` | `functional` | source's functionality call |
| `imgt` | IMGT-specific fields | from the gene-table crawl |
| `leader` | `exon1_len` | where the splice boundary sits, so the two exons stay recoverable |
| `rss` | RSS metadata | |
| `history` | `first_seen`, `last_seen` | IMGT releases an archived allele appeared in |
| `xref` | `pmid:<id>`, `imgt_accession:<n>`, `ncbi_gene`, `hgnc` | verified cross-references |
| `literature` | `pmid:<id>`, `suppmatch:<pmid>:<uid>` | **auto-mined**, unverified. `suppmatch` records *how* it matched (`sequence_exact` / `sequence_similar` / `name`) and from where. |

> **Gotcha:** `loader.add_annotation` is idempotent by *skipping* an existing natural
> key. It never updates the value. To change a value you must delete the row first.

## `sequence_links` — typed edges between sequences

Natural key `(from_sequence_id, to_sequence_id, relation, source_id)`.
Relations in use: `has_leader`, and one per signal — `has_rss_vrs`,
`has_rss_drs5`, `has_rss_drs3`, `has_rss_jrs` (coding → leader / rss). A D
allele carries two RSS links, one per flank.

## `sequence_neighbors` — a derived cache, safe to drop and rebuild

Precomputed similarity edges, `kind` = `nt` or `aa`. Keyed on **`uid`**, not `id`, with
an FK to `sequences.uid`. Denormalised (`label`, `to_species`, `to_locus`, …) so a
"find similar" is one indexed read.

Two consequences:

* A `TRUNCATE ... CASCADE` of `sequences` **empties this table**. `reload_all` checks
  whether it was populated *before* truncating and rebuilds it if so, otherwise every
  similarity query would silently fall back to the slow live path.
* It is the largest table by far (millions of rows) and the slowest thing to rebuild.

Design detail in `docs/similarity-precompute.md`.

---

## Migrations

Every schema change is an Alembic migration; no ad-hoc DDL.

`0001_initial_schema.py` is a **squashed baseline of frozen literal DDL**, generated
once from the models and never re-read from them. `alembic upgrade head` therefore
builds this exact schema on an empty database, and `downgrade base` removes it.

> **Never use `Base.metadata.create_all()` in a migration.** The pre-squash `0001` did,
> so a fresh database already had `archived` and the newer unique key, and the later
> additive migrations then collided with themselves (`DuplicateColumn`). Nobody could
> create the database from migrations. A schema change is a new revision with its own
> explicit DDL.
