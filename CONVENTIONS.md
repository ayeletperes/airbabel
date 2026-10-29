# AIRBabel (`allele_server`): conventions

A standalone reference DB for immunoglobulin and T-cell receptor **alleles, leaders, RSS** across
species, served at airbabel.org. Reference data only (no genomic or repertoire data); built so
other databases can consume it for alias + sequence-similarity resolution.

## Locked decisions

| # | Decision | Value |
|---|----------|-------|
| 1 | Database | PostgreSQL 16 (dev: `docker-compose.dev.yml`, `127.0.0.1:5433`) |
| 2 | Stack | FastAPI (REST + Jinja2 UI), SQLAlchemy 2.x, Alembic, `uv`, Python 3.12 |
| 3 | UID | Deterministic: `SPECIES-CATEGORYCODE+base32(sha256)`; coding N=8, leader/RSS N=5 |
| 4 | Protein UID | V-specific, derived from `aa_hash`: `HUMAN-VP-<base32[:6]>` (non-unique; groups synonyms) |
| 5 | Ingestion | IMGT direct from IMGT/GENE-DB's bulk FASTA (`imgt_genedb`); OGRDB direct from its v2 API (`ingest/ogrdb.py`); stdlib HTTP, FASTA-first reader |
| 6 | Scope | **V + D + J + C coding alleles + AA + aliases**, 40+ species, IG + TR; leaders + RSS populated. AA/protein UID stay V-only. C is the **secreted** form (IG only; TR constants deferred). (Was: V-only, human/mouse/rhesus, leaders/RSS deferred.) |
| 7 | Internal sets | HUSA (human) / MUSA (**rhesus macaque** curated set); MUSA aggregates KIMDB / RhGLDB+ / VRC / Guo, each credited as a first-class source |
| 8 | Interop | AIRR field names + `AlleleDescription` view; species + sources as CURIEs |
| 9 | Export | Full-dataset API export (`export.py` `EXPORT_ENABLED=True`); IMGT is CC BY 4.0, attribution required |
| 10 | Agent | MCP server (`mcp/`, `/mcp/`, `--extra mcp`), a thin read-only wrapper over the transport-agnostic query core (`queries.py`) |
| 11 | Enrichment | Opt-in `--imgt-names` (crawl IMGT gene tables) + `--ncbi-crawl` (E-utilities); cached, idempotent |

## Identity rules (non-negotiable)

- **Canonical hash** = `sha256(normalize(seq))` (`allele_server.sequence`). Normalisation uppercases,
  strips gaps, validates IUPAC, rejects RNA, preserves ambiguity. Ported verbatim from `husa_core`
  so the two databases federate.
- A sequence's identity is `(species, category, sequence_hash)`. **Names are labels, never keys** —
  always paired with a scheme (`allele_server.naming`). What each scheme is, its grammar, and a
  verified citation for it: **[docs/naming-schemes.md](docs/naming-schemes.md)**.
- AA: `sequence_aa` + `aa_hash` on V coding rows; `aa_hash(nt)` and a sourced IMGT AA collapse to the
  same digest (`hash_aa` strips `*`/gaps).
- Reverse complement is a different sequence (`rc_hash` side attribute).

## Conventions

- **Idempotent writes.** Every ingest upsert is keyed on a natural key; re-running changes zero rows.
- **Fail loudly.** No silent fallbacks; `SequenceError`/`CurieError`/`SpeciesError` carry context.
- **Source attribution on every record.** The service aggregates; it never adjudicates.
- **Every schema change is an Alembic migration.** No ad-hoc DDL. Migration `0001` is a
  squashed baseline of **frozen literal DDL**; never use `create_all` in a migration (it
  reads live models, so a fresh DB gets tomorrow's columns and the next additive
  migration collides with itself).
- **No hardcoded paths/secrets.** Config flows through `allele_server.settings` / `.env` (git-ignored).
- **CURIEs for cross-references** (`allele_server.curie`): `GENBANK:`, `OGRDB:`, `IMGT/GENE-DB:`,
  `NCBITAXON:`, internal `HUSA:`/`MUSA:`.

## Layout

```
src/allele_server/
  sequence.py  uid.py  species.py  naming.py  curie.py   identity + standards primitives
  settings.py                                            env config
  db/ models.py session.py                               schema + engine
  ingest/ fasta.py loader.py ogrdb.py imgt_genedb.py internal.py cli.py
  ingest/ imgt_crawl.py ncbi.py _ncbi_client.py          enrichment passes (IMGT alt names, NCBI xref)
  ingest/ neighbors.py                                   precomputed similarity cache build (see docs/)
  queries.py                                             transport-agnostic core
  api/ app.py routes.py schemas.py web/templates/        REST + Jinja2 UI
  mcp/ server.py                                         MCP over queries.py (read-only)
  export.py                                              full-dataset export (enabled)
alembic/versions/0001_initial_schema.py
tests/
```

## Dev commands

```bash
docker compose -f docker-compose.dev.yml up -d   # Postgres
uv run alembic upgrade head                      # schema
uv run --extra dev pytest                        # tests (need the DB)
uv run allele-server                             # serve on 127.0.0.1:8000
uv sync --extra ingest                           # ingest deps (bs4, openpyxl, requests)
uv run allele-ingest --reload                    # download + load everything (needs --extra ingest + network)
```

## When extending

- **D / J alleles**: populated. No schema change — ingested with `segment` `D`/`J`. AA stays V-only.
  **D similarity is equal-length only, and positional**: at 8–148 nt (median 19) an edit
  distance can open gaps to flatter a score, so D alleles are compared only against D
  alleles of exactly the same length, scored with **Hamming** rather than Levenshtein.
  Different-length D alleles are not ranked against each other at all — separating those
  needs an MSA. The policy is `queries.SIMILARITY_EQUAL_LENGTH` (segment -> note) plus
  `similarity_scorer()`; both halves matter, since scoring an equal-length pool with
  Levenshtein would still admit gaps. Enforced on the **read** path, not only in
  `ingest/neighbors.py` — omitting/partitioning a block alone is not enough, because a
  cache miss falls through to the live ranker. The builder partitions D blocks by length
  (`_sub_blocks`) and `_comparable` guards the cross-species paths that score one matrix
  instead of partitioning, so the cache can never hold an edge a read would refuse.
- **C (constant) alleles**: populated for IG from IMGT (`ingest/imgt_constant.py`), no
  extra download — the domains are in the same bulk FASTA. IGH is **assembled** from
  IMGT's per-domain records (`CH1` + hinge + `CH2` + `CH3`… + `CHS`) in domain order;
  IGK/IGL/IGI are published whole as `C-REGION` and taken as-is. Stored as one
  `segment='C'`, `subtype='c_secreted'` sequence with the ordered composition kept as a
  `constant/domains` annotation, so the boundaries stay recoverable — the same shape as
  the leader decision. The rule is IMGT's own, pinned by the single assembled record
  IMGT publishes (`IGHG1*04`, 990 nt) which the assembly reproduces exactly. **Deferred**:
  the membrane form (`c_membrane`, different 3' exons — a separate transcript, not a
  longer sequence) and all TR constants (`EX4UTR` mixes CDS with UTR; TRG's
  `EX2A`/`B`/`C`/`R`/`T` duplicates need a rule IMGT does not supply). See the
  `TODO(c-membrane)` / `TODO(tr-constant)` markers.
- **Counting by segment**: `dataset_stats` is segment-driven (`segments`, `by_segment`,
  `seg_locus_totals`, `locus_all_totals`) and the overview renders one table per
  segment. Prefer those over `coding_v`/`coding_j`, which are legacy views: hardcoding
  a segment pair is what hid J, and then hid C.
- **RSS**: populated from the curation tables as the complete V-RS (`v_rs`), linked via
  `sequence_links` (`has_rss_vrs`); RSS network + similarity wired up.
- **Similarity / "find similar"**: served from a precomputed neighbour cache
  (`sequence_neighbors`, built by `ingest/neighbors.py`, read by `queries.find_similar*`) with a live
  fallback. Stores every neighbour ≥ 0.75 identity (`NEIGHBOR_FLOOR`/`NEIGHBOR_K` in `queries.py`),
  keyed on `uid`; rebuilt via `--reload --neighbors`. **See [docs/similarity-precompute.md](docs/similarity-precompute.md)** for the full design.
- **Leaders**: **decided — store the spliced leader** (`exon1 + exon2` = `l_part1 + l_part2`, the
  signal-peptide CDS) as one `leader` sequence (subtype `l_spliced`), linked `has_leader`, with the
  exon1 length kept as a `leader/exon1_len` annotation so the splice boundary is recoverable. The
  genomic form (`exon1 + intron + exon2`) is **deferred** — the intron is not in the curation tables;
  when a genomic source is wired in (e.g. IMGT L-PART genomic), add it as a second
  representation (subtype `l_genomic`). See `TODO(leader-intron)` in `ingest/curation_table.py`.
- **New species**: add one entry to `allele_server.species._SPECIES` (tag, scientific name, NCBITaxon).
- **New source**: it must be **reviewed and approved** first — see
  [docs/source-policy.md](docs/source-policy.md). Approval writes one manifest at
  `src/allele_server/sources/<slug>.toml` (licence + verifiable citation + CURIE
  prefix), which is the single source of truth; `loader.get_or_create_source` refuses a
  source without one. `curie._PREFIXES` now holds only standard namespaces.
- **Export is enabled** (`export.EXPORT_ENABLED=True`): IMGT is CC BY 4.0 and other sources are published; attribute IMGT + the sources on any reuse (`export.ATTRIBUTION`).
