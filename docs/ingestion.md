# Ingestion

How data gets in. Read `docs/data-model.md` first.

Everything enters through `src/allele_server/ingest/`, and everything routes through
`loader.py`, which owns the idempotent upserts. Every pass is safe to re-run: a second
run of the same inputs changes zero rows.

Single entry point:

```bash
uv run allele-ingest [flags]        # src/allele_server/ingest/cli.py
```

---

## The two rules every pass obeys

1. **Hash-keyed upserts.** `loader.upsert_coding` / `upsert_noncoding` key on
   `(species, category, sequence_hash)`. Two sources reporting the same DNA corroborate
   one row. A source *never* creates a duplicate sequence.
2. **Fail loudly.** No silent fallbacks. A network outage raises rather than quietly
   recording "this species has no leader". (The literature passes are the deliberate
   exception: they are best-effort discovery, and a missing lead is not an error.)

---

## Sources

A source cannot be ingested until it has been reviewed and approved and has a manifest
carrying a licence and a citation — `loader.get_or_create_source` refuses anything else.
See **[docs/source-policy.md](source-policy.md)** for how a source is found, reviewed and
approved; this table is what is loaded today.

| source | kind | what it supplies | module |
|---|---|---|---|
| **IMGT/GENE-DB** | external | the multi-species backbone: V/D/J coding for 39 species across IG + TR; plus leader + RSS (**IG loci only**) | `imgt_genedb.py` |
| **IMGT (historical)** | external | *former* / retired allele names, recovered from release snapshots; drives `archived` | `imgt_history.py` |
| **IMGT (gene tables)** | external | alternative / clone names, accessions, HGNC + NCBI-Gene xrefs | `imgt_crawl.py` |
| **OGRDB** | external | curated, evidence-backed inferred germline sets | `ogrdb.py` |
| **HUSA** | internal | curated human set. Contributes **only its IgLabel**. | `curation_table.py` (+ `internal.py` for the FASTA path) |
| **MUSA** | internal | curated rhesus set, aggregating KIMDB / RhGLDB+ / VRC / Guo, each credited as a first-class source | `curation_table.py` |
| **KIARVA** | external | Karolinska Institutet Adaptive Immune Receptor Gene Variant Atlas: human IG variation (Corcoran et al. 2026) | `ki_xlsx.py` |
| **NCBI** | external | accession + Gene xrefs, and verified PubMed ids | `ncbi.py`, `_ncbi_client.py` |
| **Europe PMC / PMC-OA** | derived | literature evidence, **unverified** | `literature.py`, `literature_supp.py` |

## Passes

| pass | flag | network | disk cache | writes |
|---|---|---|---|---|
| OGRDB | *(default)* | via `--fetch` | `ogrdb_cache/` | sequences, names, attributions |
| all-species GENE-DB | `--imgt-genedb` | bulk FASTA once + GENElect | `imgt_genedb_ref.fasta`, `imgt_genelect_cache/` | coding, leaders, RSS |
| historical IMGT | `--fetch-history` to download | no (reads snapshots) | `imgt_history/` | former names, `history` annotations, `archived` |
| IMGT gene tables | `--imgt-names` | yes | `imgt_crawl_cache/` | alt names, xrefs |
| NCBI xrefs | `--ncbi-crawl` | yes | `ncbi_cache/` | `xref` annotations (incl. verified `pmid:`) |
| literature name probe | `--literature` | yes | `literature_cache/` | `literature` `pmid:` leads |
| full-text + supp mining | `--literature-supp` | yes | `literature_supp_cache/` | `literature` `pmid:` + `suppmatch:` |
| neighbour cache | `--neighbors` / `--sync-neighbors` | no | — | `sequence_neighbors` |

All caches are gitignored. Point them outside the checkout (see the `*_CACHE_DIR`
settings) if it sits in a synced folder — a full-corpus sweep is tens of gigabytes
and there is no reason to replicate it.

---

## The full reload

```bash
uv run allele-ingest --reload              # correct, complete rebuild
uv run allele-ingest --reload --neighbors  # ... and force the similarity rebuild
```

`reload.py::reload_all` TRUNCATEs and reloads everything. **The step order is a
correctness constraint, not a preference:**

1. `_genedb` — the **all-species** GENE-DB pass, the only IMGT input: IG + TR V/D/J/C
   coding from the bulk FASTA (with IMGT's published V protein), then leaders/RSS via
   GENElect. The IMGT release is recorded beside the FASTA (`<file>.release`) at download.
2. `ingest_cache` — OGRDB, for `species_keys` only. After IMGT, so a V allele both carry
   gets IMGT's published protein rather than OGRDB's translation fallback.
3. HUSA / MUSA / KIARVA curation tables.
4. `_enrich` — historical names, then **`mark_archived`**, then the IMGT-name and NCBI
   crawls. `mark_archived` flags sequences that **no current source** attests, so it
   must run after every current source above. Run it earlier and alleles that step 2
   would have supplied get marked archived and disappear from coverage and search.
5. AA backfill, then the neighbour cache.

`genedb` and `enrich` default to **on**: a bare `--reload` must produce the full
dataset. Opting out is explicit and deliberately partial (`--no-genedb`, `--no-enrich`).
`tests/test_reload.py` pins both invariants.

**Cost.** With the caches warm, a reload downloads nothing. The expensive parts are the
GENElect queries (now memoised, so only the first ever run pays) and the neighbour
rebuild (millions of rows).

**A reload truncates `sequence_annotations`**, so literature evidence is lost and must
be re-derived by re-running the two literature passes. Their disk caches make that
replay cheap — no re-crawling — but plan the order: reload, then enrich, then literature.

---

## Enrichment and evidence tiers

The literature layer is **discovery, not truth**. It writes into the `literature`
annotation domain, and the UI tiers what it finds:

* **verified** — from the curated tables / GenBank (`xref` domain).
* **sequence in paper** — the allele's exact sequence was found in the article body or a
  supplementary file (hash match). Cross-species safe: a paper about fish cannot contain
  a sequence that hashes to the human allele.
* **name mentioned** — the allele *name* appeared. A lead, not evidence: the same string
  denotes different sequences in different species. Rendered with an asterisk.

Mining is bounded on purpose: supplementary bundles over 64 MB and full text over 32 MB
are skipped (one real bundle was 323 MB of raw data), with an overall per-fetch deadline,
because `urlopen`'s timeout bounds socket *inactivity*, not total transfer, and a
trickling server once hung a crawl for nine hours.

---

## Gotchas a maintainer will hit

* **`add_annotation` never updates.** It is idempotent by *skipping* an existing
  `(sequence, source, domain, key)`. To change a value, delete the row first.
* **IMGT is many source rows.** One per species, one per historical release. Group by
  `sources.name` for display; do not assume a single IMGT row. See `_species_source`.
* **Internal sets are iglabel-only.** If HUSA or MUSA ever appears as the source of an
  IMGT name, that is stale data from an older ingest, not current behaviour.
* **The neighbour cache dies with a TRUNCATE** (FK on `sequences.uid`). `reload_all`
  detects this and rebuilds; nothing else does.
* **`detect_scheme("IGHV1-2*02")` returns `IUIS`, not `IMGT`.** Both internal ingest
  paths bucket IUIS-grammar names to `iglabel` deliberately.
