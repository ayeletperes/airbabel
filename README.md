# AIRBabel

A reference database for immunoglobulin and T-cell receptor **alleles**, **leaders** and
**RSS** across species, with deterministic per-sequence identifiers, alias resolution, and
similarity search by name or sequence. It runs at **[airbabel.org](https://airbabel.org)**;
this repository is its source (Python package `allele_server`).

Scope today: **V, D, J and C coding alleles + AA + aliases** for 40 species, with
**leaders** (spliced signal peptide) and **RSS** populated. AA and the protein UID remain
V-only by design. Constant regions cover the immunoglobulin loci in their **secreted** form:
heavy chains assembled from the domains IMGT publishes separately, light chains taken whole.
The membrane-bound form and the TR constants are not loaded. D alleles are loaded, but
because they are very short their similarity is restricted to equal-length comparison.

## What it does

- Mints a **deterministic UID** for every unique sequence: `SPECIES-CATEGORY+base32(sha256)`,
  reproducible across databases (coding UIDs are 8 chars, leader/RSS 5). V coding alleles also get
  a V-indicating **protein UID** (`HUMAN-VP-…`) derived from the AA hash.
- Stores a **network of tables** (`sequences`, `names`, `sequence_links`, `sequence_source`,
  `sequence_annotations`, `sources`) for easy pull and comparison.
- Carries the **amino-acid sequence** (`sequence_aa`) + `aa_hash` for V coding alleles, sourced from
  IMGT's V-AA or translated, so synonymous nt alleles group by protein.
- **Resolves** a name, UID, nt sequence, or AA sequence to a canonical record with all aliases and
  source tags.
- **Similarity search** by name (fuzzy), sequence (exact, sub/superstring, edit-distance), or AA.
- **AIRR-aligned**: every record carries both a `native` block and an `airr` `AlleleDescription`
  view; species is an `Ontology{NCBITAXON:…}` and sources/cross-refs are CURIEs
  (`GENBANK:X02897`, `OGRDB:…`). Use `?view=airr` for the AIRR block alone.
- **Enriched aliases**: optional passes crawl IMGT's gene tables for **name iterations**
  (clone / former / literature names like `VH26`, `DP-75`, plus reference accessions + HGNC /
  NCBI-Gene xrefs) and **NCBI E-utilities** for the GenBank accession + NCBI Gene id per allele.
- A **web UI**, a **REST API** (`/docs`) and an **MCP server** (`/mcp/`, with `--extra mcp`),
  all over one transport-agnostic query core (`queries.py`).

## Data sources

IMGT is downloaded directly from IMGT/GENE-DB (the bulk ReferenceSequences FASTA, all species);
OGRDB directly from its v2 API. Curated sets are layered on top: **HUSA** (human), **MUSA**
(rhesus macaque, aggregating KIMDB, RhGLDB+, VRC and Guo, each credited as a first-class source on
the alleles it attests) and **KIARVA**. Every stored record carries its source; the service
aggregates and never adjudicates.

> **Data terms:** the code is AGPL-3.0 (`LICENSE`); the data has its own terms, per source.
> IMGT and OGRDB are CC BY 4.0, so reuse requires attribution. See
> [`DATA-README.md`](DATA-README.md) and [airbabel.org/sources](https://airbabel.org/sources)
> for each source's licence and citation.

## Quickstart (dev)

```bash
cp .env.example .env                                   # set HUSA_TABLE / MUSA_TABLE / KI_XLSX to load the curated sets
docker compose -f docker-compose.dev.yml up -d        # Postgres on 127.0.0.1:5433
uv sync --extra ingest                                 # core + ingestion deps
uv run alembic upgrade head                            # create schema
uv run allele-ingest --fetch                           # download OGRDB into ./ogrdb_cache
uv run allele-ingest --reload                          # download IMGT, load every source, enrich
uv run allele-server                                   # serve API + UI at http://127.0.0.1:8000
```

A curated-set path left unset is skipped; one set to a missing file stops the reload.
Run the tests (needs the dev Postgres up): `uv run --extra dev pytest`.

For the design: [`docs/data-model.md`](docs/data-model.md) for the schema and how identity
is keyed, [`docs/naming-schemes.md`](docs/naming-schemes.md) for what each naming scheme is
and its citation, [`docs/similarity-precompute.md`](docs/similarity-precompute.md) for how
identity is scored and cached, [`docs/source-policy.md`](docs/source-policy.md) for what it
takes to add a source, and [`docs/ingestion.md`](docs/ingestion.md) for every ingest pass.
