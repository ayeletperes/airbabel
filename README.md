# allele_server

A standalone reference database for immunoglobulin **alleles**, **leaders**, and **RSS** across
species, with deterministic per-sequence identifiers ("iglabels"), alias resolution, and fast
similarity search by name or sequence. A simpler, focused sibling of `immunogenomic_db` (no genomic
or repertoire data); designed so `immunogenomic_db` and future databases can consume it.

Scope today: **V, D, J and C coding alleles + AA + aliases**, across 30 species, with **leaders**
(spliced signal peptide) and **RSS** populated. AA and the protein UID remain V-only by
design. Constant regions cover the immunoglobulin loci in their **secreted** form — heavy
chains assembled from the domains IMGT publishes separately, light chains taken whole; the
membrane-bound form and the TR constants are not loaded. D alleles are loaded, but because they
are very short their similarity is restricted to equal-length comparison (see below).

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
- A **web UI** plus a **REST API** (`/docs`); a transport-agnostic query core (`queries.py`) leaves
  room for an agent/MCP layer later (none in v1).

## Data sources

IMGT is downloaded directly from IMGT/GENE-DB (the bulk ReferenceSequences FASTA, all species);
OGRDB directly from its v2 API.
Internal curated sets (**HUSA** human, **MUSA** rhesus macaque) are layered on top. MUSA aggregates
several underlying rhesus databases (KIMDB, RhGLDB+, VRC, Guo), each credited as a first-class
source on the alleles it attests.

> **Redistribution:** IMGT data may not be redistributed. v1 ships with **no bulk/by-source export**;
> only per-record lookups and capped search. `export.py` keeps `EXPORT_ENABLED = False` and no route
> imports it. **This must be legally verified before flipping that flag or any public release.**

## Quickstart (dev)

```bash
cp .env.example .env
docker compose -f docker-compose.dev.yml up -d        # Postgres on 127.0.0.1:5433
uv sync                                                # core deps
uv run alembic upgrade head                            # create schema
uv sync --extra ingest                                 # ingestion deps
uv run allele-ingest --fetch                           # download OGRDB, then load
uv run allele-ingest --imgt-genedb                     # download + load IMGT (all species)
#   ... or load a pre-fetched cache, plus internal sets:
uv run allele-ingest --cache-dir ./ogrdb_cache \
    --husa /path/to/husa.fasta --musa /path/to/musa.fasta
# optional enrichment (network; cached; idempotent):
uv run allele-ingest --enrich-only --imgt-names        # IMGT clone/alt names + accessions + xrefs
uv run allele-ingest --enrich-only --ncbi-crawl        # NCBI accession + NCBI Gene per allele
uv run allele-server                                   # serve API + UI at http://127.0.0.1:8000
```

Run the tests (needs the dev Postgres up): `uv run --extra dev pytest`.

See `CONVENTIONS.md` for the design conventions and locked decisions.
