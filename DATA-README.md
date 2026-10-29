# Data provenance and redistribution terms

This repository aggregates immunoglobulin **allele / leader / RSS** reference data
from several sources. The service *aggregates and attributes* — it never adjudicates,
and every stored record carries its source. This file documents where the data comes
from and the attribution obligations when it is reused.

> **The CODE is licensed separately** (see `LICENSE`). This file governs the **DATA**.
> Code license and data terms are independent.

## Sources

**The authoritative, always-current list is [`/sources`](https://airbabel.org/sources)**, and
behind it one reviewed manifest per source at `src/allele_server/sources/<slug>.toml`
carrying that source's licence and the exact citation to use. The table below is a
summary; where the two disagree, the manifest is right. How a source gets added at all
is governed by `docs/source-policy.md`.

| Source | CURIE prefix | Content | Terms |
|--------|--------------|---------|-------|
| IMGT / IMGT/GENE-DB | `IMGT/GENE-DB:` | Official gene/allele names, reference sequences, gene tables | **CC BY 4.0 — redistributable with attribution** (2026 terms of use). Cite IMGT®. |
| OGRDB (AIRR) | `OGRDB:` | Curated germline sets | CC BY 4.0 (attribution). |
| GenBank / ENA / INSDC | `GENBANK:` `ENA:` | Nucleotide accessions referenced as cross-references | Public-domain records; we store accession IDs. |
| NCBI Gene / Taxonomy / PubMed | `NCBIGENE:` `NCBITAXON:` `PMID:` | Cross-references, species ids, literature leads | US Gov public data; attribute NCBI. |
| KIMDB | `KIMDB:` | Rhesus/cynomolgus macaque IGH alleles | No declared licence; cite Vázquez Bernat et al., *Immunity* 2021. |
| RhGLDB+ | `RHGLDB:` | Rhesus germline database (extended) | No declared licence; cite Ramesh et al., *Front Immunol* 2017. |
| VRC | `VRC:` | Rhesus macaque allele set | No declared licence; cite Kong et al., *Cell* 2019. |
| Guo et al. | `GUO:` | Rhesus macaque allele set | CC BY 4.0; cite Guo et al., *Front Immunol* 2024. |
| KIARVA | `DOI:` | Human IG gene variation atlas | No declared licence; cite Corcoran et al., *Immunity* 2026. |
| HUSA | `HUSA:` | Human novel-allele set (maintainer) | Maintainer-owned; manuscript in preparation. |
| MUSA | `MUSA:` | Rhesus curated set (aggregates KIMDB / RhGLDB+ / VRC / Guo) | No declared licence; cite Peres et al., *Immunity* 2026. Attribute the components individually as well. |

Full citation text for every row is on `/sources` and in the manifests.

## Redistribution: allowed, with attribution

As of the 2026 IMGT terms of use (https://www.imgt.org/about/termsofuse.php), **IMGT®
data is licensed CC BY 4.0** — free to share and adapt for any purpose, the only
condition being **attribution**. The other sources are drawn from **published
manuscripts** (or are maintainer-owned). The aggregated dataset is therefore
**redistributable with attribution**.

- **Cite IMGT®** — the international ImMunoGeneTics information system®,
  https://www.imgt.org — for any IMGT-derived content.
- **Attribute the originating source** of each record (the service already records
  and displays this per record; the export endpoints carry it in-band).

### Export is enabled

The full dataset is downloadable through the API (`allele_server.export.EXPORT_ENABLED
= True`): `GET /api/export` (manifest + attribution), `/api/export/alleles.fasta`,
`/api/export/alleles.json`. Every export carries source attribution, and abuse is
bounded by **rate-limiting the endpoints**, not by gating the data. The attribution
string is served with each export (`allele_server.export.ATTRIBUTION`).

## What ships in the public code repository

The public repository contains **code only** — no curated
data files, no ingestion caches, no downloaded FASTA. This is for maintainability and to
keep the ingestion/credentials private, **not** a licensing restriction: the data itself
is redistributable (above) and is served/exported by the running service.
