# Naming schemes

What each allele naming scheme is, who defines it, what a name looks like under it, and
where to read the definition. `allele_server.naming` implements the parsing; this
document is the reference for what it is parsing.

**Names are labels, never keys.** A sequence's identity is
`(species, category, sha256(normalize(seq)))`. Every name we store is paired with the
scheme it follows and the source that attributes it, and no scheme is privileged as
canonical. That is why the UID leads on an allele page and the names sit in a table
below it: the same DNA carries different names in different schemes, and the same string
denotes different DNA in different species.

**Every citation below was resolved before it was written down** (checked 2026-08-04).
Each DOI or PMID was opened and confirmed to name the article claimed for it, with the
author list, journal, volume and pages read off the record rather than recalled. This
mirrors the rule the source manifests already enforce, where `status = "verified"` means
a human resolved the identifier — and it exists because a citation once written from
memory credited the wrong authors under a DOI that led to an unrelated paper. If you add
a scheme here, resolve its citation first.

---

## At a glance

| scheme | pattern | real example | defined by |
|---|---|---|---|
| [IUIS](#iuis) | `<locus><segment><subgroup>-<gene>*<allele>` | `IGHV1-69*01` | IUIS/IMGT nomenclature subcommittee |
| [IMGT](#iuis) | as IUIS, as published in IMGT/GENE-DB | `IGHJ10*01` | IMGT |
| [VDJbase](#vdjbase) | IUIS name + `_<mutations>` | `IGLV1-36*01_c330g` | Omer et al. 2020 |
| [IgDiscover](#igdiscover) | IUIS name + `_S<n>` | `IGHV1-18*01_S0898` | Corcoran et al. 2016 |
| [OGRDB](#ogrdb) | IUIS name with `*i<n>` inferred marker | `IGHV1-69*i04` | Lees et al. 2020 |
| [IgLabel](#iglabel) | `<locus><segment>0-<4-char>*00` | `IGHV0-2TS3*00` | Lees et al. 2023 |
| [HUSA](#husa) | IgLabel grammar, our label allocation | `IGHV0-22WP*00` | internal |
| [ASC](#asc) | `IGHVF<family>-G<cluster>*<allele>` | `IGHVF2-G15*02` | Peres et al. 2023 |
| [NCBI](#ncbi) | GenBank accession | `X02897` | NCBI |

The IUIS, IMGT, VDJbase, IgDiscover, OGRDB, IgLabel and HUSA examples are real names
taken from the database. The ASC one is from the paper (see that section), and the NCBI
one illustrates accession format rather than a stored name.

**Four of the eight are implemented but carry no stored names today**: `iuis`, `husa`,
`asc` and `ncbi` are all at zero. Incoming IMGT names are stored as `imgt` rather than
`iuis`, our internal labels as `iglabel` rather than `husa`, and accessions live as
cross-references rather than as names. The parsers exist and are tested; the data has
simply not needed them yet. Do not read the enum as a description of what is in the
database. Live counts are on `/overview`.

---

## IUIS

<a id="iuis"></a>

The canonical scheme, and what almost every published allele name follows. A name is
built from the locus (`IGH`/`IGK`/`IGL`, `TRA`/`TRB`/`TRG`/`TRD`), the segment
(`V`/`D`/`J`/`C`), the subgroup number, an optional gene designation, and the allele
number after `*`:

```
IGHV1-69*01
│││ │ │  └── allele, zero-padded to two digits
│││ │ └───── gene
│││ └─────── subgroup (IMGT "subgroup", AlleleName.family)
││└───────── segment
└─────────── locus
```

We zero-pad the allele to two digits on parse and keep the raw string verbatim, so
`IGHV1-69*1` and `IGHV1-69*01` bucket together without either being rewritten in place.

IMGT also uses a dual designation for the V genes rearranged in both the TRA and TRD
loci, e.g. `TRAV14/DV4`. The primary locus stays `TRA` and the second is recovered
separately, so one gene is not split into two records.

A IUIS name carries no `_` suffix and no `i` allele marker; those belong to the schemes
below, and `parse_name` rejects them under `Scheme.IUIS` rather than quietly accepting a
name from a different grammar.

**Nomenclature authority.** The IUIS Nomenclature Subcommittee for immunoglobulins and T
cell receptors validates and assigns names for newly described IG/TR germline genes,
following IMGT rules.

> Lefranc MP. WHO-IUIS Nomenclature Subcommittee for immunoglobulins and T cell receptors
> report. *Immunogenetics*. 2007;59(12):899-902.
> [doi:10.1007/s00251-007-0260-4](https://doi.org/10.1007/s00251-007-0260-4) ·
> PMID [18046549](https://pubmed.ncbi.nlm.nih.gov/18046549/)

Subcommittee: <https://iuis.org/committees/nom/nomenclature-sub-committees/t-cell-receptor-and-immunoglobulin-nomenclature-sc/>

**IMGT** is stored as its own scheme label rather than folded into `iuis`, because it
records that IMGT specifically published the name — which matters for the release
history, where a name's IMGT span is what makes a former name recoverable. IMGT requires
a specific attribution:

> IMGT®, the international ImMunoGeneTics information system®, <https://www.imgt.org>
> (founder and director: Marie-Paule Lefranc, Montpellier, France).

---

## VDJbase

<a id="vdjbase"></a>

An IUIS name with an underscore suffix naming the nucleotide differences from the base
allele, one `<from><position><to>` triple per substitution:

```
IGLV1-36*01_c330g                  one substitution
IGLV2-14*01_g168t_t198g_c337t      three
```

The suffix makes an allele that has been observed in repertoire data but not formally
named expressible without inventing a new allele number. The grammar is what identifies
the scheme, not the source: a name carrying `_` and no other marker is bucketed
`vdjbase` wherever it came from.

> Omer A, Shemesh O, Peres A, Polak P, Shepherd AJ, Watson CT, Boyd SD, Collins AM,
> Lees W, Yaari G. VDJbase: an adaptive immune receptor genotype and haplotype database.
> *Nucleic Acids Research*. 2020;48(D1):D1051-D1056.
> [doi:10.1093/nar/gkz872](https://doi.org/10.1093/nar/gkz872) ·
> PMID [31602484](https://pubmed.ncbi.nlm.nih.gov/31602484/)

Database: <https://www.vdjbase.org/>

---

## IgDiscover

<a id="igdiscover"></a>

An IUIS name with an `_S<n>` suffix, applied to a sequence inferred from an expressed
repertoire by the IgDiscover pipeline. `_S` is checked before the general `_` rule, so an
IgDiscover name is not mistaken for a VDJbase one:

```
IGHV1-18*01_S0898
```

> Corcoran MM, Phad GE, Vázquez Bernat N, Stahl-Hennig C, Sumida N, Persson MAA,
> Martin M, Karlsson Hedestam GB. Production of individualized V gene databases reveals
> high levels of immunoglobulin genetic diversity. *Nature Communications*. 2016;7:13642.
> [doi:10.1038/ncomms13642](https://doi.org/10.1038/ncomms13642) ·
> PMID [27995928](https://pubmed.ncbi.nlm.nih.gov/27995928/)

---

## OGRDB

<a id="ogrdb"></a>

The AIRR Community's Open Germline Receptor Database marks an **inferred** allele — one
established from repertoire data rather than from genomic sequence — with an `i` before
the allele number:

```
IGHV1-69*i04
```

The marker is part of the allele field, not a suffix, so it survives zero-padding: `*i4`
and `*i04` are the same allele. `parse_name` rejects the `i` marker under IUIS, IgLabel
and HUSA, because in those schemes an allele number asserts a formally assigned name.

> Lees W, Busse CE, Corcoran M, Ohlin M, Scheepers C, Matsen FA, et al. OGRDB: a
> reference database of inferred immune receptor genes. *Nucleic Acids Research*.
> 2020;48(D1):D964-D970.
> [doi:10.1093/nar/gkz822](https://doi.org/10.1093/nar/gkz822)

Database: <https://ogrdb.airr-community.org>

---

## IgLabel

<a id="iglabel"></a>

A temporary, collision-resistant label for a sequence that has no formally assigned name
yet. It keeps the IUIS grammar but puts `0` where the subgroup goes, a four-character
code where the gene goes, and `*00` for the allele — so it is recognisably *not* a
IUIS name while still parsing under the same rules:

```
IGHV0-2TS3*00
    │ │    └── allele fixed at 00
    │ └─────── 4-character label (20 random bits, base32)
    └───────── subgroup 0 marks the name as a label, not an assignment
```

The point is coordination: a group can name a newly found allele without waiting for a
formal assignment and without colliding with another group doing the same, and the label
can later be superseded by a IUIS name without the sequence changing identity.

> Lees WD, Christley S, Peres A, Kos JT, Corrie B, Ralph D, Breden F, Cowell LG,
> Yaari G, Corcoran M, Karlsson Hedestam GB, Ohlin M, Collins AM, Watson CT, Busse CE.
> AIRR community curation and standardised representation for immunoglobulin and T cell
> receptor germline sets. *Immunoinformatics (Amsterdam)*. 2023;10:100025.
> [doi:10.1016/j.immuno.2023.100025](https://doi.org/10.1016/j.immuno.2023.100025) ·
> PMID [37388275](https://pubmed.ncbi.nlm.nih.gov/37388275/)

---

## HUSA

<a id="husa"></a>

Our internal curated sets (HUSA for human, MUSA for rhesus macaque) allocate labels under
the IgLabel grammar, so a HUSA name is an IgLabel name whose allocation we did:

```
IGHV0-22WP*00        allocated by MUSA
```

It is a distinct scheme label only to record who allocated it. Parsing rules are
identical to IgLabel, including the refusal of `_` suffixes and the `i` marker.

In practice the distinction is not currently drawn in the data: names from HUSA and MUSA
are stored under the `iglabel` scheme with the source recording who allocated them, and
the `husa` scheme label has no rows. The source is doing the work the scheme label was
meant to do, which is arguably the right split — the grammar really is IgLabel's.

**Citation.** The scheme itself is IgLabel's (cited above); what follows is the set whose
labels these are. MUSA is published:

> Peres A, Upadhyay AA, Klein V, Saha S, Rodriguez OL, Vanwinkle ZM, et al.
> Population-level genomic analysis of immunoglobulin loci variation in rhesus macaques
> reveals extensive germline diversity. *Immunity*. 2026;59(1):213-228.e6.
> [doi:10.1016/j.immuni.2025.12.002](https://doi.org/10.1016/j.immuni.2025.12.002) ·
> PMID [41494536](https://pubmed.ncbi.nlm.nih.gov/41494536/)

Reference book: <https://vdjbase.org/reference_book/Rhesus_Macaque>

The HUSA manuscript is still in preparation, so the human set has no citation yet. In
both cases the aggregated sub-sources (KIMDB, RhGLDB+, VRC, Guo) carry their own
citations and must be credited individually — see `src/allele_server/sources/*.toml`.

---

## ASC

<a id="asc"></a>

An **allele similarity cluster** groups alleles that are close enough in germline
sequence that repertoire data cannot reliably tell them apart. Clusters are built by
hierarchical clustering of pairwise Levenshtein distance over the coding region, cut at
**95% nucleotide similarity**; a looser 75% cut defines the higher-level family. The
published name replaces chromosomal position with cluster membership, and its two
components are exactly those two cuts:

```
IGHVF2-G15*02            (the allele IUIS calls IGHV3-64*02)
│   │  │   └── allele number
│   │  └────── G<n>: the allele similarity cluster (the 95% cut)
│   └───────── F<n>: allele family (the 75% cut)
└───────────── locus + segment, as in IUIS
```

There is no positional component: `F` and `G` are the whole name. That is the point —
an allele's IUIS name encodes where it sits on the chromosome, and an ASC name encodes
only what it is similar to.

It answers a different question from every other scheme here. The rest name *an allele*;
an ASC name names *the group an allele belongs to*, because genotype inference is more
robust at cluster level than at allele level when two alleles differ by a base or two.
Names convert back to IUIS, so it is an analysis-time representation rather than a
competing nomenclature.

> Peres A, Lees WD, Rodriguez OL, Lee NY, Polak P, Hope R, Kedmi M, Collins AM, Ohlin M,
> Kleinstein SH, Watson CT, Yaari G. IGHV allele similarity clustering improves genotype
> inference from adaptive immune receptor repertoire sequencing data. *Nucleic Acids
> Research*. 2023;51(16):e86.
> [doi:10.1093/nar/gkad603](https://doi.org/10.1093/nar/gkad603) ·
> PMID [37548401](https://pubmed.ncbi.nlm.nih.gov/37548401/)

> [!WARNING]
> **`Scheme.ASC` in `naming.py` does not parse this.** It accepts a bare four-character
> alphanumeric code (`_HUSA_LABEL_RE`, `^[A-Z0-9]{4}$`) — the same pattern used for the
> IgLabel gene label — and the module docstring calls that label "ASC-style". A
> published ASC name is `IGHVF2-G15*02`, which that parser rejects.
>
> Nothing is currently mislabelled in the data, because no name is stored under the
> `asc` scheme. But the enum member claims a published scheme it does not implement, so
> either it should implement the grammar above or it should be renamed to describe the
> 4-character label it actually accepts. Writing this page is what surfaced it; it is
> recorded here rather than quietly fixed, because which way to resolve it is a
> semantic decision, not an editorial one.

---

## NCBI

<a id="ncbi"></a>

Not an allele name at all: a GenBank accession, which identifies a *record* rather than a
label for a sequence. It is stored alongside the allele names because it is the one
identifier that reliably links out to the original public deposit, and it is the only
alias kind rendered as a link on an allele page.

```
X02897        legacy format
NG_001019     RefSeq format, with underscore
```

> Sayers EW, Cavanaugh M, Frisse L, Pruitt KD, Schneider VA, Underwood BA, Yankie L,
> Karsch-Mizrachi I. GenBank 2025 update. *Nucleic Acids Research*. 2025;53(D1):D56-D61.
> [doi:10.1093/nar/gkae1114](https://doi.org/10.1093/nar/gkae1114) ·
> PMID [39558184](https://pubmed.ncbi.nlm.nih.gov/39558184/)

---

## Source-specific scheme labels

`naming.Scheme` has the eight entries above, but the `names.naming_scheme` column also
holds labels naming a **source's own grammar** where that grammar is not one of the
eight: `kimdb`, `rhgldb`, `vrc`, `guo`. These are not a ninth and tenth scheme in any
published sense — they record that a name follows its source's local conventions, so it
is never silently reported as IUIS.

Some of them do carry recognisable markers (`IGHV1-105*01_W6840` from Guo, or
`IGHV1-NL_2*01_S3477` from KIMDB), and `scheme_label` will bucket a name by grammar where
the grammar is unambiguous. Where it is not, the source label is the honest answer.

This is worth knowing before writing a query over `naming_scheme`: the set of values is
not closed by the `Scheme` enum, and a filter written against the enum alone will silently
miss names.

---

## Where this is used

- `allele_server.naming` — `parse_name`, `detect_scheme`, `scheme_label`, `canonicalise`.
  `detect_scheme` guesses from grammar alone and is deliberately ordered most-specific
  first (IgLabel, then `_S`, then `*i`, then `_`, then IUIS).
- `names.naming_scheme` — the stored label, one row per (name, scheme, source, sequence).
- The allele page's **Allele names** table shows every name with its scheme and the
  sources attributing it; `detail.html` maps the stored value to a display label.
- `docs/source-policy.md` governs the citation discipline this document follows.
