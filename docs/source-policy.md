# Source policy

How a new data source is found, reviewed, approved and included. Read `docs/ingestion.md`
first — this document governs what may enter the pipeline that one describes.

**The rule in one line:** a source cannot be ingested until a human has reviewed it and
approved a manifest that carries a licence and a citation we can verify.

The service aggregates; it never adjudicates. That is only honest if every sequence we
serve can be traced back to whoever published it, under terms we actually read.

---

## The manifest is the single source of truth

One TOML file per source, at `src/allele_server/sources/<slug>.toml`. It is the only
place a source's identity, licence and citation are written down. Everything else is
derived from it:

| derived thing | how |
|---|---|
| the `sources` table row | projected by `ingest/loader.py::get_or_create_source` on every ingest |
| the CURIE prefix (`KIMDB:`, `OGRDB:` …) | registered into `curie._PREFIXES` at registry load |
| the `/sources` page | `queries.source_catalog` reads the manifest beside live counts |
| the export attribution line | `export.attribution_text()` builds it from the published manifests |

Before this policy those four lived in four hand-maintained places and had already
drifted (`curie` said `RHGLDB`, the metadata table said `RhGLDB+`, the attribution string
named neither). Editing a manifest now changes all of them at once.

Anything the manifest does not carry — how many alleles a source contributed, when it was
last loaded — stays in the database, because it is a fact about ingestion, not about the
source.

---

## Lifecycle

Every candidate source moves through these states. They are the literal values of
`source_candidates.status`, constrained in the schema so no code can invent a new one.

| state | set by | meaning |
|---|---|---|
| `discovered` | crawler | a channel returned it. Nothing has been judged. |
| `triaged` | crawler | passed the automatic filters: scored above the floor, not something we already hold, not suppressed by an earlier decision. |
| `dossier` | crawler | a report has been written. **This is the first state a human is asked to look at.** |
| `in_review` | human | a named curator has claimed it. |
| `approved` | **human only** | reviewed and accepted; a validated manifest exists. |
| `rejected` | **human only** | reviewed and declined. `suppressed = true`; it never resurfaces on any channel. |
| `deferred` | **human only** | not now. Returns to `dossier` after `deferred_until`. |
| `ingested` | `get_or_create_source` | an ingest pass actually created the `sources` row. Approval is a decision; this is delivery. |

**The automation may set any state up to `dossier`, and nothing beyond it.** The three
decision states are unreachable without a named person and a timestamp, enforced in the
schema itself:

```sql
CONSTRAINT ck_source_candidates_decided CHECK (
    status NOT IN ('approved','rejected','deferred')
    OR (decided_by IS NOT NULL AND decided_at IS NOT NULL))
```

That constraint *is* the human-review requirement. It is in the database rather than in
Python because a rule that only exists in the code that is supposed to obey it is a
convention, not a gate.

`published` is deliberately **not** a lifecycle state. It is a manifest field
(`[ingest] published`), because a source can be ingested but not published — embargoed
pre-publication data is the case — and that has to outlive the candidate row.

---

## 1. Discovery — how sources are found

> **Status: not built yet.** The `manual` channel works today — you record what you
> found and the rest of the workflow runs unchanged. The automated channels below are
> the next phase; everything downstream of them (candidates, dossiers, review,
> approval, the gate) is already in place and does not change when they land.

`allele-sources discover` will poll a set of channels, each a small cache-backed client
following the same rails as the other crawlers in `ingest/` (disk cache first, polite
rate limit, identifying User-Agent, injectable fetcher for offline tests, best-effort on
failure — a missed candidate is not an error).

| channel | what it watches | precision |
|---|---|---|
| `ogrdb` | new germline sets whose `(species, locus, version)` we do not hold | highest |
| `imgt` | new GENE-DB releases and newly covered species | highest |
| `europepmc` | papers announcing a germline set / allele database / reference directory | good |
| `zenodo`, `figshare`, `github` | DOI-bearing data deposits and repositories | mixed — but the **only** channels that hand us a machine-readable licence |
| `ncbi_bioproject`, `ncbi_sra` | germline-inference projects | lowest; tightest score floor |
| `manual` | a source you found yourself — **available now** | n/a |

Discovery never writes to `sequences`. It writes candidates, and nothing else.

### Not resurfacing what has already been judged

A crawler that re-proposes a source you rejected last month is worse than no crawler.
Two indexes are built per run and every raw hit is filtered through them:

- **already held** — every manifest's slug, name **and aliases**, every manifest DOI,
  every homepage host, every `sources.germline_set_ref` in the database. The alias list
  is what makes `RhGLDB+` and `RHGLDB` the same thing to the deduper.
- **suppressed** — anything `rejected`, `approved`, `ingested`, or `deferred` with a
  future `deferred_until`.

Cross-channel identity is the `dedupe_key`: the normalised DOI when one exists, else
`<channel>:<external_id>`. Two channels finding the same deposit both get a row — that
they agree is itself evidence — but they collapse into one dossier, and **a rejection on
either suppresses both**. This is what stops a rejected Zenodo deposit from reappearing
three months later as the paper announcing it.

Re-running a crawl over unchanged data changes zero rows, like every other pass.

---

## 2. Review — the dossier and the checklist

For each triaged candidate the automation writes a markdown dossier and records its path.
The dossier pre-fills the machine's answer to each checklist item; the human confirms or
corrects. The order matters — items 1–3 can end the review on their own.

1. **Provenance.** Who published this, where, when. Is it a primary germline set or a
   repackaging of one we already hold? Mirrors of IMGT and OGRDB are the common false
   positive.
2. **Licence.** The exact declared licence string, the URL it was read from, and the date
   it was read. If none is declared that is recorded as `none_declared` — **never guessed
   and never inferred from "it's public"**.
3. **Redistributability.** May these sequences appear in `/api/export`? Derived from (2),
   not assumed. An undeclared licence means `redistributable = false`: we may ingest and
   link to it, we may not republish it.
4. **Species and loci covered.** Mapped through `species.normalize_species`. An
   unrecognised species is flagged — it needs a `species._SPECIES` entry before the source
   can load (see `docs/maintaining.md`).
5. **Overlap with what we already hold.** How many sequences the deposit contains, how
   many are already present by `sequence.hash_sequence`, the novel remainder, and the
   nearest held neighbour for a sample of the novel ones. A source that is 100 % overlap
   is a corroboration rather than a new source — possibly still worth adding, but the
   reviewer should know that before deciding, not after.
6. **Sequence quality.** Do the sequences pass `sequence.normalize` (IUPAC, gaps
   stripped, RNA rejected)? Is the length distribution sane for the declared segment? Are
   there duplicates inside the deposit itself?
7. **Naming scheme.** What `naming.detect_scheme` reports for a sample of the names. A
   source asserting IMGT nomenclature it does not own is the failure `docs/data-model.md`
   warns about; such names are bucketed to the source's own scheme, never to `imgt`.
8. **Citation.** See below. This is the item that most often blocks an approval, and it
   should.
9. **Ingest cost.** Which module will load it, and whether it needs a new species, a new
   naming scheme, or a new parser.

Where the automation could not answer an item, the dossier says so explicitly. **An
absent number must never read as a zero** — "no sequence file was reachable" and "zero
overlap" are opposite conclusions.

---

## 3. Approval — the human step

Two commands, deliberately. The gap between them is where the review actually happens.

```bash
allele-sources draft   <channel>:<id> --manifest <slug>    # scaffold the .toml
# edit src/allele_server/sources/<slug>.toml — fill in the licence and the citation
allele-sources approve <channel>:<id> --manifest <slug> --by "Your Name" --note "..."
```

`draft` writes what the channel declared and leaves the citation empty. `approve`
validates the finished manifest and **refuses** it if the citation or licence rules below
are not met, then stamps an `[approval]` block into the TOML — approver, timestamp,
candidate id, dossier reference, decision note — and records the same decision on the
candidate row. The git diff of that block is the durable audit record; the dossier itself
is regenerable and gitignored, like every other cache.

The approver must be named, via `--by` or `$AIRBABEL_CURATOR`. There is no fallback to
`git config user.name`: an approval attributed to a guess is worse than no approval.

`reject` and `defer` take the same `--by` and require a note. A rejection is not a
failure of the crawler — it is the reviewer telling it something, and the note is what
makes the same candidate cheap to dismiss next time.

Approval happens at a terminal, not in the web UI. There is no authentication layer in
the app, and building one for a workflow with a single curator would be disproportionate.
If a browser view is wanted it should be read-only, with the decision still typed.

---

## 4. Citation and licence requirements

The policy is the one already stated in the code this replaces: *we only assert a citation
we can verify*. What changes is that "not yet filled in" stops being an indefinite state.

Every manifest must carry exactly one of:

- `citation.status = "verified"` with a non-empty `citation.text`, and a `doi` or `pmid`
  wherever one exists; **or**
- `citation.status = "none_available"` with a non-empty `citation.note` explaining why —
  for example, an internal set whose manuscript is still in preparation.

Licence takes the same shape: `licence.status` is `declared` or `none_declared`, and
`none_declared` requires a note and forces `redistributable = false`.

`none_available` is an honest answer, not a loophole. `allele-sources check` reports how
many sources are using it, so it stays visible rather than quietly becoming the default.
An invented citation is the one outcome this whole policy exists to prevent: it is worse
than a blank field, because a blank field is obviously missing and a plausible wrong
reference is not.

---

## 5. Inclusion — the gate

The gate is inside `loader.get_or_create_source`. Every ingest path already routes through
it, so there is no way to add a source that bypasses the review:

```
SourceError: no approved manifest for source 'NEWDB'. Every source must be reviewed and
approved before it can be ingested (see docs/source-policy.md).
```

The error names the remedy. Nothing silently degrades: an unapproved source is a loud
failure, in keeping with the fail-loudly rule.

The same call **projects** the manifest's licence, citation, DOI, homepage and
redistributable flag onto the `sources` row. This is the one documented exception to
"re-running an ingest changes zero rows": the first ingest after a manifest edit does
update those columns, on purpose. The alternative — a separate sync command someone has
to remember to run — is exactly the drift this design removes.

`settings.source_gate` (`strict` | `warn` | `off`) exists as an operator escape hatch for
a broken deployment. It is not a development convenience: `tests/test_sources_completeness.py`
always runs strict, so `warn` cannot quietly become the permanent state.

---

## What the tests hold in place

- `test_sources_manifest.py` — every manifest parses and validates; slugs, names, aliases
  and CURIE prefixes are unique across files; every declared species and ingest module
  actually exists.
- `test_sources_gate.py` — an unknown source raises `SourceError`; a known one lands with
  its licence and citation projected onto the row.
- `test_sources_completeness.py` — **the drift test.** Every source in the database has a
  manifest, a licence resolution and a citation resolution, and the export attribution
  line names every published redistributable source. This is the test that fails when
  someone adds a source the long way round.
- `test_source_approval.py` — approval refuses an unverifiable citation, records the
  approver, and the database itself rejects a hand-written `approved` row with no
  `decided_by`.

---

## Adding a source you already know about

The full loop, for a source that needs no discovery:

```bash
allele-sources add --channel manual --id newdb --title "NewDB germline set" --url https://…
allele-sources show manual:newdb                     # the dossier: read the checklist
allele-sources draft manual:newdb --manifest newdb   # scaffold src/allele_server/sources/newdb.toml
$EDITOR src/allele_server/sources/newdb.toml         # licence + citation — the actual work
allele-sources approve manual:newdb --manifest newdb --by "Your Name" --note "why"
# write the ingest module, wire it into reload.py in the right position, then:
uv run allele-ingest --reload
```

Retro-approving something already in the database (the twelve sources that predate this
policy) is the same flow with `--channel manual` and an `[approval] note` saying so.
