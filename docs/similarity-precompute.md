# Precomputed similarity (the neighbour cache)

How "find similar" is served. The expensive part — comparing a sequence against
every other sequence — is done **once, offline**, and stored, so a query at
request time is an indexed table read instead of a live all-vs-all computation.

## Why

Live similarity is ~1 s per call (a trigram candidate filter + `partial_ratio`
scoring over the candidates), CPU-bound, so it collapses to ~1 request/s under
concurrency. Precomputing each sequence's neighbours turns the request into a
~2.4 ms indexed lookup — ~440× faster per call — and the cost moves to a one-time
build.

## The distance

Identity is `rapidfuzz.partial_ratio` — an edit/Indel distance normalised over
the **best-matching overlap**, so an exact or fully-contained sequence scores
100% (boundary trims are not penalised) and internal substitutions/indels cost in
proportion. Computed per **nucleotide** and per **amino acid**, always **within a
single category** (coding / leader / rss). This is the same scorer the live path
uses, so cached and live results agree.

## Schema — `sequence_neighbors`

One row per directed edge `from_uid → to_uid`:

| column | meaning |
|---|---|
| `from_uid`, `to_uid` | endpoints, FK to `sequences.uid` (`ON DELETE CASCADE`) |
| `kind` | `nt` or `aa` |
| `identity` | 0..1 |
| `label`, `to_species`, `to_locus`, `to_segment`, `to_category` | **denormalised** attributes of the *to* sequence |

- Keyed on the **stable `uid`** (deterministic, content-derived), so edges
  survive id-renumbering on reload.
- `UNIQUE (from_uid, kind, to_uid)`; index `(from_uid, kind, identity)` for the
  top-k read.
- The five denormalised columns mean a read builds its result rows **with no
  join** back to `sequences` or `names` — it reads only this table.
- Created by migrations `0002` (table) and `0003` (denormalised columns).

## The build — `ingest/neighbors.py::build_neighbors`

Brute-force all-vs-all, per `(species, category, segment, kind)` block, using
`rapidfuzz.cdist` (vectorised, parallel C++):

1. Load the block's sequences; compute the full N×N identity matrix with `cdist`.
2. Per-block, compute the **denormalised attributes** once: a single batched
   query pulls every alias in the block and labels are built into a per-uid dict
   (so a block of N sequences costs one alias query, not N).
3. For each sequence, keep every neighbour with `identity >= NEIGHBOR_FLOOR`
   (**0.75**), sorted by identity, capped at `NEIGHBOR_K` (**500**) as a safety
   bound on dense families. Copy the denormalised attributes onto each edge.
4. Replace the block (per-block `DELETE` + bulk `INSERT`), so a rebuild of one
   block is idempotent and doesn't touch others.

**Segment is part of the block key.** A V is never scored against a J: they are
different features, not variants of one another, so an edge between them would be
a category error, not a weak hit. On today's data this changes nothing — real V/J
pairs top out near 0.60, well under the 0.75 floor, so no cross-segment edge was
ever written — but the builder no longer *depends* on that accident, which is what
would break first if D alleles (much shorter) or a lower floor arrived. It also
shrinks the O(n²) blocks. `segment IS NULL` (leaders, rss) is itself one block,
not a wildcard. `build_cross_species` pairs same-segment blocks only, for the same
reason: the ortholog of a V is a V. The cached read path filters on `to_segment`
as well, so a cache built before this cannot disagree with the live path.

**Floor, not top-k.** We store down to a fixed **0.75 identity floor** rather
than a flat top-k. That guarantees the cache is *complete* for any query
threshold ≥ 0.75 — every UI threshold (95/90/85/80%) is fully served, nothing in
that range falls back to live. (A flat top-200 would stop *above* 90% for dense
families, dropping those queries to the slow path.)

Full build: ~50 min, a few million edges. `NEIGHBOR_FLOOR` / `NEIGHBOR_K` are
defined once in `queries.py` and imported by the build, so the build and the read
can never drift apart.

## The read — `queries.py::_neighbors_from_table`

`find_similar` / `find_similar_aa` call it first; it returns the cached result or
`None` (→ caller computes live):

- **Gate on the build floor (0.75), not the per-uid lowest edge.** Because the
  cache holds *everything* ≥ 0.75, it's complete for any `min_identity ≥ 0.75`,
  even for a sparse allele whose only neighbour sits at 99.7% (there simply are no
  neighbours in between — so min=0.9 is served from cache, not live).
- **Capped-uid guard:** if a uid hit the 500 cap, the cache is only complete down
  to its k-th edge; a query below that falls to live.
- Reads the denormalised columns straight into `SearchHit` (no `sequences` /
  `names` join), `matched_by="cached"`.

### When it falls back to live (by design)

- `min_identity < 0.75` (below what was precomputed),
- **cross-species** search — the cache is built per-species; an "any species"
  query runs live,
- pasted sequences and reverse-complement (the query isn't a stored sequence),
- a uid with no stored edges (treated as "not built" → live, to be safe).

## Reload & incremental maintenance

- **Full rebuild** on demand via `allele-ingest --reload --neighbors`. The reload
  `TRUNCATE`s `sequences`, which **cascades** into `sequence_neighbors` (FK), so
  `reload_all` records whether the cache was non-empty **before** the truncate
  and, if so, **rebuilds it automatically** — a plain `--reload` can't silently
  leave an empty cache that degrades every search to the live path.
- **Incremental, two-sided update** for newly added alleles (no full rebuild) —
  `update_neighbors(session, new_uids)` / `sync_neighbors(session)`:
  - `O(M × N)` per touched `(species, category, segment, kind)` block, not `O(N²)`: each
    new sequence is scored against its block once, then written **both
    directions** — its own neighbour list (**outgoing**) *and* inserted into every
    existing neighbour's list (**incoming**, `ON CONFLICT DO NOTHING`, with cap
    trim).
  - `sync_neighbors` finds sequences with no cached neighbours (uids are
    content-derived, so a changed sequence is a new uid) and patches them; wired
    as `allele-ingest --sync-neighbors` after a non-reload ingest. Idempotent /
    re-runnable, and `partial_ratio` is symmetric so an incremental edge gets the
    exact identity a full rebuild would (verified: a delete-then-incremental
    reproduces the full-build edges bit-for-bit, both directions).
  - Handles **additions** only; removing an allele still needs a block rebuild.

## Performance

| | live | cached (denormalised) |
|---|---|---|
| per-call | ~1,064 ms | **~2.4 ms** (~440×) |
| concurrency | ~1 req/s | higher; now bounded by per-request Python serialisation, not similarity |

The cache removes the similarity CPU cost. The remaining concurrency ceiling is
single-process result serialisation — addressed by running multiple uvicorn
workers, not by more similarity work.

## Decisions and open items

- **Brute-force, not k-mer-prefiltered.** A k-mer prefilter was validated (it
  reproduces 99% of the brute-force top-k) but only cuts ~4–7% of comparisons on
  these sequences — V alleles share too much framework for k-mers to discriminate
  — so it isn't used. Brute-force is exact and fine at this scale; revisit a
  content-based prefilter only for a much larger or more heterogeneous pool.
- **Cross-species not cached** (per-species build) — the any-species default runs
  live. A single all-species build would cache it, at higher build cost.
- **Clustering as a build shortcut: evaluated, not adopted.** Greedy CD-HIT-style
  clustering only cut comparisons ~13.5× on these sequences (one big subgroup
  dominates the within-cluster `M²`) — a constant factor, not an asymptotic win —
  because V alleles share too much framework to separate cleanly. Incremental
  updates (above) are the real scaling lever; clustering would only be worth
  dropping in for an occasional full rebuild at much larger N.
- **RSS-similar** (`find_similar_rss`) is not wired to the cache; it stays live
  (rss is small).
