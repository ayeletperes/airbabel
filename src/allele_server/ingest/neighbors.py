"""Build the precomputed neighbour cache (``sequence_neighbors``).

All-vs-all within each (species, category, segment) via rapidfuzz ``cdist`` (vectorised
C, parallel), keeping each sequence's top-k. This is the slow similarity work
done **once, offline** so reads are an indexed lookup. Same-species only (the
default UI path); cross-species falls back to live.
"""

from __future__ import annotations

from collections import defaultdict

import numpy as np
from rapidfuzz import process
from rapidfuzz.distance import Levenshtein
from sqlalchemy import delete, func, insert, select, text
from sqlalchemy.dialects.postgresql import insert as pg_insert
from sqlalchemy.orm import Session

from ..db.models import Name, Sequence, SequenceNeighbor, Source
from .. import queries

# Rows per INSERT statement. 9 columns × this must stay under Postgres' 65535
# bind-parameter ceiling; 5000 (= 45k params) is a safe, fast batch size.
_INSERT_CHUNK = 5000


def _edge(from_id: int, to_id: int, kind: str, ident: float, m: dict) -> dict:
    """One neighbour-cache row: integer keys + the *to* sequence's denormalised
    attributes. ``m`` is the *to* sequence's :func:`_block_meta` entry (it carries
    that sequence's ``uid`` for the denormalised ``to_uid``)."""
    return {
        "from_id": from_id, "to_id": to_id, "kind": kind, "identity": ident,
        "to_uid": m.get("uid"), "label": m.get("label"),
        "to_species": m.get("to_species"), "to_locus": m.get("to_locus"),
        "to_segment": m.get("to_segment"), "to_category": m.get("to_category"),
    }


def _seg_cond(segment: str | None):
    """Predicate selecting one comparison-class block.

    The class is the segment for a coding allele and the subtype for a leader or an
    RSS (``queries.comparison_class``), because those carry no segment. Keying on
    the raw segment put every RSS in one block regardless of kind, which was
    invisible while ``v_rs`` was the only kind and wrong the moment D and J RSS
    landed: a 5'D-RS and an IGKV V-RS are both 28 nt.

    ``None`` remains a real block rather than "any", for rows that have neither.
    """
    return (
        queries.COMPARISON_CLASS.is_(None) if segment is None
        else queries.COMPARISON_CLASS == segment
    )


def _comparable(seg_a, len_a, seg_b, len_b) -> bool:
    """May these two sequences be scored against each other at all?

    Same segment always (an ortholog of a V is a V), and additionally the same length
    when that segment is only compared at equal length (D). Used where the builder
    scores one big matrix instead of partitioning the input first — the two must agree,
    or an incremental update would write edges a full rebuild refuses to.
    """
    if seg_a != seg_b:
        return False
    return not queries.equal_length_only(seg_a) or len_a == len_b


def _sub_blocks(segment: str | None, rows):
    """Split one (species, category, segment) block into the units actually compared.

    Normally that is the whole block. For a segment in
    ``queries.SIMILARITY_EQUAL_LENGTH`` (D) it is one sub-block **per sequence
    length**, because those segments are only ever compared at identical length —
    scoring across lengths and filtering afterwards would still let an edit distance
    open gaps, which over a 19 nt sequence is the metric inventing a result.

    Splitting also makes the build cheaper, not more expensive: O(n²) over one block
    becomes O(Σ nᵢ²) over the length classes, and D's length classes are small.

    Yields ``(uids, strings)`` pairs. The sub-blocks partition the block, so the
    per-sub-block delete-then-insert stays idempotent.
    """
    if not queries.equal_length_only(segment):
        yield [r[0] for r in rows], [r[1] for r in rows]
        return
    by_len: dict[int, tuple[list, list]] = {}
    for uid, s in rows:
        uids, strs = by_len.setdefault(len(s), ([], []))
        uids.append(uid)
        strs.append(s)
    for length in sorted(by_len):
        yield by_len[length]


def _segments(session: Session, category: str, sp: str | None = None) -> list[str | None]:
    """Distinct comparison classes present in a category (optionally one species).

    Blocks are keyed on the class as well as (category, species) because comparing a
    V to a J is meaningless: they are different loci-spanning features, not variants
    of one another. Nothing cross-segment clears the 0.75 floor today — V vs J tops
    out around 0.60 — so this wrote no different a cache for coding, but it stops the
    builder relying on that accident, and it shrinks the O(n^2) blocks.

    For the non-coding categories the class is the subtype, and there the separation
    is doing real work rather than guarding an accident: RSS kinds of equal length
    (a 5'D-RS and an IGKV V-RS are both 28 nt) score against each other freely.
    """
    where = [Sequence.category == category]
    if sp is not None:
        where.append(Sequence.species == sp)
    return [
        r[0] for r in session.execute(
            select(queries.COMPARISON_CLASS).where(*where).distinct()
        )
    ]


def _block_meta(
    session: Session, category: str, sp: str | None = None,
    segment: str | None = None, *, by_segment: bool = False,
) -> dict[str, dict]:
    """Per-uid denormalised attributes for a category block (one species, or all
    species when ``sp`` is None).

    Labels are computed in a single batched alias query (all names for the block
    at once) rather than one ``_aliases`` call per sequence, so a block of N
    sequences costs one query, not N.

    ``by_segment`` narrows the block to one ``segment`` (``None`` selects the rows
    whose segment IS NULL — leaders and rss — not "any segment").
    """
    where = [Sequence.category == category]
    if sp is not None:
        where.append(Sequence.species == sp)
    if by_segment:
        where.append(_seg_cond(segment))
    seq_rows = session.execute(
        select(
            Sequence.id, Sequence.uid, Sequence.species, Sequence.locus,
            Sequence.segment, Sequence.category,
        ).where(*where)
    ).all()
    id_to_uid = {r.id: r.uid for r in seq_rows}
    # One query for every alias in the block; group by sequence id.
    al_by: dict[int, list] = {}
    for sid, name, scheme, src in session.execute(
        select(Name.sequence_id, Name.name, Name.naming_scheme, Source.name)
        .join(Source, Source.id == Name.source_id)
        .where(Name.sequence_id.in_(list(id_to_uid)))
    ).all():
        al_by.setdefault(sid, []).append(
            queries.Alias(name=name, scheme=scheme, source=src)
        )
    return {
        r.uid: {
            "id": r.id, "uid": r.uid,
            "label": queries._label(al_by.get(r.id, [])),
            "to_species": r.species, "to_locus": r.locus,
            "to_segment": r.segment, "to_category": r.category,
        }
        for r in seq_rows
    }


def build_neighbors(
    session: Session, *,
    floor: float = queries.NEIGHBOR_FLOOR, k: int = queries.NEIGHBOR_K,
    species: str | None = None,
    categories=("coding", "rss", "leader"), kinds=("nt", "aa"),
) -> int:
    """(Re)build the neighbour cache. Returns the number of edges written.

    Stores every neighbour at/above ``floor`` identity (default 0.75), capped at
    ``k`` per sequence as a safety bound on dense families. Storing to a fixed
    floor — not a flat top-k — means any query whose
    threshold is ``>= floor`` is fully served from the cache: the read's
    floor-check passes for every UI threshold (95/90/85/80%), so nothing in that
    range drops to the slow live path.
    """
    thr = floor  # Levenshtein.normalized_similarity is 0..1
    total = 0
    for category in categories:
        for kind in kinds:
            if kind == "aa" and category != "coding":
                continue  # AA only on coding
            col = Sequence.sequence_aa if kind == "aa" else Sequence.sequence
            sp_list = (
                [species] if species else
                [r[0] for r in session.execute(
                    select(Sequence.species).where(Sequence.category == category).distinct()
                )]
            )
            # One block per (species, segment): a V is never scored against a J.
            for sp, seg in [
                (sp, seg) for sp in sp_list for seg in _segments(session, category, sp)
            ]:
                rows = session.execute(
                    select(Sequence.uid, col).where(
                        Sequence.category == category, Sequence.species == sp,
                        _seg_cond(seg), col.is_not(None),
                    )
                ).all()
                if len(rows) < 2:
                    continue
                # Denormalised attributes of the *to* sequence, computed once per
                # block (single batched alias query) and copied onto each edge.
                meta = _block_meta(session, category, sp, seg, by_segment=True)
                scorer = queries.similarity_scorer(seg)
                for uids, strs in _sub_blocks(seg, rows):
                    if len(uids) < 2:
                        continue
                    mat = process.cdist(
                        strs, strs, scorer=scorer, workers=-1,
                        dtype=np.float32, score_cutoff=thr,
                    )
                    id_of = {u: meta[u]["id"] for u in uids}
                    session.execute(
                        delete(SequenceNeighbor).where(
                            SequenceNeighbor.kind == kind,
                            SequenceNeighbor.from_id.in_(list(id_of.values())),
                        )
                    )
                    n = len(uids)
                    edges: list[dict] = []
                    for i in range(n):
                        row = mat[i]
                        # every neighbour at/above the floor, identity desc, capped at k
                        cand = np.where(row >= thr)[0]
                        cand = cand[np.argsort(-row[cand])]
                        cnt = 0
                        for j in cand:
                            j = int(j)
                            if j == i:
                                continue
                            edges.append(_edge(
                                id_of[uids[i]], id_of[uids[j]], kind, float(row[j]),
                                meta.get(uids[j], {}),
                            ))
                            cnt += 1
                            if cnt >= k:
                                break
                    if edges:
                        session.execute(insert(SequenceNeighbor), edges)
                        total += len(edges)
            session.commit()
    return total


def refresh_neighbor_labels(session: Session) -> int:
    """Repair the denormalised ``label`` on cached edges whose *to* allele had no
    (or a different) name when the edge was built — e.g. a name attached after the
    fact. Recomputes the label from current names for every ``to_uid`` with a NULL
    cached label and updates those edges in one set-based statement.

    (Leaders/RSS legitimately have no allele name; their NULL labels are left as is.)
    """
    uids = [u for (u,) in session.execute(
        select(SequenceNeighbor.to_uid).where(SequenceNeighbor.label.is_(None)).distinct()
    )]
    if not uids:
        return 0
    id_uid = dict(session.execute(
        select(Sequence.id, Sequence.uid).where(Sequence.uid.in_(uids))
    ).all())
    al_by: dict[int, list] = {}
    for sid, name, scheme, src in session.execute(
        select(Name.sequence_id, Name.name, Name.naming_scheme, Source.name)
        .join(Source, Source.id == Name.source_id)
        .where(Name.sequence_id.in_(list(id_uid)))
    ):
        al_by.setdefault(sid, []).append(queries.Alias(name=name, scheme=scheme, source=src))
    pairs = []
    for sid, uid in id_uid.items():
        label = queries._label(al_by.get(sid, []))
        if label:
            pairs.append({"uid": uid, "label": label})
    if not pairs:
        return 0
    session.execute(text(
        "CREATE TEMP TABLE _nl (uid varchar(32) PRIMARY KEY, label text) ON COMMIT DROP"
    ))
    session.execute(text("INSERT INTO _nl (uid, label) VALUES (:uid, :label)"), pairs)
    r = session.execute(text(
        "UPDATE sequence_neighbors sn SET label = n.label FROM _nl n "
        "WHERE sn.to_uid = n.uid AND sn.label IS NULL"
    ))
    session.commit()
    return r.rowcount


def sync_neighbors(session: Session, *, kinds=("nt", "aa")) -> int:
    """Incrementally patch the cache for sequences that have no cached neighbours
    yet — i.e. newly ingested alleles (uids are content-derived, so a changed
    sequence is a new uid). Patches BOTH same-species (``update_neighbors``) and
    cross-species (``update_cross_neighbors``) edges. Cheap to run after an
    idempotent re-ingest; a no-op when the cache is current. Genuinely isolated
    alleles (no neighbour >= floor) stay uncovered and are re-checked each run — a
    small, rare cost.
    """
    all_uids = {u for (u,) in session.execute(select(Sequence.uid))}
    covered = {u for (u,) in session.execute(
        select(Sequence.uid)
        .join(SequenceNeighbor, SequenceNeighbor.from_id == Sequence.id)
        .distinct()
    )}
    uncovered = list(all_uids - covered)
    if not uncovered:
        return 0
    n = update_neighbors(session, uncovered, kinds=kinds)
    n += update_cross_neighbors(session, uncovered, kinds=kinds)
    return n


def update_neighbors(
    session: Session, new_uids,
    *, floor: float = queries.NEIGHBOR_FLOOR, k: int = queries.NEIGHBOR_K,
    kinds=("nt", "aa"),
) -> int:
    """Two-sided incremental update for newly added sequences (additions only).

    Instead of an O(N^2) full rebuild, this is O(M*N) per touched block: each new
    sequence is scored against the whole (species, category, segment) block once,
    then we write BOTH directions —

      * **outgoing**: the new sequence's own neighbour list (every existing
        sequence at/above ``floor``, capped at ``k``);
      * **incoming**: the new sequence is inserted into every existing neighbour's
        list (``ON CONFLICT DO NOTHING``), and any list that overflows ``k`` is
        trimmed back to its top-k.

    Re-runnable and idempotent (it first clears the new uids' outgoing edges).
    Handles additions only — removing an allele still needs a rebuild of the
    affected block. ``partial_ratio`` is symmetric, so an incremental edge gets
    the same identity a full rebuild would.
    """
    new_uids = set(new_uids)
    if not new_uids:
        return 0
    thr = floor  # Levenshtein.normalized_similarity is 0..1
    # Block key carries a length class: None normally, the sequence length for a
    # segment that is only compared at equal length (D). That keeps the incremental
    # path partitioned exactly as the full rebuild partitions it -- otherwise a
    # re-ingest would quietly write the cross-length edges a rebuild refuses to.
    blocks: dict[tuple[str, str, str | None, int | None], set[str]] = {}
    for uid, sp, cat, seg, length in session.execute(
        select(Sequence.uid, Sequence.species, Sequence.category,
               queries.COMPARISON_CLASS,
               Sequence.length)
        .where(Sequence.uid.in_(new_uids))
    ):
        klass = length if queries.equal_length_only(seg) else None
        blocks.setdefault((sp, cat, seg, klass), set()).add(uid)

    total = 0
    for (sp, cat, seg, klass), block_new in blocks.items():
        meta = _block_meta(session, cat, sp, seg, by_segment=True)
        id_of = {u: meta[u]["id"] for u in meta}
        for kind in kinds:
            if kind == "aa" and cat != "coding":
                continue
            col = Sequence.sequence_aa if kind == "aa" else Sequence.sequence
            where = [Sequence.category == cat, Sequence.species == sp,
                     _seg_cond(seg), col.is_not(None)]
            if klass is not None:
                where.append(func.length(col) == klass)
            allrows = session.execute(select(Sequence.uid, col).where(*where)).all()
            all_uids = [r[0] for r in allrows]
            all_strs = [r[1] for r in allrows]
            new_pos = [i for i, u in enumerate(all_uids) if u in block_new]
            if not new_pos or len(all_uids) < 2:
                continue
            mat = process.cdist(
                [all_strs[i] for i in new_pos], all_strs,
                scorer=queries.similarity_scorer(seg), workers=-1, dtype=np.float32,
            )
            # idempotent: drop the new uids' existing outgoing edges first
            session.execute(
                delete(SequenceNeighbor).where(
                    SequenceNeighbor.kind == kind,
                    SequenceNeighbor.from_id.in_([id_of[all_uids[i]] for i in new_pos]),
                )
            )
            out_edges: list[dict] = []
            in_edges: list[dict] = []
            affected: set[str] = set()
            for r, ai in enumerate(new_pos):
                nu = all_uids[ai]
                row = mat[r]
                mnu = meta.get(nu, {})
                cand = np.where(row >= thr)[0]
                cand = cand[np.argsort(-row[cand])]
                cnt = 0
                for j in cand:
                    j = int(j)
                    ju = all_uids[j]
                    if ju == nu:
                        continue
                    ident = float(row[j])
                    if cnt < k:  # the new sequence's own list (capped)
                        out_edges.append(
                            _edge(id_of[nu], id_of[ju], kind, ident, meta.get(ju, {}))
                        )
                        cnt += 1
                    if ju not in block_new:  # existing list gains the new sequence
                        in_edges.append(_edge(id_of[ju], id_of[nu], kind, ident, mnu))
                        affected.add(ju)
            # Chunk inserts: a single statement is bounded by Postgres' 65535-param
            # limit (10 cols → ~6k rows), and an incremental batch can far exceed it.
            for c in range(0, len(out_edges), _INSERT_CHUNK):
                session.execute(insert(SequenceNeighbor), out_edges[c:c + _INSERT_CHUNK])
            total += len(out_edges)
            for c in range(0, len(in_edges), _INSERT_CHUNK):
                session.execute(
                    pg_insert(SequenceNeighbor)
                    .on_conflict_do_nothing(index_elements=["from_id", "kind", "to_id"]),
                    in_edges[c:c + _INSERT_CHUNK],
                )
            total += len(in_edges)
            # trim only the existing lists that overflowed the cap (rare)
            if affected:
                affected_ids = [id_of[u] for u in affected]
                over = session.execute(
                    select(SequenceNeighbor.from_id)
                    .where(
                        SequenceNeighbor.kind == kind,
                        SequenceNeighbor.from_id.in_(affected_ids),
                    )
                    .group_by(SequenceNeighbor.from_id)
                    .having(func.count() > k)
                ).all()
                for (fid,) in over:
                    keep = list(session.scalars(
                        select(SequenceNeighbor.to_id)
                        .where(
                            SequenceNeighbor.from_id == fid,
                            SequenceNeighbor.kind == kind,
                        )
                        .order_by(SequenceNeighbor.identity.desc())
                        .limit(k)
                    ))
                    session.execute(
                        delete(SequenceNeighbor).where(
                            SequenceNeighbor.from_id == fid,
                            SequenceNeighbor.kind == kind,
                            SequenceNeighbor.to_id.notin_(keep),
                        )
                    )
        session.commit()
    return total


def _trim_cross(session: Session, affected: set[str], kind: str, k: int) -> None:
    """Trim each affected allele's CROSS-species edges back to top-k by identity,
    leaving its same-species edges untouched (separate cap — cross never evicts
    same-species, matching ``build_cross_species``)."""
    if not affected:
        return
    rows = session.execute(
        select(Sequence.id, Sequence.species).where(Sequence.uid.in_(affected))
    ).all()
    for fid, jsp in rows:
        cond = (
            SequenceNeighbor.from_id == fid, SequenceNeighbor.kind == kind,
            SequenceNeighbor.to_species != jsp,
        )
        cnt = session.scalar(
            select(func.count()).select_from(SequenceNeighbor).where(*cond)
        )
        if not cnt or cnt <= k:
            continue
        keep = list(session.scalars(
            select(SequenceNeighbor.to_id).where(*cond)
            .order_by(SequenceNeighbor.identity.desc()).limit(k)
        ))
        session.execute(
            delete(SequenceNeighbor).where(*cond, SequenceNeighbor.to_id.notin_(keep))
        )


def update_cross_neighbors(
    session: Session, new_uids,
    *, floor: float = queries.CROSS_FLOOR, k: int = queries.NEIGHBOR_K,
    kinds=("nt", "aa"),
) -> int:
    """Two-sided incremental CROSS-species update for newly added coding sequences.

    Each new coding allele is scored against every coding sequence; its
    different-species neighbours at/above ``CROSS_FLOOR`` (0.85) become edges in
    BOTH directions. Cross edges are capped separately (``_trim_cross``) so they
    never evict same-species edges. Coding only (cross-species orthology is a
    coding concern). Idempotent: clears the new uids' existing cross edges first.
    """
    new_uids = set(new_uids)
    if not new_uids:
        return 0
    thr = floor
    total = 0
    for kind in kinds:
        col = Sequence.sequence_aa if kind == "aa" else Sequence.sequence
        rows = session.execute(
            select(Sequence.uid, col, Sequence.species, queries.COMPARISON_CLASS).where(
                Sequence.category == "coding", col.is_not(None)
            )
        ).all()
        uids = [r[0] for r in rows]
        strs = [r[1] for r in rows]
        sps = [r[2] for r in rows]
        segs = [r[3] for r in rows]
        lens = [len(s) for s in strs]
        new_pos = [i for i, u in enumerate(uids) if u in new_uids]
        if not new_pos or len(set(sps)) < 2:
            continue
        meta = _block_meta(session, "coding")  # all species
        id_of = {u: meta[u]["id"] for u in meta}
        # Unlike the per-species builder this path scores one matrix over every coding
        # row rather than per block, so comparability is enforced when picking
        # candidates (``_comparable`` below) instead of by partitioning the input.
        mat = process.cdist(
            [strs[i] for i in new_pos], strs,
            scorer=Levenshtein.normalized_similarity, workers=-1,
            dtype=np.float32, score_cutoff=thr,
        )
        # idempotent: drop these new uids' existing CROSS edges (keep same-species)
        by_sp: dict[str, list[str]] = {}
        for i in new_pos:
            by_sp.setdefault(sps[i], []).append(uids[i])
        for sp, ul in by_sp.items():
            session.execute(delete(SequenceNeighbor).where(
                SequenceNeighbor.kind == kind,
                SequenceNeighbor.from_id.in_([id_of[u] for u in ul]),
                SequenceNeighbor.to_species != sp,
            ))
        out_edges: list[dict] = []
        in_edges: list[dict] = []
        affected: set[str] = set()
        for r, ai in enumerate(new_pos):
            nu = uids[ai]
            nsp = sps[ai]
            row = mat[r]
            mnu = meta.get(nu, {})
            cand = [
                int(j) for j in np.where(row >= thr)[0]
                if sps[int(j)] != nsp
                and _comparable(segs[ai], lens[ai], segs[int(j)], lens[int(j)])
            ]
            cand.sort(key=lambda j: -row[j])
            for j in cand[:k]:
                ju = uids[j]
                ident = float(row[j])
                out_edges.append(_edge(id_of[nu], id_of[ju], kind, ident, meta.get(ju, {})))
                in_edges.append(_edge(id_of[ju], id_of[nu], kind, ident, mnu))
                affected.add(ju)
        for c in range(0, len(out_edges), _INSERT_CHUNK):
            session.execute(insert(SequenceNeighbor), out_edges[c:c + _INSERT_CHUNK])
        total += len(out_edges)
        for c in range(0, len(in_edges), _INSERT_CHUNK):
            session.execute(
                pg_insert(SequenceNeighbor)
                .on_conflict_do_nothing(index_elements=["from_id", "kind", "to_id"]),
                in_edges[c:c + _INSERT_CHUNK],
            )
        total += len(in_edges)
        _trim_cross(session, affected, kind, k)
        session.commit()
    return total


def build_cross_species(
    session: Session, *,
    floor: float = queries.CROSS_FLOOR, k: int = queries.NEIGHBOR_K,
    categories=("coding",), kinds=("nt", "aa"),
) -> int:
    """Add CROSS-species neighbour edges (>= ``floor``, the **0.85** cross floor)
    on top of the existing per-species cache, so an "any species" search is served
    from cache too.

    Cross-species <0.85 is mostly framework similarity, not orthology, so we store
    only >=0.85 and cap each allele's cross edges at ``k`` **separately** (top-k by
    identity, in-loop) — cross edges never evict same-species edges. Additive and
    idempotent: per category/kind it (re)computes, for every species A, A's edges
    to all OTHER species, replacing only A's existing cross-species edges. Reuses
    the per-species build — only the ~N² − Σ|species|² cross pairs. Coding only by
    default (cross-species orthology is a coding concern; rss/leaders stay same).
    """
    thr = floor  # Levenshtein.normalized_similarity is 0..1
    total = 0
    for category in categories:
        for kind in kinds:
            if kind == "aa" and category != "coding":
                continue
            col = Sequence.sequence_aa if kind == "aa" else Sequence.sequence
            rows = session.execute(
                select(Sequence.uid, col, Sequence.species, queries.COMPARISON_CLASS).where(
                    Sequence.category == category, col.is_not(None)
                )
            ).all()
            uids = [r[0] for r in rows]
            strs = [r[1] for r in rows]
            sps = [r[2] for r in rows]
            segs = [r[3] for r in rows]
            species = sorted(set(sps))
            if len(species) < 2:
                continue
            meta = _block_meta(session, category)  # all species
            id_of = {u: meta[u]["id"] for u in meta}
            # Blocks are (species, segment, length-class): an ortholog of a V is a V,
            # and for a segment only compared at equal length (D) an ortholog must
            # also be the same length. The class is None for every other segment, so
            # nothing else is partitioned more finely than before.
            lens = [len(s) for s in strs]
            pos: dict[tuple, list[int]] = {}
            for i, (x, g) in enumerate(zip(sps, segs)):
                klass = lens[i] if queries.equal_length_only(g) else None
                pos.setdefault((x, g, klass), []).append(i)
            # dict.fromkeys rather than sorted(): segment is None for leader/rss and
            # None does not order against str.
            blocks = list(dict.fromkeys((g, kl) for (_, g, kl) in pos))
            # Compute each unordered species PAIR once (cross identity is symmetric)
            # and accumulate both directions per from-uid — no double computation.
            cross: dict[int, list] = defaultdict(list)
            pairs = [
                (species[a], species[b], seg, kl)
                for a in range(len(species))
                for b in range(a + 1, len(species))
                for seg, kl in blocks
            ]
            for spa, spb, seg, kl in pairs:
                    A, B = pos.get((spa, seg, kl), []), pos.get((spb, seg, kl), [])
                    if not A or not B:
                        continue
                    mat = process.cdist(
                        [strs[i] for i in A], [strs[i] for i in B],
                        scorer=queries.similarity_scorer(seg), workers=-1,
                        dtype=np.float32, score_cutoff=thr,
                    )
                    for r, ai in enumerate(A):  # A -> B
                        row = mat[r]
                        for j in np.where(row >= thr)[0]:
                            cross[ai].append((float(row[j]), B[int(j)]))
                    matT = mat.T
                    for r, bi in enumerate(B):  # B -> A (the transpose)
                        row = matT[r]
                        for j in np.where(row >= thr)[0]:
                            cross[bi].append((float(row[j]), A[int(j)]))
            # idempotent: drop existing cross edges (per species; keep same-species)
            for sp in species:
                sp_ids = [
                    id_of[uids[i]]
                    for seg, kl in blocks for i in pos.get((sp, seg, kl), [])
                ]
                if sp_ids:
                    session.execute(
                        delete(SequenceNeighbor).where(
                            SequenceNeighbor.kind == kind,
                            SequenceNeighbor.from_id.in_(sp_ids),
                            SequenceNeighbor.to_species != sp,
                        )
                    )
            edges: list[dict] = []
            for fi, cands in cross.items():
                cands.sort(key=lambda x: -x[0])  # top-k cross per uid, by identity
                for ident, ti in cands[:k]:
                    edges.append(
                        _edge(id_of[uids[fi]], id_of[uids[ti]], kind, ident,
                              meta.get(uids[ti], {}))
                    )
            for c in range(0, len(edges), _INSERT_CHUNK):
                session.execute(
                    pg_insert(SequenceNeighbor)
                    .on_conflict_do_nothing(index_elements=["from_id", "kind", "to_id"]),
                    edges[c:c + _INSERT_CHUNK],
                )
            total += len(edges)
            session.commit()
    return total
