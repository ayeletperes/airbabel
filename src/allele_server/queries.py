"""Transport-agnostic query core: resolve, aliases, similarity search.

Pure functions over a SQLAlchemy ``Session`` returning Pydantic models, so the
REST API now (and an MCP layer later) share one implementation. Results are
capped; there is no bulk dump.
"""

from __future__ import annotations

import hashlib
import html as _html
import json
import re as _re
import urllib.parse
from collections import Counter
from collections.abc import Iterable
from functools import lru_cache
from pathlib import Path

from rapidfuzz import process
from rapidfuzz.distance import Hamming, Levenshtein
from rapidfuzz.fuzz import partial_ratio_alignment
from sqlalchemy import func, not_, or_, select
from sqlalchemy.orm import Session

from .api.schemas import (
    Alias,
    AliasGroup,
    AlleleDescription,
    AlignmentRow,
    AlignmentView,
    ArchivedMatch,
    DiffSummary,
    NativeRecord,
    NeighbourDiff,
    NtSub,
    Ontology,
    OrthologHint,
    PmidHit,
    PmidResult,
    RecordVersion,
    Reference,
    ResolveResult,
    RssAllele,
    RssHit,
    RssLink,
    RssSimilar,
    SearchHit,
    SearchResults,
    SequenceRecord,
    SourceTag,
    Xref,
)
from .db.models import (
    DatasetRelease,
    Name,
    Sequence,
    SequenceAnnotation,
    SequenceLink,
    SequenceNeighbor,
    SequenceSource,
    Source,
)
from .sequence import (
    SequenceError,
    _CODON_TABLE,
    hash_aa,
    hash_sequence,
    reverse_complement_hash,
)
from .settings import get_settings
from .naming import as_subgroup
from .species import SpeciesError, common_name, normalize_species
from .species import ontology as _species_ontology
from .species import scientific_name, taxon_curie
from . import curie as curie_mod
from .examples_meta import EXAMPLES
from .sources import resolve as resolve_source

# Sources shelved from the served site (unpublished internal sets). Read once at
# import — process-static config, like export's attribution constants. Lowercased for
# case-insensitive matching against Source.name. Empty set => every gate below is a
# no-op, so nothing is filtered and no extra SQL is emitted.
_SHELVED = frozenset(s.lower() for s in get_settings().shelved_sources)


def _not_shelved_src():
    """SQL predicate: this (joined) ``Source`` row is not a shelved source."""
    return func.lower(Source.name).notin_(_SHELVED)


def _attested_by(pred):
    """Correlated EXISTS: the outer ``Sequence`` is attested — via a :class:`Name` or
    a :class:`SequenceSource` — by a ``Source`` matching ``pred``.

    ``correlate(Sequence)`` keeps only the outer Sequence correlated, so the
    subquery's own ``Source`` join is not auto-correlated away when the enclosing
    query itself joins ``Source`` (rss_detail / leader_detail do).
    """
    by_name = (
        select(Name.id).join(Source, Source.id == Name.source_id)
        .where(Name.sequence_id == Sequence.id, pred).correlate(Sequence)
    )
    by_rec = (
        select(SequenceSource.id).join(Source, Source.id == SequenceSource.source_id)
        .where(SequenceSource.sequence_id == Sequence.id, pred).correlate(Sequence)
    )
    return or_(by_name.exists(), by_rec.exists())


def _visible_cond():
    """Correlated WHERE selecting only sequences the served site may show.

    A sequence is hidden iff it is attested by a shelved source **and** by no
    non-shelved one — i.e. a HUSA-only novel allele. A sequence attested by some
    non-shelved source stays (its shelved names are dropped elsewhere), and one with
    no shelved attestation at all is unaffected. ``None`` when nothing is shelved.
    """
    if not _SHELVED:
        return None
    shelved = func.lower(Source.name).in_(_SHELVED)
    return or_(not_(_attested_by(shelved)), _attested_by(_not_shelved_src()))


def _visible_seqs(stmt):
    """Add the sequence-visibility gate to a ``select(Sequence…)`` statement."""
    vis = _visible_cond()
    return stmt.where(vis) if vis is not None else stmt


def sitemap_uids(session: Session) -> list[tuple[str, str]]:
    """``(uid, category)`` for every page the sitemap lists: visible, current rows.

    Archived alleles stay reachable (and crawlable through links) but are left out,
    the same way they are left out of coverage and search.
    """
    stmt = _visible_seqs(
        select(Sequence.uid, Sequence.category)
        .where(Sequence.archived.is_(False))
        .order_by(Sequence.category, Sequence.uid)
    )
    return [(uid, cat) for uid, cat in session.execute(stmt)]


def _unshelved_names(stmt):
    """Add the shelved-name exclusion to a stmt that joins :class:`Name`.

    Makes a shelved source's labels both non-searchable and non-displayed: the join
    to ``Source`` drops name rows whose source is shelved, so a HUSA IgLabel neither
    matches a name search nor appears as an alias.
    """
    if not _SHELVED:
        return stmt
    return stmt.join(Source, Source.id == Name.source_id).where(_not_shelved_src())


def _ontology(tag: str) -> Ontology:
    o = _species_ontology(tag)
    return Ontology(id=o.id, label=o.label)

DEFAULT_CAP = 100

# Minimum candidate pool for the fuzzy ranker (:func:`_similarity_search`), kept
# independent of the display cap: with a small cap (e.g. cap=1) a ``cap * N`` pool
# would starve the ranker and let the correct germline fall outside the trigram
# prefilter. The rescore over the pool is cheap, so a generous floor costs little.
MIN_SIMILARITY_POOL = 200

# The neighbour cache (sequence_neighbors) stores, per sequence, its same-species
# neighbours down to NEIGHBOR_FLOOR and its cross-species neighbours down to the
# higher CROSS_FLOOR (cross-species <0.85 is mostly framework noise, not orthology,
# and would otherwise flood the cap). Each side is capped at NEIGHBOR_K
# *separately* (cross never evicts same-species), and NEIGHBOR_K is set above the
# observed max so the caps don't bind — so the read is complete for any threshold
# >= the relevant floor and never needs a per-uid cap check.
NEIGHBOR_FLOOR = 0.75   # same-species floor
CROSS_FLOOR = 0.85      # cross-species floor (and the gate for any-species reads)
NEIGHBOR_K = 1200       # per-side cap; > observed max (~960) so it doesn't bind

# Segments in biological order along the locus, for display. Anything not listed
# sorts after, alphabetically — an unknown segment must still appear.
SEGMENT_ORDER = ("V", "D", "J", "C")
SEGMENT_LABELS = {
    "V": "V (variable)", "D": "D (diversity)",
    "J": "J (joining)", "C": "C (constant)",
}


#: Segments compared only against sequences of **exactly the same length**, and then
#: positionally (Hamming) rather than by edit distance — with the note to display.
#:
#: D is the case. A D gene is 8–148 nt, median 19, and an edit distance is free to
#: open gaps to improve a score. Over sequences that short that is not an alignment
#: finding, it is the metric rearranging noise into a confident number: at 19 nt the
#: 0.75 floor is about four mismatches, which nearly every pair of D genes clears.
#:
#: Restricting to equal length removes the question indels would have answered, and
#: Hamming then compares position by position — which is the comparison two alleles
#: of one D gene actually support. Alleles that differ in length are simply not
#: ranked against each other; separating those needs a multiple-sequence alignment,
#: not a pairwise score.
#:
#: Enforced on the **read** path as well as in the builder. A cache miss falls
#: through to the live ranker, so a builder-only restriction would quietly serve the
#: unrestricted answer instead of the restricted one.
SIMILARITY_EQUAL_LENGTH: dict[str, str] = {
    "D": (
        "D genes are too short (8–148 nt, median 19 nt) for edit-distance scores to be "
        "interpreted reliably; the metric can open gaps that inflate identity, and at "
        "19 nt the 0.75 identity floor is about four mismatches. A D allele is therefore compared only against D "
        "alleles of exactly the same length, scored position by position. Alleles "
        "of the same gene that differ in length are not ranked here; separating "
        "those requires a multiple-sequence alignment, not a pairwise score."
    ),
}


#: The recombination signals we store: subtype -> (label, the order its parts sit in).
#: A signal is always heptamer + spacer + nonamer, but which end each is on depends
#: on which side of the coding sequence the signal lies, so the stored string reads
#: the way the locus does. Getting this wrong is not cosmetic: the page used to slice
#: the first 7 bases as "heptamer" for every kind, which mislabels the nonamer of a
#: J-RS and of a 5'D-RS as a heptamer, and vice versa.
RSS_KINDS: dict[str, tuple[str, tuple[str, str, str]]] = {
    "v_rs": ("V-RS · complete signal, 3' of V", ("heptamer", "spacer", "nonamer")),
    "d_rs_5": ("5'D-RS · complete signal, 5' of D", ("nonamer", "spacer", "heptamer")),
    "d_rs_3": ("3'D-RS · complete signal, 3' of D", ("heptamer", "spacer", "nonamer")),
    "j_rs": ("J-RS · complete signal, 5' of J", ("nonamer", "spacer", "heptamer")),
}

def rss_parts(subtype: str | None, sequence: str | None):
    """Split a stored RSS into ``((name, seq), (name, seq), (name, seq))``, or None.

    The flanking lengths are the canonical 7 and 9. A curated table may publish a
    heptamer of 6 or a nonamer of 10, and where the exact composition was recorded
    at ingest it is on the RSS row as an ``rss/composition`` annotation; this is the
    display split, and it is deliberately the same one for every source so that an
    IMGT signal -- published whole, with no components to record -- renders like the
    rest. ``None`` when the subtype is unknown or the sequence too short to split.
    """
    kind = RSS_KINDS.get(subtype or "")
    if not kind or not sequence:
        return None
    head, _mid, tail = kind[1]
    sizes = {"heptamer": 7, "nonamer": 9}
    a, b = sizes[head], sizes[tail]
    if len(sequence) < a + b + 1:
        return None
    return (
        (head, sequence[:a]),
        ("spacer", sequence[a:len(sequence) - b]),
        (tail, sequence[len(sequence) - b:]),
    )


def comparison_class(segment: str | None, subtype: str | None) -> str | None:
    """The pool a sequence may be ranked inside.

    For a coding allele that is its segment: a V is only ever compared with a V.
    Non-coding rows carry no segment -- every leader and every RSS has
    ``segment IS NULL`` -- so for those it is the subtype, which is the real
    distinction between them. That matters as soon as more than one kind of RSS
    exists: a 5'D-RS (nonamer + 12 nt spacer + heptamer) and an IGKV V-RS
    (heptamer + 12 nt spacer + nonamer) are both 28 nt and would otherwise sit in
    one pool and score against each other, which is the same mistake as ranking a
    D against a V.

    ``segment or subtype`` rather than a branch on category: a C allele carries
    both (``segment='C'``, ``subtype='c_secreted'``) and must class by segment, so
    that the membrane form, when it lands, is still compared with the secreted one.
    """
    return segment or subtype


#: SQL form of :func:`comparison_class`, for predicates and GROUP BY.
COMPARISON_CLASS = func.coalesce(Sequence.segment, Sequence.subtype)


def similarity_reason(segment: str | None) -> str | None:
    """The note explaining a restricted similarity for ``segment``, else None."""
    return SIMILARITY_EQUAL_LENGTH.get(segment) if segment else None


def equal_length_only(segment: str | None) -> bool:
    """True when ``segment`` may only be compared at identical length."""
    return bool(segment) and segment in SIMILARITY_EQUAL_LENGTH


def similarity_profile() -> dict:
    """Every parameter that determines a similarity score, as data.

    Derived from the live constants rather than written out beside them, so a change
    to ``NEIGHBOR_FLOOR`` or to the equal-length policy moves the digest a release is
    stamped with. A hand-maintained copy would drift, and a drifted provenance record
    is worse than none: it asserts reproducibility that no longer holds.

    Ordering is not incidental — the digest is taken over the JSON of this dict with
    sorted keys, so the same parameters always give the same string.
    """
    return {
        # rapidfuzz normalized_similarity: 1 - edits/length over the ungapped
        # sequences, for both nt and aa. See queries.similarity_scorer.
        "metric": "levenshtein_normalized_similarity",
        # Segments that may only be compared at identical length, and are then scored
        # positionally instead. D is the case; see SIMILARITY_EQUAL_LENGTH.
        "equal_length_metric": "hamming_normalized_similarity",
        "equal_length_segments": sorted(SIMILARITY_EQUAL_LENGTH),
        "cross_segment_comparison": False,
        "kinds": ["nt", "aa"],
        "aa_segments": ["V"],
        "neighbor_floor": NEIGHBOR_FLOOR,
        "cross_species_floor": CROSS_FLOOR,
        "neighbor_k": NEIGHBOR_K,
        # What "the sequence" means before any of the above is applied.
        "normalisation": "uppercase, gaps stripped, IUPAC validated, RNA rejected",
        "identity_key": "sha256(normalize(sequence))",
    }


def similarity_digest() -> str:
    """sha256 of :func:`similarity_profile`, the version of the distance computation."""
    blob = json.dumps(similarity_profile(), sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(blob.encode()).hexdigest()


def dataset_release(session: Session) -> dict | None:
    """The most recent stamped build, or ``None`` on a database never stamped.

    ``None`` is a real state, not an error: every instance built before
    ``dataset_releases`` existed has no row until ``allele-ingest --stamp-release``
    runs. Callers must say "unstamped" rather than invent a version.

    ``similarity_current`` compares the stamped distance parameters against the ones
    this process would use *now*. False means the code has moved since the stamp, so
    the served scores are no longer the ones the release describes -- exactly the
    drift the digest exists to catch, and the reason it is reported rather than
    silently corrected.
    """
    row = session.scalars(
        select(DatasetRelease).order_by(DatasetRelease.built_at.desc()).limit(1)
    ).first()
    if row is None:
        return None
    try:
        details = json.loads(row.details) if row.details else {}
    except ValueError:  # a hand-edited row should not take the page down
        details = {}
    return {
        "release": row.release,
        "built_at": row.built_at,
        "content_digest": row.content_digest,
        "similarity_digest": row.similarity_digest,
        "code_version": row.code_version,
        "neighbors_built": row.neighbors_built,
        "neighbor_edges": row.neighbor_edges,
        "similarity_current": row.similarity_digest == similarity_digest(),
        "details": details,
    }


def similarity_scorer(segment: str | None):
    """The distance to rank ``segment`` by.

    Hamming for the equal-length segments (positional, no gaps — the comparison is
    only defined because the lengths already match), normalised Levenshtein for
    everything else. Both are normalised to 0..1 identity, so a threshold means the
    same thing either way.
    """
    return (
        Hamming.normalized_similarity if equal_length_only(segment)
        else Levenshtein.normalized_similarity
    )


# Official allele-designation schemes vs alternative names / identifiers.
_ALLELE_SCHEMES = frozenset(
    {"iuis", "imgt", "ogrdb", "vdjbase", "igdiscover", "iglabel", "husa", "asc",
     "kimdb", "rhgldb", "vrc", "guo"}
)
def _imgt_gene_url(gene: str, species_sci: str) -> str:
    return (
        "https://www.imgt.org/genedb/GENElect?query=2+"
        f"{urllib.parse.quote(gene)}&species={urllib.parse.quote(species_sci)}"
    )


_OGRDB_SEQUENCE = "https://ogrdb.airr-community.org/sequence"


def ogrdb_sequence_url(species_tag: str | None, name: str | None) -> str | None:
    """Deep-link a gene or allele into OGRDB's per-sequence route:
    ``/sequence/<ncbi_taxon_id>/<gene-or-allele>``.

    OGRDB resolves it three ways: the **allele** record if it holds that allele; the
    **gene's** alignments if only the gene is known (or the allele is absent); else a
    not-found message. Because it answers for *any* string, callers must only offer the
    link where OGRDB actually attests the sequence (see ``build_record``): a chip on an
    allele OGRDB has never seen would advertise a curated record that does not exist.
    """
    if not name or not species_tag:
        return None
    try:
        taxon = taxon_curie(species_tag)
    except Exception:  # noqa: BLE001 - unknown species -> no link, never an error
        return None
    if not taxon or ":" not in taxon:
        return None
    tid = taxon.split(":", 1)[1]
    # '*' is a legal path char and OGRDB's examples use it unescaped.
    return f"{_OGRDB_SEQUENCE}/{tid}/{urllib.parse.quote(name, safe='*')}"


def _alias_url(alias: Alias) -> str | None:
    """Only *identifier* aliases link out. An allele designation is a label, not a
    record, so it is rendered as plain text: the outbound links for IMGT and OGRDB
    belong in **External records**, where a chip asserts "this sequence has a record
    there". An NCBI accession is itself an identifier for a record, so it keeps a link.
    """
    if alias.scheme == "ncbi" or alias.source == "ncbi_nucleotide":
        return f"https://www.ncbi.nlm.nih.gov/nuccore/{alias.name}"
    return None


# --- record assembly ------------------------------------------------------

def _aliases(session: Session, seq_id: int) -> list[Alias]:
    stmt = (
        select(Name.name, Name.naming_scheme, Source.name, Source.version, Source.note)
        .join(Source, Source.id == Name.source_id)
        .where(Name.sequence_id == seq_id)
        .order_by(Name.naming_scheme, Name.name)
    )
    if _SHELVED:  # drop names contributed by a shelved source (e.g. HUSA IgLabel)
        stmt = stmt.where(_not_shelved_src())
    rows = session.execute(stmt).all()
    out: list[Alias] = []
    for n, s, src, ver, note in rows:
        hist = bool(note and note.startswith("historical IMGT"))
        out.append(Alias(
            name=n, scheme=s, source=src,
            historical=hist, release=ver if hist else None,
        ))
    return out


def _group_aliases(aliases: list[Alias]) -> list["AliasGroup"]:
    """Collapse aliases to one entry per distinct ``(name, scheme)``, gathering
    every source that attributes it. Keeps the names table to one row per name
    (IMGT alone spans ~14 release-version rows for the same name)."""
    groups: dict[tuple, AliasGroup] = {}
    for a in aliases:
        key = (a.name, a.scheme, a.historical, a.release)
        g = groups.get(key)
        if g is None:
            groups[key] = AliasGroup(
                name=a.name, scheme=a.scheme, url=a.url,
                historical=a.historical, release=a.release, sources=[a.source],
            )
        elif a.source not in g.sources:
            g.sources.append(a.source)
    for g in groups.values():
        g.sources.sort()
    return list(groups.values())


def _source_tags(session: Session, seq_id: int) -> list[SourceTag]:
    """Every source that contributes this sequence — via a name OR an attribution —
    shown **once per source**. A source spanning many IMGT release versions is one
    source (we show the current, not every version); and HUSA / KIARVA contribute
    names without a separate attribution row, so they must be picked up from names.
    """
    best: dict[str, tuple] = {}

    def consider(n, v, cp, r, oid, note):
        hist = bool(note and note.startswith("historical IMGT"))
        gdb = (v or "").startswith("GENE-DB")
        ver = "GENE-DB" if gdb else v
        # prefer current (non-historical), then a real release over the all-species
        # GENE-DB marker, then a row carrying a link, then the latest version
        rank = (0 if hist else 1, 0 if gdb else 1, 1 if oid else 0, v or "")
        cur = best.get(n)
        if cur is None or rank > cur[0]:
            best[n] = (
                rank,
                SourceTag(name=n, version=ver, curie_prefix=cp, redistributable=r,
                          original_id=oid, url=_source_url(n, oid)),
            )

    # A shelved source (HUSA) is not shown even on a sequence it shares with a
    # visible source, so its chip never appears anywhere on the site.
    _shelved_skip = [_not_shelved_src()] if _SHELVED else []
    for n, v, cp, r, oid, note in session.execute(
        select(Source.name, Source.version, Source.curie_prefix,
               Source.redistributable, SequenceSource.original_id, Source.note)
        .join(SequenceSource, SequenceSource.source_id == Source.id)
        .where(SequenceSource.sequence_id == seq_id, *_shelved_skip)
    ):
        consider(n, v, cp, r, oid, note)
    for n, v, cp, r, note in session.execute(
        select(Source.name, Source.version, Source.curie_prefix,
               Source.redistributable, Source.note)
        .join(Name, Name.source_id == Source.id)
        .where(Name.sequence_id == seq_id, *_shelved_skip)
    ):
        consider(n, v, cp, r, None, note)
    return [tag for _, (_, tag) in sorted(best.items(), key=lambda kv: kv[0])]


def _source_url(name: str, original_id: str | None) -> str | None:
    """The source's own site, else the CURIE xref's resolved URL.

    Both halves now come from the manifest (``homepage_url`` / ``curie_url_template``)
    rather than the hardcoded per-source branches this replaces.
    """
    manifest = resolve_source(name)
    if manifest is not None and manifest.homepage_url:
        return manifest.homepage_url
    return curie_mod.resolve_url(original_id) if original_id else None


_XREF_LABELS = {"ncbi_gene": "NCBI Gene", "hgnc": "HGNC", "ncbi_url": "NCBI nucleotide"}


def _xrefs(session: Session, seq_id: int) -> list[Xref]:
    """Linkable external cross-references from the ``xref`` annotation domain."""
    rows = session.scalars(
        select(SequenceAnnotation).where(
            SequenceAnnotation.sequence_id == seq_id,
            SequenceAnnotation.domain == "xref",
        )
    )
    seen: set[tuple[str, str]] = set()
    out: list[Xref] = []
    for a in rows:
        val = a.value_text or ""
        if a.key == "match_confidence" or not val:
            continue
        if a.key.startswith("pmid:"):
            continue  # PubMed references have their own section (see _references)
        if a.key == "ncbi_url":
            label, curie, url = "NCBI nucleotide", "", val
        elif a.key.startswith("imgt_accession:"):
            label, curie, url = f"GenBank {a.key.split(':', 1)[1]}", val, curie_mod.resolve_url(val)
        else:
            label, curie, url = _XREF_LABELS.get(a.key, a.key), val, curie_mod.resolve_url(val)
        dedup = (label, curie or url)
        if dedup in seen:
            continue
        seen.add(dedup)
        out.append(Xref(label=label, curie=curie, url=url))
    return out


def _example_holds(ex: dict, rec: NativeRecord) -> bool:
    """Does the record still support what the example *says* about it?

    An example asserts a status ("archived", "retired", "IMGT removed it"), and a
    status is derived — it changes under reload. Checking only that the UID resolves
    caught a dead link but let the prose go stale: /examples described an allele as
    archived and gone from IMGT while the live record was neither, because a later
    ingest pass attributed it to a current source. Each entry now carries the claim
    in a form that can be checked, and a mismatch drops the example rather than
    publishing a false statement about our own data.
    """
    want = ex.get("expect") or {}
    checks = {
        "archived": lambda v: rec.archived is v,
        "retired": lambda v: rec.retired is v,
        "has_newer_version": lambda v: (rec.newer_version is not None) is v,
        "prev_version_count": lambda v: len(rec.prev_versions) == v,
        "former_names": lambda v: set(v) <= {a.name for a in rec.former_names},
        "names": lambda v: set(v) <= {a.name for a in rec.allele_names},
        "length": lambda v: rec.length == v,
    }
    return all(checks[k](v) for k, v in want.items() if k in checks)


def example_showcase(session: Session) -> list[dict]:
    """Resolve the curated ``/examples`` entries against live data.

    Allele examples are looked up by UID and enriched with what the page needs to
    make its point (names, former names, archived/retired flags, protein UID). An
    example is **skipped** when its UID no longer resolves *or* when the record stops
    supporting the claim it makes (``expect``, see ``_example_holds``), so a reload
    can leave the page shorter but never wrong. Search examples become a query URL.
    """
    out: list[dict] = []
    for ex in EXAMPLES:
        if ex.get("kind") == "search":
            params = {"q": ex["query"], "mode": ex.get("mode", "sequence")}
            if ex.get("partial"):
                params["partial"] = "1"
            out.append({**ex, "url": "/?" + urllib.parse.urlencode(params)})
            continue
        if ex.get("kind") == "pmid":
            # Resolved against live data for the same reason allele examples are: the
            # claim "this paper attests alleles here" is derived from annotations that
            # a reload can drop. A paper that stopped attesting anything must take its
            # example with it rather than leave the page promising a hit.
            res = alleles_for_pmid(session, ex["pmid"], cap=DEFAULT_CAP)
            want = ex.get("expect") or {}
            if res is None or res.count < want.get("min_hits", 1):
                continue
            tier = want.get("kind")
            if tier and not any(h.kind == tier for h in res.hits):
                continue
            out.append({
                **ex,
                "url": "/?" + urllib.parse.urlencode(
                    {"q": ex["pmid"], "mode": "pmid"}
                ),
                "pmid_url": res.url,
                "n_hits": res.count,
                "species": ", ".join(sorted({
                    common_name(h.species) for h in res.hits
                })),
                "names": [h.label for h in res.hits[:4] if h.label],
            })
            continue
        seq = session.scalar(select(Sequence).where(Sequence.uid == ex["uid"]))
        if seq is None:  # data changed under us -> drop the example, never 404
            continue
        rec = build_record(session, seq).native
        if not _example_holds(ex, rec):  # the claim no longer matches the record
            continue
        # Other nt alleles encoding the identical V-REGION protein — what makes the
        # "different DNA, same protein" example legible.
        siblings: list[dict] = []
        if seq.protein_uid:
            for other in session.scalars(
                select(Sequence).where(
                    Sequence.protein_uid == seq.protein_uid,
                    Sequence.uid != seq.uid,
                    Sequence.archived.is_(False),
                )
            ):
                nm = session.scalar(
                    select(Name.name).where(
                        Name.sequence_id == other.id, Name.naming_scheme == "imgt"
                    ).limit(1)
                )
                siblings.append({"uid": other.uid, "name": nm or other.uid})
        out.append({
            **ex,
            "url": f"/seq/{seq.uid}",
            "species": common_name(seq.species),
            "names": sorted({a.name for a in rec.allele_names}),
            "former": [(a.name, a.release) for a in rec.former_names],
            "archived": rec.archived,
            "retired": rec.retired,
            "protein_uid": seq.protein_uid,
            "protein_siblings": siblings,
            "n_refs": len(rec.references),
        })
    return out


def source_catalog(session: Session) -> list[dict]:
    """Every source in the database with its live counts + curated metadata.

    Powers ``/sources``: what AIRBabel aggregates, under what licence, and how to
    cite it — enough for someone to rebuild the dataset from the primary sources.
    Counts are distinct sequences *named* by that source.
    """
    rows = session.execute(
        select(Source.name, Source.kind, Source.curie_prefix,
               func.count(func.distinct(Source.version)),
               func.max(Source.version), func.max(Source.note))
        .group_by(Source.name, Source.kind, Source.curie_prefix)
    ).all()
    counts = dict(session.execute(
        select(Source.name, func.count(func.distinct(Name.sequence_id)))
        .join(Name, Name.source_id == Source.id)
        .group_by(Source.name)
    ).all())

    out: list[dict] = []
    for name, kind, cp, n_versions, version, note in rows:
        if name.lower() in _SHELVED:  # unpublished internal set: never on /sources
            continue
        manifest = resolve_source(name)
        # A source may carry its DOI in `version` (KIARVA does); the manifest wins.
        is_doi = bool(version) and version.startswith("10.")
        doi = (manifest.citation.doi if manifest else None) or (version if is_doi else None)
        licence = manifest.licence if manifest else None
        citation = manifest.citation if manifest else None
        out.append({
            "name": name,
            "label": manifest.label if manifest else name,
            "kind": kind,
            "description": (manifest.description if manifest else "") or (note or ""),
            "curie_prefix": cp,
            "versions": n_versions,
            "version": None if is_doi else version,
            "alleles": counts.get(name, 0),
            "url": manifest.homepage_url if manifest else None,
            "licence": licence.name if licence else None,
            "citation": citation.text or None if citation else None,
            "doi": doi,
            "doi_url": f"https://doi.org/{doi}" if doi else None,
            # New: lets /sources say "no citation available: <why>" instead of leaving a
            # blank cell that reads as an oversight.
            "citation_status": citation.status if citation else None,
            "citation_note": citation.note if citation else None,
            "licence_status": licence.status if licence else None,
            "redistributable": bool(licence and licence.redistributable),
            "manifest_id": manifest.slug if manifest else None,
        })
    out.sort(key=lambda d: (-d["alleles"], d["name"]))
    return out


def _accession_xrefs(session: Session, seq_id: int) -> list[Xref]:
    """GenBank accessions recorded for this sequence by any source.

    The accession is the sequence's **primary public record** — deposited in GenBank
    by the original submitters, not owned by any curator database — so it is the one
    external identifier we surface for provenance, linked straight to NCBI.
    """
    out: list[Xref] = []
    rows = session.scalars(
        select(SequenceSource.original_id)
        .where(
            SequenceSource.sequence_id == seq_id,
            SequenceSource.original_id.like("GENBANK:%"),
        )
        .distinct()
    )
    for oid in rows:
        out.append(Xref(label="GenBank", curie=oid, url=curie_mod.resolve_url(oid)))
    return out


def _dedup_xrefs(xrefs: list[Xref]) -> list[Xref]:
    seen: set[tuple[str, str]] = set()
    out: list[Xref] = []
    for x in xrefs:
        key = (x.label, x.curie or x.url or "")
        if key in seen:
            continue
        seen.add(key)
        out.append(x)
    return out


_TAG_RE = _re.compile(r"<[^>]+>")


def _clean_title(val: str | None) -> str | None:
    """Europe PMC titles carry inline markup (``<i>``, ``<sup>``) and sometimes
    double-encoded entities (``&lt;sup&gt;``, ``&amp;``). The templates auto-escape,
    so unescape (twice, for the double-encoded case) then strip tags to plain text."""
    if not val:
        return None
    txt = _TAG_RE.sub("", _html.unescape(_html.unescape(val)))
    txt = _re.sub(r"\s+", " ", txt).strip()
    return txt or None


_KIND_RANK = {"verified": 0, "sequence": 1, "name": 2}

# A paper's *species* is whatever its verified/sequence links say it is. One identical
# sequence shared with a second species must not make that species the paper's too
# (a human paper with 47 human hits and one rhesus hit is not a rhesus paper), so a
# species anchors only if it holds at least this share of the top species' evidence.
ANCHOR_MIN_SHARE = 0.2


def anchor_pairs(evidence: Iterable[tuple[str, str | None]]) -> set[tuple[str, str | None]]:
    """``(species, locus)`` pairs a paper is *about*, from its evidence-tier links.

    A bare allele **name** is a string many species' nomenclatures share, so a name
    match says nothing about which species a paper studies; only a sequence match (or
    a curator citation) does. A name link is consistent with the paper only when its
    (species, locus) is one of these pairs. No evidence -> no anchor -> every name
    link is unanchored: we do not guess a species the paper never indicated.
    """
    pairs = Counter(evidence)
    by_species: Counter[str] = Counter()
    for (sp, _locus), n in pairs.items():
        by_species[sp] += n
    top = max(by_species.values(), default=0)
    keep = {sp for sp, n in by_species.items() if n >= ANCHOR_MIN_SHARE * top}
    return {pair for pair in pairs if pair[0] in keep}


def _references(session: Session, seq_id: int) -> list[Reference]:
    """PubMed references for an allele, tiered by evidence strength.

    Three provenances collapse by PMID (strongest kind wins):

    * curator references (domain ``xref``) -> ``verified``;
    * mined hits whose companion ``suppmatch:`` annotation records a **sequence**
      match (full-text / supplement hash) -> ``sequence`` (cross-species-safe);
    * mined hits by **name** only -> ``name`` (a lead — same string, other species).

    See ``ingest/literature.py`` (name probe) and ``ingest/literature_supp.py``
    (full-text + supplement sequence mining). Auto hits are not manually checked.
    """
    rows = session.execute(
        select(SequenceAnnotation.key, SequenceAnnotation.value_text,
               SequenceAnnotation.domain, Source.name)
        .join(Source, Source.id == SequenceAnnotation.source_id)
        .where(
            SequenceAnnotation.sequence_id == seq_id,
            or_(SequenceAnnotation.key.like("pmid:%"),
                SequenceAnnotation.key.like("suppmatch:%")),
        )
    ).all()

    # First pass: the mining companions tell us how each PMID was matched.
    seq_detail: dict[str, str] = {}   # pmid -> provenance of the best sequence match
    name_pmids: set[str] = set()
    for key, val, _domain, _src in rows:
        if not key.startswith("suppmatch:"):
            continue
        pmid = key.split(":", 2)[1].strip()
        if (val or "").startswith("sequence"):
            # keep the strongest (exact > similar) provenance string
            if pmid not in seq_detail or (val or "").startswith("sequence_exact"):
                seq_detail[pmid] = val
        else:
            name_pmids.add(pmid)

    by_pmid: dict[str, Reference] = {}
    for key, val, domain, src in rows:
        if not key.startswith("pmid:"):
            continue
        pmid = key.split(":", 1)[1].strip()
        if not pmid:
            continue
        verified = domain != "literature"  # curator table = verified; mined = not
        title = _clean_title(val) if (domain == "literature" and val and not val.startswith("PMID")) else None
        if verified:
            kind, detail = "verified", None
        elif pmid in seq_detail:
            kind, detail = "sequence", seq_detail[pmid]
        else:
            kind, detail = "name", None
        cur = by_pmid.get(pmid)
        if cur is None:
            by_pmid[pmid] = Reference(
                pmid=pmid, url=f"https://pubmed.ncbi.nlm.nih.gov/{pmid}",
                title=title, verified=verified, kind=kind, detail=detail, source=src,
            )
        else:
            cur.title = cur.title or title
            cur.source = cur.source or src
            if _KIND_RANK[kind] < _KIND_RANK[cur.kind]:  # promote to the stronger tier
                cur.kind, cur.detail, cur.verified = kind, detail, verified
    return sorted(
        by_pmid.values(),
        key=lambda r: (_KIND_RANK[r.kind], int(r.pmid) if r.pmid.isdigit() else 0),
    )


def normalize_pmid(raw: str) -> str | None:
    """Pull a bare PMID out of what a user is likely to paste.

    Accepts ``29163486``, ``PMID:29163486``, ``pmid: 29163486`` and a PubMed URL.
    Returns None when there is no plausible PMID, so the caller can fall through to
    another search mode instead of running a lookup that cannot match.
    """
    if not raw:
        return None
    s = raw.strip()
    if not s:
        return None
    low = s.lower()
    if "ncbi.nlm.nih.gov" in low or low.startswith(("http://", "https://")):
        # .../pubmed/29163486, possibly with a trailing slash or query string
        parts = [p for p in _re.split(r"[/?#]", s) if p]
        for part in reversed(parts):
            if part.isdigit() and 4 <= len(part) <= 9:
                return part
        return None
    if low.startswith("pmid"):
        s = s[4:].lstrip(": ").strip()
    if s.isdigit() and 4 <= len(s) <= 9:
        return s
    # A bare number embedded in something else is deliberately NOT accepted: an
    # allele name contains digits, and guessing here would hijack name search.
    return None


def alleles_for_pmid(
    session: Session, pmid: str, species: str | None = None, cap: int = DEFAULT_CAP,
    include_name: bool = False,
) -> PmidResult | None:
    """Every allele a PubMed record attests, tiered by how it attests them.

    The inverse of ``_references``, and it must agree with it: the same PMID on the
    same allele has to report the same tier from either direction. So the tier is
    derived the same way — ``xref`` domain is curator-verified; a mined hit is
    ``sequence`` when a companion ``suppmatch:`` annotation records a full-text or
    supplement hash match for *that* sequence, and ``name`` otherwise.

    Archived alleles are included. Excluding them would defeat the reason they are
    kept: a paper cited an allele, and the whole point is that the citation still
    leads somewhere. They are flagged, and a superseded record links to its successor.

    **Name-only links are not evidence** and are withheld by default: a supplementary
    table that lists ``IGHV1-69*01`` matches that string in every species carrying it
    (one spatial-transcriptomics paper "attested" 1,627 alleles this way). With
    ``include_name`` they are returned, labelled ``name`` (unverified), but only when
    their (species, locus) agrees with the paper's own verified/sequence links
    (:func:`anchor_pairs`). ``n_name`` counts those consistent leads and
    ``n_name_unanchored`` the ones dropped for disagreeing, so nothing disappears
    unreported.

    Returns None when ``pmid`` is not a PMID at all, so callers can fall through.
    """
    bare = normalize_pmid(pmid)
    if bare is None:
        return None
    tags = _species_tags(species)
    key = f"pmid:{bare}"
    rows = session.execute(
        select(SequenceAnnotation.sequence_id, SequenceAnnotation.domain,
               SequenceAnnotation.value_text, Source.name)
        .join(Source, Source.id == SequenceAnnotation.source_id)
        .where(SequenceAnnotation.key == key)
    ).all()
    if not rows:
        return PmidResult(
            pmid=bare, url=f"https://pubmed.ncbi.nlm.nih.gov/{bare}", count=0, hits=[]
        )
    # Companion mining annotations: which sequences this paper matched *by sequence*.
    seq_matched: dict[int, str] = {}
    for sid, val in session.execute(
        select(SequenceAnnotation.sequence_id, SequenceAnnotation.value_text).where(
            or_(SequenceAnnotation.key == f"suppmatch:{bare}",
                SequenceAnnotation.key.like(f"suppmatch:{bare}:%")),
        )
    ).all():
        if (val or "").startswith("sequence"):
            if sid not in seq_matched or (val or "").startswith("sequence_exact"):
                seq_matched[sid] = val
    title: str | None = None
    best: dict[int, tuple[str, str | None, str | None]] = {}
    for sid, domain, val, src in rows:
        if domain == "literature" and val and not val.startswith("PMID"):
            title = title or _clean_title(val)
        if domain != "literature":
            kind, detail = "verified", None
        elif sid in seq_matched:
            kind, detail = "sequence", seq_matched[sid]
        else:
            kind, detail = "name", None
        cur = best.get(sid)
        # One sequence can be annotated by several sources; the strongest tier wins,
        # exactly as _references promotes on the allele page.
        if cur is None or _KIND_RANK[kind] < _KIND_RANK[cur[0]]:
            best[sid] = (kind, detail, src)
    seqs = list(session.execute(
        _visible_seqs(select(Sequence).where(Sequence.id.in_(list(best))))
    ).scalars())
    # Anchors come from the whole paper, before any species filter narrows the view.
    anchors = anchor_pairs(
        (s.species, s.locus) for s in seqs if best[s.id][0] != "name"
    )
    if tags is not None:
        seqs = [s for s in seqs if s.species in tags]
    n = {"verified": 0, "sequence": 0, "name": 0}
    n_unanchored = 0
    keep: list[Sequence] = []
    for s in seqs:
        kind = best[s.id][0]
        if kind == "name" and (s.species, s.locus) not in anchors:
            n_unanchored += 1
            continue
        n[kind] += 1
        if kind != "name" or include_name:
            keep.append(s)
    seqs = keep
    seqs.sort(key=lambda s: (_KIND_RANK[best[s.id][0]], s.species, s.uid))
    capped = len(seqs) > cap
    page = seqs[:cap]
    versions = _versions_for(session, page)
    hits: list[PmidHit] = []
    for s in page:
        kind, detail, src = best[s.id]
        newer = versions.get(s.id, ([], None))[1]
        hits.append(PmidHit(
            uid=s.uid, label=_label(_aliases(session, s.id)), species=s.species,
            locus=s.locus, sequence_type=s.segment, kind=kind, detail=detail,
            source=src, archived=s.archived,
            superseded_by=newer.uid if newer else None,
        ))
    return PmidResult(
        pmid=bare, url=f"https://pubmed.ncbi.nlm.nih.gov/{bare}", title=title,
        count=len(seqs), capped=capped,
        n_verified=n["verified"], n_sequence=n["sequence"], n_name=n["name"],
        n_name_unanchored=n_unanchored, include_name=include_name,
        hits=hits,
    )


def _germline_set_ref(session: Session, seq_id: int) -> str | None:
    return session.scalar(
        select(Source.germline_set_ref)
        .join(SequenceSource, SequenceSource.source_id == Source.id)
        .where(
            SequenceSource.sequence_id == seq_id,
            Source.germline_set_ref.is_not(None),
        )
        .limit(1)
    )


def _rss_links(session: Session, seq_id: int) -> list[RssLink]:
    """RSS sequences linked to a coding allele, with the source(s) of each link."""
    rows = session.execute(
        select(SequenceLink.relation, Sequence.uid, Sequence.subtype,
               Sequence.sequence, Source.name)
        .join(Sequence, Sequence.id == SequenceLink.to_sequence_id)
        .join(Source, Source.id == SequenceLink.source_id)
        .where(
            SequenceLink.from_sequence_id == seq_id,
            SequenceLink.relation.like("has_rss_%"),
            *([_not_shelved_src()] if _SHELVED else []),
        )
        .distinct()
        .order_by(SequenceLink.relation)
    ).all()
    by_uid: dict[tuple, RssLink] = {}
    for rel, uid, sub, seq, src in rows:  # same RSS may be linked by several sources
        link = by_uid.setdefault(
            (rel, uid), RssLink(uid=uid, subtype=sub, sequence=seq, relation=rel, sources=[])
        )
        if src and src not in link.sources:
            link.sources.append(src)
    return list(by_uid.values())


def _composition(session: Session, seq_id: int) -> str | None:
    """Domain composition of an assembled constant region, e.g.
    ``CH1:294+H:45+CH2:330+CH3-CHS:321``, or None when the source published the
    sequence whole.

    Surfaced on the record because a sequence we assembled and one a source
    published must not look alike: the reader is entitled to know which domains
    went in, in what order, and where the boundaries fall. Written by
    ``ingest/imgt_constant.py``.
    """
    return session.scalar(
        select(SequenceAnnotation.value_text).where(
            SequenceAnnotation.sequence_id == seq_id,
            SequenceAnnotation.domain == "constant",
            SequenceAnnotation.key == "domains",
        )
    )


def _leader_links(session: Session, seq_id: int) -> list[RssLink]:
    """Spliced leader (exon1+exon2) sequences linked to a coding allele."""
    rows = session.execute(
        select(
            SequenceLink.relation, Sequence.uid, Sequence.subtype,
            Sequence.sequence, Sequence.id, Source.name,
        )
        .join(Sequence, Sequence.id == SequenceLink.to_sequence_id)
        .join(Source, Source.id == SequenceLink.source_id)
        .where(
            SequenceLink.from_sequence_id == seq_id,
            SequenceLink.relation == "has_leader",
            *([_not_shelved_src()] if _SHELVED else []),
        )
        .distinct()
    ).all()
    by_uid: dict[str, RssLink] = {}
    for rel, uid, sub, seq, leader_id, src in rows:
        if uid not in by_uid:
            exon1 = session.scalar(
                select(SequenceAnnotation.value_int).where(
                    SequenceAnnotation.sequence_id == leader_id,
                    SequenceAnnotation.domain == "leader",
                    SequenceAnnotation.key == "exon1_len",
                ).limit(1)
            )
            by_uid[uid] = RssLink(
                uid=uid, subtype=sub, sequence=seq, relation=rel,
                exon1_len=exon1, sources=[],
            )
        if src and src not in by_uid[uid].sources:
            by_uid[uid].sources.append(src)
    return list(by_uid.values())


# Preference order for the representative allele label when a sequence carries
# names from several schemes. IMGT/IUIS nomenclature first, then VDJbase,
# IgDiscover, the remaining curated schemes, and the internal IgLabel last.
_LABEL_PREF = (
    "imgt", "iuis", "vdjbase", "igdiscover", "ogrdb",
    "rhgldb", "kimdb", "vrc", "guo", "husa", "asc", "iglabel",
)


def _label(aliases: list[Alias], prefer: str | None = None) -> str | None:
    """Pick a representative allele name by scheme preference (see _LABEL_PREF).

    ``prefer`` is the text the caller searched on. A sequence carries a name from
    every source that attests it, and scheme preference alone picks the same one
    whatever was asked -- so a search for ``IGKV2`` could answer with rows reading
    ``IGKV10-QHYS*05`` and ``lib4kappa_39_S0151``. Both are correct records (each
    also carries an ``IGKV2-...`` name) and both look like a broken search. Where
    a name containing the query exists, it is the one shown; scheme preference
    then decides among those, so the choice is still deterministic and still
    prefers a designation over a lab identifier.
    """
    if not aliases:
        return None
    pool = aliases
    if prefer:
        needle = prefer.strip().lower()
        matched = [a for a in aliases if needle and needle in a.name.lower()]
        if matched:
            pool = matched
    for pref in _LABEL_PREF:
        for a in pool:
            if a.scheme == pref:
                return a.name
    return pool[0].name


def _name_spans(
    session: Session, names: list[str], species: list[str]
) -> dict[tuple[int, str], tuple[str | None, str | None, bool]]:
    """For every sequence in ``species`` carrying one of ``names`` as an IMGT
    designation, the release range that name was seen in and whether a **current**
    (non-historical) source still attributes it.

    Keyed ``(sequence_id, name)``: the unit of a version chain is a name on a
    sequence, not a sequence, since one record can carry several designations.
    """
    if not names or not species:
        return {}
    rows = session.execute(
        select(Name.sequence_id, Name.name, Source.version, Source.note)
        .join(Source, Source.id == Name.source_id)
        .join(Sequence, Sequence.id == Name.sequence_id)
        .where(
            Name.name.in_(names),
            Name.naming_scheme == "imgt",
            Sequence.species.in_(species),
        )
    ).all()
    out: dict[tuple[int, str], tuple[str | None, str | None, bool]] = {}
    for sid, nm, ver, note in rows:
        first, last, current = out.get((sid, nm), (None, None, False))
        if note and note.startswith("historical IMGT"):
            v = ver or ""
            first = v if first is None else min(first, v)
            last = v if last is None else max(last, v)
        else:
            current = True
        out[(sid, nm)] = (first, last, current)
    return out


def _overlaps(a: tuple[str | None, str | None], b: tuple[str | None, str | None]) -> bool:
    """Do two release ranges intersect? Two records that held one name in the same
    release are concurrent — an ambiguity, not a version chain."""
    if a[0] is None or b[0] is None:
        return False  # a record with no historical span never conflicts on dates
    return a[0] <= b[1] and b[0] <= a[1]  # type: ignore[operator]


def _versions_for(
    session: Session, seqs: list[Sequence]
) -> dict[int, tuple[list[RecordVersion], RecordVersion | None]]:
    """Other records of the **same allele** for a batch of sequences: same IMGT
    designation, disjoint release range, different sequence.

    IMGT re-cuts a record without renaming it (a partial V-REGION later completed, a
    long record trimmed). Identity is the hash, so each extent is its own row and the
    older one is left with no current name. Deriving the chain from the release
    ranges is exact, unlike inferring it from similarity: an extent change is charged
    like a mutation, so a record can score below the neighbour floor against its own
    successor.

    Batched deliberately — three queries whatever the input size, so a page of search
    hits costs the same as one allele page.
    """
    out: dict[int, tuple[list[RecordVersion], RecordVersion | None]] = {
        s.id: ([], None) for s in seqs
    }
    if not seqs:
        return out
    mine_names = session.execute(
        select(Name.sequence_id, Name.name).where(
            Name.sequence_id.in_([s.id for s in seqs]), Name.naming_scheme == "imgt"
        ).distinct()
    ).all()
    if not mine_names:
        return out
    names = sorted({nm for _sid, nm in mine_names})
    species = sorted({s.species for s in seqs})
    spans = _name_spans(session, names, species)
    wanted = {sid for s in seqs for (sid, _nm) in spans if sid != s.id}
    if not wanted:
        return out
    rows = {
        s.id: s for s in session.execute(
            select(Sequence).where(Sequence.id.in_(sorted(wanted)))
        ).scalars()
    }
    # A record IMGT still designates has not been superseded by anything: if one of
    # its old names now sits on another hash, that name was *reassigned* — the item-6
    # rename, one sequence changing names — not this record being re-cut. Mouse
    # IGHV2-9-1*01 moved to a different sequence in 2015 while its original record
    # carried on as IGHV2-9*02, and without this it reads as superseded by the
    # sequence that took its old name.
    involved = sorted({s.id for s in seqs} | set(rows))
    current_named = {
        sid for (sid,) in session.execute(
            select(Name.sequence_id)
            .join(Source, Source.id == Name.source_id)
            .where(
                Name.sequence_id.in_(involved),
                Name.naming_scheme == "imgt",
                or_(Source.note.is_(None), Source.note.notlike("historical IMGT%")),
            ).distinct()
        )
    }
    by_name: dict[str, list[int]] = {}
    for (sid, nm) in spans:
        by_name.setdefault(nm, []).append(sid)
    for seq in seqs:
        older: list[tuple[str, RecordVersion]] = []
        newer: list[tuple[str, RecordVersion]] = []
        for sid0, nm in mine_names:
            if sid0 != seq.id:
                continue
            mine = spans.get((seq.id, nm))
            if mine is None:
                continue
            for sid in by_name.get(nm, ()):
                if sid == seq.id or sid not in rows:
                    continue
                other = rows[sid]
                if other.species != seq.species:
                    continue  # one name, two species is ambiguity, never a version
                first, last, current = spans[(sid, nm)]
                if _overlaps((first, last), (mine[0], mine[1])):
                    continue  # concurrent holders of one name: not a version chain
                if current and mine[2]:
                    continue  # two current records for a name: ambiguity, not a chain
                if other.sequence in seq.sequence:
                    relation = "trims"
                elif seq.sequence in other.sequence:
                    relation = "extends"
                else:
                    relation = "differs"
                ver = RecordVersion(
                    uid=other.uid, name=nm, length=other.length, relation=relation,
                    current=current, first_seen=first, last_seen=last,
                )
                forward = (
                    (current and not mine[2])
                    or (not mine[2] and not current and (last or "") > (mine[1] or ""))
                )
                if forward:
                    # Only a record with no current designation of its own can have
                    # been superseded; otherwise IMGT still publishes it.
                    if seq.id not in current_named:
                        newer.append((last or "9999", ver))
                elif sid not in current_named:
                    older.append((last or "", ver))
        older.sort(key=lambda t: (t[0], t[1].uid))
        newer.sort(key=lambda t: (t[0], t[1].uid))
        # Several successors can only happen if IMGT split the name again; the last
        # one wins, and the rest stay visible through the loser's own page.
        out[seq.id] = ([v for _k, v in older], newer[-1][1] if newer else None)
    return out


def _record_versions(
    session: Session, seq: Sequence
) -> tuple[list[RecordVersion], RecordVersion | None]:
    """``_versions_for`` for a single record — see there."""
    return _versions_for(session, [seq]).get(seq.id, ([], None))


def _archived_matches(
    session: Session, q: str, species: str | None, cap: int = 5
) -> list[ArchivedMatch]:
    """Records an archived-excluding name search deliberately skipped.

    Called only when a name search found nothing. Silence is indistinguishable from
    "no such allele", which is wrong twice over: the record exists, and where it was
    merely re-cut the searcher is one link from the record they wanted.
    """
    q = q.strip()
    if not q:
        return []
    try:
        seqs = [s for s in resolve_exact_name(session, q, species) if s.archived]
    except SpeciesError:
        return []
    if not seqs:
        return []
    seqs = sorted(seqs, key=lambda s: s.uid)[:cap]
    vers = _versions_for(session, seqs)
    # Derive the release range from the snapshot sources that name each record rather
    # than from the ``history`` annotations: those are written by ``mark_archived``,
    # and the hint should not go blank because an enrichment pass has not run.
    spans: dict[int, tuple[str, str]] = {}
    for sid, ver in session.execute(
        select(Name.sequence_id, Source.version)
        .join(Source, Source.id == Name.source_id)
        .where(
            Name.sequence_id.in_([s.id for s in seqs]),
            Source.note.like("historical IMGT%"),
        )
    ).all():
        v = ver or ""
        first, last = spans.get(sid, (v, v))
        spans[sid] = (min(first, v), max(last, v))
    out: list[ArchivedMatch] = []
    for seq in seqs:
        first, last = spans.get(seq.id, (None, None))
        newer = vers.get(seq.id, ([], None))[1]
        out.append(ArchivedMatch(
            uid=seq.uid, name=q,
            first_seen=first, last_seen=last,
            superseded_by=newer.uid if newer else None,
            superseded_by_name=newer.name if newer else None,
            superseded_by_length=newer.length if newer else None,
            superseded_relation=newer.relation if newer else None,
        ))
    return out


def _mark_superseded(session: Session, results: SearchResults) -> SearchResults:
    """Flag hits IMGT has since re-cut at a different extent, so a superseded record
    never sits in a result list looking like the current one.

    Keyed on the hits themselves rather than the caller's rows, so every search entry
    point gets it the same way. Four queries, independent of hit count.
    """
    if not results.hits:
        return results
    rows = list(session.execute(
        select(Sequence).where(Sequence.uid.in_([h.uid for h in results.hits]))
    ).scalars())
    vers = _versions_for(session, rows)
    by_uid = {s.uid: vers.get(s.id, ([], None))[1] for s in rows}
    for h in results.hits:
        newer = by_uid.get(h.uid)
        if newer is not None:
            h.superseded_by = newer.uid
    return results


def build_record(session: Session, seq: Sequence) -> SequenceRecord:
    """Assemble the native + AIRR views of a sequence row."""
    aliases = _aliases(session, seq.id)
    current = [a for a in aliases if not a.historical]
    # Prefer a current designation as the label; a retired allele (current names
    # all gone) falls back to its historical name.
    label = _label(current) if current else _label(aliases)
    gene = label.split("*")[0] if label and "*" in label else None
    species_sci = scientific_name(seq.species)
    for a in aliases:
        a.url = _alias_url(a)
    allele_names = [a for a in current if a.scheme in _ALLELE_SCHEMES]
    other_names = [a for a in current if a.scheme not in _ALLELE_SCHEMES]
    # Former IMGT names = historical names no longer in use (dedup by name).
    current_names = {a.name for a in current}
    former_names: list[Alias] = []
    seen_former: set[str] = set()
    for a in aliases:
        if a.historical and a.name not in current_names and a.name not in seen_former:
            former_names.append(a)
            seen_former.add(a.name)
    prev_versions, newer_version = _record_versions(session, seq)
    # "Retired" must mean the *allele* is gone, not that this record lost its name to
    # a re-cut of the same allele. A superseded record has a live successor, so it is
    # not retired — saying so contradicts the successor's own page.
    retired = not current and newer_version is None
    # Archived alleles (DB flag from mark_archived) carry first/last-seen release
    # dates as history annotations — surface them on the record.
    first_seen = last_seen = None
    if seq.archived:
        hist = dict(session.execute(
            select(SequenceAnnotation.key, SequenceAnnotation.value_text).where(
                SequenceAnnotation.sequence_id == seq.id,
                SequenceAnnotation.domain == "history",
            )
        ).all())
        first_seen, last_seen = hist.get("first_seen"), hist.get("last_seen")
    distinct_names = list(dict.fromkeys(a.name for a in (current or aliases)))
    # External records = records that actually exist elsewhere: the GenBank accession
    # (primary public deposit) and, *only when OGRDB attests this sequence*, a deep link
    # to OGRDB's entry using OGRDB's own name for it. We deliberately do not link alleles
    # OGRDB has never seen: the route answers for any string, so a chip on an allele OGRDB
    # lacks would imply a curated record that isn't there.
    xrefs = _xrefs(session, seq.id) + _accession_xrefs(session, seq.id)
    imgt_names = [a.name for a in current if a.source == "IMGT"]
    if imgt_names and gene and not seq.archived:
        # IMGT's page is per gene, not per allele; name the allele on the chip anyway.
        # Never for an archived allele: IMGT removed it, so a chip would advertise a
        # record that is no longer there (the former names are listed separately).
        xrefs.append(Xref(label="IMGT", curie=sorted(imgt_names)[0],
                          url=_imgt_gene_url(gene, species_sci)))
    ogrdb_names = [a.name for a in current if a.source == "OGRDB"]
    if ogrdb_names:
        name = sorted(ogrdb_names)[0]
        ogrdb_link = ogrdb_sequence_url(seq.species, name)
        if ogrdb_link:
            xrefs.append(Xref(label="OGRDB", curie=name, url=ogrdb_link))
    native = NativeRecord(
        uid=seq.uid,
        species=seq.species,
        category=seq.category,
        segment=seq.segment,
        locus=seq.locus,
        secondary_locus=seq.secondary_locus,
        subtype=seq.subtype,
        sequence=seq.sequence,
        sequence_hash=seq.sequence_hash,
        sequence_aa=seq.sequence_aa,
        aa_hash=seq.aa_hash,
        aa_origin=seq.aa_origin,
        protein_uid=seq.protein_uid,
        functional=seq.functional,
        length=seq.length,
        has_ambiguity=seq.has_ambiguity,
        allele_names=allele_names,
        other_names=other_names,
        allele_name_groups=_group_aliases(allele_names),
        other_name_groups=_group_aliases(other_names),
        former_names=former_names,
        retired=retired,
        prev_versions=prev_versions,
        newer_version=newer_version,
        archived=seq.archived,
        similarity_note=similarity_reason(seq.segment),
        composition=_composition(session, seq.id) if seq.segment == "C" else None,
        first_seen=first_seen,
        last_seen=last_seen,
        xrefs=_dedup_xrefs(xrefs),
        references=_references(session, seq.id),
        rss=_rss_links(session, seq.id),
        leaders=_leader_links(session, seq.id),
        sources=_source_tags(session, seq.id),
    )
    airr = AlleleDescription(
        label=label,
        aliases=[n for n in distinct_names if n != label],
        coding_sequence=seq.sequence,
        sequence_aa=seq.sequence_aa,
        functional=seq.functional,
        locus=seq.locus,
        sequence_type=seq.segment,
        species=_ontology(seq.species),
        germline_set_ref=_germline_set_ref(session, seq.id),
    )
    return SequenceRecord(native=native, airr=airr)


def _hit(
    session: Session, seq: Sequence, score: float, matched_by: str,
    prefer: str | None = None,
) -> SearchHit:
    aliases = _aliases(session, seq.id)
    return SearchHit(
        uid=seq.uid,
        label=_label(aliases, prefer),
        species=seq.species,
        category=seq.category,
        locus=seq.locus,
        sequence_type=seq.segment,
        length=seq.length,
        score=round(score, 4),
        matched_by=matched_by,
        aliases=list(dict.fromkeys(a.name for a in aliases)),
    )


# --- resolve --------------------------------------------------------------

def resolve(
    session: Session, query: str, species: str | None = None
) -> ResolveResult | None:
    """Resolve a name, UID, nt sequence, or AA sequence to one record.

    A name or sequence can be reused across species (identity is per-species, so
    the *same* nucleotide sequence in two species is two records, and the *same*
    allele name is regularly shared between species). When ``species`` is given,
    resolution is scoped to it, which is the reliable way to pin an ambiguous
    query to the intended species. Without it, the pick is deterministic and
    prefers a current (non-archived) record over an archived one, but a
    cross-species name still resolves to only one of its species; use
    :func:`resolve_exact_name` (exposed as ``/api/search/name?exact=1``) to list
    every species that carries a name.
    """
    q = query.strip()
    if not q:
        return None
    cond = _species_cond(_species_tags(species))

    def _pick(stmt):
        if cond is not None:
            stmt = stmt.where(cond)
        # A HUSA-only sequence is not on the site: resolve by any key returns nothing.
        stmt = _visible_seqs(stmt)
        # Prefer a current record over an archived one, then order stably so the
        # choice is deterministic when the key is shared across species.
        stmt = stmt.order_by(
            Sequence.archived.asc(), Sequence.species.asc(), Sequence.uid.asc()
        )
        return session.scalar(stmt)

    seq = _pick(select(Sequence).where(Sequence.uid == q))
    if seq is not None:
        return ResolveResult(matched_by="uid", record=build_record(session, seq))

    seq = _pick(select(Sequence).where(Sequence.protein_uid == q))
    if seq is not None:
        return ResolveResult(matched_by="protein_uid", record=build_record(session, seq))

    seq = _pick(_unshelved_names(
        select(Sequence)
        .join(Name, Name.sequence_id == Sequence.id)
        .where(Name.name == q)
    ))
    if seq is not None:
        return ResolveResult(matched_by="name", record=build_record(session, seq))

    # As a nucleotide sequence (hash, then reverse-complement).
    try:
        digest = hash_sequence(q)
    except SequenceError:
        digest = None
    if digest is not None:
        seq = _pick(select(Sequence).where(Sequence.sequence_hash == digest))
        if seq is not None:
            return ResolveResult(matched_by="sequence_hash", record=build_record(session, seq))
        seq = _pick(select(Sequence).where(Sequence.rc_hash == digest))
        if seq is not None:
            return ResolveResult(matched_by="rc_hash", record=build_record(session, seq))

    # As an amino-acid sequence.
    try:
        aah = hash_aa(q)
    except SequenceError:
        aah = None
    if aah is not None:
        seq = _pick(select(Sequence).where(Sequence.aa_hash == aah))
        if seq is not None:
            return ResolveResult(matched_by="aa_hash", record=build_record(session, seq))

    return None


def list_aliases(session: Session, uid: str) -> list[Alias]:
    """All aliases for the sequence with the given UID."""
    seq = session.scalar(select(Sequence).where(Sequence.uid == uid))
    if seq is None:
        return []
    return _aliases(session, seq.id)


# --- search ---------------------------------------------------------------

def _species_tag(species: str | None) -> str | None:
    return normalize_species(species) if species else None


def _species_tags(species: str | None) -> list[str] | None:
    """Normalise a species selection to a list of tags, or None for "any".

    Accepts one species or several comma-separated (``human,rhesus``). Each is
    normalised through :func:`normalize_species`, so aliases and common names work
    exactly as they do for a single value, and an unknown species still raises
    ``SpeciesError`` rather than being silently dropped — a filter that quietly
    ignores what you asked for is worse than one that refuses.

    Order is preserved and duplicates removed. Empty or whitespace-only entries are
    skipped, so a trailing comma is harmless; a selection that reduces to nothing is
    None, i.e. unfiltered.
    """
    if not species:
        return None
    tags = list(dict.fromkeys(
        normalize_species(part) for part in species.split(",") if part.strip()
    ))
    return tags or None


def _species_cond(tags: list[str] | None):
    """WHERE condition restricting ``Sequence.species`` to ``tags`` (or None)."""
    if not tags:
        return None
    return Sequence.species == tags[0] if len(tags) == 1 else Sequence.species.in_(tags)


def _source_exists(source: str):
    """A correlated ``EXISTS`` keeping sequences contributed by ``source``.

    A sequence is kept if it carries **at least one** name (:class:`Name`) OR one
    record (:class:`SequenceSource`) from a source whose ``name`` matches ``source``
    case-insensitively. An unknown source name matches nothing, so the result is
    simply empty -- the same fail-soft behaviour as an unknown ``locus`` (no error).
    """
    name_match = func.lower(Source.name) == source.strip().lower()
    by_name = (
        select(Name.id)
        .join(Source, Source.id == Name.source_id)
        .where(Name.sequence_id == Sequence.id, name_match)
    )
    by_rec = (
        select(SequenceSource.id)
        .join(Source, Source.id == SequenceSource.source_id)
        .where(SequenceSource.sequence_id == Sequence.id, name_match)
    )
    return or_(by_name.exists(), by_rec.exists())


def _advanced_conds(functional: bool | None, source: str | None) -> list:
    """Combined WHERE conditions for the functionality + source filters.

    ``Sequence.functional`` is a tri-state boolean (``True``=functional/AIRR ``F``,
    ``False``=non-functional, ``NULL``=unknown). ``functional`` unset -> no filter;
    an explicit ``True``/``False`` matches only that value, so NULL/unknown rows are
    excluded by either. The raw AIRR token is preserved as a ``functionality``
    annotation; this filter is on the boolean only.
    """
    conds: list = []
    if functional is not None:
        conds.append(Sequence.functional.is_(functional))
    if source and source.strip():
        conds.append(_source_exists(source))
    return conds


def _apply_advanced(
    session: Session, seqs: list[Sequence],
    functional: bool | None, source: str | None,
) -> list[Sequence]:
    """Restrict an already-materialised list to rows passing the advanced filters.

    Used by name search, where results come back as ORM objects rather than a live
    statement; the DB does the source ``EXISTS`` in one query, order preserved.
    """
    conds = _advanced_conds(functional, source)
    if not conds:
        return seqs
    ids = [s.id for s in seqs]
    if not ids:
        return seqs
    keep = set(session.scalars(
        select(Sequence.id).where(Sequence.id.in_(ids), *conds)
    ))
    return [s for s in seqs if s.id in keep]


def list_sources(session: Session) -> list[str]:
    """Distinct source names, sorted -- populates the advanced-search source picker.

    Names collapse across versions (a source's historical snapshots share its name),
    so this is the set a user can filter on. Cheap; safe per-request like
    :func:`species_options`.
    """
    stmt = select(Source.name).distinct().order_by(Source.name)
    if _SHELVED:
        stmt = stmt.where(_not_shelved_src())
    return list(session.scalars(stmt))


def _partial_identity(a: str, b: str) -> float:
    """Local, containment-tolerant % identity: the identity of the shorter sequence
    against its best-matching **window** of the longer one.

    rapidfuzz finds the optimal window (Indel alignment); we then report the exact
    edit-distance (Levenshtein) identity over that window — so a partial/truncated
    query fully contained in a longer allele scores ~1.0 instead of being penalised
    for the uncovered flanks. Contrast the default (global) identity, which counts
    the missing flanks as edits. Used only in the opt-in "partial / contained" mode.
    """
    if not a or not b:
        return 0.0
    al = partial_ratio_alignment(a, b)
    if al is None:
        return 0.0
    if len(a) <= len(b):
        short, window = a, b[al.dest_start:al.dest_end]
    else:
        short, window = b, a[al.src_start:al.src_end]
    if not window:
        return 0.0
    return Levenshtein.normalized_similarity(short, window)


def species_options(session: Session) -> list[dict]:
    """Distinct current coding-allele species, each with its common-name label —
    for the search-bar species picker. Cheap; safe to call on every request (unlike
    the full :func:`dataset_stats`, which only runs on the landing page)."""
    rows = session.execute(
        select(Sequence.species)
        .where(Sequence.category == "coding", Sequence.archived.is_(False))
        .distinct()
    ).scalars().all()
    return sorted(
        ({"tag": sp, "label": common_name(sp)} for sp in rows),
        key=lambda s: s["label"],
    )


def _like_escape(s: str) -> str:
    """Escape ILIKE wildcards so a name fragment is matched literally."""
    return s.replace("\\", "\\\\").replace("%", "\\%").replace("_", "\\_")


def _name_ilike(
    session: Session, pattern: str, tags: list[str] | None
) -> list[Sequence]:
    """Distinct non-archived sequences whose name matches an ILIKE ``pattern``."""
    stmt = _unshelved_names(
        select(Sequence)
        .join(Name, Name.sequence_id == Sequence.id)
        .where(Name.name.ilike(pattern), Sequence.archived.is_(False))
        .distinct()
    )
    cond = _species_cond(tags)
    if cond is not None:
        stmt = stmt.where(cond)
    return list(session.scalars(stmt))


def _subgroup_ilike(
    session: Session, subgroup: str, tags: list[str] | None
) -> list[Sequence]:
    """Every non-archived sequence whose name is in ``subgroup``.

    Two patterns, both anchored, never a bare prefix: ``IGHV1-%`` for the ordinary
    ``IGHV1-69*01`` shape and ``IGHV1*%`` for genes named without a gene part
    (``IGHV6*01``). ``IGHV1%`` would be wrong — it swallows ``IGHV10-1*01``, which is
    the same boundary error the 2026-07 name-search rework removed from gene search.
    """
    stmt = _unshelved_names(
        select(Sequence)
        .join(Name, Name.sequence_id == Sequence.id)
        .where(
            or_(
                Name.name.ilike(f"{_like_escape(subgroup)}-%"),
                Name.name.ilike(f"{_like_escape(subgroup)}*%"),
            ),
            Sequence.archived.is_(False),
        )
        .distinct()
    )
    cond = _species_cond(tags)
    if cond is not None:
        stmt = stmt.where(cond)
    return list(session.scalars(stmt))


def search_by_name(
    session: Session,
    q: str,
    cap: int = DEFAULT_CAP,
    species: str | None = None,
    exact: bool = False,
    partial: bool = False,
    functional: bool | None = None,
    source: str | None = None,
    subgroup: bool = False,
) -> SearchResults:
    """Look up records by allele or gene name.

    Default (precise): an exact allele/name match, or -- for a bare gene name --
    the alleles of exactly that gene, matched to the allele boundary (``*``) so
    that ``IGHV1-69`` returns only its own alleles and ``IGHV1-6`` does not
    silently pull ``IGHV1-69``/``IGHV1-68``/... The query resolves to a real
    allele or gene, or to nothing: no fuzzy guess.

    ``exact=True``: only records carrying the query as an exact allele name
    (across species unless scoped) -- lists every species that reuses a name.

    ``partial=True``: browse mode. Every allele whose name *contains* the query,
    so a fragment (``IGHV1-6``, ``1-69``) returns all variations to refine among.

    ``subgroup=True``: treat the query as an IMGT subgroup (``IGHV1``) and return
    every allele in it. A bare subgroup is not a parseable name, so this cannot be
    reached by the ordinary paths. Returns nothing for a query that is not a
    subgroup — ``IGHV1-69`` is a gene, and reinterpreting it here would undo the
    precision the default path exists to provide.

    The default path also *falls back* to a subgroup lookup, but only when the query
    resolved to nothing at all and is itself a bare subgroup. It can therefore never
    widen a query that already matched something.

    Name matches are not sequence identity: the score is 1.0 and the caller must
    not present it as identity. Real % identity comes from the sequence/neighbour
    queries.
    """
    q = q.strip()
    tags = _species_tags(species)

    if subgroup:
        sg = as_subgroup(q)
        seqs = _subgroup_ilike(session, sg, tags) if sg else []
    elif exact:
        # Archived alleles stay reachable by UID / resolve, not through search.
        seqs = [s for s in resolve_exact_name(session, q, species) if not s.archived]
    elif partial:
        seqs = _name_ilike(session, f"%{_like_escape(q)}%", tags)
    else:
        seqs = [s for s in resolve_exact_name(session, q, species) if not s.archived]
        if not seqs and "*" not in q:
            # Bare gene name: its alleles, bounded at the allele separator.
            seqs = _name_ilike(session, f"{_like_escape(q)}*%", tags)
        if not seqs:
            # Still nothing, and the query is a bare subgroup: its alleles.
            sg = as_subgroup(q)
            if sg:
                seqs = _subgroup_ilike(session, sg, tags)

    # Advanced functionality / source filters restrict the resolved records.
    seqs = _apply_advanced(session, seqs, functional, source)
    seqs = sorted(seqs, key=lambda s: s.uid)
    capped = len(seqs) > cap
    hits = [_hit(session, s, 1.0, "name", prefer=q) for s in seqs[:cap]]
    # Nothing current matched: say what was skipped rather than return silence.
    skipped = _archived_matches(session, q, species) if not hits else []
    return _mark_superseded(session, SearchResults(
        query=q, count=len(hits), capped=capped, hits=hits, archived_matches=skipped,
    ))


def _similarity_search(
    session: Session, query: str, column, value_of, cap: int,
    tags: list[str] | None = None, category: str = "coding", score_cutoff: float = 0.0,
    partial: bool = False, locus: str | None = None, segment: str | None = None,
    functional: bool | None = None, source: str | None = None,
) -> list[tuple[Sequence, float, str]]:
    """Shared exact + trigram-candidate + rapidfuzz ranking over a text column.

    ``partial=False`` (default): **global** normalised-Levenshtein identity — a
    partial/contained sequence is penalised for the uncovered length (a 200 nt
    fragment of a 300 nt allele scores ~0.67). ``partial=True``: **local**
    identity over the best-matching window (:func:`_partial_identity`) — a
    contained fragment scores ~1.0. Opt-in via the search UI's match mode.

    ``segment`` (V/D/J/C) scopes the pool to one segment. It is passed only when the
    query's own segment is known — i.e. the search starts from a stored allele — and
    stays None both for a pasted raw sequence (segment unknown) and for leader/RSS,
    whose rows carry ``segment`` NULL and are scoped by ``category`` instead.

    For a segment in ``SIMILARITY_EQUAL_LENGTH`` (D) the pool is further restricted to
    candidates of the **query's exact length** and scored with Hamming. Both halves
    matter: the length predicate is what makes a positional score well-defined, and
    scoring the restricted pool with Levenshtein would still let gaps in.
    """
    # Gather candidates with a trigram pre-filter, ordered before the cap so a loose
    # threshold can't let an unordered LIMIT drop genuinely high-identity candidates.
    # The pool is decoupled from ``cap`` (MIN_SIMILARITY_POOL floor): a small cap
    # must not shrink the candidate set and miss the correct germline.
    pool = max(cap * 5, MIN_SIMILARITY_POOL)
    conds = [Sequence.category == category, column.is_not(None),
             Sequence.archived.is_(False)]
    vis = _visible_cond()  # exclude shelved-only sequences from the candidate pool
    if vis is not None:
        conds.append(vis)
    if segment:
        # DNA has only 64 possible trigrams, so the `%` prefilter barely separates
        # segments at all: a J query trigram-matches 11,276 of the 11,279 V alleles
        # (V-vs-J similarity ~0.47, well over the 0.3 threshold). V outnumbers J ~5:1
        # and so both crowds the bounded pool and lands in the results — a J allele's
        # "neighbours" came back V-only. They are not comparable sequences anyway
        # (the highest V-vs-J identity over all 21.3M pairs is 0.60), so scope the
        # pool whenever the query's segment is known.
        #
        # Scoped on the comparison class, not the raw segment, so the same guard
        # covers the non-coding categories: every RSS carries segment NULL, so
        # without this a 5'D-RS would be ranked against the equal-length V-RS of an
        # IGKV. The two are the same length and different elements.
        conds.append(COMPARISON_CLASS == segment)
        if equal_length_only(segment):
            # Same length, in the database, before the pool cap — not filtered after
            # ranking. Otherwise differently-sized candidates would occupy the bounded
            # pool and push out the same-length alleles that are the whole point.
            conds.append(func.length(column) == len(query))
    if partial:
        # Containment search: WHOLE-STRING trigram similarity is length-sensitive —
        # a short fragment scores low against a long allele (a 57 nt piece of a
        # 296 nt allele ~0.19) and never enters the pool. Rank by WORD similarity
        # instead (how well the query matches a *portion* of each sequence).
        #
        # We deliberately do NOT gate on the `<%` boolean operator. Its threshold
        # (pg_trgm.word_similarity_threshold, default 0.6) drops a *divergent* query
        # wholesale: a somatically mutated antibody window can word-match its own
        # germline at only ~0.4, so `<%` returns an empty pool and the search yields
        # no hit at all. Ordering without a threshold always surfaces the best
        # candidates; CPU stays bounded by the pool LIMIT (and, at the API layer, by
        # a concurrency limiter). `_partial_identity` below then rescores the pool by
        # exact local identity, so the loose trigram order only has to pull the
        # parent allele *into* the pool, not rank it first.
        order = func.word_similarity(query, column)
    else:
        conds.append(column.op("%")(query))
        order = func.similarity(column, query)
    stmt = select(Sequence).where(*conds)
    sp_cond = _species_cond(tags)
    if sp_cond is not None:
        stmt = stmt.where(sp_cond)
    if locus:
        # Scope to a known chain; a shared TRA/TRD gene matches under either locus.
        stmt = stmt.where(or_(Sequence.locus == locus, Sequence.secondary_locus == locus))
    # Advanced filters (functionality + source) narrow the candidate pool itself, so
    # a filtered-out allele never displaces a genuine hit within the bounded LIMIT.
    for cond in _advanced_conds(functional, source):
        stmt = stmt.where(cond)
    # A bounded ranked pool keeps scoring cheap; the displayed top is well within it.
    cands = list(session.scalars(stmt.order_by(order.desc()).limit(pool)))
    strings = [value_of(s) or "" for s in cands]
    if partial:
        # Local containment-tolerant identity: best-window edit-distance identity
        # per candidate (each pair is two C calls — fast enough over the pool).
        cut = max(0.0, score_cutoff)
        scored = sorted(
            ((s, _partial_identity(query, s), i) for i, s in enumerate(strings) if s),
            key=lambda t: t[1], reverse=True,
        )
        scored = [t for t in scored if t[1] >= cut]
    else:
        # Score the whole pool in one vectorised C call instead of a Python loop.
        # Global normalised Levenshtein (true % identity). score_cutoff lets
        # rapidfuzz prune the edit-distance band and skip pairs that cannot reach
        # the threshold — a large speedup for high-identity search.
        scored = process.extract(
            query, strings, scorer=similarity_scorer(segment), limit=None,
            score_cutoff=max(0.0, score_cutoff),
        )  # -> [(string, score, index), ...] sorted by score desc, ≥ cutoff only
    ranked: list[tuple[Sequence, float, str]] = []
    for stored, score, idx in scored:
        if not stored:
            continue
        seq = cands[idx]
        if stored == query:
            ranked.append((seq, 1.0, "exact"))
        elif query in stored or stored in query:
            ranked.append((seq, score, "substring"))
        else:
            ranked.append((seq, score, "similar"))
    return ranked


def search_by_sequence(
    session: Session, seq_str: str, cap: int = DEFAULT_CAP,
    species: str | None = None, category: str = "coding", partial: bool = False,
    locus: str | None = None, functional: bool | None = None,
    source: str | None = None,
) -> SearchResults:
    """Nucleotide search within ``category`` (coding / leader / rss): exact hash,
    sub/superstring, and edit-distance ranking. ``partial=True`` scores by local
    (best-window) identity so a truncated query matches its allele near 100%.
    ``locus`` restricts the search to one chain (IGH/IGK/IGL/TRA/...)."""
    from .sequence import normalize

    try:
        norm = normalize(seq_str)
    except SequenceError as exc:
        raise ValueError(str(exc)) from exc

    tags = _species_tags(species)
    loc = locus.upper() if locus else None
    # Hash-first: a pasted sequence that *is* a stored allele resolves instantly via
    # the indexed sequence_hash — no edit-distance scan. Only fall through to fuzzy
    # identity ranking when the sequence is genuinely not in the database.
    digest = hash_sequence(seq_str)
    exact_stmt = _visible_seqs(select(Sequence).where(
        Sequence.sequence_hash == digest, Sequence.category == category,
        Sequence.archived.is_(False),
    ))
    sp_cond = _species_cond(tags)
    if sp_cond is not None:
        exact_stmt = exact_stmt.where(sp_cond)
    if loc:
        exact_stmt = exact_stmt.where(
            or_(Sequence.locus == loc, Sequence.secondary_locus == loc))
    for cond in _advanced_conds(functional, source):
        exact_stmt = exact_stmt.where(cond)
    exact = list(session.scalars(exact_stmt))  # ≥1: same sequence can span species
    if exact:
        hits = [_hit(session, s, 1.0, "exact") for s in exact[:cap]]
        return _mark_superseded(session, SearchResults(
            query=f"<{len(norm)} nt>", count=len(hits),
            capped=len(exact) > cap, hits=hits,
        ))
    ranked = _similarity_search(
        session, norm, Sequence.sequence, lambda s: s.sequence, cap, tags,
        category=category, partial=partial, locus=loc,
        functional=functional, source=source,
    )
    hits = [_hit(session, s, score, mb) for s, score, mb in ranked[:cap]]
    return _mark_superseded(session, SearchResults(
        query=f"<{len(norm)} nt>", count=len(hits), capped=len(ranked) > cap, hits=hits
    ))


def search_by_aa(
    session: Session, seq_aa: str, cap: int = DEFAULT_CAP, species: str | None = None,
    partial: bool = False, locus: str | None = None,
    functional: bool | None = None, source: str | None = None,
) -> SearchResults:
    """Amino-acid search: exact aa_hash, then fuzzy on sequence_aa. ``partial=True``
    scores by local (best-window) identity for truncated AA queries. ``locus``
    restricts the search to one chain (IGH/IGK/IGL/TRA/...)."""
    from .sequence import normalize_aa

    try:
        norm = normalize_aa(seq_aa)
    except SequenceError as exc:
        raise ValueError(str(exc)) from exc

    tags = _species_tags(species)
    sp_cond = _species_cond(tags)
    loc = locus.upper() if locus else None
    ranked = _similarity_search(
        session, norm, Sequence.sequence_aa, lambda s: s.sequence_aa, cap, tags,
        partial=partial, locus=loc, functional=functional, source=source,
    )
    digest = hash_aa(seq_aa)
    if not any(s.aa_hash == digest for s, _, _ in ranked):
        stmt = _visible_seqs(select(Sequence).where(
            Sequence.aa_hash == digest, Sequence.archived.is_(False)
        ))
        if sp_cond is not None:
            stmt = stmt.where(sp_cond)
        if loc:
            stmt = stmt.where(or_(Sequence.locus == loc, Sequence.secondary_locus == loc))
        for cond in _advanced_conds(functional, source):
            stmt = stmt.where(cond)
        exact = session.scalar(stmt.limit(1))
        if exact is not None:
            ranked.insert(0, (exact, 1.0, "exact"))
    hits = [_hit(session, s, score, mb) for s, score, mb in ranked[:cap]]
    return _mark_superseded(session, SearchResults(
        query=f"<{len(norm)} aa>", count=len(hits), capped=len(ranked) > cap, hits=hits
    ))


def length_distribution(
    session: Session, n_buckets: int = 24, species: str | None = None,
    segment: str | None = None,
) -> dict | None:
    """Histogram of coding-allele length → count, split per species (overview chart).

    With ``species`` the histogram is restricted to that one species (own scale);
    otherwise it spans every species, split by colour. ``segment`` restricts it to
    one segment: mixing V/D/J/C on one axis puts unrelated length scales together.
    """
    tag = _species_tag(species)
    stmt = (
        select(Sequence.length, Sequence.species, func.count())
        .where(Sequence.category == "coding")
        .group_by(Sequence.length, Sequence.species)
    )
    if tag is not None:
        stmt = stmt.where(Sequence.species == tag)
    if segment is not None:
        stmt = stmt.where(Sequence.segment == segment)
    rows = session.execute(stmt).all()
    if not rows:
        return None
    lengths = [r[0] for r in rows]
    lo, hi = min(lengths), max(lengths)
    total = sum(r[2] for r in rows)
    species = sorted({r[1] for r in rows})
    span = max(hi - lo, 1)
    nb = min(n_buckets, span + 1)
    step = max(1, -(-span // nb))  # ceil division → bucket width in nt
    nb = -(-span // step) + 1
    buckets = [
        {"lo": lo + k * step, "hi": min(lo + (k + 1) * step - 1, hi),
         "count": 0, "by_species": {}}
        for k in range(nb)
    ]
    for length, sp, c in rows:
        b = buckets[(length - lo) // step]
        b["count"] += c
        b["by_species"][sp] = b["by_species"].get(sp, 0) + c
    max_count = max((b["count"] for b in buckets), default=1)
    return {
        "buckets": buckets, "max_count": max_count, "min_len": lo, "max_len": hi,
        "total": total, "species": species,
    }


def source_concordance(session: Session, species: str | None = None) -> dict | None:
    """Pairwise source co-occurrence over coding sequences (overview matrix): for
    each pair of sources, how many coding sequences carry a name from both (the
    diagonal is each source's own sequence count).

    With ``species`` the matrix is restricted to one species — the only setting in
    which concordance is meaningful, since sources from different species never
    share a sequence.
    """
    tag = _species_tag(species)
    stmt = (
        select(Source.name, Name.sequence_id)
        .join(Name, Name.source_id == Source.id)
        .join(Sequence, Sequence.id == Name.sequence_id)
        .where(
            Sequence.category == "coding",
            Sequence.archived.is_(False),  # archived alleles aren't current coverage
            # count current names only — historical release snapshots are aliases,
            # not current coverage, so they must not inflate a source's totals.
            or_(Source.note.is_(None), Source.note.notlike("historical IMGT%")),
            *((_not_shelved_src(),) if _SHELVED else ()),  # no shelved source row/col
        )
        .distinct()
    )
    if tag is not None:
        stmt = stmt.where(Sequence.species == tag)
    rows = session.execute(stmt).all()
    if not rows:
        return None
    seqs_by_src: dict[str, set[int]] = {}
    for name, seq_id in rows:
        seqs_by_src.setdefault(name, set()).add(seq_id)
    sources = sorted(seqs_by_src)
    totals = {s: len(seqs_by_src[s]) for s in sources}
    matrix = [
        [len(seqs_by_src[a] & seqs_by_src[b]) for b in sources]
        for a in sources
    ]
    off_diag = [
        matrix[i][j]
        for i in range(len(sources)) for j in range(len(sources)) if i != j
    ]
    return {
        "sources": sources, "totals": totals, "matrix": matrix,
        "max": max(off_diag, default=0) or 1,
    }


def conservation_profile(
    session: Session, uids, kind: str = "nt", max_cols: int = 400
) -> dict | None:
    """Per-position base/residue composition across a set of sequences (by UID).

    Trivial alignment: positions are compared directly over the **modal length**
    (sequences of a different length are excluded and reported), so no MSA is
    needed — fine for the high-identity sets the similarity pages produce.
    """
    column = Sequence.sequence_aa if kind == "aa" else Sequence.sequence
    seqs = [
        s for s in session.scalars(
            select(column).where(Sequence.uid.in_(list(uids)))
        ) if s
    ]
    if len(seqs) < 2:
        return None
    modal = Counter(len(s) for s in seqs).most_common(1)[0][0]
    aligned = [s for s in seqs if len(s) == modal]
    n = len(aligned)
    cols = []
    for i in range(min(modal, max_cols)):
        counts = Counter(s[i] for s in aligned)
        # (residue, fraction, count) — count feeds the per-column hover tooltip.
        cols.append([(ch, round(c / n, 4), c) for ch, c in counts.most_common()])
    return {
        "kind": kind, "n_aligned": n, "n_total": len(seqs),
        "length": modal, "cols": cols,
    }


def dataset_stats(session: Session) -> dict:
    """Per-species counts + contributing sources for the landing-page overview.

    Counts **current** alleles only; archived (historical-only) alleles and the
    historical release snapshots are excluded from coverage and reported separately.
    """
    # Shelved-only sequences (HUSA novel alleles) are not on the site, so they are
    # not in the coverage counts either. `*_vis` threads the same gate through every
    # per-species / per-segment / per-locus tally below (all of which use `_current`).
    vis = _visible_cond()
    _vis = (vis,) if vis is not None else ()
    _current = (Sequence.category == "coding", Sequence.archived.is_(False), *_vis)
    coding = session.scalar(select(func.count()).select_from(Sequence).where(*_current))
    archived = session.scalar(
        select(func.count()).select_from(Sequence)
        .where(Sequence.category == "coding", Sequence.archived.is_(True), *_vis)
    )
    rss = session.scalar(
        select(func.count()).select_from(Sequence)
        .where(Sequence.category == "rss", *_vis)
    )
    leader = session.scalar(
        select(func.count()).select_from(Sequence)
        .where(Sequence.category == "leader", *_vis)
    )

    def _by(category: str) -> dict[str, int]:
        return dict(
            session.execute(
                select(Sequence.species, func.count())
                .where(Sequence.category == category, Sequence.archived.is_(False), *_vis)
                .group_by(Sequence.species)
            ).all()
        )

    coding_by, rss_by, leader_by = _by("coding"), _by("rss"), _by("leader")
    archived_by = dict(session.execute(
        select(Sequence.species, func.count())
        .where(Sequence.category == "coding", Sequence.archived.is_(True), *_vis)
        .group_by(Sequence.species)
    ).all())

    # Coding alleles are no longer V-only: J is loaded for 29 species. Count per
    # segment and keep the two apart. The landing page used to print the *combined*
    # coding total under the label "V alleles", which both overstated V and hid J;
    # and the by-locus table mixed a V-only matrix with V+J row totals, so its rows
    # did not add up to their own totals.
    by_segment = dict(session.execute(
        select(Sequence.segment, func.count())
        .where(*_current, Sequence.segment.is_not(None))
        .group_by(Sequence.segment)
    ).all())

    def _seg_by(segment: str) -> dict[str, int]:
        return dict(session.execute(
            select(Sequence.species, func.count())
            .where(*_current, Sequence.segment == segment)
            .group_by(Sequence.species)
        ).all())

    # Every segment present, in biological order, so adding one (C now, D later)
    # extends the page instead of quietly falling outside it. The V/J-specific names
    # below are views onto this, kept because templates and callers already use them.
    segments = [s for s in SEGMENT_ORDER if s in by_segment] + sorted(
        s for s in by_segment if s not in SEGMENT_ORDER
    )
    seg_species: dict[str, dict[str, int]] = {s: _seg_by(s) for s in segments}
    v_by, j_by = seg_species.get("V", {}), seg_species.get("J", {})

    # Per (species, locus) counts, per segment: powers the by-locus overview table and
    # the per-species browse links (IGH/IGK/IGL/TRA…).
    seg_loci_counts: dict[str, dict[str, dict[str, int]]] = {}
    seg_locus_totals: dict[str, dict[str, int]] = {}

    def _tally(rows, counts, totals) -> None:
        for sp, loc, c in rows:
            counts.setdefault(sp, {})[loc] = counts.setdefault(sp, {}).get(loc, 0) + c
            totals[loc] = totals.get(loc, 0) + c

    def _locus_rows(segment: str, column):
        return session.execute(
            select(Sequence.species, column, func.count())
            .where(*_current, Sequence.segment == segment, column.is_not(None))
            .group_by(Sequence.species, column)
        ).all()

    for seg in segments:
        counts: dict[str, dict[str, int]] = {}
        totals: dict[str, int] = {}
        _tally(_locus_rows(seg, Sequence.locus), counts, totals)
        # Shared TRA/TRD genes count under their second locus too (TRAV14/DV4 -> TRD).
        _tally(_locus_rows(seg, Sequence.secondary_locus), counts, totals)
        seg_loci_counts[seg] = counts
        seg_locus_totals[seg] = totals

    loci_counts = seg_loci_counts.get("V", {})
    locus_totals = seg_locus_totals.get("V", {})
    _order = ["IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD"]
    _seen_loci = {loc for t in seg_locus_totals.values() for loc in t}
    all_loci = [x for x in _order if x in _seen_loci] + sorted(
        x for x in _seen_loci if x not in _order
    )

    # current sources (exclude historical release snapshots + shelved sources) per
    # species
    _not_hist = or_(Source.note.is_(None), Source.note.notlike("historical IMGT%"))
    _src_conds = (_not_hist, *((_not_shelved_src(),) if _SHELVED else ()))
    src_by: dict[str, list[str]] = {}
    for sp, name in session.execute(
        select(Sequence.species, Source.name)
        .join(Name, Name.sequence_id == Sequence.id)
        .join(Source, Source.id == Name.source_id)
        .where(*_current, *_src_conds)
        .group_by(Sequence.species, Source.name)
    ).all():
        src_by.setdefault(sp, []).append(name)

    species = sorted(
        (
            {
                "tag": sp,
                "label": common_name(sp),
                "coding": coding_by[sp],
                "v": v_by.get(sp, 0),
                "j": j_by.get(sp, 0),
                "rss": rss_by.get(sp, 0),
                "leader": leader_by.get(sp, 0),
                "archived": archived_by.get(sp, 0),
                "loci": [x for x in all_loci if x in loci_counts.get(sp, {})],
                "loci_counts": loci_counts.get(sp, {}),
                # Segment-keyed views, so a page can render every segment present
                # without naming any of them.
                "by_segment": {s: seg_species[s].get(sp, 0) for s in segments},
                "seg_loci_counts": {
                    s: seg_loci_counts[s].get(sp, {}) for s in segments
                },
                "sources": sorted(src_by.get(sp, [])),
            }
            for sp in coding_by
        ),
        key=lambda s: -s["coding"],
    )
    sources = list(session.scalars(
        select(Source.name).where(*_src_conds).distinct().order_by(Source.name)
    ))
    # When the served data was last written. The only build date we actually hold:
    # sources carry per-release versions (see /sources) but the *current* IMGT
    # reference is ingested without one, so anything called a "release date" here
    # would be invented. This is checkable — it is the load timestamp.
    last_loaded = session.scalar(select(func.max(Sequence.created_at)))
    return {
        "last_loaded": last_loaded,
        # "coding" is the sum over segments; per-segment counts are `by_segment`
        # (they gained C without a code change, and will gain D the same way).
        "coding": coding,
        "by_segment": by_segment,
        "segments": segments,
        "segment_labels": {s: SEGMENT_LABELS.get(s, s) for s in segments},
        "seg_locus_totals": seg_locus_totals,
        "seg_loci_counts": seg_loci_counts,
        # Alleles per locus summed over every segment — the locus chips on the landing
        # page. Summed here, not in the template, so a new segment is included by
        # construction rather than by someone remembering to add a term.
        "locus_all_totals": {
            loc: sum(t.get(loc, 0) for t in seg_locus_totals.values())
            for loc in all_loci
        },
        "archived": archived,
        "rss": rss,
        "leader": leader,
        "sources": sources,
        "species": species,
        "all_loci": all_loci,
        "locus_totals": locus_totals,
    }


def resolve_exact_name(
    session: Session, name: str, species: str | None = None
) -> list[Sequence]:
    """Return every sequence carrying ``name`` as an exact allele name.

    Across species unless ``species`` is given. >1 result means the allele name
    is shared between species (or designations) without a single sequence.
    """
    cond = _species_cond(_species_tags(species))
    stmt = _unshelved_names(
        select(Sequence)
        .join(Name, Name.sequence_id == Sequence.id)
        .where(Name.name == name.strip())
        .distinct()
    )
    if cond is not None:
        stmt = stmt.where(cond)
    return list(session.scalars(stmt))


def _find_similar(
    session: Session, uid: str, column, value_of, kind_label: str,
    min_identity: float, cap: int, same_species: bool,
    target_species: str | None = None,
    functional: bool | None = None, source: str | None = None,
) -> SearchResults | None:
    seq = session.scalar(select(Sequence).where(Sequence.uid == uid))
    if seq is None:
        return None
    qval = value_of(seq)
    if not qval:
        return SearchResults(
            query=f"no {kind_label} sequence for {uid}", count=0, capped=False, hits=[]
        )
    # An explicit target species wins over same_species: restrict the scan to those
    # species (e.g. a rhesus allele -> only its closest *human* orthologs). More than
    # one may be named, which is the point of the selection: "compare me against
    # human and mouse, nothing else".
    if target_species:
        tags = _species_tags(target_species)
    else:
        tags = [seq.species] if same_species else None
    # Search within the target's own category, so a leader finds similar leaders
    # (and RSS finds RSS), not coding alleles — and, for coding, within its own
    # segment, so a J allele finds J rather than the far more numerous V. ``segment``
    # is NULL on leader/RSS rows, which leaves those searches category-scoped only.
    # Pass the threshold so the scorer can prune everything below it.
    ranked = _similarity_search(
        session, qval, column, value_of, cap, tags, category=seq.category,
        score_cutoff=min_identity, segment=seq.segment,
        functional=functional, source=source,
    )
    hits: list[SearchHit] = []
    for s, score, mb in ranked:
        if s.id == seq.id or score < min_identity:
            continue
        hits.append(_hit(session, s, score, mb))
    return SearchResults(
        query=f"{kind_label} similar to {uid}", count=len(hits[:cap]),
        capped=len(hits) > cap, hits=hits[:cap],
    )


def _linked_allele_rollup(
    session: Session, rows, cap: int
) -> tuple[list[RssAllele], list[str]]:
    """Roll (Sequence, source-name) link rows up per allele into the capped allele
    list + the overall source set. Shared by rss_detail and leader_detail."""
    by_allele: dict[int, tuple[Sequence, set[str]]] = {}
    for seq, src_name in rows:
        by_allele.setdefault(seq.id, (seq, set()))[1].add(src_name)
    alleles = [
        RssAllele(
            uid=seq.uid, label=_label(_aliases(session, seq.id)),
            species=seq.species, locus=seq.locus, sources=sorted(srcs),
        )
        for seq, srcs in list(by_allele.values())[:cap]
    ]
    overall = sorted({s for _seq, srcs in by_allele.values() for s in srcs})
    return alleles, overall


def rss_detail(
    session: Session, uid: str, cap: int = DEFAULT_CAP
) -> tuple[Sequence, list[RssAllele], list[str]] | None:
    """For an RSS UID: the RSS sequence, the coding alleles linked to it (each with
    the source(s) of *that* allele->RSS link), and the overall source set."""
    rss = session.scalar(
        select(Sequence).where(Sequence.uid == uid, Sequence.category == "rss")
    )
    if rss is None:
        return None
    rows = session.execute(
        _visible_seqs(
            select(Sequence, Source.name)
            .join(SequenceLink, SequenceLink.from_sequence_id == Sequence.id)
            .join(Source, Source.id == SequenceLink.source_id)
            .where(
                SequenceLink.to_sequence_id == rss.id,
                *([_not_shelved_src()] if _SHELVED else []),
            )
        )
    ).all()
    alleles, overall = _linked_allele_rollup(session, rows, cap)
    return rss, alleles, overall


def leader_detail(
    session: Session, uid: str, cap: int = DEFAULT_CAP
) -> tuple[Sequence, list[RssAllele], list[str], int | None] | None:
    """For a leader UID: the leader sequence, the coding alleles that carry it
    (each with the source(s) of that allele->leader link), the source set, and
    the exon1 length (splice boundary)."""
    leader = session.scalar(
        select(Sequence).where(Sequence.uid == uid, Sequence.category == "leader")
    )
    if leader is None:
        return None
    exon1_len = session.scalar(
        select(SequenceAnnotation.value_int).where(
            SequenceAnnotation.sequence_id == leader.id,
            SequenceAnnotation.domain == "leader",
            SequenceAnnotation.key == "exon1_len",
        ).limit(1)
    )
    rows = session.execute(
        _visible_seqs(
            select(Sequence, Source.name)
            .join(SequenceLink, SequenceLink.from_sequence_id == Sequence.id)
            .join(Source, Source.id == SequenceLink.source_id)
            .where(
                SequenceLink.to_sequence_id == leader.id,
                SequenceLink.relation == "has_leader",
                *([_not_shelved_src()] if _SHELVED else []),
            )
        )
    ).all()
    alleles, overall = _linked_allele_rollup(session, rows, cap)
    return leader, alleles, overall, exon1_len


def find_similar_rss(
    session: Session, uid: str, min_identity: float = 0.0,
    cap: int = DEFAULT_CAP, same_species: bool = True,
) -> RssSimilar | None:
    """RSS sequences ranked by identity to the RSS named by ``uid`` (excl. self)."""
    rss = session.scalar(
        select(Sequence).where(Sequence.uid == uid, Sequence.category == "rss")
    )
    if rss is None:
        return None
    tags = [rss.species] if same_species else None
    ranked = _similarity_search(
        session, rss.sequence, Sequence.sequence, lambda s: s.sequence, cap,
        tags, category="rss",
        # Rank a V-RS against V-RSs only, a 5'D-RS against 5'D-RSs only. Both are
        # 28 nt where the spacer is 12, so without this they are indistinguishable
        # to a pairwise score and would come back as each other's best neighbours.
        segment=comparison_class(rss.segment, rss.subtype),
    )
    hits: list[RssHit] = []
    for s, score, _mb in ranked:
        if s.id == rss.id or score < min_identity:
            continue
        n = session.scalar(
            select(func.count(func.distinct(SequenceLink.from_sequence_id)))
            .where(SequenceLink.to_sequence_id == s.id)
        ) or 0
        hits.append(RssHit(
            uid=s.uid, sequence=s.sequence, species=s.species,
            subtype=s.subtype, score=round(score, 4), allele_count=n,
        ))
    return RssSimilar(
        uid=rss.uid, sequence=rss.sequence, count=len(hits[:cap]),
        capped=len(hits) > cap, hits=hits[:cap],
    )


def _neighbors_from_table(
    session: Session, uid: str, kind: str, min_identity: float, cap: int,
    species: list[str] | None = None, floor: float = NEIGHBOR_FLOOR,
    functional: bool | None = None, source: str | None = None,
    segment: str | None = None,
) -> SearchResults | None:
    """Serve from the precomputed cache, or None (→ caller computes live).

    The cache stores same-species neighbours >= ``NEIGHBOR_FLOOR`` and cross-species
    >= ``CROSS_FLOOR``, each capped above the observed max so the cap never binds.
    So for a threshold >= the relevant ``floor`` (0.75 same-species, 0.85 for an
    any-species read) the cache is complete — filter + limit straight in SQL (no
    full-list scan, no per-uid cap check). ``species`` restricts to that species
    (same-species search); None returns all species (cross-species). Falls to live
    only below ``floor`` or for an uid that isn't in the cache at all.

    ``segment`` restricts to neighbours of the query's own segment, matching the
    live path. The builder blocks on segment so it never writes a cross-segment
    edge, but a cache built before that (or a floor lowered below the ~0.60 V-vs-J
    ceiling) could contain one, and the two paths must not disagree.
    """
    if min_identity < floor - 1e-9:
        return None  # below this read's precomputed floor → live
    # Cache is keyed on integer sequences.id: resolve the query uid once, then filter
    # + order + limit in SQL on the (from_id, kind, identity) index. The denormalised
    # to_uid/label/species/locus/type come straight off the row (no per-hit join).
    from_id = session.scalar(select(Sequence.id).where(Sequence.uid == uid))
    if from_id is None:
        return None  # unknown uid → live (to be safe)
    q = (
        select(
            SequenceNeighbor.to_uid, SequenceNeighbor.identity,
            SequenceNeighbor.label, SequenceNeighbor.to_species,
            SequenceNeighbor.to_category, SequenceNeighbor.to_locus,
            SequenceNeighbor.to_segment, Sequence.length,
        )
        # exclude archived (historical-only) alleles from the surfaced neighbours;
        # an archived allele's own page still shows its similar *current* alleles.
        .join(Sequence, Sequence.id == SequenceNeighbor.to_id)
        .where(
            SequenceNeighbor.from_id == from_id, SequenceNeighbor.kind == kind,
            SequenceNeighbor.identity >= min_identity,
            Sequence.archived.is_(False),
        )
    )
    if species:
        q = q.where(
            SequenceNeighbor.to_species == species[0] if len(species) == 1
            else SequenceNeighbor.to_species.in_(species)
        )
    if segment is not None:
        q = q.where(SequenceNeighbor.to_segment == segment)
    # A cached edge may point at a shelved-only neighbour (the builder does not
    # exclude them); drop it on read so the cache can never surface HUSA. Correlates
    # on the ``to`` Sequence joined above, exactly like the advanced source filter.
    q = _visible_seqs(q)
    # Advanced filters apply to the neighbour (``to``) sequence, which is already
    # joined above; the source EXISTS correlates on that same ``Sequence``.
    for cond in _advanced_conds(functional, source):
        q = q.where(cond)
    rows = session.execute(
        q.order_by(SequenceNeighbor.identity.desc()).limit(cap)
    ).all()
    if not rows and session.scalar(
        select(SequenceNeighbor.from_id)
        .where(SequenceNeighbor.from_id == from_id, SequenceNeighbor.kind == kind)
        .limit(1)
    ) is None:
        return None  # uid not in the cache at all → live (to be safe)
    hits = [
        SearchHit(
            uid=r[0], label=r[2], species=r[3], category=r[4], locus=r[5],
            sequence_type=r[6], length=r[7], score=round(r[1], 4),
            matched_by="cached", aliases=[],
        )
        for r in rows
    ]
    return SearchResults(
        query=f"{kind} similar to {uid}", count=len(hits), capped=False, hits=hits
    )


def _neighbor_scope(
    sp_q: str | None, cat_q: str | None, kind: str,
    same_species: bool, target_species: str | None,
) -> tuple[list[str] | None, float, bool]:
    """Resolve a neighbour lookup to ``(species_tags, floor, cache_ok)``.

    The cache holds same-species edges >= ``NEIGHBOR_FLOOR`` (0.75) and cross-species
    edges >= ``CROSS_FLOOR`` (0.85); cross is built for coding + all aa only. The
    floor returned here is the read's completeness gate: below it the caller must go
    live, because the cache provably does not hold everything the query asks for.

    ``target_species`` may name several. The floor is then the **stricter** of the
    two, unless the selection is exactly the query's own species:

    * ``{own}``            -> 0.75, the same-species cache is complete
    * ``{other}``          -> 0.85
    * ``{own, other, ...}``-> 0.85

    The mixed case is the one worth stating. Serving it at 0.75 would return every
    same-species neighbour down to 0.75 but only those cross-species neighbours at
    0.85 and above, because that is all the builder stored — a result set that looks
    complete and silently favours the query's own species. Taking the stricter floor
    means a mixed selection below 0.85 falls to the live path, which computes both
    sides at the same threshold.
    """
    cross_cached = kind == "aa" or cat_q == "coding"
    if target_species:
        tgt = _species_tags(target_species)
        if tgt == [sp_q]:
            return tgt, NEIGHBOR_FLOOR, True
        return tgt, CROSS_FLOOR, cross_cached
    if same_species:
        return ([sp_q] if sp_q else None), NEIGHBOR_FLOOR, True
    return None, CROSS_FLOOR, cross_cached


def find_similar(
    session: Session, uid: str, min_identity: float = 0.0,
    cap: int = DEFAULT_CAP, same_species: bool = True,
    target_species: str | None = None,
    functional: bool | None = None, source: str | None = None,
) -> SearchResults | None:
    """Coding alleles ranked by **nucleotide** identity to ``uid`` (excl. self).

    ``target_species`` restricts the neighbours to a chosen species or several,
    comma-separated (``"human"``, ``"human,mouse"``) — the cross-species "closest
    ortholog" lookup — overriding ``same_species``. See :func:`_neighbor_scope` for
    which identity floor a mixed selection uses, and why it is the stricter one.
    ``functional`` / ``source`` are the advanced filters (see :func:`search_by_sequence`),
    applied to the neighbour on both the cached and the live path.
    """
    row = session.execute(
        select(Sequence.species, Sequence.category, Sequence.segment)
        .where(Sequence.uid == uid)
    ).first()
    sp_q, cat_q, seg_q = (row[0], row[1], row[2]) if row else (None, None, None)
    tag, floor, cache_ok = _neighbor_scope(
        sp_q, cat_q, "nt", same_species, target_species
    )
    if cache_ok:
        cached = _neighbors_from_table(
            session, uid, "nt", min_identity, cap, species=tag, floor=floor,
            functional=functional, source=source, segment=seg_q,
        )
        if cached is not None:
            return cached
    return _find_similar(
        session, uid, Sequence.sequence, lambda s: s.sequence, "nt",
        min_identity, cap, same_species, target_species,
        functional=functional, source=source,
    )


def find_similar_aa(
    session: Session, uid: str, min_identity: float = 0.0,
    cap: int = DEFAULT_CAP, same_species: bool = True,
    target_species: str | None = None,
    functional: bool | None = None, source: str | None = None,
) -> SearchResults | None:
    """Coding alleles ranked by **amino-acid** identity to ``uid`` (excl. self).

    ``target_species`` restricts to a chosen species or several, comma-separated
    (cross-species AA is fully cached).
    ``functional`` / ``source`` are the advanced filters applied to the neighbour.
    """
    row = session.execute(
        select(Sequence.species, Sequence.category, Sequence.segment)
        .where(Sequence.uid == uid)
    ).first()
    sp_q, cat_q, seg_q = (row[0], row[1], row[2]) if row else (None, None, None)
    tag, floor, cache_ok = _neighbor_scope(
        sp_q, cat_q, "aa", same_species, target_species
    )
    if cache_ok:
        cached = _neighbors_from_table(
            session, uid, "aa", min_identity, cap, species=tag, floor=floor,
            functional=functional, source=source, segment=seg_q,
        )
        if cached is not None:
            return cached
    return _find_similar(
        session, uid, Sequence.sequence_aa, lambda s: s.sequence_aa, "aa",
        min_identity, cap, same_species, target_species,
        functional=functional, source=source,
    )


# --- diff with consequence (INT-1) + IMGT alignment (SIM-3) ----------------
#
# One receptor_utils-backed helper, shared by the diff (region labels) and the
# alignment view. V-only: region labels and IMGT numbering are defined for
# V-REGION. D/J/C fall back to positional (equal-length) diffs, matching the
# Hamming-only policy in CLAUDE.md.

# IMGT V-REGION region boundaries, by IMGT unique-numbering codon (1-based).
_V_REGION_BOUNDS = ((26, "FR1"), (38, "CDR1"), (55, "FR2"), (65, "CDR2"), (104, "FR3"))


def _region_of(imgt_codon: int | None) -> str | None:
    if imgt_codon is None:
        return None
    for hi, name in _V_REGION_BOUNDS:
        if imgt_codon <= hi:
            return name
    return "CDR3"


@lru_cache(maxsize=1)
def _imgt_gapped_refs() -> dict[tuple[str, str], tuple[tuple[str, str], ...]]:
    """Every IMGT-gapped V-REGION in IMGT's gapped bulk FASTA
    (``settings.imgt_genedb_fasta_gapped``), keyed ``(species tag, locus)``. Empty when
    the file is absent -- the caller then reports the alignment as unavailable rather
    than guessing a numbering. Never silently fabricated.

    Parsed inline (a few lines) rather than via ``ingest.fasta`` on purpose: the query
    server is independent of the ``ingest`` package, which the public build strips
    (see ``tests/test_public_build.py``)."""
    path = Path(get_settings().imgt_genedb_fasta_gapped)
    if not path.exists():
        return {}
    refs: dict[tuple[str, str], list[tuple[str, str]]] = {}
    key = name = None
    chunks: list[str] = []

    def flush():
        if key is not None:
            refs.setdefault(key, []).append((name, "".join(chunks).upper()))

    for line in path.read_text().splitlines():
        if line.startswith(">"):
            flush()
            key, chunks = None, []
            # accession|allele|species[_strain]|functionality|label|...
            p = line[1:].split("|")
            if len(p) < 5 or p[4].strip() != "V-REGION":
                continue
            try:
                tag = normalize_species(" ".join(p[2].split("_")[0].split()[:2]))
            except SpeciesError:
                continue
            key, name = (tag, p[1][:3]), p[1]
        elif key is not None:
            chunks.append(line.strip())
    flush()
    return {k: tuple(v) for k, v in refs.items()}


def _imgt_gapped_ref(species: str, locus: str) -> tuple[tuple[str, str], ...] | None:
    """The IMGT-gapped V reference for ``(species, locus)``: ``(name, gapped_seq)``
    pairs, or None when IMGT publishes none (or the file is not deployed)."""
    return _imgt_gapped_refs().get((species, locus))


def _imgt_number(species: str, locus: str | None, nt_seq: str) -> tuple[str, str] | None:
    """IMGT-number one V-REGION nucleotide sequence: return ``(gapped_nt,
    gapped_aa)`` or None. Delegates the numbering to
    ``receptor_utils.number_v.gap_sequence`` (gap against the closest gapped
    reference) — we do not reinvent IMGT numbering."""
    if not locus:
        return None
    ref = _imgt_gapped_ref(species, locus)
    if not ref:
        return None
    from receptor_utils import number_v

    gapped_ref = {k: v for k, v in ref}
    ungapped_ref = {k: v.replace(".", "") for k, v in ref}
    try:
        gnt, gaa, _notes = number_v.gap_sequence(nt_seq, gapped_ref, ungapped_ref)
    except Exception:  # noqa: BLE001 - a sequence number can't be assigned -> no labels
        return None
    return gnt, gaa


def _imgt_codon_map(gapped_nt: str) -> dict[int, int]:
    """Map ungapped-codon index (0-based) -> IMGT codon number (1-based) from a
    gapped nt sequence. IMGT V gaps are whole codons on codon boundaries, so a
    codon containing '.' is a gap position that advances the IMGT number but not
    the ungapped index.
    # ponytail: assumes whole-codon gaps (true for IMGT V numbering); a partial-
    # codon gap would mismap, but number_v never emits one for V-REGION."""
    out: dict[int, int] = {}
    ungapped = 0
    for ci in range(0, len(gapped_nt) - len(gapped_nt) % 3, 3):
        if "." in gapped_nt[ci:ci + 3]:
            continue
        out[ungapped] = ci // 3 + 1
        ungapped += 1
    return out


def _aa_of(codon: str) -> str | None:
    return _CODON_TABLE.get(codon, "X") if len(codon) == 3 else None


def _pair_subs(
    q_seq: str, o_seq: str, imgt_map: dict[int, int] | None
) -> tuple[list[NtSub], int, int]:
    """Per-position substitutions between two EQUAL-LENGTH sequences, with the
    per-codon protein consequence. Returns (subs, n_nonsyn_codons, n_syn_subs)."""
    subs: list[NtSub] = []
    changed_codons: set[int] = set()
    n_syn = 0
    for p in range(len(q_seq)):
        if q_seq[p] == o_seq[p]:
            continue
        ci = p // 3
        cf, ct = q_seq[ci * 3:ci * 3 + 3], o_seq[ci * 3:ci * 3 + 3]
        af, at = _aa_of(cf), _aa_of(ct)
        syn = (af is not None and af == at)
        if syn:
            n_syn += 1
        elif af is not None:
            changed_codons.add(ci)
        imgt = imgt_map.get(ci) if imgt_map else None
        subs.append(NtSub(
            pos=p + 1, from_nt=q_seq[p], to_nt=o_seq[p],
            imgt_codon=imgt, region=_region_of(imgt),
            codon_from=cf, codon_to=ct, aa_from=af, aa_to=at, synonymous=syn,
        ))
    return subs, len(changed_codons), n_syn


def _query_imgt_map(seq: Sequence) -> tuple[dict[int, int] | None, bool]:
    """``(ungapped-codon -> IMGT-codon map, region_labeled)`` for a V query, else
    ``(None, False)`` — non-V, or the IMGT-gapped reference is not cached."""
    if seq.segment != "V":
        return None, False
    numbered = _imgt_number(seq.species, seq.locus, seq.sequence)
    if numbered is None:
        return None, False
    return _imgt_codon_map(numbered[0]), True


def _neighbour_diff(
    seq: Sequence, other: Sequence, identity: float, label: str | None,
    imgt_map: dict[int, int] | None,
) -> NeighbourDiff:
    """The diff from ``seq`` to one ``other`` sequence. Equal length -> per-position
    substitutions with consequence; a length difference (indel) reports the delta
    only, since an edit distance there would open gaps the equal-length policy
    refuses."""
    delta = len(other.sequence) - len(seq.sequence)
    if delta == 0:
        subs, n_nonsyn, n_syn = _pair_subs(seq.sequence, other.sequence, imgt_map)
        note = None
    else:
        subs, n_nonsyn, n_syn = [], 0, 0
        note = (f"differs in length by {abs(delta)} nt (indel); "
                "per-position substitutions not computed")
    return NeighbourDiff(
        uid=other.uid, label=label, identity=identity,
        same_length=(delta == 0), length_delta=delta,
        n_subs=len(subs), n_nonsyn=n_nonsyn, n_syn=n_syn, subs=subs, note=note,
    )


def diff_to_neighbours(
    session: Session, uid: str, cap: int = 12, min_identity: float = 0.90,
) -> DiffSummary | None:
    """The exact nt substitutions from a record to its nearest allele(s), each
    with its protein-level synonymous / non-synonymous consequence (INT-1).

    The nearest alleles come from :func:`find_similar` (the same neighbour cache
    the rest of the site uses); we then join on the **sequence** (never the name)
    to compute the diff. V substitutions carry an IMGT FR/CDR region label when
    the IMGT-gapped reference is cached; D/J/C report positional diffs only."""
    seq = session.scalar(select(Sequence).where(Sequence.uid == uid))
    if seq is None:
        return None
    label = _label(_aliases(session, seq.id))
    sim = find_similar(session, uid, min_identity=min_identity, cap=cap)
    hits = sim.hits if sim else []
    by_uid = {
        s.uid: s for s in session.scalars(
            select(Sequence).where(Sequence.uid.in_([h.uid for h in hits]))
        )
    }
    imgt_map, region_labeled = _query_imgt_map(seq)
    out = [
        _neighbour_diff(seq, by_uid[h.uid], h.score, h.label, imgt_map)
        for h in hits if h.uid in by_uid
    ]
    # Stable tiebreak: the cache orders only by identity, so equal-identity ties come
    # back in arbitrary DB order. Sort by (identity desc, uid) for determinism.
    out.sort(key=lambda nb: (-nb.identity, nb.uid))
    return DiffSummary(
        uid=seq.uid, label=label, segment=seq.segment,
        region_labeled=region_labeled,
        region_source="IMGT/GENE-DB gapped V reference" if region_labeled else None,
        neighbours=out,
    )


def _build_alignment(
    seq: Sequence, ref_seq: Sequence, seq_label: str | None, ref_label: str | None,
) -> AlignmentView | None:
    """IMGT-gapped alignment of ``seq`` against ``ref_seq``, carrying both the
    nucleotide and amino-acid gapped rows (nt is the primary view). Both are
    numbered independently to IMGT coordinates, so equal positions compare
    directly. Returns None when either can't be numbered."""
    q_num = _imgt_number(seq.species, seq.locus, seq.sequence)
    r_num = _imgt_number(ref_seq.species, ref_seq.locus, ref_seq.sequence)
    if q_num is None or r_num is None:
        return None
    q_nt, q_aa = q_num
    r_nt, r_aa = r_num
    nt_len, aa_len = max(len(q_nt), len(r_nt)), max(len(q_aa), len(r_aa))
    q_nt, r_nt = q_nt.ljust(nt_len, "."), r_nt.ljust(nt_len, ".")
    q_aa, r_aa = q_aa.ljust(aa_len, "."), r_aa.ljust(aa_len, ".")
    return AlignmentView(
        uid=seq.uid, label=seq_label,
        against_uid=ref_seq.uid, against_label=ref_label,
        nt_length=nt_len, aa_length=aa_len,
        nt_regions=[_region_of(i // 3 + 1) or "" for i in range(nt_len)],
        aa_regions=[_region_of(i + 1) or "" for i in range(aa_len)],
        rows=[
            AlignmentRow(uid=ref_seq.uid, label=ref_label, is_reference=True,
                         gapped_nt=r_nt, gapped_aa=r_aa),
            AlignmentRow(uid=seq.uid, label=seq_label, is_reference=False,
                         gapped_nt=q_nt, gapped_aa=q_aa),
        ],
    )


def imgt_alignment(
    session: Session, uid: str, against: str | None = None,
) -> AlignmentView | None:
    """IMGT-gapped alignment (nt primary, aa secondary) of a V record against a
    reference allele (SIM-3). The reference is ``against`` if given, else the
    gene's ``*01``, else the nearest neighbour. None for non-V or when the
    IMGT-gapped reference is not cached."""
    seq = session.scalar(select(Sequence).where(Sequence.uid == uid))
    if seq is None or seq.segment != "V":
        return None
    label = _label(_aliases(session, seq.id))
    ref_seq = None
    if against:
        ref_seq = session.scalar(select(Sequence).where(Sequence.uid == against))
    if ref_seq is None and label and "*" in label:
        star01 = f"{label.split('*')[0]}*01"
        for cand in resolve_exact_name(session, star01, seq.species):
            if cand.uid != seq.uid and cand.segment == "V":
                ref_seq = cand
                break
    if ref_seq is None:
        sim = find_similar(session, uid, min_identity=0.75, cap=1)
        if sim and sim.hits:
            ref_seq = session.scalar(
                select(Sequence).where(Sequence.uid == sim.hits[0].uid)
            )
    if ref_seq is None:
        return None
    return _build_alignment(
        seq, ref_seq, label, _label(_aliases(session, ref_seq.id))
    )


def _alignment_text(query: Sequence, hit: Sequence) -> str | None:
    """Formatted IMGT alignment of the query vs one hit via
    ``receptor_utils.sequence_alignment.create_alignment`` — the shared tool VDJbase
    uses, not a hand-rolled renderer. The query anchors it (reference row, shown with
    its AA); the hit shows only its differences (identical codons ``---``, changed
    bases in lower case, non-silent residue changes above). None if either sequence
    can't be IMGT-numbered."""
    q = _imgt_number(query.species, query.locus, query.sequence)
    h = _imgt_number(hit.species, hit.locus, hit.sequence)
    if q is None or h is None:
        return None
    from receptor_utils.sequence_alignment import create_alignment

    # create_alignment takes the alphabetically-first key as the reference; force the
    # query to sort first so it anchors the alignment. Rows are labelled by UID.
    qkey, hkey = query.uid, hit.uid
    if not qkey < hkey:
        qkey = " " + qkey
    return create_alignment({qkey: q[0], hkey: h[0]}, sequence_type="V", codon_wrap=30) or None


def compare_pair(
    session: Session, uid: str, other_uid: str,
) -> tuple[NeighbourDiff, str | None] | None:
    """The on-demand comparison surfaced per similarity hit: the nt diff (with
    consequence + region) and the formatted IMGT alignment of ``uid`` against
    ``other_uid``. Returns None if either UID is unknown; the alignment is None for
    non-V."""
    seq = session.scalar(select(Sequence).where(Sequence.uid == uid))
    other = session.scalar(select(Sequence).where(Sequence.uid == other_uid))
    if seq is None or other is None:
        return None
    other_label = _label(_aliases(session, other.id))
    imgt_map, _ = _query_imgt_map(seq)
    identity = round(Levenshtein.normalized_similarity(seq.sequence, other.sequence), 4)
    diff = _neighbour_diff(seq, other, identity, other_label, imgt_map)
    align = _alignment_text(seq, other) if seq.segment == "V" else None
    return diff, align


# --- protein group (INT-2) -------------------------------------------------

def protein_group(session: Session, seq: Sequence) -> list[SearchHit]:
    """The other coding alleles whose V-REGION encodes the identical protein
    (same ``protein_uid``, different DNA). Keyed on the protein hash, never a
    name. Empty for a sequence with no protein UID (non-V, or no AA)."""
    if not seq.protein_uid:
        return []
    others = list(session.scalars(
        select(Sequence).where(
            Sequence.protein_uid == seq.protein_uid,
            Sequence.uid != seq.uid,
            Sequence.archived.is_(False),
        )
    ))
    hits = [
        SearchHit(
            uid=o.uid, label=_label(_aliases(session, o.id)), species=o.species,
            category=o.category, locus=o.locus, sequence_type=o.segment,
            length=o.length, score=1.0, matched_by="protein_uid",
            aliases=list(dict.fromkeys(a.name for a in _aliases(session, o.id))),
        )
        for o in others
    ]
    return sorted(hits, key=lambda h: (h.species, h.uid))


# --- reciprocal cross-species nearest (SIM-2) ------------------------------

def reciprocal_orthologs(
    session: Session, uid: str, min_identity: float = CROSS_FLOOR,
) -> list[OrthologHint]:
    """For each other species, this record's nearest allele there, flagged when
    the relationship is reciprocal (the record is also that allele's nearest in
    the record's own species). A **sequence-based** hint, explicitly not a
    curated orthology call. Joins on sequence identity via the neighbour cache,
    never on names."""
    seq = session.scalar(select(Sequence).where(Sequence.uid == uid))
    if seq is None:
        return []
    cross = find_similar(
        session, uid, min_identity=min_identity, cap=400, same_species=False,
    )
    if cross is None:
        return []
    best: dict[str, SearchHit] = {}  # hits are identity-desc, so first per species is best
    for h in cross.hits:
        if h.species == seq.species:
            continue
        best.setdefault(h.species, h)
    out: list[OrthologHint] = []
    for sp, h in best.items():
        back = find_similar(
            session, h.uid, min_identity=min_identity, cap=1,
            same_species=False, target_species=seq.species,
        )
        top = back.hits[0] if back and back.hits else None
        out.append(OrthologHint(
            species=sp, uid=h.uid, label=h.label, identity=h.score,
            reciprocal=(top is not None and top.uid == seq.uid),
            back_uid=top.uid if top else None,
            back_label=top.label if top else None,
        ))
    out.sort(key=lambda o: (-o.identity, o.species))
    return out
