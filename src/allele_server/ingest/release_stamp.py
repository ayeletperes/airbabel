"""Stamp a build of the served dataset: what was loaded, and how it is scored.

A load timestamp says when rows were written. It does not say *which* data, and it
says nothing about the **distance computation** — yet every identity score, every
neighbour list and every threshold the site shows is a function of that. Two instances
can hold the same alleles and rank them differently, and nothing on either page would
reveal it. So a release pins both halves and digests each; see
``alembic/versions/0009_dataset_releases.py`` for why the table sits outside the
reload's TRUNCATE.

Digests are computed **from the database**, never from the ingest's own bookkeeping.
Counting what you believe you wrote is how a load that silently dropped rows still
reports a full dataset; hashing what is actually stored cannot.
"""

from __future__ import annotations

import hashlib
import json
from datetime import UTC, datetime

from sqlalchemy import func, select
from sqlalchemy.orm import Session

from .. import __version__, queries
from ..db.models import DatasetRelease, Name, Sequence, SequenceNeighbor, Source

#: Rows per fetch when streaming the digest inputs. The digest is order-sensitive, so
#: the ORDER BY is not optional — it is what makes two databases comparable.
_CHUNK = 10_000


def content_digest(session: Session) -> tuple[str, dict]:
    """sha256 over the sequence identity set and the name set, plus the counts.

    Two streams, each sorted and fed into one hash with a labelled separator:

    * ``species | category | sequence_hash`` — the identity triple this database is
      keyed on (CLAUDE.md: identity is the sequence, names are labels).
    * ``name | scheme | source`` — the label layer, joined to the source *name* rather
      than ``source_id``, because a reload restarts that sequence and an id-based
      digest would change on a rebuild that loaded identical data.

    Both, because an enrichment pass that adds official names changes what the service
    answers without touching a single sequence. A digest over sequences alone would
    call that the same release.
    """
    h = hashlib.sha256()

    h.update(b"sequences\n")
    n_seq = 0
    for sp, cat, sh in session.execute(
        select(Sequence.species, Sequence.category, Sequence.sequence_hash)
        .order_by(Sequence.species, Sequence.category, Sequence.sequence_hash)
        .execution_options(yield_per=_CHUNK)
    ):
        h.update(f"{sp}|{cat}|{sh}\n".encode())
        n_seq += 1

    h.update(b"names\n")
    n_name = 0
    for nm, scheme, src in session.execute(
        select(Name.name, Name.naming_scheme, Source.name)
        .join(Source, Source.id == Name.source_id)
        .order_by(Name.name, Name.naming_scheme, Source.name)
        .execution_options(yield_per=_CHUNK)
    ):
        h.update(f"{nm}|{scheme}|{src}\n".encode())
        n_name += 1

    return h.hexdigest(), {"sequences": n_seq, "names": n_name}


def _inventory(session: Session) -> dict:
    """Counts and source versions, for the record's ``details`` blob."""
    by_category = dict(session.execute(
        select(Sequence.category, func.count()).group_by(Sequence.category)
    ).all())
    by_segment = dict(session.execute(
        select(Sequence.segment, func.count())
        .where(Sequence.category == "coding", Sequence.archived.is_(False),
               Sequence.segment.is_not(None))
        .group_by(Sequence.segment)
    ).all())
    archived = session.scalar(
        select(func.count()).select_from(Sequence)
        .where(Sequence.category == "coding", Sequence.archived.is_(True))
    )
    species = session.scalar(
        select(func.count(func.distinct(Sequence.species)))
        .where(Sequence.category == "coding", Sequence.archived.is_(False))
    )
    # The same tallies behind the served-site gate (queries._visible_cond), which is
    # what /overview, dataset_stats and the paper report. `by_category` / `coding_by_
    # segment` above stay RAW (they match the digested set); these hide shelved-only
    # (HUSA-only) sequences. Use served_* when quoting what the site shows.
    vis = queries._visible_cond()
    gate = (vis,) if vis is not None else ()
    served_by_category = dict(session.execute(
        select(Sequence.category, func.count()).where(*gate).group_by(Sequence.category)
    ).all())
    served_by_segment = dict(session.execute(
        select(Sequence.segment, func.count())
        .where(Sequence.category == "coding", Sequence.archived.is_(False),
               Sequence.segment.is_not(None), *gate)
        .group_by(Sequence.segment)
    ).all())
    served_archived = session.scalar(
        select(func.count()).select_from(Sequence)
        .where(Sequence.category == "coding", Sequence.archived.is_(True), *gate)
    )
    # Every source with the versions actually ingested -- the reference releases this
    # build stands on. A release that named no source versions would not let anyone
    # rebuild it.
    sources: dict[str, list[str]] = {}
    for name, version in session.execute(
        select(Source.name, Source.version).distinct().order_by(Source.name)
    ).all():
        sources.setdefault(name, [])
        if version:
            sources[name].append(version)
    return {
        "by_category": by_category,
        "coding_by_segment": by_segment,
        "archived": archived,
        "served_by_category": served_by_category,
        "served_coding_by_segment": served_by_segment,
        "served_archived": served_archived,
        "species": species,
        "sources": {k: sorted(v) for k, v in sources.items()},
    }


def _next_release(session: Session, now: datetime) -> str:
    """``YYYY.MM.DD``, with ``.N`` for the Nth build of the same day.

    Date-based rather than a counter: the release string is read by people comparing
    two deployments, and a date answers "is mine older?" without a lookup. The suffix
    exists because rebuilding twice in a day is routine during a data fix.
    """
    base = now.strftime("%Y.%m.%d")
    taken = set(session.scalars(
        select(DatasetRelease.release).where(DatasetRelease.release.like(f"{base}%"))
    ))
    if base not in taken:
        return base
    n = 2
    while f"{base}.{n}" in taken:
        n += 1
    return f"{base}.{n}"


def stamp_release(session: Session) -> dict:
    """Write one ``dataset_releases`` row describing the database as it stands now.

    Safe to run against a database this process did not build -- that is the backfill
    path for instances that predate the table (``allele-ingest --stamp-release``).
    Everything is read back from the database, so a stamp cannot claim more than is
    actually there.
    """
    now = datetime.now(UTC)
    digest, counted = content_digest(session)
    edges = session.scalar(select(func.count()).select_from(SequenceNeighbor))
    details = _inventory(session)
    details["digested"] = counted
    # The full profile, not just its digest: a reader comparing two releases needs to
    # see *which* parameter moved, not only that one did.
    details["similarity"] = queries.similarity_profile()

    row = DatasetRelease(
        release=_next_release(session, now),
        built_at=now,
        content_digest=digest,
        similarity_digest=queries.similarity_digest(),
        code_version=__version__,
        neighbors_built=bool(edges),
        neighbor_edges=edges,
        details=json.dumps(details, sort_keys=True, separators=(",", ":")),
    )
    session.add(row)
    session.flush()
    return {
        "release": row.release,
        "content_digest": digest,
        "similarity_digest": row.similarity_digest,
        "neighbor_edges": edges,
        **details,
    }
