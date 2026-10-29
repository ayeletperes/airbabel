"""The candidate queue: recording what was found, and who decided what about it.

Transport-agnostic, like ``queries.py`` — the CLI in this package drives it, and the
discovery channels will drive the same functions. Nothing here may set a decision
status without a named curator; the database enforces that too
(``ck_source_candidates_decided``), so this is belt and braces on purpose.

See ``docs/source-policy.md``.
"""

from __future__ import annotations

import json
import os
from datetime import UTC, datetime

from sqlalchemy import select
from sqlalchemy.orm import Session

from ..db.models import SourceCandidate
from . import SourceError, registry

__all__ = [
    "add",
    "curator",
    "decide",
    "dedupe_key_for",
    "get",
    "listing",
    "normalise_doi",
    "parse_ref",
    "suppressed",
]

_CURATOR_ENV = "AIRBABEL_CURATOR"


def parse_ref(ref: str) -> tuple[str, str]:
    """Split a ``<channel>:<external_id>`` CLI handle.

    Raises:
        SourceError: If the handle is malformed or names an unknown channel.
    """
    channel, sep, external_id = ref.partition(":")
    if not sep or not external_id.strip():
        raise SourceError(
            f"{ref!r} is not a candidate reference; expected <channel>:<external_id>, "
            f"e.g. manual:newdb"
        )
    if channel not in SourceCandidate.CHANNELS:
        raise SourceError(
            f"unknown channel {channel!r}; known: {sorted(SourceCandidate.CHANNELS)}"
        )
    return channel, external_id.strip()


def normalise_doi(doi: str | None) -> str | None:
    """Strip the many ways a DOI is written down, so two channels agree on one key."""
    if not doi:
        return None
    value = doi.strip().casefold()
    for prefix in ("https://doi.org/", "http://doi.org/", "doi.org/", "doi:"):
        if value.startswith(prefix):
            value = value[len(prefix):]
    return value.strip("/") or None


def dedupe_key_for(channel: str, external_id: str, doi: str | None) -> str:
    """Cross-channel identity for a candidate.

    The DOI when there is one, so the Zenodo deposit and the paper announcing it collapse
    to one decision; otherwise the channel-local id, which at least dedupes re-crawls.
    """
    return normalise_doi(doi) or f"{channel}:{external_id}"


def curator(explicit: str | None = None) -> str:
    """The name to record against a decision.

    Raises:
        SourceError: if neither ``--by`` nor ``$AIRBABEL_CURATOR`` is set. There is no
            fallback to ``git config user.name``: an approval attributed to a guess is
            worse than no approval: fail loudly, never fall back silently.
    """
    name = (explicit or os.environ.get(_CURATOR_ENV) or "").strip()
    if not name:
        raise SourceError(
            "a decision must name its curator: pass --by 'Your Name' or set "
            f"${_CURATOR_ENV}. Approvals are attributed, never anonymous."
        )
    return name


def get(session: Session, ref: str) -> SourceCandidate | None:
    """The candidate for a ``<channel>:<external_id>`` handle, if we hold one."""
    channel, external_id = parse_ref(ref)
    return session.scalar(
        select(SourceCandidate).where(
            SourceCandidate.channel == channel,
            SourceCandidate.external_id == external_id,
        )
    )


def require(session: Session, ref: str) -> SourceCandidate:
    """As ``get``, but raises when there is no such candidate."""
    candidate = get(session, ref)
    if candidate is None:
        raise SourceError(
            f"no candidate {ref!r}. List what is queued with `allele-sources list`, or "
            f"record a new one with `allele-sources add --channel manual --id ...`."
        )
    return candidate


def add(
    session: Session,
    *,
    channel: str,
    external_id: str,
    title: str | None = None,
    url: str | None = None,
    doi: str | None = None,
    licence: str | None = None,
    species_hint: str | None = None,
    loci_hint: str | None = None,
    score: float | None = None,
    evidence: dict | None = None,
    status: str = "discovered",
) -> tuple[SourceCandidate, bool]:
    """Record a candidate, or refresh one we already hold. Returns ``(row, created)``.

    Idempotent in the way every other pass is: re-seeing an unchanged candidate bumps
    ``last_seen_at`` and fills in anything newly learned, and **never** touches the
    status of a row a human has decided. A crawler cannot un-reject a source.
    """
    if channel not in SourceCandidate.CHANNELS:
        raise SourceError(
            f"unknown channel {channel!r}; known: {sorted(SourceCandidate.CHANNELS)}"
        )
    if status not in SourceCandidate.STATUSES:
        raise SourceError(f"unknown status {status!r}")
    if status in SourceCandidate.DECIDED:
        raise SourceError(
            f"add() cannot set {status!r}: a decision needs a curator — use decide()"
        )

    existing = session.scalar(
        select(SourceCandidate).where(
            SourceCandidate.channel == channel,
            SourceCandidate.external_id == external_id,
        )
    )
    now = datetime.now(UTC)
    if existing is not None:
        existing.last_seen_at = now
        # Fill gaps only: a channel that learns the DOI on a later pass should record
        # it, but must not overwrite something a curator corrected by hand.
        for field, value in (
            ("title", title), ("url", url), ("doi", doi), ("licence", licence),
            ("species_hint", species_hint), ("loci_hint", loci_hint),
        ):
            if value and not getattr(existing, field):
                setattr(existing, field, value)
        if existing.doi and not existing.dedupe_key:
            existing.dedupe_key = dedupe_key_for(channel, external_id, existing.doi)
        session.flush()
        return existing, False

    candidate = SourceCandidate(
        channel=channel,
        external_id=external_id,
        dedupe_key=dedupe_key_for(channel, external_id, doi),
        title=title,
        url=url,
        doi=normalise_doi(doi),
        licence=licence,
        species_hint=species_hint,
        loci_hint=loci_hint,
        status=status,
        score=score,
        evidence=json.dumps(evidence) if evidence else None,
        first_seen_at=now,
        last_seen_at=now,
    )
    session.add(candidate)
    session.flush()
    return candidate, True


def decide(
    session: Session,
    candidate: SourceCandidate,
    *,
    status: str,
    by: str,
    note: str | None = None,
    manifest_id: str | None = None,
    until: datetime | None = None,
) -> SourceCandidate:
    """Record a human decision on a candidate.

    ``rejected`` also sets ``suppressed``, which is what keeps the candidate — and
    anything sharing its ``dedupe_key`` — from being re-proposed by any channel.

    Raises:
        SourceError: for an unknown status, a missing curator, a ``deferred`` with no
            date, or an ``approved`` whose manifest is not itself approved.
    """
    if status not in SourceCandidate.DECIDED:
        raise SourceError(
            f"{status!r} is not a decision; expected one of {SourceCandidate.DECIDED}"
        )
    who = curator(by)
    if status == "deferred" and until is None:
        raise SourceError("defer needs --until: an open-ended deferral is a rejection")
    if status == "approved":
        if not manifest_id:
            raise SourceError("approve needs --manifest: approval is of a manifest")
        manifest = registry().get(manifest_id)
        if manifest is None:
            raise SourceError(
                f"no manifest {manifest_id!r}; scaffold one with "
                f"`allele-sources draft {candidate.ref} --manifest {manifest_id}`"
            )
        if not manifest.approved:
            raise SourceError(
                f"manifest {manifest_id!r} has no [approval] block yet — the approve "
                f"command writes it; this candidate cannot be marked approved first"
            )
        candidate.manifest_id = manifest_id

    candidate.status = status
    candidate.decided_by = who
    candidate.decided_at = datetime.now(UTC)
    candidate.decision_note = note
    candidate.deferred_until = until if status == "deferred" else None
    if status == "rejected":
        candidate.suppressed = True
    session.flush()
    return candidate


def claim(session: Session, candidate: SourceCandidate, *, by: str) -> SourceCandidate:
    """Mark a candidate as being reviewed by a named curator."""
    candidate.decided_by = curator(by)
    candidate.status = "in_review"
    session.flush()
    return candidate


def listing(
    session: Session, *, status: str | None = None, channel: str | None = None
) -> list[SourceCandidate]:
    """Candidates, newest first, optionally filtered."""
    stmt = select(SourceCandidate)
    if status:
        stmt = stmt.where(SourceCandidate.status == status)
    if channel:
        stmt = stmt.where(SourceCandidate.channel == channel)
    return list(session.scalars(stmt.order_by(SourceCandidate.last_seen_at.desc())))


def suppressed(session: Session) -> tuple[set[tuple[str, str]], set[str]]:
    """What discovery must not re-propose: ``(channel, external_id)`` pairs and keys.

    Covers anything rejected, already approved or ingested, and anything deferred whose
    date has not yet arrived. Returning the ``dedupe_key`` set as well is what stops a
    rejected deposit from coming back through a different channel.
    """
    now = datetime.now(UTC)
    pairs: set[tuple[str, str]] = set()
    keys: set[str] = set()
    for row in session.scalars(select(SourceCandidate)):
        blocked = (
            row.suppressed
            or row.status in ("approved", "rejected", "ingested")
            or (
                row.status == "deferred"
                and row.deferred_until is not None
                and row.deferred_until > now
            )
        )
        if blocked:
            pairs.add((row.channel, row.external_id))
            if row.dedupe_key:
                keys.add(row.dedupe_key)
    return pairs, keys
