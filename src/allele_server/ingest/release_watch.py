"""Notice when IMGT publishes a release, and say how fast releases arrive.

**Report only.** Nothing here writes to the database, ingests, or registers a species.
It answers "what has IMGT published that we have not taken, and how often does that
happen", and leaves acting on the answer to a person. Moving to auto-ingest later is a
deliberate step, not a flag flip: see ``docs/source-policy.md`` for why admission stays
human.

**Detect from the archive, ingest from imgt.org.** Detection wants something cheap and
history-bearing, so it reads the directory listing of Jamie Heather's ``genedb-releases``
archive (https://github.com/JamieHeather/genedb-releases) — one HTTP call instead of a
7.6 MB FASTA download and a full diff. Ingest keeps coming direct from IMGT
(``imgt_genedb._BULK_URL``), because serving a mirror of a mirror is not what this
database is for. ``imgt.toml`` records the channel and the credit.

**A dead tracker looks exactly like a quiet week, and that is the failure that matters.**
If the archive stops being updated — maintainer moves on, repo renamed, rate limit, a
silently changed path — a checker that only reports "no new releases" reports it forever
and nobody learns anything. So staleness is a first-class finding here, not the absence
of one: ``ReleaseReport.stale`` is set from the observed cadence rather than a guessed
constant, because the right threshold for a weekly archive is not the right one for an
annual one.

Detection lag is bounded by the archive's cadence, not IMGT's: we learn when Jamie's
snapshot lands, not when IMGT publishes. That is a property to state, not a defect.
"""

from __future__ import annotations

import re
import statistics
from collections.abc import Callable
from dataclasses import dataclass, field
from datetime import date, datetime, timezone

from sqlalchemy import select
from sqlalchemy.orm import Session

_API = "https://api.github.com/repos/JamieHeather/genedb-releases/contents/releases"

# Directory names look like ``2026-01-04_GENEDB_20261-7``: ISO date, the channel the
# snapshot came through, then IMGT's own ``YYYYWW-V`` release id.
_DIR_RE = re.compile(r"^(\d{4}-\d{2}-\d{2})_([A-Za-z-]+)_(\d{4,6}-\d+)$")

# How many observed intervals to average when deciding what "overdue" means here.
_CADENCE_WINDOW = 12
# A release cadence is never perfectly regular, so allow a generous multiple of the
# recent median before calling the tracker stale. Floor it so a very fast cadence does
# not produce a hair-trigger, and cap it so an archive that used to be annual cannot
# hide a year of silence.
_STALE_FACTOR = 4
_STALE_MIN_DAYS = 45
_STALE_MAX_DAYS = 400


@dataclass(frozen=True)
class ArchiveRelease:
    """One snapshot directory in the archive."""

    dirname: str
    released: date
    channel: str        # WBM (Wayback) / JH-GeneDL / GENEDB
    imgt_release: str   # IMGT's own YYYYWW-V id

    @property
    def label(self) -> str:
        return self.released.isoformat()


@dataclass
class ReleaseReport:
    """What the checker found. Every field is an observation, not an action."""

    checked_at: datetime
    available: list[ArchiveRelease] = field(default_factory=list)
    ingested_labels: list[str] = field(default_factory=list)
    #: Archive releases newer than the newest one we hold.
    new_releases: list[ArchiveRelease] = field(default_factory=list)
    #: Releases older than our newest that we never took. Not "new", but not held —
    #: worth seeing, because it is the difference between a gap and a backlog.
    missing_older: list[ArchiveRelease] = field(default_factory=list)
    #: Days between consecutive releases over the recent window.
    median_interval_days: float | None = None
    mean_interval_days: float | None = None
    releases_per_year: float | None = None
    days_since_latest: int | None = None
    stale_after_days: int | None = None
    stale: bool = False
    #: Set when the listing could not be read at all. An empty result and a failed
    #: fetch must never look the same.
    error: str | None = None

    @property
    def ok(self) -> bool:
        return self.error is None and not self.stale

    def lines(self) -> list[str]:
        """A human-readable report. The loud things go first."""
        out: list[str] = []
        if self.error:
            out.append(f"ERROR  could not read the release archive: {self.error}")
            out.append("       Treat this as 'unknown', never as 'no new releases'.")
            return out
        if self.stale:
            out.append(
                f"STALE  no new snapshot for {self.days_since_latest} days "
                f"(overdue past {self.stale_after_days}). The archive may have stopped "
                f"updating; a dead tracker reports 'no changes' forever."
            )
        out.append(f"archive     {len(self.available)} releases, "
                   f"{self.available[0].label} .. {self.available[-1].label}"
                   if self.available else "archive     none found")
        tail = ", ".join(self.ingested_labels[-3:])
        more = "..." if len(self.ingested_labels) > 3 else ""
        out.append(f"ingested    {len(self.ingested_labels)} ({tail}{more})")
        if self.median_interval_days is not None:
            out.append(
                f"rate        every {self.median_interval_days:.0f} days median "
                f"({self.mean_interval_days:.0f} mean) over the last "
                f"{min(_CADENCE_WINDOW, max(len(self.available) - 1, 0))} intervals "
                f"= {self.releases_per_year:.0f}/year"
            )
        if self.days_since_latest is not None:
            out.append(f"latest      {self.available[-1].label} "
                       f"({self.days_since_latest} days ago), IMGT "
                       f"{self.available[-1].imgt_release}")
        if self.new_releases:
            out.append(f"NEW         {len(self.new_releases)} release(s) we do not hold: "
                       + ", ".join(r.label for r in self.new_releases[:6])
                       + ("..." if len(self.new_releases) > 6 else ""))
        else:
            out.append("NEW         none — nothing newer than what we hold")
        if self.missing_older:
            out.append(f"not taken   {len(self.missing_older)} older release(s) exist in "
                       f"the archive that we never ingested (by design: we sample, "
                       f"we do not mirror)")
        return out


def parse_release_dir(name: str) -> ArchiveRelease | None:
    """Parse one archive directory name, or None if it is not a release directory."""
    m = _DIR_RE.match(name.strip())
    if m is None:
        return None
    try:
        released = date.fromisoformat(m.group(1))
    except ValueError:
        return None
    return ArchiveRelease(
        dirname=name.strip(), released=released,
        channel=m.group(2), imgt_release=m.group(3),
    )


def _default_fetcher(url: str) -> str:
    import urllib.request

    req = urllib.request.Request(
        url, headers={"Accept": "application/vnd.github+json",
                      "User-Agent": "allele-server-release-watch"},
    )
    with urllib.request.urlopen(req, timeout=60) as resp:  # noqa: S310 - fixed https URL
        return resp.read().decode("utf-8", "replace")


def list_archive_releases(
    fetcher: Callable[[str], str] | None = None, *, url: str = _API,
) -> list[ArchiveRelease]:
    """Every release snapshot the archive offers, oldest first.

    One directory listing, no downloads. Entries that do not parse are skipped rather
    than raised on: the archive may hold README or tooling directories, and a checker
    that dies on an unexpected sibling is a checker that stops noticing.
    """
    import json

    raw = (fetcher or _default_fetcher)(url)
    entries = json.loads(raw)
    if isinstance(entries, dict):  # GitHub error payloads are objects, not arrays
        raise RuntimeError(str(entries.get("message") or entries)[:200])
    out = [
        rel for e in entries
        if e.get("type") == "dir" and (rel := parse_release_dir(e.get("name", "")))
    ]
    return sorted(out, key=lambda r: (r.released, r.dirname))


def ingested_release_labels(session: Session) -> list[str]:
    """Release labels we already hold, from the historical IMGT source rows.

    Read off ``sources`` rather than a table of our own: the release history is already
    recorded there by ``imgt_history.ingest_imgt_release``, one source row per release,
    and inventing a second place to write it down would create two answers to one
    question. ``reload`` truncates ``sources``, but it also re-ingests the snapshots, so
    the two move together.
    """
    from ..db.models import Source

    rows = session.execute(
        select(Source.version).where(Source.note.like("historical IMGT%")).distinct()
    ).scalars()
    return sorted(v for v in rows if v)


def check_releases(
    session: Session,
    *,
    fetcher: Callable[[str], str] | None = None,
    now: datetime | None = None,
) -> ReleaseReport:
    """Compare the archive against what we hold, and measure the release rate.

    Never raises on a fetch failure: it records the failure in ``error`` so a caller
    cannot mistake "could not look" for "nothing new".
    """
    checked_at = now or datetime.now(timezone.utc)
    report = ReleaseReport(checked_at=checked_at)
    try:
        report.available = list_archive_releases(fetcher)
    except Exception as exc:  # noqa: BLE001 - any failure must surface as 'unknown'
        report.error = f"{type(exc).__name__}: {exc}"
        return report
    report.ingested_labels = ingested_release_labels(session)
    if not report.available:
        report.error = "the archive listing parsed to zero releases"
        return report

    held = set(report.ingested_labels)
    newest_held = max(held) if held else ""
    report.new_releases = [r for r in report.available if r.label > newest_held]
    report.missing_older = [
        r for r in report.available if r.label <= newest_held and r.label not in held
    ]

    days = [
        (b.released - a.released).days
        for a, b in zip(report.available, report.available[1:])
    ]
    recent = days[-_CADENCE_WINDOW:]
    if recent:
        report.median_interval_days = float(statistics.median(recent))
        report.mean_interval_days = float(statistics.fmean(recent))
        if report.median_interval_days > 0:
            report.releases_per_year = 365.25 / report.median_interval_days
    latest = report.available[-1].released
    report.days_since_latest = (checked_at.date() - latest).days
    # Threshold from the observed cadence, clamped: a weekly archive and an annual one
    # should not share a constant, but neither should be able to hide indefinitely.
    base = report.median_interval_days or _STALE_MIN_DAYS
    report.stale_after_days = int(
        min(max(base * _STALE_FACTOR, _STALE_MIN_DAYS), _STALE_MAX_DAYS)
    )
    report.stale = report.days_since_latest > report.stale_after_days
    return report
