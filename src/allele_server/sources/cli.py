"""``allele-sources`` — the curation CLI: review candidates, approve sources.

A separate entry point from ``allele-ingest`` because it is a different verb. Ingest
loads data; this decides what may be loaded at all. Plain ``argparse`` with subparsers,
matching the rest of the project (no click, no typer).

    allele-sources add     --channel manual --id newdb --title "NewDB" --url https://…
    allele-sources list    [--status dossier] [--channel manual] [--json]
    allele-sources show    manual:newdb                 # print the dossier
    allele-sources draft   manual:newdb --manifest newdb
    allele-sources claim   manual:newdb --by "Your Name"
    allele-sources approve manual:newdb --manifest newdb --by "Your Name" --note "why"
    allele-sources reject  manual:newdb --by "Your Name" --note "why"
    allele-sources defer   manual:newdb --by "Your Name" --until 2027-01-01 --note "why"
    allele-sources check   [--strict]

Approval is two commands on purpose. ``draft`` scaffolds a manifest with an empty
citation; the human fills it in; ``approve`` validates and stamps it. The gap between
them is where the review actually happens — see ``docs/source-policy.md``.
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import UTC, datetime
from pathlib import Path

from sqlalchemy import select

from ..db.models import Source, SourceCandidate
from ..db.session import session_scope
from ..settings import get_settings
from . import SourceError, registry, resolve, slugify
from . import candidates as cands
from .dossier import render, write_dossier

_MANIFEST_DIR = Path(__file__).resolve().parent


def _manifest_path(slug: str) -> Path:
    """Where a manifest is written: the overlay dir if set, else this package."""
    overlay = get_settings().source_manifest_dir
    return (Path(overlay) if overlay else _MANIFEST_DIR) / f"{slug}.toml"


def _toml_str(value: str) -> str:
    """Quote a string for TOML. Multi-line values use a literal block."""
    if "\n" in value:
        body = value.replace("\\", "\\\\").replace('"""', '\\"\\"\\"')
        return f'"""\\\n{body}"""'
    return '"' + value.replace("\\", "\\\\").replace('"', '\\"') + '"'


def _scaffold(candidate: SourceCandidate, slug: str) -> str:
    """A manifest with everything the channel told us, and nothing it did not.

    The citation is left blank deliberately: filling it in is the review.
    """
    licence = candidate.licence
    name = (candidate.title or slug).strip()
    lines = [
        "schema_version = 1",
        "",
        f"slug = {_toml_str(slug)}",
        f"name = {_toml_str(name)}",
        f"label = {_toml_str(name)}",
        'aliases = []          # every other spelling this source is written as',
        'role = "external"     # internal | external | derived',
        f"description = {_toml_str(candidate.title or 'TODO: one line, what this is')}",
    ]
    if candidate.url:
        lines.append(f"homepage_url = {_toml_str(candidate.url)}")
    lines += [
        f'curie_prefix = "{slug.upper()}"   # the namespace this source\'s ids use',
        "",
        "[licence]",
    ]
    if licence:
        lines += [
            'status = "declared"',
            f"name = {_toml_str(licence)}",
            'url = ""       # TODO: the page you read the licence on',
            "redistributable = true",
            'note = ""      # TODO: where and when you read it',
        ]
    else:
        lines += [
            'status = "none_declared"',
            "redistributable = false",
            'note = ""      # REQUIRED: where you looked and what you found',
        ]
    lines += [
        "",
        "[citation]",
        '# Fill this in. "verified" needs text; "none_available" needs a justification.',
        'status = "verified"',
        'text = ""',
    ]
    if candidate.doi:
        lines.append(f"doi = {_toml_str(candidate.doi)}")
    lines += [
        "",
        "[scope]",
        "species = []      # species._SPECIES tags"
        + (f"  # hint: {candidate.species_hint}" if candidate.species_hint else ""),
        "loci = []"
        + (f"         # hint: {candidate.loci_hint}" if candidate.loci_hint else ""),
        "segments = []",
        "categories = []",
        "",
        "[ingest]",
        'module = ""       # e.g. allele_server.ingest.curation_table',
        "published = true",
        "",
        "# [approval] is written by `allele-sources approve` — do not add it by hand.",
        "",
    ]
    return "\n".join(lines)


def _stamp_approval(path: Path, candidate: SourceCandidate, by: str, note: str) -> None:
    """Append the ``[approval]`` block. The resulting git diff is the audit record."""
    text = path.read_text(encoding="utf-8").rstrip("\n")
    marker = "# [approval] is written by"
    text = "\n".join(ln for ln in text.split("\n") if not ln.startswith(marker))
    block = [
        "",
        "[approval]",
        f"candidate = {_toml_str(candidate.ref)}",
        f"approved_by = {_toml_str(by)}",
        f"approved_at = {datetime.now(UTC).strftime('%Y-%m-%dT%H:%M:%SZ')}",
    ]
    if candidate.dossier_path:
        block.append(f"dossier = {_toml_str(candidate.dossier_path)}")
    if note:
        block.append(f"note = {_toml_str(note)}")
    path.write_text(text.rstrip("\n") + "\n" + "\n".join(block) + "\n", encoding="utf-8")


# --- subcommands ----------------------------------------------------------------------

def _cmd_add(args) -> int:
    with session_scope() as session:
        candidate, created = cands.add(
            session, channel=args.channel, external_id=args.id, title=args.title,
            url=args.url, doi=args.doi, licence=args.licence,
            species_hint=args.species, loci_hint=args.loci,
        )
        path = write_dossier(session, candidate)
        print(f"{'recorded' if created else 'refreshed'} {candidate.ref} "
              f"[{candidate.status}]")
        print(f"dossier: {path}")
    return 0


def _cmd_list(args) -> int:
    with session_scope() as session:
        rows = cands.listing(session, status=args.status, channel=args.channel)
        if args.json:
            print(json.dumps([{
                "ref": r.ref, "status": r.status, "title": r.title, "url": r.url,
                "doi": r.doi, "score": r.score, "manifest_id": r.manifest_id,
                "decided_by": r.decided_by,
                "decided_at": r.decided_at.isoformat() if r.decided_at else None,
                "decision_note": r.decision_note,
            } for r in rows], indent=2))
            return 0
        if not rows:
            print("no candidates" + (f" with status {args.status!r}" if args.status else ""))
            return 0
        for r in rows:
            who = f"  ({r.decided_by})" if r.decided_by else ""
            print(f"{r.status:<11} {r.ref:<46} {(r.title or '')[:44]}{who}")
    return 0


def _cmd_show(args) -> int:
    with session_scope() as session:
        candidate = cands.require(session, args.candidate)
        path = Path(candidate.dossier_path) if candidate.dossier_path else None
        # Prefer the stored file, but never fail because a cache was cleaned up: the
        # dossier is regenerable from the evidence column by design.
        print(path.read_text(encoding="utf-8") if path and path.is_file()
              else render(candidate))
    return 0


def _cmd_draft(args) -> int:
    with session_scope() as session:
        candidate = cands.require(session, args.candidate)
        slug = args.manifest or slugify(candidate.title or candidate.external_id)
        path = _manifest_path(slug)
        if path.exists() and not args.force:
            print(f"{path} already exists; pass --force to overwrite", file=sys.stderr)
            return 1
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(_scaffold(candidate, slug), encoding="utf-8")
        print(f"scaffolded {path}")
        print("Fill in the licence and the citation, then:")
        print(f"  allele-sources approve {candidate.ref} --manifest {slug} "
              f"--by 'Your Name' --note '...'")
    return 0


def _cmd_claim(args) -> int:
    with session_scope() as session:
        candidate = cands.claim(session, cands.require(session, args.candidate), by=args.by)
        print(f"{candidate.ref} claimed by {candidate.decided_by}")
    return 0


def _cmd_approve(args) -> int:
    with session_scope() as session:
        candidate = cands.require(session, args.candidate)
        who = cands.curator(args.by)
        slug = args.manifest
        path = _manifest_path(slug)
        if not path.is_file():
            raise SourceError(
                f"no manifest at {path}. Scaffold it first:\n"
                f"  allele-sources draft {candidate.ref} --manifest {slug}"
            )
        # Validate BEFORE stamping, so an unverifiable citation fails early and
        # legibly rather than leaving a half-approved manifest behind.
        registry(refresh=True)
        manifest = resolve(slug)
        if manifest is None:
            raise SourceError(f"{path} did not load as manifest {slug!r}")
        if manifest.approved:
            # The manifest and the candidate row are two halves of one decision, and
            # they can end up disagreeing: a manifest written by hand, or a run that
            # stamped the file and then died before `decide`. Bailing here left the
            # candidate stuck at `dossier` forever with no way to fix it through the
            # CLI -- the one command that should reconcile them refused to.
            #
            # So: if the candidate agrees, this is a genuine no-op. If it does not,
            # bring it into line with the manifest, which is the source of truth (its
            # [approval] block is the audit record, and it is what is committed).
            if candidate.status == "approved" and candidate.manifest_id == slug:
                print(f"{slug} is already approved by {manifest.approval.approved_by}",
                      file=sys.stderr)
                return 1
            was = candidate.status  # read before decide() mutates it
            cands.decide(
                session, candidate, status="approved",
                by=manifest.approval.approved_by or who,
                note=args.note or "reconciled with an already-approved manifest",
                manifest_id=slug,
            )
            print(f"{slug} was already approved by {manifest.approval.approved_by}; "
                  f"brought candidate {candidate.ref} into line "
                  f"({was!r} -> 'approved').")
            return 0

        _stamp_approval(path, candidate, who, args.note or "")
        registry(refresh=True)
        cands.decide(session, candidate, status="approved", by=who,
                     note=args.note, manifest_id=slug)
        print(f"approved {candidate.ref} as manifest {slug!r} (by {who})")
        print(f"  {path}")
        print("Commit that file: its [approval] block is the audit record.")
    return 0


def _cmd_reject(args) -> int:
    with session_scope() as session:
        candidate = cands.require(session, args.candidate)
        cands.decide(session, candidate, status="rejected", by=args.by, note=args.note)
        print(f"rejected {candidate.ref} (by {candidate.decided_by}). "
              f"It will not be proposed again on any channel.")
    return 0


def _cmd_defer(args) -> int:
    until = datetime.strptime(args.until, "%Y-%m-%d").replace(tzinfo=UTC)
    with session_scope() as session:
        candidate = cands.require(session, args.candidate)
        cands.decide(session, candidate, status="deferred", by=args.by,
                     note=args.note, until=until)
        print(f"deferred {candidate.ref} until {args.until} (by {candidate.decided_by})")
    return 0


def _cmd_check(args) -> int:
    """The gate report: what would fail, and what is merely pending."""
    problems: list[str] = []
    reg = registry(refresh=True)
    print(f"{len(reg)} manifests loaded.\n")

    unapproved = sorted(m.slug for m in reg.values() if not m.approved)
    if unapproved:
        problems.append(f"manifests with no [approval] block: {', '.join(unapproved)}")

    pending = sorted(
        m.slug for m in reg.values() if m.citation.status == "none_available"
    )
    undeclared = sorted(
        m.slug for m in reg.values() if m.licence.status == "none_declared"
    )
    print(f"citation none_available ({len(pending)}): {', '.join(pending) or '—'}")
    print(f"licence none_declared  ({len(undeclared)}): {', '.join(undeclared) or '—'}")

    with session_scope() as session:
        for name, manifest_id in session.execute(
            select(Source.name, Source.manifest_id).distinct()
        ):
            if resolve(name) is None:
                problems.append(f"source {name!r} in the database has no manifest")
            elif manifest_id is None:
                problems.append(
                    f"source {name!r} has no manifest_id — it predates the gate; "
                    f"re-run the ingest that creates it"
                )
            elif manifest_id not in reg:
                problems.append(
                    f"source {name!r} points at manifest {manifest_id!r}, which is gone"
                )
        orphans = session.scalars(
            select(SourceCandidate.manifest_id)
            .where(SourceCandidate.manifest_id.is_not(None))
            .distinct()
        ).all()
        for slug in orphans:
            if slug not in reg:
                problems.append(
                    f"candidate approved against manifest {slug!r}, which no longer exists"
                )
        # The other direction, which nothing checked: an approved manifest whose
        # candidate never reached a decided state. That is exactly the shape a
        # hand-written manifest leaves behind, and it reads on `list` as though the
        # source is still awaiting review when it has in fact been ingested.
        by_ref = {
            m.approval.candidate: m.slug
            for m in reg.values()
            if m.approved and m.approval and m.approval.candidate
        }
        retro: list[str] = []
        for ref, slug in sorted(by_ref.items()):
            row = cands.get(session, ref)
            if row is None:
                # Not drift: the sources that predate the review queue were
                # retro-approved during the source-policy migration and never had a
                # candidate row. Counting those as problems would make this report
                # cry wolf on nine of thirteen manifests, and a gate report nobody
                # trusts is worse than none.
                retro.append(slug)
            elif row.status not in ("approved", "ingested"):
                problems.append(
                    f"manifest {slug!r} is approved but candidate {ref!r} is still "
                    f"{row.status!r} — run: allele-sources approve {ref} "
                    f"--manifest {slug}"
                )
        if retro:
            print(f"retro-approved, no candidate row ({len(retro)}): "
                  f"{', '.join(sorted(retro))}")

    if problems:
        print("\nPROBLEMS:")
        for p in problems:
            print(f"  - {p}")
        return 1 if args.strict else 0
    print("\nNo problems.")
    return 0


# --- parser ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="allele-sources",
        description="Review and approve data sources. See docs/source-policy.md.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    def _cand(p):
        p.add_argument("candidate", help="<channel>:<external_id>, e.g. manual:newdb")

    def _by(p, required=True):
        p.add_argument("--by", default=None, required=False,
                       help="the curator recording this decision "
                            "(or set $AIRBABEL_CURATOR). Never inferred.")
        p.add_argument("--note", default=None,
                       help="why. Required in spirit; a rejection without a reason "
                            "costs the next reviewer the same work.")

    p = sub.add_parser("add", help="record a candidate you found yourself")
    p.add_argument("--channel", default="manual", choices=SourceCandidate.CHANNELS)
    p.add_argument("--id", required=True, help="external id, unique within the channel")
    p.add_argument("--title")
    p.add_argument("--url")
    p.add_argument("--doi")
    p.add_argument("--licence")
    p.add_argument("--species")
    p.add_argument("--loci")
    p.set_defaults(func=_cmd_add)

    p = sub.add_parser("list", help="list candidates")
    p.add_argument("--status", choices=SourceCandidate.STATUSES)
    p.add_argument("--channel", choices=SourceCandidate.CHANNELS)
    p.add_argument("--json", action="store_true")
    p.set_defaults(func=_cmd_list)

    p = sub.add_parser("show", help="print a candidate's dossier")
    _cand(p)
    p.set_defaults(func=_cmd_show)

    p = sub.add_parser("draft", help="scaffold a manifest from a candidate")
    _cand(p)
    p.add_argument("--manifest", default=None, help="manifest slug (default: from title)")
    p.add_argument("--force", action="store_true", help="overwrite an existing manifest")
    p.set_defaults(func=_cmd_draft)

    p = sub.add_parser("claim", help="mark a candidate as yours to review")
    _cand(p)
    _by(p)
    p.set_defaults(func=_cmd_claim)

    p = sub.add_parser("approve", help="approve a source (validates, then stamps)")
    _cand(p)
    _by(p)
    p.add_argument("--manifest", required=True, help="manifest slug to approve")
    p.set_defaults(func=_cmd_approve)

    p = sub.add_parser("reject", help="reject a candidate; it never resurfaces")
    _cand(p)
    _by(p)
    p.set_defaults(func=_cmd_reject)

    p = sub.add_parser("defer", help="revisit a candidate later")
    _cand(p)
    _by(p)
    p.add_argument("--until", required=True, metavar="YYYY-MM-DD")
    p.set_defaults(func=_cmd_defer)

    p = sub.add_parser("check", help="report manifest/database drift")
    p.add_argument("--strict", action="store_true", help="exit 1 on any problem")
    p.set_defaults(func=_cmd_check)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    try:
        return args.func(args)
    except SourceError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
