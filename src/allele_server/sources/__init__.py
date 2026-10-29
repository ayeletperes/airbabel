"""The source registry: one reviewed, approved manifest per data source.

This package is the single source of truth for what AIRBabel aggregates, under what
terms, and how to credit each contributor. It replaces four hand-maintained registries
that nothing kept in sync — ``sources_meta.SOURCE_META``, the source half of
``curie._PREFIXES``, ``export.ATTRIBUTION``, and the metadata columns of the ``sources``
table — with ``<slug>.toml`` files validated on load.

Read ``docs/source-policy.md`` for the policy those manifests encode. In short: a source
cannot be ingested until a human has approved a manifest carrying a licence and a
citation we can verify, and ``ingest/loader.get_or_create_source`` refuses anything else.

TOML rather than YAML because ``tomllib`` is stdlib in 3.12 and the query server's
dependency set stays frozen (``pyproject.toml`` keeps ingest-only libraries in an extra);
data files rather than a Python dict because approval is a diff a human reads.
"""

from __future__ import annotations

import tomllib
from importlib import resources
from pathlib import Path

from pydantic import ValidationError

from ..settings import get_settings
from ._model import Approval, Citation, Ingest, Licence, Scope, SourceManifest

__all__ = [
    "Approval",
    "Citation",
    "Ingest",
    "Licence",
    "Scope",
    "SourceError",
    "SourceManifest",
    "published",
    "registry",
    "resolve",
    "resolve_strict",
    "slugify",
]


class SourceError(ValueError):
    """Raised when a source has no approved manifest, or a manifest is invalid.

    Sibling of ``SequenceError`` / ``CurieError`` / ``SpeciesError``: carries enough
    context to act on, never a silent fallback.
    """


# slug -> manifest, and every alias casefolded -> slug. Built once; `refresh=True`
# rebuilds (the tests install a temporary manifest directory between cases).
_MANIFESTS: dict[str, SourceManifest] | None = None
_INDEX: dict[str, str] = {}


def slugify(name: str) -> str:
    """A manifest slug from a source name: ``"RhGLDB+"`` -> ``"rhgldb"``.

    Used to suggest a slug in error messages and by ``allele-sources draft``; the
    authoritative slug is always the one written in the manifest.
    """
    out = [c if c.isalnum() else "_" for c in name.casefold()]
    slug = "".join(out).strip("_")
    while "__" in slug:
        slug = slug.replace("__", "_")
    return slug or "source"


def _manifest_files() -> list[Path]:
    """Every ``*.toml`` to load: the ones shipped in this package, then the overlay.

    The overlay (``settings.source_manifest_dir``) lets a deployment hold a source the
    repository does not — an unpublished collaborator set, say. A duplicate slug across
    the two is an error, not a silent override: which one won would be invisible.
    """
    here = resources.files(__package__)
    paths = sorted(Path(str(p)) for p in here.iterdir() if str(p).endswith(".toml"))
    overlay = get_settings().source_manifest_dir
    if overlay:
        extra = Path(overlay)
        if not extra.is_dir():
            raise SourceError(f"source_manifest_dir {overlay!r} is not a directory")
        paths.extend(sorted(extra.glob("*.toml")))
    return paths


def _load_one(path: Path) -> SourceManifest:
    try:
        raw = tomllib.loads(path.read_text(encoding="utf-8"))
    except (OSError, tomllib.TOMLDecodeError) as exc:
        raise SourceError(f"{path.name}: cannot read manifest: {exc}") from exc
    try:
        manifest = SourceManifest.model_validate(raw)
    except ValidationError as exc:
        raise SourceError(f"{path.name}: invalid manifest:\n{exc}") from exc
    if manifest.slug != path.stem:
        raise SourceError(
            f"{path.name}: slug {manifest.slug!r} does not match the filename; "
            f"rename the file to {manifest.slug}.toml or fix the slug"
        )
    return manifest


def registry(*, refresh: bool = False) -> dict[str, SourceManifest]:
    """Every manifest, keyed by slug. Memoised.

    Loads and validates each file, checks that no two manifests claim the same slug,
    name, alias or private CURIE prefix, and registers each source's prefix with
    ``allele_server.curie`` so ``KIMDB:...`` resolves without a second registry.

    A manifest may name a *standard* namespace as its ``curie_prefix`` — KIARVA
    identifies its rows by ``DOI:``, the mined layers by ``PMID:``, the NCBI pass by
    ``GENBANK:``. Those are shared public namespaces, so they are neither owned nor
    exclusive; only a source's own prefix is registered and checked for collisions.

    Raises:
        SourceError: on an unreadable, invalid or conflicting manifest.
    """
    global _MANIFESTS
    if _MANIFESTS is not None and not refresh:
        return _MANIFESTS

    # Deferred: curie imports nothing from here at module level (a top-level import
    # would be a cycle), and this is what makes manifest prefixes resolvable.
    from .. import curie as _curie

    manifests: dict[str, SourceManifest] = {}
    index: dict[str, str] = {}
    owned: dict[str, str] = {}
    for path in _manifest_files():
        manifest = _load_one(path)
        if manifest.slug in manifests:
            raise SourceError(f"duplicate source slug {manifest.slug!r} ({path})")
        manifests[manifest.slug] = manifest
        for key in manifest.keys:
            owner = index.get(key)
            if owner is not None and owner != manifest.slug:
                raise SourceError(
                    f"{manifest.slug!r} and {owner!r} both claim the name/alias {key!r}; "
                    f"a source name must resolve to exactly one manifest"
                )
            index[key] = manifest.slug
        prefix = manifest.curie_prefix
        if prefix and prefix not in _curie.STANDARD_PREFIXES:
            owner = owned.get(prefix)
            if owner is not None:
                raise SourceError(
                    f"{manifest.slug!r} and {owner!r} both claim CURIE prefix {prefix!r}"
                )
            owned[prefix] = manifest.slug

    _MANIFESTS = manifests
    _INDEX.clear()
    _INDEX.update(index)

    # Rebuild rather than add to: a refresh after a manifest is renamed or removed must
    # not leave its namespace behind for ``curie.split`` to keep matching.
    _curie.reset_source_prefixes()
    for slug, prefix in ((v, k) for k, v in owned.items()):
        manifest = manifests[slug]
        _curie.register_prefix(
            prefix, manifest.description or manifest.label, manifest.curie_url_template
        )
    return _MANIFESTS


def resolve(name_or_slug: str | None) -> SourceManifest | None:
    """The manifest for a source name, slug or alias; None if we hold none.

    Matching is case-insensitive and covers every spelling the manifest declares, which
    is how ``RhGLDB+`` and ``RHGLDB`` reach the same record.
    """
    if not name_or_slug:
        return None
    mans = registry()
    slug = _INDEX.get(name_or_slug.strip().casefold())
    return mans.get(slug) if slug else None


def resolve_strict(name_or_slug: str) -> SourceManifest:
    """As ``resolve``, but raises with the remedy when there is no manifest.

    Raises:
        SourceError: if no manifest matches, or the one that does is unapproved.
    """
    manifest = resolve(name_or_slug)
    if manifest is None:
        slug = slugify(name_or_slug)
        raise SourceError(
            f"no manifest for source {name_or_slug!r}. Every source must be reviewed and "
            f"approved before it can be ingested (see docs/source-policy.md).\n"
            f"  known sources: {sorted(registry())}\n"
            f"  to add it:\n"
            f"    allele-sources add --channel manual --id {slug} --title {name_or_slug!r}\n"
            f"    allele-sources draft manual:{slug} --manifest {slug}\n"
            f"    $EDITOR src/allele_server/sources/{slug}.toml   # licence + citation\n"
            f"    allele-sources approve manual:{slug} --manifest {slug} --by 'Your Name'"
        )
    if not manifest.approved:
        raise SourceError(
            f"source {manifest.name!r} has a manifest but no [approval] block: it has not "
            f"been reviewed (see docs/source-policy.md).\n"
            f"  allele-sources approve manual:{manifest.slug} "
            f"--manifest {manifest.slug} --by 'Your Name' --note '...'"
        )
    return manifest


def published() -> list[SourceManifest]:
    """Approved, published manifests, ordered by label — the credited contributors.

    Excludes ``[ingest] published = false``: the mined literature layers aggregate other
    people's articles rather than contributing a germline set, and crediting them
    alongside IMGT would misrepresent both.
    """
    return sorted(
        (m for m in registry().values() if m.ingest.published and m.approved),
        key=lambda m: m.label.casefold(),
    )
