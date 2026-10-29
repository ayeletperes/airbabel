"""CURIE (compact URI) helpers for source and cross-reference identifiers.

AIRR Community tools (OGRDB, VDJbase, the user's ``digby_backend``) identify
sources and accessions with CURIEs like ``NCBITAXON:9606``, ``GENBANK:X02897``,
``OGRDB:Human_IGH:2021.11``. We store cross-references in this form so records
are interoperable without a remapping step.

The registry has two halves. ``_STANDARD`` are public namespaces that belong to
nobody in particular — a taxonomy id, an accession, a DOI — and are listed here.
**Source** prefixes (``OGRDB``, ``KIMDB``, ``HUSA`` …) are *not* listed here: each
one is declared by its manifest in ``allele_server.sources`` and registered on load
via :func:`register_prefix`. Two hand-maintained lists is how ``RHGLDB`` here and
``RhGLDB+`` in the metadata table drifted apart; there is now one.
"""

from __future__ import annotations

__all__ = [
    "STANDARD_PREFIXES",
    "CurieError",
    "build",
    "is_curie",
    "known_prefixes",
    "register_prefix",
    "resolve_url",
    "split",
    "validate",
]


class CurieError(ValueError):
    """Raised when a CURIE is malformed or uses an unknown prefix."""


# Standard namespaces -> short description. The local id follows the *last* colon, so
# prefixes themselves may contain neither leading/trailing spaces nor a colon.
_STANDARD: dict[str, str] = {
    "NCBITAXON": "NCBI Taxonomy id (species)",
    "GENBANK": "GenBank/INSDC nucleotide accession",
    "ENA": "European Nucleotide Archive accession",
    "DOI": "Digital Object Identifier",
    "ORCID": "researcher identifier",
    "NCBIGENE": "NCBI Gene id",
    "HGNC": "HUGO Gene Nomenclature Committee id",
    "PMID": "PubMed identifier",
}

# Public namespaces nobody owns. A source manifest may declare one of these as its
# ``curie_prefix`` (KIARVA identifies its rows by DOI, the mined layers by PMID), which
# is why the registry does not treat them as exclusive.
STANDARD_PREFIXES: frozenset[str] = frozenset(_STANDARD)

# Resolvable prefixes -> URL template. A template without ``{id}`` is a constant landing
# page (OGRDB has no per-set permalink), which is still better than no link at all.
_URL_TEMPLATES: dict[str, str] = {
    "GENBANK": "https://www.ncbi.nlm.nih.gov/nuccore/{id}",
    "ENA": "https://www.ncbi.nlm.nih.gov/nuccore/{id}",
    "NCBIGENE": "https://www.ncbi.nlm.nih.gov/gene/{id}",
    "HGNC": "https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:{id}",
    "DOI": "https://doi.org/{id}",
    "ORCID": "https://orcid.org/{id}",
    "PMID": "https://pubmed.ncbi.nlm.nih.gov/{id}",
}

# Live registry: the standard namespaces plus whatever the manifests register.
_PREFIXES: dict[str, str] = dict(_STANDARD)

# Which of the above came from a manifest, so a registry reload can drop the ones whose
# manifest is gone. Without this a renamed or deleted source leaves a ghost namespace
# that ``split()`` would still match.
_FROM_MANIFEST: set[str] = set()

# Guards the deferred manifest import below. Set *before* loading, so the
# register_prefix() calls that loading makes cannot recurse back into it.
_sources_loaded = False


def _ensure_sources_loaded() -> None:
    """Load the source manifests so their CURIE prefixes are known.

    Imported inside the function on purpose: ``allele_server.sources`` imports this
    module, so a top-level import here would be a cycle.
    """
    global _sources_loaded
    if _sources_loaded:
        return
    _sources_loaded = True
    try:
        from . import sources as _sources

        _sources.registry()
    except Exception:
        _sources_loaded = False
        raise


def reset_source_prefixes() -> None:
    """Drop every manifest-supplied prefix. Called by the registry before reloading.

    The standard namespaces are untouched — they are not the manifests' to remove.
    """
    for key in _FROM_MANIFEST:
        _PREFIXES.pop(key, None)
        _URL_TEMPLATES.pop(key, None)
    _FROM_MANIFEST.clear()


def register_prefix(prefix: str, description: str, url_template: str | None = None) -> None:
    """Register a source-supplied CURIE prefix. Called by the manifest registry.

    Idempotent: re-registering the same prefix with the same template is a no-op, so a
    registry refresh is free.

    Raises:
        CurieError: If the prefix is blank, contains a colon, or collides with a
            standard namespace (a source cannot claim ``DOI:`` as its own).
    """
    key = prefix.strip()
    if not key or ":" in key:
        raise CurieError(f"invalid CURIE prefix {prefix!r}")
    if key in _STANDARD:
        raise CurieError(
            f"CURIE prefix {key!r} is a standard namespace; a source may cite it but "
            f"cannot own it"
        )
    _PREFIXES[key] = description
    _FROM_MANIFEST.add(key)
    if url_template:
        _URL_TEMPLATES[key] = url_template


def known_prefixes() -> tuple[str, ...]:
    """Return the registered CURIE prefixes (standard namespaces + source manifests)."""
    _ensure_sources_loaded()
    return tuple(_PREFIXES)


def build(prefix: str, local_id: str) -> str:
    """Construct a CURIE from a known prefix and a local id.

    Raises:
        CurieError: If the prefix is unknown or the local id is blank.
    """
    _ensure_sources_loaded()
    if prefix not in _PREFIXES:
        raise CurieError(f"unknown CURIE prefix {prefix!r}; known: {sorted(_PREFIXES)}")
    local = local_id.strip()
    if not local:
        raise CurieError(f"empty local id for prefix {prefix!r}")
    return f"{prefix}:{local}"


def split(curie: str) -> tuple[str, str]:
    """Split a CURIE into ``(prefix, local_id)``.

    The prefix is matched greedily against the known registry so that values
    containing colons (``OGRDB:Human_IGH:2021.11``) keep the remainder intact.

    Raises:
        CurieError: If no known prefix matches.
    """
    if not isinstance(curie, str) or ":" not in curie:
        raise CurieError(f"not a CURIE: {curie!r}")
    _ensure_sources_loaded()
    # Longest first: a manifest may register a prefix that extends another
    # (``IMGT`` vs ``IMGT/GENE-DB``), and insertion order must not decide the winner.
    for prefix in sorted(_PREFIXES, key=len, reverse=True):
        if curie.startswith(prefix + ":"):
            return prefix, curie[len(prefix) + 1:]
    raise CurieError(f"unknown CURIE prefix in {curie!r}; known: {sorted(_PREFIXES)}")


def validate(curie: str) -> str:
    """Return ``curie`` unchanged if well-formed, else raise.

    Raises:
        CurieError: If the prefix is unknown or the local id is blank.
    """
    prefix, local = split(curie)
    if not local.strip():
        raise CurieError(f"empty local id in CURIE {curie!r}")
    return curie


def is_curie(value: str) -> bool:
    """Return whether ``value`` is a well-formed CURIE with a known prefix."""
    try:
        split(value)
    except CurieError:
        return False
    return True


def resolve_url(curie: str) -> str | None:
    """Return an external web URL for a CURIE, or None if not linkable.

    Source prefixes resolve through the ``curie_url_template`` their manifest declares;
    standard namespaces through ``_URL_TEMPLATES``. A template with no ``{id}`` is a
    landing page rather than a permalink.
    """
    try:
        prefix, local = split(curie)
    except CurieError:
        return None
    template = _URL_TEMPLATES.get(prefix)
    if not template:
        return None
    return template.format(id=local) if "{id}" in template else template
