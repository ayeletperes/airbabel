"""The source manifest schema.

One manifest per data source, validated on load. The manifest is the single place a
source's identity, licence and citation are written down; the ``sources`` table row,
the CURIE prefix, the ``/sources`` page and the export attribution line are all derived
from it (see ``docs/source-policy.md``).

The validators here are the policy, expressed as code. Two rules do the real work:

* a citation is either **verified** with text, or **none_available** with a written
  justification. There is no third state, because "not filled in yet" is how nine
  ``TODO(citation)`` entries survived for a year.
* an undeclared licence forces ``redistributable = false``. We may ingest and link to a
  source whose terms we could not read; we may not republish it.

``extra="forbid"`` everywhere: a typo'd key is a loud failure, not a silently ignored one.
"""

from __future__ import annotations

import re
from datetime import datetime
from typing import Literal

from pydantic import BaseModel, ConfigDict, model_validator

__all__ = [
    "Approval",
    "Citation",
    "Ingest",
    "Licence",
    "Scope",
    "SourceManifest",
]

# Slugs are filenames, CURIE-adjacent identifiers and CLI arguments at once, so keep
# them boring: lowercase, no spaces, no punctuation beyond the underscore.
_SLUG_RE = re.compile(r"^[a-z0-9][a-z0-9_]*$")


class _Base(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)


class Licence(_Base):
    """What the source publishes for *its* data, and whether we may republish it.

    Note that a sequence's underlying GenBank accession is the primary public deposit
    and is not owned by any curator database (see ``queries._accession_xrefs``) — this
    field is about the curated *set*, not the sequences' original deposit.
    """

    status: Literal["declared", "none_declared"]
    name: str | None = None
    url: str | None = None
    redistributable: bool = False
    note: str = ""

    @model_validator(mode="after")
    def _check(self) -> Licence:
        if self.status == "declared" and not (self.name or "").strip():
            raise ValueError("licence.status='declared' requires licence.name")
        if self.status == "none_declared":
            if not self.note.strip():
                raise ValueError(
                    "licence.status='none_declared' requires a note saying where you "
                    "looked and what you found"
                )
            if self.redistributable:
                raise ValueError(
                    "cannot set redistributable=true with no declared licence; "
                    "an undeclared licence is not permission"
                )
        return self


class Citation(_Base):
    """How to credit the source. Only ever a citation we can verify.

    ``none_available`` is an honest answer for an internal set whose manuscript is still
    in preparation. It is not a way to skip the work: the note is mandatory, and
    ``allele-sources check`` reports how many sources are using it.
    """

    status: Literal["verified", "none_available"]
    text: str = ""
    doi: str | None = None
    pmid: str | None = None
    note: str = ""

    @model_validator(mode="after")
    def _check(self) -> Citation:
        if self.status == "verified" and not self.text.strip():
            raise ValueError("citation.status='verified' requires citation.text")
        if self.status == "none_available":
            if not self.note.strip():
                raise ValueError(
                    "citation.status='none_available' requires a note explaining why no "
                    "citation exists"
                )
            if self.text.strip():
                raise ValueError(
                    "citation.status='none_available' but citation.text is set; if you "
                    "have a citation, verify it and set status='verified'"
                )
        return self


class Scope(_Base):
    """What the source covers. Empty tuple means "unrestricted", not "nothing"."""

    species: tuple[str, ...] = ()
    loci: tuple[str, ...] = ()
    segments: tuple[str, ...] = ()
    categories: tuple[str, ...] = ()


class Ingest(_Base):
    """How the source is loaded, and whether it is part of the published dataset.

    ``published=False`` covers a source we hold but do not advertise as a contributor —
    the mined literature layers are aggregation over someone else's articles, not a
    germline set, and crediting them alongside IMGT would misrepresent both.
    """

    module: str
    published: bool = True
    note: str = ""
    #: May new data from this source be ingested without a per-release review?
    #:
    #: Scope is deliberately narrow: **additions only** — new sequences and new names.
    #: Three event classes always queue for a human regardless of this flag, because
    #: none of them can be applied safely by a machine:
    #:
    #: * a **new species** needs a code change (``species._SPECIES`` carries the tag,
    #:   scientific name and NCBITaxon), so it cannot be auto-applied at all;
    #: * a **removal** is the open semantic question of what "removed" means versus
    #:   ``archived``, and auto-applying one silently drops data we serve;
    #: * a **rename** needs a decision about which name is canonical and whether the
    #:   superseded one stays searchable.
    #:
    #: Default ``False``: a source is reviewed per release unless someone deliberately
    #: decided otherwise, and that decision is reviewed here alongside the licence and
    #: the citation rather than in a separate list that would drift from this file.
    auto_update: bool = False

    @model_validator(mode="after")
    def _check(self) -> Ingest:
        if not self.module.strip():
            raise ValueError(
                "ingest.module must name the module that loads this source "
                "(e.g. allele_server.ingest.curation_table)"
            )
        return self


class Approval(_Base):
    """The human decision. Its presence is what unlocks ingestion.

    Written by ``allele-sources approve``; the git diff of this block is the audit record.
    """

    candidate: str
    approved_by: str
    approved_at: datetime
    dossier: str | None = None
    note: str = ""

    @model_validator(mode="after")
    def _check(self) -> Approval:
        if not self.approved_by.strip():
            raise ValueError("approval.approved_by must name a person")
        return self


class SourceManifest(_Base):
    """A reviewed, approved data source."""

    schema_version: int
    slug: str
    name: str
    label: str
    aliases: tuple[str, ...] = ()
    role: Literal["internal", "external", "derived"]
    description: str
    homepage_url: str | None = None
    curie_prefix: str | None = None
    curie_url_template: str | None = None
    licence: Licence
    citation: Citation
    scope: Scope = Scope()
    ingest: Ingest
    approval: Approval | None = None

    @model_validator(mode="after")
    def _check(self) -> SourceManifest:
        if self.schema_version != 1:
            raise ValueError(f"unsupported schema_version {self.schema_version!r}; expected 1")
        if not _SLUG_RE.match(self.slug):
            raise ValueError(
                f"slug {self.slug!r} must be lowercase alphanumeric with underscores"
            )
        if not self.name.strip():
            raise ValueError("name must not be blank (it is the sources.name written to the DB)")
        if self.curie_url_template and not self.curie_prefix:
            raise ValueError("curie_url_template needs a curie_prefix")
        if self.curie_prefix and self.curie_prefix != self.curie_prefix.strip():
            raise ValueError(f"curie_prefix {self.curie_prefix!r} has surrounding whitespace")
        if ":" in (self.curie_prefix or ""):
            raise ValueError(f"curie_prefix {self.curie_prefix!r} must not contain a colon")
        return self

    @property
    def approved(self) -> bool:
        return self.approval is not None

    @property
    def keys(self) -> tuple[str, ...]:
        """Every string that should resolve to this manifest, casefolded.

        The alias list is what makes ``RhGLDB+`` and ``RHGLDB`` the same source to the
        resolver and to the discovery deduper.
        """
        return tuple(dict.fromkeys(
            k.casefold() for k in (self.slug, self.name, self.label, *self.aliases) if k
        ))
