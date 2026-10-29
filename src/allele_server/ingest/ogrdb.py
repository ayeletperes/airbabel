"""OGRDB: fetch the germline sets from its v2 API, and ingest them.

:func:`fetch_to_cache` writes each set, split by chain, into a directory tree::

    {cache}/{species}/vdj/airrc_{species}_{CHAIN}.fasta        OGRDB (ungapped)

and :func:`ingest_cache` reads that tree (no network) and upserts into the DB.
CHAIN spans V, D, J and IGHC (see ``_CHAINS``). IMGT is not loaded here: it comes
from ``imgt_genedb``.
"""

from __future__ import annotations

import logging
from pathlib import Path

from sqlalchemy.orm import Session

from ..curie import build as curie_build
from ..naming import locus_of, scheme_label, segment_of
from ..species import SpeciesError, normalize_species
from . import imgt_constant, loader
from .fasta import parse_record, read_fasta

_LOG = logging.getLogger(__name__)

# Every chain we fetch and ingest. D exists only in the heavy chain; the light
# chains have V and J alone. Widening this map is why `segment` below is derived
# from the allele name rather than hardcoded -- the two changes have to move
# together, or a J record lands in the table labelled V.
_CHAINS = {
    "IGHV": "IGH", "IGHD": "IGH", "IGHJ": "IGH", "IGHC": "IGH",
    "IGKV": "IGK", "IGKJ": "IGK",
    "IGLV": "IGL", "IGLJ": "IGL",
}

# OGRDB publishes constants as a separate set (germline_set_id 9606.IGHC) carrying two
# forms per allele, distinguished only by the name: a bare name is the **membrane**
# form, and a ``_SC`` suffix is the **secreted** one. They share their first ~1000 nt
# and then diverge at the 3' end (IGHA1*01 is 1211 nt, IGHA1*01_SC is 1058).
#
# This database stores C as the secreted region (``imgt_constant.SUBTYPE``), with the
# membrane form explicitly deferred -- see TODO(c-membrane). So only the ``_SC`` records
# are admitted. Taking the bare ones would file a membrane transcript under the same
# segment as a secreted one with nothing to tell them apart, which is worse than not
# having them.
#
# They do not hash-match our IMGT-assembled C, and the reason is systematic rather than
# biological: ours carries one extra 5' nucleotide, so OGRDB's secreted form is exactly
# our sequence minus its first base (1059/1058, 1131/1130, 1359/1358, ... across every
# shared allele). They are therefore distinct sequences under this database's own
# identity rule and are stored as such. Re-cutting ours to match would change 753
# existing UIDs to make a join look tidy, and identity is the sequence.
_OGRDB_SECRETED_SUFFIX = "_SC"
_SPECIES_KEYS = ("human", "mouse", "rhesus")
_API = "https://ogrdb.airr-community.org/api_v2"


def ingest_ogrdb(session: Session, cache_dir: str | Path, species_keys=None) -> int:
    """Ingest OGRDB V/D/J/C alleles as corroborating aliases/attribution. Rows touched.

    ``species_keys=None`` means *every species present in the cache*, discovered from
    the directory tree rather than from a hardcoded tuple. The tuple is what kept
    trout and Peromyscus out even after ``fetch_ogrdb_v`` had downloaded them: two
    lists of species had to agree, and they silently did not.
    """
    cache_dir = Path(cache_dir)
    # A missing cache is a wrong path, not "OGRDB has nothing": say so rather than
    # load zero rows and let a reload look complete.
    if not cache_dir.is_dir():
        raise FileNotFoundError(f"OGRDB cache not found: {cache_dir} (run --fetch)")
    if species_keys is None:
        species_keys = sorted(
            d.name for d in cache_dir.iterdir()
            if d.is_dir() and (d / "vdj").is_dir()
        )
    source = loader.get_or_create_source(
        session, name="OGRDB", version=None, curie_prefix="OGRDB",
        note="OGRDB / AIRR-C curated sets",
    )
    n = 0
    for sp in species_keys:
        tag = normalize_species(sp)
        for chain, locus in _CHAINS.items():
            path = cache_dir / sp / "vdj" / f"airrc_{sp}_{chain}.fasta"
            for header, seq in read_fasta(path):
                rec = parse_record(header, seq)
                if not rec.allele or "*" not in rec.allele or not rec.sequence:
                    continue
                segment = segment_of(rec.allele) or chain[3:4]
                subtype = None
                if segment == "C":
                    # Secreted only; see _OGRDB_SECRETED_SUFFIX.
                    if not rec.allele.endswith(_OGRDB_SECRETED_SUFFIX):
                        continue
                    subtype = imgt_constant.SUBTYPE
                seq_row = loader.upsert_coding(
                    session, species=tag, segment=segment,
                    locus=locus_of(rec.allele) or locus, nt=rec.sequence,
                    subtype=subtype,
                )
                loader.add_name(session, seq_row, source, rec.allele, scheme_label(rec.allele))
                loader.add_attribution(
                    session, seq_row, source,
                    original_id=curie_build("OGRDB", rec.allele),
                    original_sequence=rec.sequence,
                )
                # OGRDB provides no AA; fall back to translation if still unset.
                # V only -- there is no protein for a D or J segment.
                if segment == "V":
                    loader.set_v_aa(session, seq_row)
                n += 1
    return n


def ogrdb_species(timeout: int = 60) -> dict[str, tuple[str, str]]:
    """``{our tag: (taxon id, scientific name)}`` for every species OGRDB serves.

    **Asked, not assumed.** This used to be a three-entry dict (human, mouse, rhesus)
    written when those were the whole database, so the two sets OGRDB added since --
    *Oncorhynchus mykiss* and *Peromyscus leucopus*, 440 sequences between them --
    were never fetched. Nothing reported that: the loop simply had no key for them.
    Reading ``/germline/species`` means the next species OGRDB publishes arrives on
    the next run instead of waiting for someone to notice.

    A species OGRDB serves that ``allele_server.species`` does not know is skipped
    with a warning rather than guessed at: minting a tag here would create an
    identity, and identities are not the fetcher's to invent.
    """
    import json
    import urllib.request

    req = urllib.request.Request(
        _API + "/germline/species",
        headers={"User-Agent": "allele-server/0.1"},
    )
    payload = json.loads(urllib.request.urlopen(req, timeout=timeout).read().decode())
    out: dict[str, tuple[str, str]] = {}
    for entry in payload.get("species", []):
        label = (entry.get("label") or "").strip()
        taxon = str(entry.get("id") or "").strip()
        if not (label and taxon):
            continue
        try:
            out[normalize_species(label)] = (taxon, label)
        except SpeciesError:
            _LOG.warning(
                "OGRDB serves %r (taxon %s), which allele_server.species does not "
                "know; skipping. Add it to _SPECIES to ingest it.", label, taxon,
            )
    return out


def fetch_ogrdb_v(cache_dir: str | Path, species_keys=None) -> None:
    """Fetch OGRDB V/D/J/C FASTAs for every species OGRDB serves, via its v2 API.

    Discovers the species list and each species' germline sets, takes each set's
    latest release, and writes it split by chain as ``airrc_{sp}_{locus}{segment}.fasta``.

    ``species_keys=None`` (the default) means *everything OGRDB has*. Pass a tuple to
    restrict it. The old default was the three-species tuple, which silently capped
    the fetch -- see :func:`ogrdb_species`.

    A set that cannot be fetched is logged and skipped, so one bad set does not cost
    the rest; its previous files (if any) stay on disk.
    """
    import json
    import urllib.parse
    import urllib.request

    from .imgt_genedb import _read_fasta_text

    discovered = ogrdb_species()
    if species_keys is not None:
        wanted = {normalize_species(k) for k in species_keys}
        discovered = {k: v for k, v in discovered.items() if k in wanted}
    # Cache paths are keyed on the LOWERCASE species key ("human", "trout"), which is
    # the layout every other pass reads and writes. Canonical tags are upper-case, so
    # keying paths on them silently created a second, parallel tree.
    sci = {tag.lower(): label for tag, (_, label) in discovered.items()}
    tax = {tag.lower(): taxon for tag, (taxon, _) in discovered.items()}

    def get_text(path: str) -> str:
        req = urllib.request.Request(_API + path, headers={"User-Agent": "allele-server/0.1"})
        with urllib.request.urlopen(req, timeout=60) as resp:  # noqa: S310 - fixed host
            return resp.read().decode("utf-8")

    cache_dir = Path(cache_dir)
    for sp in sorted(sci):
        try:
            sets = json.loads(get_text(f"/germline/sets/{tax[sp]}"))
        except Exception as e:  # noqa: BLE001
            _LOG.warning("OGRDB: no set list for %s (taxon %s): %s", sp, tax[sp], e)
            continue
        key = next((k for k in sets if isinstance(sets[k], list)), None)
        for entry in sets.get(key, []) if key else []:
            locus = entry.get("locus")
            if locus not in ("IGH", "IGK", "IGL"):
                continue
            set_id = urllib.parse.quote(entry["germline_set_id"], safe=".")
            try:
                latest = json.loads(get_text(f"/germline/set/{set_id}/latest"))
                version = str(latest["GermlineSet"][0]["release_version"])
                # OGRDB serves the human sets' extended form under an _ex suffix.
                ex = "_ex" if sci[sp] == "Homo sapiens" else ""
                text = get_text(f"/germline/set/{set_id}/{version}/ungapped{ex}")
            except Exception as e:  # noqa: BLE001 - per set; one bad set must not cost the rest
                _LOG.warning("OGRDB: could not fetch %s %s (%s): %s",
                             sp, entry.get("germline_set_name"), set_id, e)
                continue
            # The germline set arrives complete, V/D/J/C together: split it into the
            # per-chain filenames ingest_ogrdb reads.
            by_seg: dict[str, dict[str, str]] = {}
            for header, seq in _read_fasta_text(text):
                name = header.split()[0]
                seg = segment_of(name) or (name[3:4].upper() if len(name) > 3 else "")
                if f"{locus}{seg}" in _CHAINS:
                    by_seg.setdefault(seg, {})[name] = seq.upper()
            vdj = cache_dir / sp / "vdj"
            vdj.mkdir(parents=True, exist_ok=True)
            for seg, recs in by_seg.items():
                (vdj / f"airrc_{sp}_{locus}{seg}.fasta").write_text(
                    "".join(f">{n}\n{q}\n" for n, q in recs.items()))


def ingest_cache(session: Session, cache_dir: str | Path, species_keys=None) -> dict[str, int]:
    """Ingest OGRDB from a cache dir. Returns row counts.

    IMGT is not read from here any more: ``imgt_genedb`` loads it for every species
    from IMGT's bulk FASTA, so this tree holds only what OGRDB serves.
    """
    return {"ogrdb": ingest_ogrdb(session, cache_dir, species_keys)}


def fetch_to_cache(cache_dir: str | Path, species_keys=None) -> None:
    """Populate the cache with OGRDB V/D/J/C FASTAs (needs network; stdlib only)."""
    fetch_ogrdb_v(cache_dir, species_keys)
