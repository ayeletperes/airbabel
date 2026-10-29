"""Species normalisation: scientific name / key -> short UID tag.

Kept independent of the ingest dependencies so the query server
never pulls them in. The tag is what appears in a UID, e.g. ``HUMAN`` in
``HUMAN-V-A3F2C10B``. Adding a species is a one-line edit to ``_SPECIES``.
"""

from __future__ import annotations

from dataclasses import dataclass

__all__ = [
    "Ontology",
    "SpeciesError",
    "common_name",
    "known_tags",
    "normalize_species",
    "ontology",
    "scientific_name",
    "taxon_curie",
]

# Human-friendly common names for display (spelled out to avoid ambiguity, e.g.
# "Rhesus macaque" rather than a bare "Rhesus").
_COMMON: dict[str, str] = {
    "HUMAN": "Human", "MOUSE": "Mouse", "RHESUS": "Rhesus macaque",
    "CYNO": "Cynomolgus macaque", "PIGTAIL": "Pig-tailed macaque",
    "CHIMP": "Chimpanzee", "GORILLA": "Gorilla", "ORANGSUM": "Sumatran orangutan",
    "ORANGBOR": "Bornean orangutan", "SMANGABEY": "Sooty mangabey",
    "BABOON": "Olive baboon", "LEMUR": "Ring-tailed lemur", "COW": "Cow",
    "SHEEP": "Sheep", "GOAT": "Goat", "PIG": "Pig", "HORSE": "Horse",
    "DOG": "Dog", "CAT": "Cat", "FERRET": "Ferret", "MINK": "American mink",
    "RABBIT": "Rabbit", "RAT": "Rat", "RATBLACK": "Black rat",
    "HAMSTER": "Golden hamster", "NMOLERAT": "Naked mole-rat", "ALPACA": "Alpaca",
    "CAMEL": "Dromedary camel", "DOLPHIN": "Bottlenose dolphin",
    "PLATYPUS": "Platypus", "CHICKEN": "Chicken", "ZEBRAFISH": "Zebrafish",
    "SALMON": "Atlantic salmon", "TROUT": "Rainbow trout", "GRASSCARP": "Grass carp",
    "PEROMYSCUS": "White-footed mouse",
    "SPRETUS": "Algerian mouse",
    "MUSPAHARI": "Shrew mouse", "MUSMINUTOIDES": "Pygmy mouse",
    "MUSCOOKII": "Cook's mouse", "MUSSAXICOLA": "Spiny mouse",
}


def common_name(tag: str) -> str:
    """Return a human-friendly common name for a canonical tag."""
    return _COMMON.get(tag.upper(), tag.capitalize())


class SpeciesError(ValueError):
    """Raised when a species value cannot be mapped to a known tag."""


@dataclass(frozen=True, slots=True)
class Ontology:
    """AIRR ``Ontology`` value: a CURIE id plus a human label."""

    id: str
    label: str


# tag -> (scientific name, NCBITaxon CURIE, accepted aliases). Aliases match
# case-insensitively with spaces/underscores normalised.
_SPECIES: dict[str, tuple[str, str, frozenset[str]]] = {
    "HUMAN": (
        "Homo sapiens",
        "NCBITAXON:9606",
        frozenset({"homo sapiens", "human", "hs"}),
    ),
    "MOUSE": (
        "Mus musculus",
        "NCBITAXON:10090",
        frozenset({"mus musculus", "mouse", "mm"}),
    ),
    "RHESUS": (
        "Macaca mulatta",
        "NCBITAXON:9544",
        frozenset({"macaca mulatta", "rhesus", "rhesus macaque", "rhmac"}),
    ),
    # IMGT/GENE-DB additional species (all-species ingest). NCBITaxon ids fetched
    # from NCBI taxonomy; aliases include the IMGT scientific name.
    "CYNO": ("Macaca fascicularis", "NCBITAXON:9541", frozenset({"macaca fascicularis", "cynomolgus", "cyno"})),
    "PIGTAIL": ("Macaca nemestrina", "NCBITAXON:9545", frozenset({"macaca nemestrina", "pig-tailed macaque"})),
    "CHIMP": ("Pan troglodytes", "NCBITAXON:9598", frozenset({"pan troglodytes", "chimpanzee", "chimp"})),
    "GORILLA": ("Gorilla gorilla", "NCBITAXON:9593", frozenset({"gorilla gorilla", "gorilla"})),
    "ORANGSUM": ("Pongo abelii", "NCBITAXON:9601", frozenset({"pongo abelii", "sumatran orangutan"})),
    "ORANGBOR": ("Pongo pygmaeus", "NCBITAXON:9600", frozenset({"pongo pygmaeus", "bornean orangutan"})),
    "SMANGABEY": ("Cercocebus atys", "NCBITAXON:9531", frozenset({"cercocebus atys", "sooty mangabey"})),
    "BABOON": ("Papio anubis", "NCBITAXON:9555", frozenset({"papio anubis", "olive baboon", "baboon"})),
    "LEMUR": ("Lemur catta", "NCBITAXON:9447", frozenset({"lemur catta", "ring-tailed lemur", "lemur"})),
    "COW": ("Bos taurus", "NCBITAXON:9913", frozenset({"bos taurus", "cattle", "bovine", "cow"})),
    "SHEEP": ("Ovis aries", "NCBITAXON:9940", frozenset({"ovis aries", "sheep"})),
    "GOAT": ("Capra hircus", "NCBITAXON:9925", frozenset({"capra hircus", "goat"})),
    "PIG": ("Sus scrofa", "NCBITAXON:9823", frozenset({"sus scrofa", "pig", "swine"})),
    "HORSE": ("Equus caballus", "NCBITAXON:9796", frozenset({"equus caballus", "horse"})),
    "DOG": ("Canis lupus", "NCBITAXON:9612", frozenset({"canis lupus", "canis lupus familiaris", "dog"})),
    "CAT": ("Felis catus", "NCBITAXON:9685", frozenset({"felis catus", "cat"})),
    "FERRET": ("Mustela putorius", "NCBITAXON:9668", frozenset({"mustela putorius", "mustela putorius furo", "ferret"})),
    "MINK": ("Neogale vison", "NCBITAXON:452646", frozenset({"neogale vison", "mustela vison", "american mink", "mink"})),
    "RABBIT": ("Oryctolagus cuniculus", "NCBITAXON:9986", frozenset({"oryctolagus cuniculus", "rabbit"})),
    "RAT": ("Rattus norvegicus", "NCBITAXON:10116", frozenset({"rattus norvegicus", "rat"})),
    "RATBLACK": ("Rattus rattus", "NCBITAXON:10117", frozenset({"rattus rattus", "black rat"})),
    "HAMSTER": ("Mesocricetus auratus", "NCBITAXON:10036", frozenset({"mesocricetus auratus", "golden hamster", "hamster"})),
    "NMOLERAT": ("Heterocephalus glaber", "NCBITAXON:10181", frozenset({"heterocephalus glaber", "naked mole-rat"})),
    "ALPACA": ("Vicugna pacos", "NCBITAXON:30538", frozenset({"vicugna pacos", "alpaca"})),
    "CAMEL": ("Camelus dromedarius", "NCBITAXON:9838", frozenset({"camelus dromedarius", "dromedary", "camel"})),
    "DOLPHIN": ("Tursiops truncatus", "NCBITAXON:9739", frozenset({"tursiops truncatus", "bottlenose dolphin", "dolphin"})),
    "PLATYPUS": ("Ornithorhynchus anatinus", "NCBITAXON:9258", frozenset({"ornithorhynchus anatinus", "platypus"})),
    "CHICKEN": ("Gallus gallus", "NCBITAXON:9031", frozenset({"gallus gallus", "chicken"})),
    "ZEBRAFISH": ("Danio rerio", "NCBITAXON:7955", frozenset({"danio rerio", "zebrafish"})),
    "SALMON": ("Salmo salar", "NCBITAXON:8030", frozenset({"salmo salar", "atlantic salmon", "salmon"})),
    "TROUT": ("Oncorhynchus mykiss", "NCBITAXON:8022", frozenset({"oncorhynchus mykiss", "rainbow trout", "trout"})),
    # OGRDB publishes an IGH/IGK/IGL set for the white-footed mouse; it is a distinct
    # genus from Mus and must not fold into MOUSE.
    "PEROMYSCUS": ("Peromyscus leucopus", "NCBITAXON:10041", frozenset({"peromyscus leucopus", "white-footed mouse", "peromyscus"})),
    # Added 2026-08-14: IMGT publishes 240 V/D/J records for it, which every reload
    # was counting under "unmapped:" and skipping. Taxon verified against NCBI
    # esummary (7959 -> Ctenopharyngodon idella, grass carp), not taken from a guess
    # at the pattern.
    "GRASSCARP": ("Ctenopharyngodon idella", "NCBITAXON:7959", frozenset({"ctenopharyngodon idella", "grass carp", "grasscarp"})),
    "SPRETUS": ("Mus spretus", "NCBITAXON:10096", frozenset({"mus spretus", "algerian mouse"})),
    "MUSPAHARI": ("Mus pahari", "NCBITAXON:10093", frozenset({"mus pahari", "shrew mouse"})),
    "MUSMINUTOIDES": ("Mus minutoides", "NCBITAXON:10105", frozenset({"mus minutoides", "pygmy mouse"})),
    "MUSCOOKII": ("Mus cookii", "NCBITAXON:10098", frozenset({"mus cookii", "cook's mouse"})),
    "MUSSAXICOLA": ("Mus saxicola", "NCBITAXON:10094", frozenset({"mus saxicola", "spiny mouse"})),
}

# Reverse lookup built once: every accepted spelling -> tag.
_LOOKUP: dict[str, str] = {}
for _tag, (_sci, _curie, _aliases) in _SPECIES.items():
    _LOOKUP[_tag.lower()] = _tag
    _LOOKUP[_sci.lower()] = _tag
    for _a in _aliases:
        _LOOKUP[_a] = _tag


def _key(value: str) -> str:
    return value.strip().lower().replace("_", " ")


def normalize_species(value: str) -> str:
    """Return the canonical UID tag for a species value.

    Accepts the scientific name (``"Homo sapiens"`` / ``"Homo_sapiens"``), a
    common name (``"human"``), or the tag itself (``"HUMAN"``).

    Raises:
        SpeciesError: If ``value`` is blank or not a known species. The error
            lists known tags so adding one is obvious.
    """
    if not isinstance(value, str) or not value.strip():
        raise SpeciesError("species value is empty")
    tag = _LOOKUP.get(_key(value))
    if tag is None:
        raise SpeciesError(
            f"unknown species {value!r}; known: {sorted(_SPECIES)}. "
            "Add it to allele_server.species._SPECIES."
        )
    return tag


def _entry(tag: str) -> tuple[str, str, frozenset[str]]:
    entry = _SPECIES.get(tag.upper())
    if entry is None:
        raise SpeciesError(f"unknown species tag {tag!r}; known: {sorted(_SPECIES)}")
    return entry


def scientific_name(tag: str) -> str:
    """Return the scientific name for a canonical tag."""
    return _entry(tag)[0]


def taxon_curie(tag: str) -> str:
    """Return the NCBITaxon CURIE for a canonical tag (e.g. ``NCBITAXON:9606``)."""
    return _entry(tag)[1]


def ontology(tag: str) -> Ontology:
    """Return the AIRR ``Ontology`` (NCBITaxon CURIE + scientific name) for a tag."""
    sci, curie, _ = _entry(tag)
    return Ontology(id=curie, label=sci)


def known_tags() -> tuple[str, ...]:
    """Return the known species tags, sorted."""
    return tuple(sorted(_SPECIES))
