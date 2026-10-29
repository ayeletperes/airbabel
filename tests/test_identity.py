"""Tests for the identity primitives: sequence, uid, species, naming, curie."""

from __future__ import annotations

import pytest

from allele_server import curie, naming, species
from allele_server.sequence import (
    SequenceError,
    aa_hash,
    hash_aa,
    hash_sequence,
    is_ambiguous,
    normalize,
    reverse_complement_hash,
    translate_v,
)
from allele_server.uid import (
    Category,
    UidError,
    base32_of_hash,
    mint_protein_uid,
    mint_uid,
)

# A real-ish in-frame V fragment (multiple of 3).
NT = "CAGGTGCAGCTGGTGCAGTCTGGGGCTGAGGTGAAGAAGCCTGGG"


# --- sequence -------------------------------------------------------------

def test_normalize_strips_gaps_and_uppercases():
    assert normalize("ac-g.t\n") == "ACGT"


def test_normalize_rejects_rna_and_junk():
    with pytest.raises(SequenceError):
        normalize("ACGU")
    with pytest.raises(SequenceError):
        normalize("ACGZ")
    with pytest.raises(SequenceError):
        normalize("....")


def test_hash_sequence_is_normalisation_invariant():
    assert hash_sequence("acgt") == hash_sequence("AC-G.T")
    assert len(hash_sequence("ACGT")) == 64


def test_is_ambiguous():
    assert is_ambiguous("ACGTN")
    assert not is_ambiguous("ACGT")


def test_reverse_complement_hash():
    # Non-palindrome: rc differs from the forward hash.
    assert reverse_complement_hash("AAAC") != hash_sequence("AAAC")
    # rc(AAAC) == GTTT.
    assert reverse_complement_hash("AAAC") == hash_sequence("GTTT")
    # ACGT is its own reverse complement (palindrome).
    assert reverse_complement_hash("ACGT") == hash_sequence("ACGT")


def test_aa_hash_matches_hash_aa_of_translation():
    assert aa_hash(NT) == hash_aa(translate_v(NT))


def test_hash_aa_ignores_stop_and_gaps():
    assert hash_aa("QVQL*") == hash_aa("qvql")
    assert hash_aa("Q.V-Q L") == hash_aa("QVQL")


# --- uid ------------------------------------------------------------------

def test_base32_of_hash_length_and_validation():
    token = base32_of_hash(hash_sequence(NT))
    assert token.isupper() and "=" not in token
    with pytest.raises(UidError):
        base32_of_hash("not-hex")
    with pytest.raises(UidError):
        base32_of_hash("ab")  # too short


def test_mint_uid_deterministic_and_namespaced():
    h = hash_sequence(NT)
    a = mint_uid("HUMAN", Category.CODING, "V", h)
    b = mint_uid("HUMAN", Category.CODING, "V", h)
    assert a == b
    assert a.startswith("HUMAN-V")
    assert len(a) == len("HUMAN-V") + 8  # coding default N=8
    # different species / segment / category -> different UID
    assert mint_uid("MOUSE", Category.CODING, "V", h) != a
    assert mint_uid("HUMAN", Category.CODING, "D", h) != a
    assert mint_uid("HUMAN", Category.LEADER, None, h) != a
    assert len(mint_uid("HUMAN", Category.LEADER, None, h)) == len("HUMAN-L") + 5


def test_mint_uid_requires_segment_for_coding():
    with pytest.raises(UidError):
        mint_uid("HUMAN", Category.CODING, None, hash_sequence(NT))


def test_mint_uid_extends_on_collision():
    h = hash_sequence(NT)
    other = hash_sequence("TTTTAAAACCCCGGGG")
    base = mint_uid("HUMAN", Category.CODING, "V", h)

    # Simulate that `base` is already taken by a *different* hash.
    def exists_hash(uid: str) -> str | None:
        return other if uid == base else None

    extended = mint_uid("HUMAN", Category.CODING, "V", h, exists_hash)
    assert extended != base
    assert len(extended) == len(base) + 1
    # Idempotent: if the same hash already holds it, reuse without extending.
    assert mint_uid("HUMAN", Category.CODING, "V", h, lambda u: h if u == base else None) == base


def test_mint_protein_uid_groups_by_aa():
    h_aa = aa_hash(NT)
    puid = mint_protein_uid("HUMAN", "V", h_aa)
    assert puid.startswith("HUMAN-VP-")
    assert len(puid) == len("HUMAN-VP-") + 6
    # A synonymous nt change (same protein) -> same protein UID.
    syn = NT[:-3] + "GGA"  # GGG/GGA both Gly
    if translate_v(syn) == translate_v(NT):
        assert mint_protein_uid("HUMAN", "V", aa_hash(syn)) == puid


# --- species --------------------------------------------------------------

def test_normalize_species_accepts_many_spellings():
    for v in ["Homo sapiens", "homo_sapiens", "human", "HUMAN", "hs"]:
        assert species.normalize_species(v) == "HUMAN"
    assert species.normalize_species("Mus musculus") == "MOUSE"
    assert species.normalize_species("rhesus macaque") == "RHESUS"


def test_species_ontology_and_curie():
    assert species.taxon_curie("HUMAN") == "NCBITAXON:9606"
    onto = species.ontology("MOUSE")
    assert onto.id == "NCBITAXON:10090" and onto.label == "Mus musculus"


def test_unknown_species_raises():
    with pytest.raises(species.SpeciesError):
        species.normalize_species("drosophila melanogaster")


# --- naming ---------------------------------------------------------------

def test_detect_and_parse_schemes():
    assert naming.detect_scheme("IGHV1-2*02") is naming.Scheme.IUIS
    assert naming.detect_scheme("IGHV1-2*02_a4g") is naming.Scheme.VDJBASE
    assert naming.detect_scheme("IGHV4-39*i03") is naming.Scheme.OGRDB
    assert naming.detect_scheme("IGHV1-CO5H*00") is naming.Scheme.IGLABEL
    assert naming.detect_scheme("nonsense") is None


def test_parse_name_pads_allele_and_keeps_raw():
    n = naming.parse_name("IGHV1-2*2", naming.Scheme.IUIS)
    assert n.allele == "02" and n.raw == "IGHV1-2*2"
    assert n.canonical() == "IGHV1-2*02"


def test_segment_and_locus_helpers():
    assert naming.segment_of("IGHV1-2*02") == "V"
    assert naming.locus_of("IGHV1-2*02") == "IGH"
    assert naming.locus_of("IGKV1D-13*01") == "IGK"
    assert naming.locus_of("nonsense") is None


def test_tcr_names_parse():
    # TR loci (TRA/TRB/TRG/TRD) parse the same grammar as IG; IMGT nomenclature.
    assert naming.segment_of("TRBV20-1*01") == "V"
    assert naming.locus_of("TRBV20-1*01") == "TRB"
    assert naming.locus_of("TRAV1-2*01") == "TRA"
    assert naming.locus_of("TRGV9*01") == "TRG"
    assert naming.locus_of("TRDV2*01") == "TRD"
    assert naming.detect_scheme("TRBV20-1*01") is naming.Scheme.IUIS
    assert naming.scheme_label("TRBV20-1*01") == "imgt"


# --- curie ----------------------------------------------------------------

def test_curie_build_and_split():
    assert curie.build("GENBANK", "X02897") == "GENBANK:X02897"
    assert curie.split("GENBANK:X02897") == ("GENBANK", "X02897")


def test_curie_handles_multi_colon_ogrdb_ref():
    ref = "OGRDB:Human_IGH:2021.11"
    assert curie.validate(ref) == ref
    assert curie.split(ref) == ("OGRDB", "Human_IGH:2021.11")


def test_curie_rejects_unknown_prefix():
    assert not curie.is_curie("FOO:bar")
    with pytest.raises(curie.CurieError):
        curie.build("FOO", "bar")


def test_dual_locus_trav_dv_names():
    # IMGT shared TRA/TRD V genes carry a combined name and must parse: primary
    # locus TRA, segment V, and a second locus TRD recovered separately.
    for n in ("TRAV14/DV4*01", "TRAV13-4/DV7*01", "TRAV14D-3/DV8*01",
              "TRAV15-1/DV6-1*01", "TRAV15D-1/DV6D-1*01"):
        assert naming.locus_of(n) == "TRA", n
        assert naming.segment_of(n) == "V", n
        assert naming.dual_locus_of(n) == "TRD", n
    # ordinary single-locus names carry no second locus
    for n in ("TRAV23*01", "TRDV1*01", "IGHV1-2*02", "TRBV2*01"):
        assert naming.dual_locus_of(n) is None, n
    # the combined name round-trips through parse/canonical (the /DV is kept)
    p = naming.parse_name("TRAV14/DV4*01", naming.detect_scheme("TRAV14/DV4*01"))
    assert p.canonical() == "TRAV14/DV4*01"
    assert p.dual == "/DV4"


def test_subgroup_recognition_is_strict_in_both_directions():
    """``as_subgroup`` asks "is this a subgroup?"; ``subgroup_of`` asks "which one?".

    A bare subgroup is deliberately not a parseable *name* (``_NAME_RE`` requires
    ``*<allele>``), so it needs its own recogniser. Keeping that recogniser strict is
    what stops a subgroup lookup from quietly reinterpreting a precise gene query.
    """
    # genuine subgroups, across IG and TR, including the letter/slash family forms
    for s in ("IGHV1", "IGHV10", "IGKV2", "TRBV6", "IGHD3", "TRAV14"):
        assert naming.as_subgroup(s) == s, s
    # case is normalised, whitespace tolerated
    assert naming.as_subgroup("  ighv1 ") == "IGHV1"
    # NOT subgroups: genes, alleles, junk, wrong types
    for s in ("IGHV1-69", "IGHV1-69*01", "IGHV", "V1", "", "IGXV1", "IGHV1-"):
        assert naming.as_subgroup(s) is None, s
    assert naming.as_subgroup(None) is None

    # the inverse: which subgroup a full name belongs to
    assert naming.subgroup_of("IGHV1-69*01") == "IGHV1"
    assert naming.subgroup_of("IGHV10-1*01") == "IGHV10"   # not IGHV1
    assert naming.subgroup_of("IGHV6*01") == "IGHV6"       # no gene part
    assert naming.subgroup_of("TRAV14/DV4*01") == "TRAV14"
    assert naming.subgroup_of("not-a-name") is None


def test_constant_gene_names_resolve_to_segment_c():
    """IMGT does not name heavy-chain constant genes with a 'C' segment letter -- it
    names them by isotype (IGHG1, IGHM, IGHA1), and IGKC/TRAC carry a C with no family
    number. Both shapes fall outside the V/D/J grammar, so segment_of used to return
    None for a row we store as segment 'C'."""
    for name, locus in [
        ("IGHG1*01", "IGH"), ("IGHG4*01", "IGH"), ("IGHM*01", "IGH"),
        ("IGHA1*02", "IGH"), ("IGHE*01", "IGH"), ("IGHGP*01", "IGH"),
        ("IGKC*01", "IGK"), ("IGLC2*01", "IGL"),
        ("TRAC*01", "TRA"), ("TRBC1*01", "TRB"), ("TRDC*01", "TRD"),
        # Non-human isotypes: IgY (birds), IgZ/IgT (teleosts), duplicated IgG.
        ("IGHY*01", "IGH"), ("IGHZ*01", "IGH"), ("IGHG2A*01", "IGH"),
        ("IGHG5-1*01", "IGH"), ("IGHDD1P*01", "IGH"),
    ]:
        assert naming.segment_of(name) == "C", name
        assert naming.locus_of(name) == locus, name


def test_igd_the_isotype_does_not_swallow_igd_the_diversity_segment():
    """The one real collision: IGHD*01 is the constant delta gene, IGHD3-10*01 is a
    diversity gene. The V/D/J grammar is tried first and requires a family number, so
    the presence of one decides. Verified against the bulk reference FASTA: no
    D-REGION record is named IGHD*NN."""
    assert naming.segment_of("IGHD*01") == "C"
    assert naming.segment_of("IGHD*02") == "C"
    assert naming.segment_of("IGHD3-10*01") == "D"
    assert naming.segment_of("IGHD1-1*01") == "D"
    assert naming.segment_of("IGHD7-27*01") == "D"


def test_the_constant_grammar_does_not_widen_anything_else():
    for name, seg in [("IGHV1-69*01", "V"), ("IGHJ4*02", "J"),
                      ("TRAV14/DV4*01", "V"), ("IGKV1-5*01", "V")]:
        assert naming.segment_of(name) == seg, name
    for junk in ["not-a-name", "IGH*01", "", "IGHQ"]:
        assert naming.segment_of(junk) is None, junk
