"""Curated showcase examples for ``/examples``.

Each entry points at a **real record** and names the one thing it demonstrates.
Together they are the argument for why this database exists: allele names are
unstable (renamed, retired, removed), ambiguous (the same string in two species),
and non-unique (one sequence, several designations), so the *sequence* is the key.

Entries are resolved defensively at render time (``queries.example_showcase``): if a
UID no longer exists after a reload, the example is skipped rather than 404-ing.
UIDs are deterministic from the sequence, so they survive reingest of the same data.

**``expect`` is not optional decoration.** Every entry that asserts a *status* must
state it there in checkable form, because a status is derived and moves under
reload. This page once described an allele as archived and gone from IMGT while the
live record was neither — the UID resolved, so nothing caught it. ``_example_holds``
now drops an example whose record stopped supporting its prose. Keep the two in
step: if you edit a notice, edit the expectation with it.

Do not hardcode counts ("one of 731 archived alleles") in a notice. They are derived
too, and prose cannot be checked.
"""

from __future__ import annotations

from typing import Any

EXAMPLES: list[dict[str, Any]] = [
    {
        "id": "duplicate-designations",
        "kind": "allele",
        "uid": "HUMAN-VQLNRZF6P",
        "expect": {"names": ["IGHV1-69*01", "IGHV1-69D*01"], "archived": False},
        "title": "One sequence, two official IMGT names",
        "lede": "Names are not unique in either direction.",
        "notice": "IGHV1-69 and its duplicated paralogue IGHV1-69D carry the identical sequence at "
                  "*01: this one piece of DNA is recorded as both IGHV1-69*01 and IGHV1-69D*01. Search "
                  "either name and AIRBabel resolves to this same record. The converse also occurs: the string "
                  "IGHV1-69*01 is also used for entirely different sequences in chicken, gorilla, mouse "
                  "and rhesus. One sequence, many names; one name, many sequences. Only the sequence "
                  "identifies the allele.",
    },
    {
        "id": "former-name",
        "kind": "allele",
        "uid": "HUMAN-VACLUDJCT",
        "expect": {"archived": False, "former_names": ["IGKV1D-17*02"],
                   "names": ["IGKV1-17*03"], "has_newer_version": False},
        "title": "A former IMGT name still finds its sequence",
        "lede": "Names are unstable; the sequence is not.",
        "notice": "IMGT called this sequence IGKV1D-17*02 in the 2013-01-20 release and has called it "
                  "IGKV1-17*03 from 2015-03-19 onwards, a different gene and a different allele number "
                  "for the identical DNA. Searching the old name still lands here, because the release "
                  "snapshots attach every former designation to the sequence's hash rather than to a "
                  "name. One sequence, two official names, held one after the other.",
    },
    {
        "id": "record-revision",
        "kind": "allele",
        "uid": "MOUSE-VLVYQVPRP",
        "expect": {"archived": False, "length": 321, "prev_version_count": 2,
                   "has_newer_version": False, "names": ["IGHV1-2*01"]},
        "title": "Same name, three different sequences",
        "lede": "A reference record is not fixed, even when its name is.",
        "notice": "IMGT has published mouse IGHV1-2*01 over three sequence extents: 297 nt up to "
                  "2015-03-19, 294 nt up to 2022-01-17, and the 321 nt record shown here, which contains "
                  "both earlier ones. The name never changed and the allele was never withdrawn; only how "
                  "much sequence IMGT publishes under it did. Because identity is the hash, each extent is "
                  "its own record. The chain is listed below under Earlier records, so it remains visible rather "
                  "than looking like three unrelated alleles. Note also that IGHV1-2*01 is a live name in "
                  "human too, on an entirely different sequence.",
    },
    {
        "id": "archived",
        "kind": "allele",
        "uid": "HUMAN-V3LLA7BIC",
        "expect": {"archived": True, "retired": True, "has_newer_version": False,
                   "former_names": ["IGHV1-69*07"]},
        "title": "An allele designation no longer present in IMGT",
        "lede": "Reference-set membership changes across releases.",
        "notice": "IGHV1-69*07 was in every IMGT release from 2013-01-20 to 2024-01-07 and is in none "
                  "since. No other record carries the name and no current source publishes the sequence, "
                  "so this is a change in current reference-set membership rather than a rename or a change "
                  "in published sequence extent. It is deferred from coverage counts and default search, but "
                  "remains resolvable for literature follow-up.",
    },
    {
        "id": "alias-conventions",
        "kind": "allele",
        "uid": "HUMAN-V35RLULJX",
        "expect": {"archived": False, "names": ["IGHV1-2*09"]},
        "title": "The same allele under two naming conventions",
        "lede": "Different sources use different naming grammars.",
        "notice": "IGHV1-2*09 (an IMGT designation) and IGHV1-2*04_S3434 (a source's mutation-suffix name) "
                  "denote the identical sequence. AIRBabel resolves them to one record by hashing the sequence, "
                  "not by parsing the names.",
    },
    {
        "id": "same-protein",
        "kind": "allele",
        "uid": "HUMAN-VOST4CB5R",
        "expect": {"archived": False, "names": ["IGHV3-25*05"]},
        "title": "One silent substitution, identical amino-acid sequence",
        "lede": "Nucleotide identity and amino-acid identity are different questions.",
        "notice": "IGHV3-25*05 and IGHV3-25*03 differ at a single nucleotide (position 279, C to G). The "
                  "substitution is silent, so their V-REGION amino-acid sequences are identical and both "
                  "carry the same protein UID. Nucleotide-level and amino-acid-level identity answer "
                  "different questions, so the record keeps both.",
    },
    {
        "id": "literature-ogrdb",
        "kind": "allele",
        "uid": "HUMAN-V4ZQZGA5J",
        "expect": {"archived": False, "names": ["IGHV1-2*02"]},
        "title": "Evidence: sequence-level literature support and OGRDB curation",
        "lede": "Not every citation is equal.",
        "notice": "IGHV1-2*02 is the germline behind VRC01-class HIV antibodies. Its references separate papers "
                  "whose text or supplement contains this exact sequence (a hash match, robust across species) "
                  "from papers that merely mention the name (a lead, marked with an asterisk). External records "
                  "link to OGRDB's curated, evidence-backed entry for the same allele.",
    },
    {
        "id": "pmid-reverse",
        "kind": "pmid",
        "pmid": "29163486",
        # Checked against live data, like the allele examples: the paper must still
        # attest at least this many alleles at this evidence tier, or the example is
        # dropped. The floor is deliberately below the real count (235 on dev and lab
        # as of 2026-08-06) so an ordinary reload does not silently delete the entry,
        # while a pass that lost the citations entirely still would.
        "expect": {"min_hits": 100, "kind": "verified"},
        "title": "Start from a paper, not from a name",
        "lede": "A citation can be used as a query, not only as a record annotation.",
        "notice": "Ramesh et al. 2017 (Front Immunol 8:1407) assembled the rhesus macaque "
                  "immunoglobulin loci de novo. Searching its PubMed ID returns AIRBabel records "
                  "linked to that publication, each labelled with its evidence tier; "
                  "here as a curator-verified citation carried by the MUSA rhesus set. That is the "
                  "reverse of an allele page's reference list: given a paper, which AIRBabel allele records "
                  "are linked to it? Archived records are returned on purpose, because a "
                  "paper citing a name that IMGT has since retired is exactly the case where the "
                  "citation still has to lead somewhere.",
    },
    {
        "id": "partial-sequence",
        "kind": "search",
        "query": "CAGGTTCAGCTGGTGCAGTCTGGAGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAG",
        "mode": "sequence",
        "partial": True,
        "title": "Paste a fragment, not the whole allele",
        "lede": "User-supplied repertoire sequences may be partial.",
        "notice": "This is a 57 nt substring of IGHV1-18*01. With partial matching on, it resolves at 100% "
                  "local identity. You do not need a full-length sequence, and you do not need to know any name.",
    },
]


# One-click example inputs for the search page: what to type, per mode. These are
# *inputs*, not claims about a result, so unlike EXAMPLES they carry no ``expect``
# block — the same reason the ``kind: "search"`` entries above are not validated
# either. Each is a real query verified against the loaded data; a reload that
# changed the answer would change the result page, not make the page lie.
#
# ``params`` is the query string for "/" (or the explicit ``url`` for the two that
# are not searches). Keep the list short: it sits under the search box.
SEARCH_EXAMPLES: list[dict[str, Any]] = [
    {
        "label": "human IGHV1-69*01",
        "what": "allele name",
        # Scoped to human on purpose. Unscoped, this name matches five species on
        # five different sequences and the search answers with a disambiguation
        # picker -- true, and the point of the service, but not what a first click
        # labelled "allele name" should land on. The ambiguity gets its own worked
        # example on /tutorial.
        "note": "Resolves to one record that carries two official IMGT names. The "
                "species is part of the example on purpose: unscoped, this name "
                "matches five species on five different sequences.",
        "params": {"q": "IGHV1-69*01", "mode": "auto", "species": "human"},
    },
    {
        "label": "HUMAN-VQLNRZF6P",
        "what": "server UID",
        "note": "The same record as the chip before it, reached by its identifier "
                "instead of a name. A UID is minted from the sequence hash, so it "
                "is unambiguous where a name is not.",
        "params": {"q": "HUMAN-VQLNRZF6P", "mode": "auto"},
    },
    {
        "label": "IGKV1D-17*02",
        "what": "retired name",
        "note": "A name IMGT stopped using in 2015; it still resolves to its sequence, "
                "now published as IGKV1-17*03.",
        "params": {"q": "IGKV1D-17*02", "mode": "auto"},
    },
    {
        "label": "partial sequence",
        "what": "nucleotide",
        "note": "A partial sequence, no name: which allele contains it? "
                "(partial matching on)",
        "params": {
            "q": "CAGGTTCAGCTGGTGCAGTCTGGAGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAG",
            "mode": "sequence", "partial": "1",
        },
    },
    {
        "label": "protein fragment",
        "what": "amino acid",
        "note": "A protein fragment: AIRBabel records whose translated V-REGION "
                "contains it. (partial matching on)",
        "params": {
            "q": "QVQLVQSGAEVKKPGSSVKVSCKASGGTFSSYAISWVRQAPGQGLEWMGGIIPIFGTANY",
            "mode": "aa", "partial": "1",
        },
    },
    {
        "label": "PMID 29163486",
        "what": "publication",
        "note": "AIRBabel records linked to this paper, by evidence tier.",
        "params": {"q": "29163486", "mode": "pmid"},
    },
    {
        "label": "mouse IGHV5-17*01 → human",
        "what": "cross-species",
        "note": "The closest human alleles to a mouse allele, ranked by identity.",
        "url": "/seq/MOUSE-V3GPCFPTO/similar?to_species=HUMAN&min=0.8&kind=nt",
    },
]
