#!/usr/bin/env python3
"""AIRBabel REST API — a runnable tour of the six things you'll actually do.

Standard library only: no install step, no API key, no registration.

    python airbabel_example.py                        # against the public server
    python airbabel_example.py http://127.0.0.1:8000  # against your own instance

Each function is independent — copy the one you need. Every call returns plain
JSON, so `json.load` is the whole client.

Data is aggregated from IMGT, OGRDB and published curated sets; attribute the
originating sources on reuse (see <BASE>/sources and <BASE>/api/export).
"""

from __future__ import annotations

import json
import sys
import urllib.error
import urllib.parse
import urllib.request

BASE = "https://airbabel.org"

# A 57 nt fragment of a human IGHV1-18 allele: the kind of partial sequence a
# primer-trimmed read or a figure in a paper actually gives you.
FRAGMENT = "CAGGTTCAGCTGGTGCAGTCTGGAGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAG"


def _call(path: str, params: dict | None = None, body: dict | None = None):
    """GET (or POST, when `body` is given) one endpoint and return parsed JSON.

    Raises on anything that is not 2xx, with the server's own message: the API
    answers 404 for "no such record" and 422 for "that is not a sequence / not a
    PMID", and those are answers worth reading rather than swallowing.
    """
    url = BASE.rstrip("/") + path
    if params:
        url += "?" + urllib.parse.urlencode(params)
    data = json.dumps(body).encode() if body is not None else None
    headers = {"content-type": "application/json"} if data else {}
    req = urllib.request.Request(url, data=data, headers=headers)
    try:
        with urllib.request.urlopen(req, timeout=60) as resp:
            return json.load(resp)
    except urllib.error.HTTPError as exc:
        detail = exc.read().decode("utf-8", "replace")[:400]
        raise RuntimeError(f"{exc.code} {exc.reason} for {url}\n{detail}") from exc


# --- 0. which build am I talking to? -------------------------------------------
# Record this next to any result you keep. `release` names the data; the digests make
# it checkable, and `similarity_digest` is what makes a reported % identity
# reproducible -- it pins the metric, the per-segment policy and the cache floors.

def dataset_version() -> dict:
    return _call("/version")


# --- 1. resolve a name ---------------------------------------------------------
# Any name, from any scheme or source, including names a reference set has since
# renamed or withdrawn. `matched_by` says how it matched; `species=` scopes a name
# that several species reuse.

def resolve_name(name: str, species: str | None = None) -> dict:
    params = {"query": name}
    if species:
        params["species"] = species
    return _call("/api/resolve", params)


# --- 2. resolve a nucleotide sequence ------------------------------------------
# Identity is the sequence hash, so an exact sequence match is exact regardless of
# what anyone calls it. `partial=true` scores over the best-matching window, which
# is what you want for a truncated read: a contained fragment scores ~1.0.

def resolve_sequence(seq: str, partial: bool = True, species: str | None = None) -> dict:
    params: dict = {"partial": str(partial).lower()}
    if species:
        params["species"] = species
    return _call("/api/search/sequence", params, body={"sequence": seq})


# --- 3. retrieve a sequence record ---------------------------------------------
# The full record: sequence, translated V-REGION, every name with its scheme and
# source, cross-references and literature evidence.

def get_record(uid: str) -> dict:
    return _call(f"/api/sequences/{uid}")


# --- 4. get similar alleles ----------------------------------------------------
# Ranked by % identity (normalised Levenshtein). `to_species` is the cross-species
# nearest-neighbour lookup: "the closest human alleles to this mouse one".

def find_similar(uid: str, min_identity: float = 0.9, kind: str = "nt",
                 to_species: str | None = None) -> dict:
    params: dict = {"min_identity": min_identity, "kind": kind}
    if to_species:
        params["to_species"] = to_species
        params["same_species"] = "false"
    return _call(f"/api/sequences/{uid}/similar", params)


# --- 5. AIRR AlleleDescription -------------------------------------------------
# The same record under AIRR field names, with species and sources as CURIEs —
# what you want when the output feeds an AIRR-compliant pipeline.

def airr_description(name_or_uid: str, species: str | None = None) -> dict:
    params = {"query": name_or_uid, "view": "airr"}
    if species:
        params["species"] = species
    return _call("/api/resolve", params)


# --- 6. query a PMID -----------------------------------------------------------
# The reverse of a record's reference list: given a paper, which sequences does it
# stand behind? Each hit carries its evidence tier (verified / sequence / name),
# and the per-tier totals are reported separately, never summed.

def alleles_for_pmid(pmid: str, species: str | None = None) -> dict:
    params = {"species": species} if species else None
    return _call(f"/api/pmid/{pmid}", params)


def main(base: str = BASE) -> None:
    global BASE
    BASE = base
    print(f"# AIRBabel at {BASE}")
    ver = dataset_version()
    ds = ver.get("dataset")
    if ds:
        print(f"# dataset release {ds['release']}  "
              f"content {ds['content_digest'][:12]}…  "
              f"distance {ds['similarity_digest'][:12]}…")
        if not ds["similarity_current"]:
            print("# WARNING: the server's distance parameters have changed since it "
                  "was stamped; scores are not the ones this release describes.")
    else:
        print("# dataset: unstamped -- cite the access date, not a version")
    print()

    hit = resolve_name("IGHV1-69*01", species="human")
    native = hit["record"]["native"]
    uid = native["uid"]
    print(f"1. resolve name  IGHV1-69*01 -> {uid}  "
          f"({native['species']} {native['locus']} {native['length']} nt, "
          f"matched_by={hit['matched_by']})")
    print("   names:", ", ".join(n["name"] for n in native["allele_names"]))

    # A name a reference set stopped using still resolves, because former
    # designations hang off the sequence hash rather than off the name.
    old = resolve_name("IGKV1D-17*02")
    print(f"   retired name IGKV1D-17*02 -> {old['record']['native']['uid']}")

    res = resolve_sequence(FRAGMENT, partial=True)
    top = res["hits"][0]
    print(f"\n2. resolve sequence  {len(FRAGMENT)} nt fragment -> {res['count']} hits; "
          f"best {top['label']} ({top['uid']}) at {top['score']:.3f} [{top['matched_by']}]")

    rec = get_record(uid)["native"]
    print(f"\n3. record  {uid}  protein_uid={rec.get('protein_uid')}  "
          f"functional={rec.get('functional')}")

    sim = find_similar(uid, min_identity=0.95, kind="nt")
    print(f"\n4. similar (nt, >=95%, same species): {sim['count']}")
    for h in sim["hits"][:5]:
        print(f"     {h['score']:.3f}  {h['uid']}  {h['label']}")

    # Cross-species: the closest human alleles to a mouse allele.
    cross = find_similar("MOUSE-V3GPCFPTO", min_identity=0.8, to_species="human")
    print(f"   mouse IGHV5-17*01 -> closest human: {cross['count']}")
    for h in cross["hits"][:3]:
        print(f"     {h['score']:.3f}  {h['uid']}  {h['label']}")

    # Scope to a species: IGHV1-69*01 is a live name on a *different* sequence in
    # chicken, gorilla, mouse and rhesus, so an unscoped resolve is entitled to
    # answer with any of them. This is the whole point of the service, and it bites
    # scripts first.
    airr = airr_description("IGHV1-69*01", species="human")
    print(f"\n5. AIRR AlleleDescription  label={airr['label']}  "
          f"locus={airr['locus']}  species={airr['species']['id']}  "
          f"aliases={airr['aliases']}")

    pm = alleles_for_pmid("29163486")
    print(f"\n6. PMID {pm['pmid']}: {pm['count']} alleles attested "
          f"(verified {pm['n_verified']}, sequence {pm['n_sequence']}, "
          f"name-only {pm['n_name']})")
    print(f"   {pm['url']}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else BASE)
