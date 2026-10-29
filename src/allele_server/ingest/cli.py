"""``allele-ingest`` CLI: download (optional) and load the DB."""

from __future__ import annotations

import argparse
import os

from ..db.session import session_scope
from ..settings import get_settings
from . import ogrdb


def pmid_list(raw: list[str] | None) -> list[str] | None:
    """The PMIDs named on the command line, or ``None`` for "none were named".

    Accepts both shapes a user will type -- space-separated after the flag and
    comma-separated in one token -- because the flag is how a paper that the name
    probe cannot see gets mined at all, and a silently-dropped id would look like
    the paper simply had no matches. ``None`` (not ``[]``) is the empty answer, so
    it can be handed straight to ``enrich_literature_supp(pmids=...)``, where
    ``None`` means "fall back to the discovered set".
    """
    out = [p.strip() for item in (raw or []) for p in item.split(",") if p.strip()]
    return list(dict.fromkeys(out)) or None


def main(argv: list[str] | None = None) -> int:
    settings = get_settings()
    parser = argparse.ArgumentParser(description="Ingest V alleles into allele_server.")
    parser.add_argument(
        "--cache-dir", default=settings.ogrdb_cache_dir,
        help="OGRDB download cache directory",
    )
    parser.add_argument(
        "--fetch", action="store_true",
        help="download OGRDB before ingesting (needs network); "
             "IMGT is downloaded by --imgt-genedb / --reload",
    )
    parser.add_argument(
        "--check-releases", action="store_true",
        help="report only: which IMGT releases the archive offers that we do not hold, "
             "and how fast releases arrive. Writes nothing and ingests nothing; a "
             "failed check is reported as unknown, never as 'no new releases'.",
    )
    parser.add_argument(
        "--fetch-history", action="store_true",
        help="download the curated historical IMGT/GENE-DB release snapshots into "
             "imgt_history_dir (for former/retired names; re-ingested on --reload --enrich)",
    )
    parser.add_argument(
        "--imgt-genedb", action="store_true",
        help="ingest ALL IMGT/GENE-DB species: bulk V (one download) + per-species "
             "leader (L-PART1+L-PART2) + RSS (V-RS) via GENElect, then incrementally "
             "patch the similarity cache (no full rebuild). Direct from imgt.org.",
    )
    parser.add_argument(
        "--species", nargs="+", default=list(ogrdb._SPECIES_KEYS),
        help="species keys (human mouse rhesus)",
    )
    parser.add_argument(
        "--imgt-names", action="store_true",
        help="enrich: crawl IMGT gene tables for alt/clone names + accessions + HGNC/NCBI-Gene xrefs",
    )
    parser.add_argument(
        "--ncbi-crawl", action="store_true",
        help="enrich: NCBI E-utilities xref (accession, NCBI Gene) per allele",
    )
    parser.add_argument(
        "--literature", action="store_true",
        help="enrich: auto-collect PubMed leads per allele via Europe PMC name "
             "match (stored as UNVERIFIED references; safe to run on a schedule)",
    )
    parser.add_argument(
        "--literature-scope", choices=("title_abstract", "all"), default="title_abstract",
        help="with --literature: where the allele name must appear. 'title_abstract' "
             "(default) is precise; 'all' searches full text too, yielding many more "
             "(noisier, asterisked) leads and a larger article set for --literature-supp",
    )
    parser.add_argument(
        "--literature-supp", action="store_true",
        help="enrich: mine PMC-OA supplementary files of discovered papers for "
             "allele names + sequences, matched back to the DB (UNVERIFIED; "
             "reuses the pmid: leads from --literature as the article set)",
    )
    parser.add_argument(
        "--literature-segments", nargs="+", metavar="SEG", default=["V"],
        choices=("V", "D", "J", "C"),
        help="with --literature: which segments to probe. Defaults to V alone, "
             "which is where the pass started; D/J/C carry allele names too and "
             "get no leads at all unless they are named here.",
    )
    parser.add_argument(
        "--literature-all-species", action="store_true",
        help="run the literature passes over every species in the database, "
             "ignoring --species. That flag defaults to human/mouse/rhesus, and in "
             "the supplement pass it gates which *sequences* a match may bind to — "
             "so without this, a paper about any other species matches and is then "
             "discarded.",
    )
    parser.add_argument(
        "--literature-pmid", nargs="+", metavar="PMID", default=None,
        help="mine these PubMed IDs' supplements directly, instead of only the "
             "papers the name probe found (implies --literature-supp). The probe "
             "discovers papers by allele-name string, so a genome/assembly paper "
             "that publishes sequences but no IMGT names is unreachable through it "
             "— naming the PMID is the only way in.",
    )
    parser.add_argument(
        "--literature-discover", action="store_true",
        help="enrich: find the articles to mine by sweeping the open-access corpus "
             "for gene-family names (IGHV, TRBV, ...) instead of querying Europe PMC "
             "once per allele, then mine each article once (implies "
             "--literature-supp). ~23 queries reach ~12,600 articles in seconds, "
             "where the per-allele probe needs ~17,800 queries to reach fewer.",
    )
    parser.add_argument(
        "--literature-oa-resolve", action="store_true",
        help="with --literature-discover / --literature-pmid: do not trust Europe "
             "PMC's OPEN_ACCESS flag as the only gate. Each gene-family query is "
             "also run without the flag and the extra articles are put to OpenAlex "
             "(no key needed, cached per PMID); one with a legal OA copy is kept, "
             "and OpenAlex's PMC location supplies a PMCID Europe PMC has not "
             "indexed. Recovered articles with no PMC copy are kept and reported "
             "as no_pmcid by the miner rather than dropped.",
    )
    parser.add_argument(
        "--literature-max-articles", type=int, default=None, metavar="N",
        help="with --literature-discover: stop after N articles. Use it to try a "
             "run end to end before committing a machine to the whole corpus.",
    )
    parser.add_argument(
        "--literature-keep-unanchored-names", action="store_true",
        help="with the supplement pass: also store bare *name* matches the article's "
             "own sequence evidence does not support. Off by default: a name is a "
             "string many species share, so a table that merely lists names matched "
             "over a thousand alleles in papers that carry no allele sequence at all. "
             "By default a name match is kept only in a (species, locus) the same "
             "article gave sequence evidence for, and an article with none keeps none.",
    )
    parser.add_argument(
        "--enrich-only", action="store_true",
        help="skip base ingest; run only the requested enrichment passes",
    )
    parser.add_argument(
        "--reload", action="store_true",
        help="consolidated full reload: TRUNCATE everything and load all sources "
             "(the all-species IMGT/GENE-DB pass incl. TCR, OGRDB, "
             "HUSA, MUSA, KI) then enrichment, then the similarity cache",
    )
    parser.add_argument(
        "--enrich", action="store_true",
        help="with --reload: DEPRECATED, enrichment is now on by default "
             "(use --no-enrich to skip it)",
    )
    parser.add_argument(
        "--no-enrich", action="store_true",
        help="with --reload: skip the enrichment passes. Drops historical IMGT names, "
             "the archived flags, IMGT clone-names and NCBI xrefs",
    )
    parser.add_argument(
        "--no-genedb", action="store_true",
        help="with --reload: skip the all-species IMGT/GENE-DB pass, the only IMGT "
             "input. The rebuild is then deliberately partial: no IMGT at all",
    )
    parser.add_argument(
        "--neighbors", action="store_true",
        help="with --reload: rebuild the precomputed similarity cache "
             "(sequence_neighbors) after ingest/enrichment",
    )
    parser.add_argument(
        "--sync-neighbors", action="store_true",
        help="after a (non-reload) ingest: incrementally patch the similarity "
             "cache for newly added alleles, two-sided — no full rebuild",
    )
    parser.add_argument(
        "--stamp-release", action="store_true",
        help="stamp the database as it stands now: content digest, distance-"
             "computation digest, counts and source versions (dataset_releases). "
             "--reload does this automatically; use this to version a database "
             "built before releases were recorded, or after an incremental ingest",
    )
    args = parser.parse_args(argv)

    if args.stamp_release:
        from .release_stamp import stamp_release

        with session_scope() as session:
            info = stamp_release(session)
        print(f"Stamped release {info['release']}")
        print(f"  content digest:    {info['content_digest']}")
        print(f"  similarity digest: {info['similarity_digest']}")
        print(f"  sequences: {info['digested']['sequences']}  "
              f"names: {info['digested']['names']}  "
              f"species: {info['species']}")
        print(f"  neighbour edges:   {info['neighbor_edges']:,}")
        return 0

    if args.check_releases:
        from .release_watch import check_releases

        with session_scope() as session:
            report = check_releases(session)
        for line in report.lines():
            print(line)
        # Exit non-zero on "could not look" or "the tracker may be dead", so a cron
        # entry surfaces it instead of mailing a clean report nobody reads.
        return 0 if report.ok else 1

    if args.fetch:
        print(f"Fetching OGRDB into {args.cache_dir} ...")
        ogrdb.fetch_to_cache(args.cache_dir, args.species)

    if args.fetch_history:
        from .imgt_history import fetch_imgt_history

        print(f"Fetching historical IMGT/GENE-DB snapshots into {settings.imgt_history_dir} ...")
        got = fetch_imgt_history(settings.imgt_history_dir)
        print(f"  {len(got)} release snapshots present")

    if args.imgt_genedb:
        from ..species import known_tags
        from .imgt_constant import ingest_constant
        from .imgt_genedb import (
            fetch_reference,
            fetch_reference_aa,
            fetch_reference_gapped,
            ingest_leaders_rss,
            ingest_reference,
            ingest_reference_aa,
        )
        from .neighbors import refresh_neighbor_labels, sync_neighbors

        print("All-species IMGT/GENE-DB ingest (V/J + C + leader + RSS) ...")
        fasta = fetch_reference(settings.imgt_genedb_fasta)
        fetch_reference_gapped(settings.imgt_genedb_fasta_gapped)  # diff & align
        with session_scope() as session:
            print(f"  V/D/J coding: {ingest_reference(session, fasta, segments=('V', 'D', 'J'))}")
        with session_scope() as session:
            print(f"  C coding:   {ingest_constant(session, fasta)}")
        with session_scope() as session:
            aa_fasta = fetch_reference_aa(settings.imgt_genedb_fasta_aa)
            print(f"  V AA:       {ingest_reference_aa(session, aa_fasta)}")
        with session_scope() as session:
            print(f"  leader/RSS: {ingest_leaders_rss(session, sorted(known_tags()), cache_dir=settings.imgt_genelect_cache_dir)}")
        with session_scope() as session:
            print(f"  similarity (incremental): {sync_neighbors(session)} edges")
        with session_scope() as session:
            print(f"  label refresh: {refresh_neighbor_labels(session)} edges")
        return 0  # standalone — don't fall through to the base OGRDB ingest

    if args.reload:
        from .reload import reload_all

        print("Consolidated full reload (truncate + load all sources) ...")
        with session_scope() as session:
            counts = reload_all(
                session, args.cache_dir,
                husa_path=settings.husa_table, musa_path=settings.musa_table,
                ki_path=settings.ki_xlsx, species_keys=tuple(args.species),
                genedb=not args.no_genedb, enrich=not args.no_enrich,
                neighbors=args.neighbors,
            )
        for k, v in counts.items():
            print(f"  {k}: {v}")
        return 0

    with session_scope() as session:
        if not args.enrich_only:
            counts = ogrdb.ingest_cache(session, args.cache_dir, args.species)
            print(f"OGRDB rows: {counts['ogrdb']} (IMGT: use --imgt-genedb)")
        if args.imgt_names:
            from .imgt_crawl import crawl_imgt_names

            print("Crawling IMGT gene tables for name iterations ...")
            print(f"IMGT names: {crawl_imgt_names(session, cache_dir=settings.imgt_crawl_cache_dir, species_keys=tuple(args.species))}")
        if args.ncbi_crawl:
            from ._ncbi_client import NCBIClient
            from .ncbi import enrich_ncbi

            print("Crawling NCBI E-utilities for cross-references ...")
            client = NCBIClient(cache_dir=settings.ncbi_cache_dir, api_key=settings.ncbi_api_key)
            print(f"NCBI: {enrich_ncbi(session, client, species_keys=args.species)}")
        # --species defaults to human/mouse/rhesus, which the base ingest needs. For
        # the literature passes that default is a trap: in the supplement pass it
        # gates which *sequences* a match may bind to, so a paper about any other
        # species is matched and then silently discarded. None = every species.
        lit_species = None if args.literature_all_species else args.species
        if args.literature:
            from .literature import EuropePMCClient, enrich_literature

            segs = tuple(dict.fromkeys(args.literature_segments))
            print(f"Collecting PubMed leads for {'/'.join(segs)} "
                  "(Europe PMC name match; unverified) ...")
            lc = EuropePMCClient(settings.literature_cache_dir, scope=args.literature_scope)
            print(f"Literature: {enrich_literature(session, lc, species_keys=lit_species, segments=segs)}")
        if args.literature_supp or args.literature_pmid or args.literature_discover:
            from .literature_supp import PMCSuppClient, enrich_literature_supp

            # Explicit PMIDs bypass the name probe entirely: enrich_literature_supp
            # already takes a `pmids` list, it just had no way in from the command
            # line, which made supplement mining reachable only for papers that
            # name an allele in their text.
            oa_client = None
            if args.literature_oa_resolve:
                from .openalex import OpenAlexClient

                # Lives under the literature cache: same opt-in pass, same disk.
                oa_client = OpenAlexClient(
                    os.path.join(settings.literature_cache_dir, "openalex")
                )
            pmids = pmid_list(args.literature_pmid)
            if not pmids and args.literature_discover:
                from .literature import EuropePMCClient, discover_articles

                print("Sweeping the open-access corpus for articles to mine ...")
                dc = EuropePMCClient(settings.literature_cache_dir)
                pmids, per_prefix, oa_stats = discover_articles(
                    session, dc, species_keys=lit_species,
                    max_articles=args.literature_max_articles,
                    oa_client=oa_client,
                )
                print(f"  {len(pmids)} distinct articles from {len(per_prefix)} "
                      f"gene-family queries: {per_prefix}")
                if oa_client is not None:
                    print(f"  OA resolution (OpenAlex): {oa_stats}")
            scope = f"{len(pmids)} named PMID(s)" if pmids else "discovered papers"
            print(f"Mining PMC-OA supplementary files ({scope}) for allele evidence "
                  "(unverified) ...")
            sc = PMCSuppClient(settings.literature_supp_cache_dir, oa_client=oa_client)
            print("Literature supp: " + str(enrich_literature_supp(
                session, sc, pmids=pmids,
                species_keys=None if pmids else lit_species,
                # A corpus sweep is hours of work; a run that dies partway must
                # not start over. Only for the sweep: naming PMIDs explicitly is
                # a deliberate "read this paper", and should read it.
                resume=args.literature_discover and not args.literature_pmid,
                anchor_names=not args.literature_keep_unanchored_names,
            )))
        if args.sync_neighbors:
            from .neighbors import sync_neighbors

            print("Syncing similarity cache (incremental) ...")
            print(f"  neighbour edges written: {sync_neighbors(session)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
