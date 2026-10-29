"""Crawl IMGT Repertoire gene tables for allele name iterations + xrefs.

The IMGT Repertoire "Genes table" (per species/group) lists, per IMGT allele:
its **Clone names** (alternative / former / literature names such as ``DP-75``,
``VH1-2``) and the reference + secondary **accession numbers**. We parse that
table and attach the clone names as ``imgt_alt`` aliases and the accessions as
``GENBANK:`` cross-references to the matching sequence (found by its IMGT allele
name already in ``names``).

GENElect query 4.1 additionally gives per-gene **HGNC** and **NCBI Gene** ids,
attached as ``xref`` annotations to every allele of that gene.

Parsing resolves ``rowspan``/``colspan`` into a dense grid, then maps columns by
their header text — robust to IMGT's merged subgroup/gene cells. Network-bound,
cached, polite; ``bs4`` lives in the ``ingest`` extra.
"""

from __future__ import annotations

import re
import time
import urllib.parse
import urllib.request
from collections.abc import Callable, Iterator
from dataclasses import dataclass, field
from pathlib import Path

from sqlalchemy import select
from sqlalchemy.orm import Session

from ..curie import build as curie_build
from ..db.models import Name, Sequence
from ..species import normalize_species
from . import loader

_REPERTOIRE = "https://www.imgt.org/IMGTrepertoire/index.php"
_GENEDB = "https://www.imgt.org/genedb/GENElect"
_UA = "allele-server/0.1 (contact: aperes7@gmail.com)"

# IMGT species URL form per our species key.
_SPECIES_URL = {"human": "Homo+sapiens", "mouse": "Mus+musculus", "rhesus": "Macaca+mulatta"}
_V_GROUPS = ("IGHV", "IGKV", "IGLV")

_ALLELE_RE = re.compile(r"^IG[HKL][VDJ][A-Za-z0-9()/-]*\*\d+$")
_ACCESSION_RE = re.compile(r"^[A-Z]{1,2}\d{5,6}(?:\.\d+)?$")

# IMGT Repertoire cells append presentational footnote markers to the allele
# name, e.g. ``IGLV1-47*01(42)``, keyed to notes printed under the table. They
# say nothing about the allele, so we drop them before matching. The marker is
# only stripped when a well-formed allele name is left behind, so clan names,
# which carry parentheses of their own (``IGHV(II)-15-1*01``), survive intact.
_FOOTNOTED_ALLELE_RE = re.compile(r"^(IG[HKL][VDJ][A-Za-z0-9()/-]*\*\d+)(?:\(\d+\))+$")


def _allele_cell(cell: str) -> str:
    """Normalise a gene-table cell to a bare IMGT allele name.

    Strips IMGT's layout whitespace and any trailing footnote marker; the
    result is what ``_ALLELE_RE`` is written against.
    """
    text = cell.strip().replace(" ", "")
    m = _FOOTNOTED_ALLELE_RE.match(text)
    return m.group(1) if m else text


@dataclass
class GeneTableRecord:
    allele: str
    clone_names: list[str] = field(default_factory=list)
    accessions: list[str] = field(default_factory=list)


# --- HTTP (cached, polite) -------------------------------------------------

def _fetch(url: str, cache_dir: Path, fetcher: Callable[[str], str] | None) -> str:
    import hashlib

    cache_dir.mkdir(parents=True, exist_ok=True)
    path = cache_dir / f"{hashlib.sha256(url.encode()).hexdigest()}.html"
    if path.exists():
        return path.read_text(encoding="utf-8")
    if fetcher is not None:
        body = fetcher(url)
    else:
        time.sleep(0.5)  # polite
        req = urllib.request.Request(url, headers={"User-Agent": _UA})
        with urllib.request.urlopen(req, timeout=120) as resp:
            body = resp.read().decode("utf-8", "replace")
    path.write_text(body, encoding="utf-8")
    return body


# --- HTML table grid -------------------------------------------------------

def _grid(table) -> list[list[str]]:
    """Resolve a bs4 <table> into a dense text grid, expanding row/colspans."""
    rows = table.find_all("tr")
    grid: list[list[str | None]] = []
    pending: dict[tuple[int, int], str] = {}  # (row, col) -> text from a rowspan
    for r, tr in enumerate(rows):
        if len(grid) <= r:
            grid.append([])
        row = grid[r]
        c = 0
        for cell in tr.find_all(["td", "th"], recursive=False):
            while (r, c) in pending:
                _set(grid, r, c, pending.pop((r, c)))
                c += 1
            text = " ".join(cell.get_text(" ", strip=True).split())
            try:
                rowspan = int(cell.get("rowspan", 1))
            except ValueError:
                rowspan = 1
            try:
                colspan = int(cell.get("colspan", 1))
            except ValueError:
                colspan = 1
            for dc in range(colspan):
                _set(grid, r, c, text)
                for dr in range(1, rowspan):
                    pending[(r + dr, c)] = text
                c += 1
        # flush any pending cells that begin at columns beyond this row's cells
    # normalise to strings
    width = max((len(row) for row in grid), default=0)
    return [[(row[i] if i < len(row) and row[i] is not None else "") for i in range(width)]
            for row in grid]


def _set(grid: list[list], r: int, c: int, val: str) -> None:
    while len(grid) <= r:
        grid.append([])
    row = grid[r]
    while len(row) <= c:
        row.append(None)
    row[c] = val


# --- parsing ---------------------------------------------------------------

def parse_gene_table(html: str) -> list[GeneTableRecord]:
    """Parse an IMGT Repertoire gene-table page into per-allele records."""
    from bs4 import BeautifulSoup

    soup = BeautifulSoup(html, "html.parser")
    best: list[GeneTableRecord] = []
    for table in soup.find_all("table"):
        header_text = table.get_text(" ", strip=True).lower()
        if "imgt allele name" not in header_text or "clone names" not in header_text:
            continue
        grid = _grid(table)
        if not grid:
            continue
        # Data starts at the first row carrying an allele name; everything above
        # is header. IMGT uses a two-row header (group header + sub-header), so
        # we combine all header rows per column into one label string.
        data_start = next(
            (i for i, row in enumerate(grid)
             if any(_ALLELE_RE.match(_allele_cell(c)) for c in row)),
            None,
        )
        if not data_start:  # None or 0 (need at least one header row)
            continue
        width = len(grid[data_start])
        header = [
            " ".join(
                grid[r][c] for r in range(data_start) if c < len(grid[r])
            ).lower()
            for c in range(width)
        ]
        allele_cols = [i for i, h in enumerate(header) if "imgt allele name" in h]
        clone_cols = [i for i, h in enumerate(header) if "clone name" in h]
        acc_cols = [i for i, h in enumerate(header) if "accession" in h]
        if not allele_cols:
            continue
        records = list(
            _rows_to_records(grid[data_start:], allele_cols[0], clone_cols, acc_cols)
        )
        # The page carries several matching tables; keep the richest data table.
        if len(records) > len(best):
            best = records
    return best


def _split_tokens(cell: str) -> list[str]:
    return [t for t in re.split(r"[\s,;]+", cell.strip()) if t and t not in ("-", "and")]


def _looks_like_clone(token: str) -> bool:
    """A plausible clone / former name (DP-75, VH1-2, hv1263) — not noise/accession."""
    t = token.strip()
    if len(t) < 3 or len(t) > 18:
        return False
    if not re.search(r"[A-Za-z]", t):
        return False
    if _ACCESSION_RE.match(t):
        return False
    if _ALLELE_RE.match(_allele_cell(t)):  # an allele name, footnoted or not
        return False
    if "_" in t or "." in t:  # contig / genome-assembly ids
        return False
    if re.match(r"^[A-Z]{2,}\d{4,}", t):  # accession-shaped
        return False
    return True


def _rows_to_records(rows, allele_col, clone_cols, acc_cols) -> Iterator[GeneTableRecord]:
    for row in rows:
        if allele_col >= len(row):
            continue
        allele = _allele_cell(row[allele_col])
        if not _ALLELE_RE.match(allele):
            continue
        clones: list[str] = []
        for ci in clone_cols:
            if ci < len(row):
                clones += [t for t in _split_tokens(row[ci]) if _looks_like_clone(t)]
        accs: list[str] = []
        for ci in acc_cols:
            if ci < len(row):
                accs += [t for t in _split_tokens(row[ci]) if _ACCESSION_RE.match(t)]
        yield GeneTableRecord(
            allele=allele,
            clone_names=list(dict.fromkeys(clones)),
            accessions=list(dict.fromkeys(accs)),
        )


def parse_genedb_xref(html: str) -> dict[str, dict[str, str]]:
    """Parse GENElect 4.1 cross-ref table: gene -> {hgnc, ncbi_gene}."""
    from bs4 import BeautifulSoup

    soup = BeautifulSoup(html, "html.parser")
    out: dict[str, dict[str, str]] = {}
    for table in soup.find_all("table"):
        if "ncbi gene" not in table.get_text(" ", strip=True).lower():
            continue
        for tr in table.find_all("tr"):
            text = tr.get_text(" ", strip=True)
            gene = next((t for t in re.split(r"\s+", text) if re.match(r"^IG[HKL][VDJ]", t)), None)
            if not gene:
                continue
            xref: dict[str, str] = {}
            for a in tr.find_all("a", href=True):
                m = re.search(r"hgnc_id/HGNC:(\d+)", a["href"])
                if m:
                    xref["hgnc"] = m.group(1)
                m = re.search(r"/gene/(\d+)", a["href"])
                if m:
                    xref["ncbi_gene"] = m.group(1)
            if xref:
                out.setdefault(gene.split("*")[0], {}).update(xref)
    return out


# --- attach to DB ----------------------------------------------------------

def _gene_of(allele: str) -> str:
    return allele.split("*")[0]


def crawl_imgt_names(
    session: Session,
    *,
    cache_dir: str | Path,
    species_keys: tuple[str, ...] = ("human", "mouse", "rhesus"),
    groups: tuple[str, ...] = _V_GROUPS,
    fetcher: Callable[[str], str] | None = None,
) -> dict[str, int]:
    """Crawl gene tables + GENElect xrefs and attach aliases/annotations.

    Returns counts: ``alt_names``, ``accession_xrefs``, ``gene_xrefs``, ``alleles``.
    """
    cache_dir = Path(cache_dir)
    source = loader.get_or_create_source(
        session, name="IMGT", version=None, curie_prefix="IMGT/GENE-DB",
    )
    stats = {"alt_names": 0, "accession_xrefs": 0, "gene_xrefs": 0, "alleles": 0}

    # Index our IMGT-named sequences by (species_tag, allele-name) for matching.
    for sp in species_keys:
        tag = normalize_species(sp)
        sp_url = _SPECIES_URL[sp]
        for group in groups:
            # Gene table: clone names + accessions per allele.
            gt_url = (
                f"{_REPERTOIRE}?section=LocusGenes&repertoire=genetable"
                f"&species={sp_url}&group={group}"
            )
            try:
                records = parse_gene_table(_fetch(gt_url, cache_dir, fetcher))
            except Exception:  # noqa: BLE001 - degrade; a group may be absent
                records = []
            # GENElect 4.1 cross-refs (HGNC / NCBI Gene) per gene.
            xref_url = (
                f"{_GENEDB}?query={urllib.parse.quote('4.1 ' + group)}"
                f"&species={urllib.parse.quote(sp_url.replace('+', ' '))}"
            )
            try:
                gene_xrefs = parse_genedb_xref(_fetch(xref_url, cache_dir, fetcher))
            except Exception:  # noqa: BLE001
                gene_xrefs = {}

            for rec in records:
                seqs = _sequences_for_name(session, tag, rec.allele)
                if not seqs:
                    continue
                stats["alleles"] += 1
                for seq in seqs:
                    for clone in rec.clone_names:
                        loader.add_name(session, seq, source, clone, "imgt_alt")
                        stats["alt_names"] += 1
                    for acc in rec.accessions:
                        loader.add_annotation(
                            session, seq, source, domain="xref",
                            key=f"imgt_accession:{acc}", value=curie_build("GENBANK", acc),
                        )
                        stats["accession_xrefs"] += 1
                    gx = gene_xrefs.get(_gene_of(rec.allele), {})
                    if gx.get("hgnc"):
                        loader.add_annotation(
                            session, seq, source, domain="xref", key="hgnc",
                            value=curie_build("HGNC", gx["hgnc"]),
                        )
                        stats["gene_xrefs"] += 1
                    if gx.get("ncbi_gene"):
                        loader.add_annotation(
                            session, seq, source, domain="xref", key="ncbi_gene",
                            value=curie_build("NCBIGENE", gx["ncbi_gene"]),
                        )
                        stats["gene_xrefs"] += 1
    return stats


def _sequences_for_name(session: Session, species_tag: str, allele: str) -> list[Sequence]:
    return list(
        session.scalars(
            select(Sequence)
            .join(Name, Name.sequence_id == Sequence.id)
            .where(Sequence.species == species_tag, Name.name == allele)
            .distinct()
        )
    )
