#!/usr/bin/env python3
"""
fetch_ncbi_orthologs.py

Downloads homologous coding sequences (CDS) from NCBI Datasets Orthologs by
NCBI Gene ID, Gene Symbol, or NCBI Datasets URL (e.g.,
https://www.ncbi.nlm.nih.gov/datasets/gene/7157), curates 1 canonical in-frame
ortholog per species (matching the MANE Select / RefSeq Select canonical length
to avoid alternative splice N-terminal extensions or truncated isoforms),
formats clean `>Genus_species` FASTA headers, and optionally launches
`hyphy3 characterize` (the HyPhy-3 One-Button analysis).
"""

import argparse
import json
import os
import re
import shutil
import statistics
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple


STOP_CODONS = {"TAA", "TAG", "TGA"}


def find_executable(name: str, extra_paths: List[str]) -> Optional[str]:
    found = shutil.which(name)
    if found:
        return found
    for p in extra_paths:
        cand = os.path.expanduser(p)
        if os.path.isfile(cand) and os.access(cand, os.X_OK):
            return cand
    return None


def parse_gene_query(query: str) -> Tuple[str, str]:
    """
    Returns ('gene-id', id_str) or ('symbol', symbol_str).
    Accepts:
      - '7157'
      - 'https://www.ncbi.nlm.nih.gov/datasets/gene/7157'
      - 'https://www.ncbi.nlm.nih.gov/gene/7157'
      - 'TP53'
    """
    q = query.strip()
    m_url = re.search(r"ncbi\.nlm\.nih\.gov/(?:datasets/)?gene/(\d+)", q)
    if m_url:
        return ("gene-id", m_url.group(1))
    if q.isdigit():
        return ("gene-id", q)
    return ("symbol", q)


def sanitize_species_name(name: str) -> str:
    s = re.sub(r"[^A-Za-z0-9]+", "_", name.strip()).strip("_")
    return s if s else "Unknown_species"


def is_valid_inframe_cds(seq: str, require_start_atg: bool = True) -> Tuple[bool, str]:
    """
    Checks that a nucleotide sequence is a multiple of 3, optionally starts with ATG,
    and has no internal premature stop codons (excluding the terminal codon).
    """
    s = seq.upper().replace("U", "T")
    if len(s) < 60:
        return (False, "too_short")
    if len(s) % 3 != 0:
        return (False, "not_multiple_of_3")
    if require_start_atg and not s.startswith("ATG"):
        return (False, "no_start_ATG")
    # Check internal codons (all except the very last codon) for stop codons or high N content
    n_count = s.count("N")
    if n_count > 0.05 * len(s):
        return (False, "excessive_Ns")
    for i in range(0, len(s) - 3, 3):
        codon = s[i : i + 3]
        if codon in STOP_CODONS:
            return (False, f"internal_stop_at_nt_{i+1}")
    return (True, "ok")


def kmer_set(seq: str, k: int = 12) -> set:
    return {seq[i : i + k] for i in range(0, max(0, len(seq) - k + 1), 3)}


def jaccard_similarity(k1: set, k2: set) -> float:
    if not k1 or not k2:
        return 0.0
    inter = len(k1 & k2)
    union = len(k1 | k2)
    return inter / union if union > 0 else 0.0


def parse_ncbi_zip(
    zip_path: Path,
    selection_mode: str = "canonical",
    include_accession_in_header: bool = False,
    require_start_atg: bool = True,
    max_length_dev_frac: float = 0.25,
) -> Tuple[List[Tuple[str, str, Dict]], Dict]:
    """
    Extracts and curates 1 CDS per species from an NCBI Datasets gene zip package.
    Returns (selected_records, summary_stats).
    """
    select_accessions: Dict[str, Tuple[str, str]] = {}  # species -> (acc, selectCategory)
    gene_symbol = "GENE"

    with zipfile.ZipFile(zip_path, "r") as zf:
        names = set(zf.namelist())
        if "ncbi_dataset/data/cds.fna" not in names:
            raise RuntimeError(
                f"Archive {zip_path} does not contain ncbi_dataset/data/cds.fna. "
                "Make sure `--include cds` was passed to `datasets download`."
            )

        # 1. Read product_report.jsonl if available to identify MANE_SELECT / REFSEQ_SELECT
        if "ncbi_dataset/data/product_report.jsonl" in names:
            with zf.open("ncbi_dataset/data/product_report.jsonl") as pf:
                for raw_line in pf:
                    line = raw_line.decode("utf-8", errors="replace").strip()
                    if not line:
                        continue
                    obj = json.loads(line)
                    sp = obj.get("taxname", "")
                    sym = obj.get("symbol", "")
                    if sym and gene_symbol == "GENE":
                        gene_symbol = sym
                    for tr in obj.get("transcripts", []):
                        cat = tr.get("selectCategory", "")
                        acc = tr.get("accessionVersion", "")
                        if cat in ("MANE_SELECT", "REFSEQ_SELECT") and acc and sp:
                            # Prefer MANE_SELECT over REFSEQ_SELECT
                            if sp not in select_accessions or cat == "MANE_SELECT":
                                select_accessions[sp] = (acc, cat)

        # 2. Parse all FASTA records from ncbi_dataset/data/cds.fna
        raw_records: List[Tuple[str, str]] = []
        with zf.open("ncbi_dataset/data/cds.fna") as cf:
            hdr: Optional[str] = None
            buf: List[str] = []
            for raw_line in cf:
                line = raw_line.decode("utf-8", errors="replace").strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if hdr is not None:
                        raw_records.append((hdr, "".join(buf)))
                    hdr = line
                    buf = []
                else:
                    buf.append(line.upper().replace("U", "T"))
            if hdr is not None:
                raw_records.append((hdr, "".join(buf)))

    # 3. Group valid in-frame CDS isoforms by species
    by_species: Dict[str, List[Dict]] = {}
    rejected_counts: Dict[str, int] = {}

    for hdr, seq in raw_records:
        if "low quality" in hdr.lower() or "partial" in hdr.lower():
            rejected_counts["low_quality_or_partial"] = rejected_counts.get("low_quality_or_partial", 0) + 1
            continue

        ok, reason = is_valid_inframe_cds(seq, require_start_atg=require_start_atg)
        if not ok:
            rejected_counts[reason] = rejected_counts.get(reason, 0) + 1
            continue

        first_tok = hdr[1:].split()[0]
        acc = first_tok.split(":")[0]
        m_org = re.search(r"\[organism=([^\]]+)\]", hdr)
        sp_raw = m_org.group(1).strip() if m_org else acc
        m_tr = re.search(r"\[transcript=([^\]]+)\]", hdr)
        tr_var = m_tr.group(1).strip() if m_tr else ""
        m_gid = re.search(r"\[GeneID=(\d+)\]", hdr)
        gid = m_gid.group(1) if m_gid else ""

        is_mane_or_select = (
            sp_raw in select_accessions and select_accessions[sp_raw][0] == acc
        )
        select_cat = select_accessions[sp_raw][1] if is_mane_or_select else ""

        by_species.setdefault(sp_raw, []).append(
            {
                "species": sp_raw,
                "accession": acc,
                "gene_id": gid,
                "transcript": tr_var,
                "select_category": select_cat,
                "is_curated_nm": acc.startswith("NM_"),
                "length_nt": len(seq),
                "seq": seq,
                "header": hdr,
            }
        )

    if not by_species:
        raise RuntimeError("No valid in-frame coding sequences passed quality filters.")

    # 4. Determine canonical reference sequence & target length
    #    Priority:
    #      (a) MANE_SELECT / REFSEQ_SELECT isoform (e.g. Homo sapiens NM_000546.6 = 1182 nt)
    #      (b) Mode/median of single-isoform species + curated NM_ sequences
    canonical_seq: Optional[str] = None
    canonical_len: int = 0
    canonical_source: str = ""

    for sp, items in by_species.items():
        for it in items:
            if it["select_category"] == "MANE_SELECT":
                canonical_seq = it["seq"]
                canonical_len = it["length_nt"]
                canonical_source = f"{sp} ({it['accession']}, MANE_SELECT, {canonical_len} nt)"
                break
        if canonical_seq:
            break

    if not canonical_seq:
        for sp, items in by_species.items():
            for it in items:
                if it["select_category"] == "REFSEQ_SELECT":
                    canonical_seq = it["seq"]
                    canonical_len = it["length_nt"]
                    canonical_source = f"{sp} ({it['accession']}, REFSEQ_SELECT, {canonical_len} nt)"
                    break
            if canonical_seq:
                break

    if not canonical_seq:
        # Fallback: median length of primary isoforms across species
        rep_lengths = [ max(items, key=lambda x: (x["is_curated_nm"], x["length_nt"]))["length_nt"] for items in by_species.values() ]
        canonical_len = int(statistics.median(rep_lengths))
        # Pick the sequence closest to median length, preferring NM_
        best_item = None
        for items in by_species.values():
            for it in items:
                key = (-abs(it["length_nt"] - canonical_len), it["is_curated_nm"])
                if best_item is None or key > (-abs(best_item["length_nt"] - canonical_len), best_item["is_curated_nm"]):
                    best_item = it
        if best_item:
            canonical_seq = best_item["seq"]
            canonical_len = best_item["length_nt"]
            canonical_source = f"{best_item['species']} ({best_item['accession']}, median-anchor, {canonical_len} nt)"

    ref_kmers = kmer_set(canonical_seq) if canonical_seq else set()

    # 5. Select 1 isoform per species
    selected: List[Tuple[str, str, Dict]] = []
    dropped_outliers = 0

    for sp in sorted(by_species.keys()):
        items = by_species[sp]
        if selection_mode == "longest":
            chosen = max(
                items,
                key=lambda x: (
                    x["select_category"] != "",
                    x["length_nt"],
                    x["is_curated_nm"],
                ),
            )
        else:
            # 'canonical' mode:
            # 1. Exact MANE_SELECT / REFSEQ_SELECT if annotated for this species
            # 2. Minimal length distance to canonical_len, breaking ties by:
            #    - Curated NM_ over predicted XM_
            #    - Higher in-frame 12-mer Jaccard similarity to the canonical reference
            #    - Transcript X1 / variant 1 over higher-numbered variants
            def score_isoform(it: Dict) -> Tuple:
                is_select = 1 if it["select_category"] in ("MANE_SELECT", "REFSEQ_SELECT") else 0
                len_diff = abs(it["length_nt"] - canonical_len)
                # Penalize length differences in codon steps; within identical length, rank by k-mer similarity & NM_
                k_sim = jaccard_similarity(kmer_set(it["seq"]), ref_kmers) if ref_kmers else 0.0
                is_x1 = 1 if it["transcript"] in ("", "1", "X1", "variant 1") else 0
                return (
                    is_select,
                    -len_diff,
                    round(k_sim, 4),
                    1 if it["is_curated_nm"] else 0,
                    is_x1,
                )

            chosen = max(items, key=score_isoform)

        # Filter out extreme fragmentary or chimeric length outliers (> max_length_dev_frac from canonical_len)
        if canonical_len > 0 and max_length_dev_frac > 0:
            rel_dev = abs(chosen["length_nt"] - canonical_len) / float(canonical_len)
            if rel_dev > max_length_dev_frac:
                dropped_outliers += 1
                continue

        sp_tag = sanitize_species_name(sp)
        if include_accession_in_header:
            acc_clean = re.sub(r"[^A-Za-z0-9]+", "_", chosen["accession"])
            fasta_id = f"{sp_tag}_{acc_clean}"
        else:
            fasta_id = sp_tag

        selected.append((fasta_id, chosen["seq"], chosen))

    stats = {
        "gene_symbol": gene_symbol,
        "total_raw_cds": len(raw_records),
        "total_species_found": len(by_species),
        "selected_species_count": len(selected),
        "dropped_length_outliers": dropped_outliers,
        "canonical_reference": canonical_source,
        "canonical_length_nt": canonical_len,
        "rejected_counts": rejected_counts,
    }
    return (selected, stats)


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Fetch NCBI Datasets gene orthologs (by Gene ID, Symbol, or NCBI URL), "
            "curate 1 canonical in-frame CDS per species, and optionally run "
            "`hyphy3 characterize`."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Examples:
  # 1. Fetch TP53 (Gene ID 7157 or URL) across Primates and run hyphy3 characterize:
  %(prog)s https://www.ncbi.nlm.nih.gov/datasets/gene/7157 --ortholog primates --run

  # 2. Fetch BRCA1 across Mammals (curated 1-per-species FASTA only):
  %(prog)s BRCA1 --ortholog mammals -o brca1_mammals.fna

  # 3. Fetch IFITM3 (10410) across Chiroptera (bats) and run with 8 threads:
  %(prog)s 10410 --ortholog chiroptera --run --threads 8
""",
    )
    parser.add_argument(
        "gene",
        help="NCBI Gene ID (e.g. 7157), Gene Symbol (e.g. TP53), or NCBI Datasets Gene URL",
    )
    parser.add_argument(
        "--ortholog",
        "-g",
        default="primates",
        help="Taxonomic group for orthologs (e.g. 'primates', 'mammals', 'carnivora', 'aves', 'rodentia', 'all'; default: primates)",
    )
    parser.add_argument(
        "--taxon",
        default="human",
        help="Anchor species when querying by gene symbol (default: human)",
    )
    parser.add_argument(
        "--output",
        "-o",
        default="",
        help="Output curated unaligned CDS FASTA path (default: <symbol>_<ortholog>_cds.fna)",
    )
    parser.add_argument(
        "--isoform-mode",
        choices=["canonical", "longest"],
        default="canonical",
        help="Isoform selection strategy per species: 'canonical' matches MANE/RefSeq Select length & homology (default); 'longest' picks the longest ORF",
    )
    parser.add_argument(
        "--max-length-deviation",
        type=float,
        default=0.25,
        help="Maximum fractional CDS length deviation from the canonical ortholog before filtering a sequence as truncated/chimeric (default: 0.25; 0 to disable)",
    )
    parser.add_argument(
        "--include-accession",
        action="store_true",
        help="Append RefSeq accession to species name in FASTA header (e.g. >Homo_sapiens_NM_000546_6)",
    )
    parser.add_argument(
        "--keep-zip",
        default="",
        help="Optional path to save the raw NCBI Datasets .zip package",
    )
    parser.add_argument(
        "--run",
        action="store_true",
        help="Immediately run `hyphy3 characterize --unaligned <output.fna>` after curating the ortholog FASTA",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=8,
        help="Number of OpenMP threads for `hyphy3 characterize` (default: 8)",
    )
    parser.add_argument(
        "--hyphy3-bin",
        default="",
        help="Explicit path to `hyphy3` or `hyphy3-characterize` binary",
    )
    parser.add_argument(
        "--datasets-bin",
        default="",
        help="Explicit path to NCBI `datasets` CLI binary",
    )

    args, hyphy_args = parser.parse_known_args()

    # 1. Locate NCBI `datasets` CLI
    datasets_bin = args.datasets_bin or find_executable(
        "datasets",
        [
            "~/miniconda3/bin/datasets",
            "~/anaconda3/bin/datasets",
            "/opt/homebrew/bin/datasets",
            "/usr/local/bin/datasets",
        ],
    )
    if not datasets_bin:
        print(
            "Error: Could not find NCBI `datasets` CLI. Install via `conda install -c conda-forge ncbi-datasets-cli` or pass --datasets-bin.",
            file=sys.stderr,
        )
        return 1

    mode, query_val = parse_gene_query(args.gene)
    ortholog_taxa = [t.strip() for t in args.ortholog.split(",") if t.strip()]

    print(f"── NCBI Datasets Ortholog Fetcher ──────────────────────────────────────")
    print(f"  Query           : {query_val} ({mode})")
    print(f"  Ortholog Clade  : {', '.join(ortholog_taxa)}")
    print(f"  Isoform Filter  : {args.isoform_mode} (max length dev: {args.max_length_deviation:.0%})")

    with tempfile.TemporaryDirectory(prefix="hyphy3_ncbi_") as tmpdir:
        zip_path = Path(args.keep_zip) if args.keep_zip else (Path(tmpdir) / "ncbi_dataset.zip")
        cmd = [datasets_bin, "download", "gene", mode, query_val]
        if mode == "symbol":
            cmd.extend(["--taxon", args.taxon])
        for ot in ortholog_taxa:
            cmd.extend(["--ortholog", ot])
        cmd.extend(
            [
                "--include",
                "cds,product-report",
                "--no-progressbar",
                "--filename",
                str(zip_path),
            ]
        )

        print(f"  Downloading     : {' '.join(cmd)}")
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode != 0 or not zip_path.is_file():
            print(f"Error running NCBI datasets CLI:\n{res.stderr}", file=sys.stderr)
            return 1

        selected, stats = parse_ncbi_zip(
            zip_path,
            selection_mode=args.isoform_mode,
            include_accession_in_header=args.include_accession,
            require_start_atg=True,
            max_length_dev_frac=args.max_length_deviation,
        )

    gene_sym = stats["gene_symbol"].lower()
    clade_slug = re.sub(r"[^A-Za-z0-9]+", "_", args.ortholog.lower()).strip("_")
    out_fasta = Path(args.output) if args.output else Path(f"{gene_sym}_{clade_slug}_cds.fna")
    out_meta = out_fasta.with_suffix(".manifest.json")

    with open(out_fasta, "w", encoding="utf-8") as f:
        for fasta_id, seq, _ in selected:
            f.write(f">{fasta_id}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i : i + 80] + "\n")

    manifest = {
        "query": args.gene,
        "ortholog_filter": ortholog_taxa,
        "isoform_mode": args.isoform_mode,
        "stats": stats,
        "selected_orthologs": [
            {
                "taxon_id": fid,
                "species": meta["species"],
                "accession": meta["accession"],
                "gene_id": meta["gene_id"],
                "transcript": meta["transcript"],
                "select_category": meta["select_category"],
                "length_nt": meta["length_nt"],
                "length_codons": meta["length_nt"] // 3,
            }
            for fid, _, meta in selected
        ],
    }
    with open(out_meta, "w", encoding="utf-8") as mf:
        json.dump(manifest, mf, indent=2)

    lengths_codons = [m["length_nt"] // 3 for _, _, m in selected]
    min_c = min(lengths_codons) if lengths_codons else 0
    max_c = max(lengths_codons) if lengths_codons else 0
    med_c = int(statistics.median(lengths_codons)) if lengths_codons else 0

    print(f"  Canonical Anchor: {stats['canonical_reference']}")
    print(
        f"  Curated Output  : {stats['selected_species_count']} species "
        f"(from {stats['total_raw_cds']} raw CDS isoforms; "
        f"{stats['dropped_length_outliers']} length outliers dropped)"
    )
    print(f"  Codon Lengths   : median {med_c} codons (range: {min_c}..{max_c})")
    print(f"  FASTA Written   : {out_fasta}")
    print(f"  Manifest JSON   : {out_meta}")
    print(f"────────────────────────────────────────────────────────────────────────")

    if args.run:
        script_dir = Path(__file__).resolve().parent
        repo_root = script_dir.parent
        hyphy3_bin = args.hyphy3_bin or find_executable(
            "hyphy3-characterize",
            [
                str(repo_root / "build" / "hyphy3-characterize"),
                str(repo_root / "build" / "hyphy3"),
                "hyphy3",
            ],
        )
        if not hyphy3_bin:
            print(
                "Error: Could not find `hyphy3` or `hyphy3-characterize` binary. Build it first or pass --hyphy3-bin.",
                file=sys.stderr,
            )
            return 1

        run_cmd = [hyphy3_bin]
        if Path(hyphy3_bin).name == "hyphy3":
            run_cmd.append("characterize")
        run_cmd.extend(["--unaligned", str(out_fasta), "--threads", str(args.threads)])

        extra = [a for a in hyphy_args if a != "--"]
        run_cmd.extend(extra)

        print(f"\nLaunching HyPhy-3 One-Button Characterization:\n  {' '.join(run_cmd)}\n")
        return subprocess.call(run_cmd)

    return 0


if __name__ == "__main__":
    sys.exit(main())
