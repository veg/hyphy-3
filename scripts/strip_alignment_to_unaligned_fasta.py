#!/usr/bin/env python3
"""
strip_alignment_to_unaligned_fasta.py

Reads one or more multiple sequence alignments in NEXUS (sequential or interleaved)
or FASTA (.fa, .fna, .fasta, .nex, .nexus, optionally .gz) format, strips all gap
characters ('-', '.') while preserving in-frame coding triplets (or arbitrary
nucleotide/protein sequences with --raw), and writes clean unaligned FASTA files.

Usage examples:
  # Convert a single NEXUS or FASTA alignment to an unaligned FASTA file:
  python3 scripts/strip_alignment_to_unaligned_fasta.py tests/data/camelid.nex -o camelid_unaligned.fa

  # Batch-convert multiple NEXUS/FASTA alignments into an output directory:
  python3 scripts/strip_alignment_to_unaligned_fasta.py tests/data/*.nex -d unaligned_out/
"""

import argparse
import gzip
import os
import re
import sys
from pathlib import Path
from typing import Dict, List, Tuple


def read_text_or_gz(filepath: Path) -> str:
    with open(filepath, "rb") as f:
        magic = f.read(2)
    if magic == b"\x1f\x8b" or str(filepath).endswith(".gz"):
        with gzip.open(filepath, "rt", encoding="utf-8", errors="replace") as gz:
            return gz.read()
    with open(filepath, "r", encoding="utf-8", errors="replace") as f:
        return f.read()


def parse_fasta(content: str) -> List[Tuple[str, str]]:
    records: List[Tuple[str, str]] = []
    curr_name = None
    curr_chunks: List[str] = []

    for raw_line in content.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if curr_name is not None:
                records.append((curr_name, "".join(curr_chunks)))
            header = line[1:].strip()
            curr_name = header.split()[0] if header else f"taxon_{len(records) + 1}"
            curr_chunks = []
        elif curr_name is not None:
            curr_chunks.append(re.sub(r"\s+", "", line))

    if curr_name is not None:
        records.append((curr_name, "".join(curr_chunks)))
    return records


def strip_nexus_comments(text: str) -> str:
    """Remove [...] NEXUS comments (handling nested brackets)."""
    out = []
    depth = 0
    for ch in text:
        if ch == "[":
            depth += 1
        elif ch == "]":
            if depth > 0:
                depth -= 1
        elif depth == 0:
            out.append(ch)
    return "".join(out)


def parse_nexus(content: str) -> List[Tuple[str, str]]:
    clean = strip_nexus_comments(content)
    # Locate MATRIX ... ; block inside DATA or CHARACTERS block
    m = re.search(r"\bMATRIX\b(.*?);", clean, flags=re.IGNORECASE | re.DOTALL)
    if not m:
        raise ValueError("Could not find a MATRIX ... ; block in NEXUS file.")

    matrix_body = m.group(1)
    order: List[str] = []
    seq_map: Dict[str, List[str]] = {}

    # Tokenize line by line; supports both standard single-line/interleaved and multi-line NEXUS
    current_taxon = None
    ntax_match = re.search(r"\bNTAX\s*=\s*(\d+)", clean, flags=re.IGNORECASE)
    nchar_match = re.search(r"\bNCHAR\s*=\s*(\d+)", clean, flags=re.IGNORECASE)
    expected_nchar = int(nchar_match.group(1)) if nchar_match else None

    for raw_line in matrix_body.splitlines():
        line = raw_line.strip()
        if not line:
            continue

        # Handle quoted taxon names ('Taxon Name' or "Taxon Name") or whitespace-separated name + sequence
        if line[0] in ("'", '"'):
            quote_char = line[0]
            end_q = line.find(quote_char, 1)
            if end_q != -1:
                name = line[1:end_q].strip().replace(" ", "_")
                rest = re.sub(r"\s+", "", line[end_q + 1 :])
                if name not in seq_map:
                    order.append(name)
                    seq_map[name] = []
                if rest:
                    seq_map[name].append(rest)
                current_taxon = name
                continue

        parts = line.split()
        if len(parts) >= 2:
            # Determine whether parts[0] is a taxon name or a continuation of sequence characters
            first = parts[0]
            is_seq_continuation = (
                current_taxon is not None
                and expected_nchar is not None
                and sum(len(x) for x in seq_map[current_taxon]) < expected_nchar
                and first not in seq_map
                and re.fullmatch(r"[ACGTURYKMSWBDHVNX?\-.]+", first, flags=re.IGNORECASE) is not None
            )
            if is_seq_continuation:
                seq_map[current_taxon].append("".join(parts))
            else:
                name = first
                seq_part = "".join(parts[1:])
                if name not in seq_map:
                    order.append(name)
                    seq_map[name] = []
                seq_map[name].append(seq_part)
                current_taxon = name
        elif len(parts) == 1:
            token = parts[0]
            if (
                current_taxon is not None
                and expected_nchar is not None
                and sum(len(x) for x in seq_map[current_taxon]) < expected_nchar
            ):
                seq_map[current_taxon].append(token)
            else:
                if token not in seq_map:
                    order.append(token)
                    seq_map[token] = []
                current_taxon = token

    return [(name, "".join(seq_map[name])) for name in order]


def parse_alignment(filepath: Path) -> List[Tuple[str, str]]:
    content = read_text_or_gz(filepath)
    first_non_empty = ""
    for line in content.splitlines():
        s = line.strip()
        if s:
            first_non_empty = s
            break

    if first_non_empty.startswith(">"):
        return parse_fasta(content)
    if first_non_empty.upper().startswith("#NEXUS") or re.search(r"\bBEGIN\s+(DATA|CHARACTERS)\b", content, re.I):
        return parse_nexus(content)
    # Fallback: try FASTA first, then NEXUS
    recs = parse_fasta(content)
    if recs:
        return recs
    return parse_nexus(content)


def strip_gaps_from_sequence(seq: str, codon_aware: bool = True, strip_terminal_stop: bool = False) -> str:
    seq_up = seq.upper().replace("U", "T")
    if codon_aware and len(seq_up) >= 3 and len(seq_up) % 3 == 0:
        codons: List[str] = []
        for i in range(0, len(seq_up), 3):
            trip = seq_up[i : i + 3]
            # Drop completely gapped or missing triplets
            if trip in ("---", "...", "???", "NNN"):
                continue
            # If a triplet has a partial frame-shift gap (e.g. 'A-G'), either strip '-' or replace with 'N'
            if "-" in trip or "." in trip or "?" in trip:
                cleaned_trip = "".join(c for c in trip if c not in "-.?")
                if not cleaned_trip:
                    continue
                if len(cleaned_trip) < 3:
                    # Replace partial within-codon gap with N to preserve reading frame
                    trip = "".join("N" if c in "-.?" else c for c in trip)
            codons.append(trip)
        if strip_terminal_stop and codons and codons[-1] in ("TAA", "TAG", "TGA"):
            codons.pop()
        return "".join(codons)

    # Raw character-level gap stripping (for non-divisible-by-3 nucleotide or amino acid sequences)
    return re.sub(r"[-.?]", "", seq_up)


def wrap_fasta_seq(seq: str, width: int = 80) -> str:
    if width <= 0 or len(seq) <= width:
        return seq
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def convert_file(
    in_path: Path,
    out_path: Path,
    codon_aware: bool = True,
    strip_terminal_stop: bool = False,
    wrap_width: int = 0,
) -> Tuple[int, int, int]:
    records = parse_alignment(in_path)
    if not records:
        raise ValueError(f"No sequences found in {in_path}")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    lengths: List[int] = []

    with open(out_path, "w", encoding="utf-8") as out_f:
        for name, raw_seq in records:
            unaligned = strip_gaps_from_sequence(
                raw_seq, codon_aware=codon_aware, strip_terminal_stop=strip_terminal_stop
            )
            lengths.append(len(unaligned))
            out_f.write(f">{name}\n{wrap_fasta_seq(unaligned, wrap_width)}\n")

    return len(records), min(lengths), max(lengths)


def default_output_path(in_path: Path, out_dir: Path = None) -> Path:
    name = in_path.name
    if name.endswith(".gz"):
        name = name[:-3]
    stem = Path(name).stem
    out_name = f"{stem}_unaligned.fa"
    return (out_dir / out_name) if out_dir else in_path.with_name(out_name)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Strip gaps from NEXUS or FASTA alignments and save as unaligned FASTA."
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        type=Path,
        help="Input NEXUS (.nex, .nexus) or FASTA (.fa, .fna, .fasta) alignment file(s).",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Output unaligned FASTA file path (valid when a single input file is given).",
    )
    parser.add_argument(
        "-d",
        "--outdir",
        type=Path,
        default=None,
        help="Directory to write <stem>_unaligned.fa files when converting one or more alignments.",
    )
    parser.add_argument(
        "--raw",
        action="store_true",
        help="Strip all '-', '.', '?' characters directly without enforcing in-frame codon triplets.",
    )
    parser.add_argument(
        "--strip-stop",
        action="store_true",
        help="Strip terminal stop codon (TAA/TAG/TGA) if present at the 3' end.",
    )
    parser.add_argument(
        "-w",
        "--wrap",
        type=int,
        default=0,
        help="Wrap output FASTA lines at N characters (default: 0 = single line per sequence).",
    )

    args = parser.parse_args()

    if args.output and len(args.inputs) > 1:
        print("Error: --output (-o) can only be used with a single input file. Use --outdir (-d) for multiple files.", file=sys.stderr)
        return 1

    for in_path in args.inputs:
        out_path = args.output if args.output else default_output_path(in_path, args.outdir)
        try:
            n_taxa, min_nt, max_nt = convert_file(
                in_path=in_path,
                out_path=out_path,
                codon_aware=not args.raw,
                strip_terminal_stop=args.strip_stop,
                wrap_width=args.wrap,
            )
            if not args.raw and min_nt % 3 == 0 and max_nt % 3 == 0:
                len_desc = f"{min_nt // 3}..{max_nt // 3} codons ({min_nt}..{max_nt} nt)"
            else:
                len_desc = f"{min_nt}..{max_nt} bp/aa"
            print(f"[OK] {in_path} -> {out_path} ({n_taxa} taxa, unaligned lengths: {len_desc})")
        except Exception as exc:
            print(f"[ERROR] {in_path}: {exc}", file=sys.stderr)
            return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
