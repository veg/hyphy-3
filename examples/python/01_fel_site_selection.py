#!/usr/bin/env python3
"""
Example 1: Fixed Effects Likelihood (FEL) Analysis in Python
Tests for site-by-site pervasive diversifying and purifying selection.
"""

import sys
from pathlib import Path

# Add python bindings directory to sys.path if running from repository
repo_root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(repo_root / "python"))

import hyphy3 as hp

def main():
    data_dir = repo_root / "benchmarks" / "data"
    aln_path = str(data_dir / "cd2.fna")
    tree_path = str(data_dir / "cd2.nwk")

    print(f"Loading alignment: {aln_path}")
    print(f"Loading tree:      {tree_path}")

    aln = hp.Alignment.load(aln_path)
    tree = hp.Tree.from_newick_file(tree_path)

    print(f"Dataset: {aln.num_taxa} sequences, {aln.num_codons} codons, {len(aln.patterns)} unique patterns")

    # Fit baseline GTR + MG94 and run FEL site-by-site LRTs
    print("\nRunning FEL analysis...")
    analyzer = hp.FELAnalyzer.create_and_fit(tree, aln, pvalue_threshold=0.10)
    results = analyzer.run()

    print(f"Baseline GTR ln L : {analyzer.gtr_log_l:.2f}")
    print(f"Baseline MG94 ln L: {analyzer.mg94_log_l:.2f} (omega = {analyzer.mg94_omega:.4f})")
    print(f"Analyzed {len(results)} codon sites.")

    # Filter significant sites
    pos_sites = [r for r in results if r.p_value < 0.10 and r.beta > r.alpha]
    neg_sites = [r for r in results if r.p_value < 0.10 and r.beta < r.alpha]

    print(f"\n[Summary of Selection]")
    print(f"  Pervasive Diversifying Sites (dN > dS, p < 0.10): {len(pos_sites)}")
    for s in pos_sites:
        omega_site = (s.beta / s.alpha) if s.alpha > 0 else float("inf")
        print(f"    Site {s.site + 1:3d}: alpha={s.alpha:6.3f}, beta={s.beta:6.3f} (dN/dS={omega_site:6.2f}), LRT={s.lrt:6.2f}, p={s.p_value:.4f}")

    print(f"  Pervasive Purifying Sites (dN < dS, p < 0.10): {len(neg_sites)}")
    print(f"    First 5 sites: {[s.site + 1 for s in neg_sites[:5]]}...")

    # Export Datamonkey-compatible JSON
    json_str = analyzer.to_json()
    out_file = "cd2.FEL.json"
    with open(out_file, "w") as f:
        f.write(json_str)
    print(f"\nSaved Datamonkey JSON to {out_file}")

if __name__ == "__main__":
    main()
