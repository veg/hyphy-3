#!/usr/bin/env python3
"""
Example 2: Mixed Effects Model of Evolution (MEME) Analysis in Python
Tests for site-by-site episodic diversifying selection using empirical Bayes mixtures.
"""

import sys
from pathlib import Path

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

    print(f"Running MEME analysis on {aln.num_codons} codons...")
    meme = hp.MEMEAnalyzer.create_and_fit(tree, aln, pvalue_threshold=0.10)
    results = meme.run()

    print(f"Global MG94 Baseline ln L: {meme.global_log_l:.2f}")

    # Identify sites undergoing episodic selection
    sig_sites = [r for r in results if r.p_value < 0.10 and r.beta_plus > 1.0]

    print(f"\n[MEME Results]")
    print(f"  Total sites analyzed: {len(results)}")
    print(f"  Sites with evidence of episodic positive selection (p < 0.10, beta+ > 1.0): {len(sig_sites)}")

    for s in sig_sites:
        print(f"    Codon {s.site + 1:3d}: alpha={s.alpha:5.2f}, beta-={s.beta_minus:5.2f} (p-={s.p_minus:4.2f}), "
              f"beta+={s.beta_plus:6.2f} (p+={s.p_plus:4.2f}), LRT={s.lrt:5.2f}, p={s.p_value:.4f}")

    # Export Datamonkey JSON
    out_file = "cd2.MEME.json"
    with open(out_file, "w") as f:
        f.write(meme.to_json())
    print(f"\nSaved Datamonkey JSON to {out_file}")

if __name__ == "__main__":
    main()
