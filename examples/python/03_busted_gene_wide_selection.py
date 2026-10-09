#!/usr/bin/env python3
"""
Example 3: Branch-site Unrestricted Statistical Test (BUSTED) & BUSTED-S (SRV) in Python
Tests for gene-wide episodic diversifying selection and synonymous rate variation.
"""

import sys
from pathlib import Path

repo_root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(repo_root / "python"))

import hyphy3 as hp

def main():
    data_dir = repo_root / "benchmarks" / "data"
    aln_path = str(data_dir / "adh.fna")
    tree_path = str(data_dir / "adh.nwk")

    print(f"Loading alignment: {aln_path}")
    print(f"Loading tree:      {tree_path}")

    aln = hp.Alignment.load(aln_path)
    tree = hp.Tree.from_newick_file(tree_path)

    # 1. Standard BUSTED (K=3)
    print("\n--- 1. Standard BUSTED (K=3) ---")
    analyzer = hp.BUSTEDAnalyzer.create_and_fit(tree, aln)
    res_standard = analyzer.run()

    print(f"Unconstrained ln L : {res_standard.unconstrained.log_likelihood:.2f}")
    print(f"Constrained ln L   : {res_standard.constrained.log_likelihood:.2f}")
    print(f"LRT                : {res_standard.lrt:.4f}")
    print(f"p-value            : {res_standard.p_value:.6e}")
    print(f"Alt omegas         : {[round(w, 4) for w in res_standard.unconstrained.test_distribution.omegas]}")
    print(f"Alt weights        : {[round(p, 4) for p in res_standard.unconstrained.test_distribution.weights]}")

    # 2. BUSTED-S with Synonymous Rate Variation (SRV)
    print("\n--- 2. BUSTED-S with Synonymous Rate Variation (SRV, M=3) ---")
    settings = hp.BUSTEDSettings()
    settings.srv = True
    settings.num_rate_classes = 3
    settings.num_syn_rate_classes = 3
    settings.refine_branch_lengths = True

    res_srv = analyzer.run(settings)

    print(f"SRV Unconstrained ln L : {res_srv.unconstrained.log_likelihood:.2f}")
    print(f"SRV Constrained ln L   : {res_srv.constrained.log_likelihood:.2f}")
    print(f"SRV LRT                : {res_srv.lrt:.4f}")
    print(f"SRV p-value            : {res_srv.p_value:.6e}")
    if res_srv.unconstrained.test_distribution.syn_rates:
        print(f"Inferred Syn Rates     : {[round(a, 4) for a in res_srv.unconstrained.test_distribution.syn_rates]}")
        print(f"Inferred Syn Weights   : {[round(q, 4) for q in res_srv.unconstrained.test_distribution.syn_weights]}")

    # Inspect per-site Evidence Ratios (Bayes Factors for positive selection)
    bfs = res_srv.evidence_ratios
    high_bf_sites = [(idx + 1, bf) for idx, bf in enumerate(bfs) if bf > 2.0]
    print(f"\nSites with Evidence Ratio > 2.0: {len(high_bf_sites)}")
    for site, bf in high_bf_sites[:5]:
        print(f"  Site {site:3d}: Bayes Factor = {bf:6.2f}")

    # Export JSON
    out_file = "adh.BUSTED_S.json"
    with open(out_file, "w") as f:
        f.write(analyzer.to_json(res_srv))
    print(f"\nSaved Datamonkey JSON to {out_file}")

if __name__ == "__main__":
    main()
