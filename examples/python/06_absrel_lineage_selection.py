#!/usr/bin/env python3
"""
Example 6: Adaptive Branch-Site Random Effects Likelihood (aBSREL) in Python
Infers the optimal number of rate categories per branch and tests each lineage
for episodic diversifying selection with Holm-Bonferroni multiple testing correction.
"""

import sys
from pathlib import Path

repo_root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(repo_root / "python"))

import hyphy3 as hp

def main():
    data_path = str(repo_root / "benchmarks" / "data" / "bglobin.nex")
    print(f"Loading alignment and tree from: {data_path}")

    aln = hp.Alignment.load(data_path)
    tree = hp.Tree.from_newick(aln.embedded_tree_newick)

    print(f"Dataset: {aln.num_taxa} sequences, {aln.num_codons} codons, {tree.num_nodes()} tree nodes.")

    # 1. Configure aBSREL Settings
    settings = hp.ABSRELSettings()
    settings.max_rate_classes = 3
    settings.p_threshold = 0.05
    settings.test_branches = "All"

    print("\n--- Running aBSREL Inference & Branch Testing ---")
    analyzer = hp.ABSRELAnalyzer.create(tree, aln, settings)
    res = analyzer.run()

    # 2. Inspect Model Fits
    print(f"\nModel Fits:")
    print(f"  Nucleotide GTR Log-L    : {res.gtr_fit.log_likelihood:.2f} (AICc: {res.gtr_fit.aicc:.2f})")
    print(f"  Baseline MG94xREV Log-L : {res.baseline_fit.log_likelihood:.2f} (AICc: {res.baseline_fit.aicc:.2f})")
    print(f"  Full Adaptive Log-L     : {res.full_adaptive_fit.log_likelihood:.2f} (AICc: {res.full_adaptive_fit.aicc:.2f})")
    print(f"  Tested Branches         : {len(res.tested_branches)}")
    print(f"  Positive Branches       : {len(res.positive_branches)}")
    print(f"  Execution Time          : {res.runtime_seconds:.2f}s")

    # 3. Report Significant Branches
    print("\nBranches detected under episodic positive selection (p_adj <= 0.05):")
    for bname in res.positive_branches:
        bres = res.branches[bname]
        rates_str = ", ".join(f"w_{k+1}={r:.2f} ({p*100:.1f}%)" for k, (r, p) in enumerate(zip(bres.rate_distribution.rates, bres.rate_distribution.weights)))
        print(f"  * {bname:10s} | LRT={bres.lrt:5.2f} | p_raw={bres.uncorrected_p_value:.3e} | p_adj={bres.corrected_p_value:.3e} | {bres.sites_ebf_100} sites @ EBF>=100")
        print(f"    Distribution: [{rates_str}]")

    # 4. JSON Serialization
    json_out = res.to_json(tree, aln)
    out_file = repo_root / "bglobin.ABSREL.json"
    with open(out_file, "w") as f:
        f.write(json_out)
    print(f"\nSaved full Datamonkey JSON to: {out_file}")

if __name__ == "__main__":
    main()
