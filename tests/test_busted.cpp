#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

#include "hyphy/core/types.hpp"
#include "hyphy/core/genetic_code.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/analyses/busted.hpp"

#include <nlohmann/json.hpp>
#include <fstream>
#include <iostream>

using namespace hyphy::core;
using namespace hyphy::analyses;

TEST_CASE("BUSTED: CD2 Dataset (10 taxa, 187 codons)") {
    auto aln = Alignment::load("/Users/sergei/Development/hyphy/hyphy-next/benchmarks/data/cd2.fna");
    auto tree = Tree::from_newick_file("/Users/sergei/Development/hyphy/hyphy-next/benchmarks/data/cd2.nwk");

    auto busted = BUSTEDAnalyzer::create_and_fit(tree, aln);
    auto res = busted.run();

    std::cout << "\n[BUSTED CD2 Results]\n"
              << "  Unconstrained ln L : " << res.unconstrained.log_likelihood << "\n"
              << "  Constrained ln L   : " << res.constrained.log_likelihood << "\n"
              << "  LRT                : " << res.lrt << "\n"
              << "  p-value            : " << res.p_value << "\n"
              << "  Alt omegas         : [" << res.unconstrained.test_distribution.omegas[0] << ", "
                                          << res.unconstrained.test_distribution.omegas[1] << ", "
                                          << res.unconstrained.test_distribution.omegas[2] << "]\n"
              << "  Alt weights        : [" << res.unconstrained.test_distribution.weights[0] << ", "
                                          << res.unconstrained.test_distribution.weights[1] << ", "
                                          << res.unconstrained.test_distribution.weights[2] << "]\n"
              << "  Null omegas        : [" << res.constrained.test_distribution.omegas[0] << ", "
                                          << res.constrained.test_distribution.omegas[1] << ", "
                                          << res.constrained.test_distribution.omegas[2] << "]\n"
              << "  Runtime            : " << res.runtime_seconds << " seconds\n";

    CHECK(res.lrt >= 0.0);
    CHECK(res.p_value >= 0.0);
    CHECK(res.p_value <= 1.0);
    CHECK(res.unconstrained.test_distribution.omegas[2] >= 1.0);
    CHECK(res.constrained.test_distribution.omegas[2] == 1.0);
    CHECK(res.evidence_ratios.size() == 187);

    // Test Legacy JSON export
    auto legacy_out = busted.to_legacy_json(res);
    CHECK(legacy_out.contains("test results"));
    CHECK(legacy_out.contains("fits"));
    CHECK(legacy_out["fits"].contains("Unconstrained model"));
    CHECK(legacy_out["fits"].contains("Constrained model"));

    // Test Modern JSON export
    auto modern_out = busted.to_modern_json(res);
    CHECK(modern_out.contains("statistical_tests"));
    CHECK(modern_out.contains("model_fits"));
    CHECK(modern_out["model_fits"].contains("unconstrained"));
    CHECK(modern_out["model_fits"].contains("constrained"));
    CHECK(modern_out.contains("site_results"));
}

TEST_CASE("BUSTED: ADH Dataset (23 taxa, 254 codons)") {
    auto aln = Alignment::load("/Users/sergei/Development/hyphy/tests/data/adh.nex");
    REQUIRE(!aln.embedded_tree_newick.empty());
    auto tree = Tree::from_newick(aln.embedded_tree_newick);

    auto busted = BUSTEDAnalyzer::create_and_fit(tree, aln);
    auto res = busted.run();

    std::cout << "\n[BUSTED ADH Results]\n"
              << "  Unconstrained ln L : " << res.unconstrained.log_likelihood << "\n"
              << "  Constrained ln L   : " << res.constrained.log_likelihood << "\n"
              << "  LRT                : " << res.lrt << "\n"
              << "  p-value            : " << res.p_value << "\n"
              << "  Alt omegas         : [" << res.unconstrained.test_distribution.omegas[0] << ", "
                                          << res.unconstrained.test_distribution.omegas[1] << ", "
                                          << res.unconstrained.test_distribution.omegas[2] << "]\n"
              << "  Alt weights        : [" << res.unconstrained.test_distribution.weights[0] << ", "
                                          << res.unconstrained.test_distribution.weights[1] << ", "
                                          << res.unconstrained.test_distribution.weights[2] << "]\n"
              << "  Null omegas        : [" << res.constrained.test_distribution.omegas[0] << ", "
                                          << res.constrained.test_distribution.omegas[1] << ", "
                                          << res.constrained.test_distribution.omegas[2] << "]\n"
              << "  Runtime            : " << res.runtime_seconds << " seconds\n";

    CHECK(res.lrt >= 0.0);
    CHECK(res.p_value >= 0.0);
    CHECK(res.p_value <= 1.0);
    CHECK(res.unconstrained.test_distribution.omegas[2] >= 1.0);
    CHECK(res.constrained.test_distribution.omegas[2] == 1.0);
    CHECK(res.evidence_ratios.size() == 254);
    CHECK(!res.unconstrained.branch_lengths.empty());
}

TEST_CASE("BUSTED: Automatic K Model Selection (AICc Step-Up)") {
    auto aln = Alignment::load("/Users/sergei/Development/hyphy/tests/data/adh.nex");
    auto tree = Tree::from_newick(aln.embedded_tree_newick);

    auto busted = BUSTEDAnalyzer::create_and_fit(tree, aln);
    BUSTEDSettings settings;
    settings.auto_select_k = true;
    settings.max_k = 3;
    settings.refine_branch_lengths = true;

    auto res = busted.run(settings);

    std::cout << "\n[BUSTED Auto-K Results]\n"
              << "  Optimal K          : " << res.optimal_k << "\n"
              << "  Unconstrained ln L : " << res.unconstrained.log_likelihood << "\n"
              << "  Constrained ln L   : " << res.constrained.log_likelihood << "\n"
              << "  LRT                : " << res.lrt << "\n"
              << "  p-value            : " << res.p_value << "\n";

    // For ADH, K=2 is strongly preferred over K=1 (delta AICc ~ -91), and K=3 is penalized (+15 AICc)
    CHECK(res.optimal_k == 2);
    CHECK(res.unconstrained.test_distribution.omegas.size() == 2);
    CHECK(res.constrained.test_distribution.omegas.size() == 2);
    CHECK(res.unconstrained.test_distribution.omegas[1] >= 1.0);
    CHECK(res.constrained.test_distribution.omegas[1] == 1.0);
    CHECK(res.lrt >= 0.0);
    CHECK(res.p_value <= 0.05); // Positive selection still robustly detected under optimal K=2
}

TEST_CASE("BUSTED-S: Synonymous Rate Variation (SRV) on ADH") {
    auto aln = Alignment::load("/Users/sergei/Development/hyphy/tests/data/adh.nex");
    auto tree = Tree::from_newick(aln.embedded_tree_newick);

    auto busted = BUSTEDAnalyzer::create_and_fit(tree, aln);
    BUSTEDSettings settings;
    settings.srv = true;
    settings.num_rate_classes = 3;
    settings.num_syn_rate_classes = 3;
    settings.refine_branch_lengths = true;

    auto res = busted.run(settings);

    std::cout << "\n[BUSTED-S (SRV) ADH Results]\n"
              << "  Unconstrained ln L : " << res.unconstrained.log_likelihood << "\n"
              << "  Constrained ln L   : " << res.constrained.log_likelihood << "\n"
              << "  LRT                : " << res.lrt << "\n"
              << "  p-value            : " << res.p_value << "\n";
    if (!res.unconstrained.test_distribution.syn_rates.empty()) {
        std::cout << "  Syn rates          : ["
                  << res.unconstrained.test_distribution.syn_rates[0] << ", "
                  << res.unconstrained.test_distribution.syn_rates[1] << ", "
                  << res.unconstrained.test_distribution.syn_rates[2] << "]\n"
                  << "  Syn weights        : ["
                  << res.unconstrained.test_distribution.syn_weights[0] << ", "
                  << res.unconstrained.test_distribution.syn_weights[1] << ", "
                  << res.unconstrained.test_distribution.syn_weights[2] << "]\n";
    }

    CHECK(res.lrt >= 0.0);
    CHECK(res.p_value >= 0.0);
    CHECK(res.p_value <= 1.0);
    CHECK(res.unconstrained.test_distribution.syn_rates.size() == 3);
    CHECK(res.unconstrained.test_distribution.syn_weights.size() == 3);

    // Verify mean synonymous rate is ~ 1.0
    Scalar mean_syn = 0.0;
    for (size_t m = 0; m < 3; ++m) {
        mean_syn += res.unconstrained.test_distribution.syn_rates[m] * res.unconstrained.test_distribution.syn_weights[m];
    }
    CHECK(std::abs(mean_syn - 1.0) < 0.05);

    auto j = busted.to_legacy_json(res);
    CHECK(j["analysis"]["settings"]["srv"] == "Yes");
    CHECK(j["fits"]["Unconstrained model"]["Rate Distributions"].contains("Synonymous site-to-site rates"));

    auto j_mod = busted.to_modern_json(res);
    CHECK(j_mod["analysis"]["settings"]["srv"] == true);
    CHECK(j_mod["model_fits"]["unconstrained"]["rate_distributions"].contains("synonymous"));
}

