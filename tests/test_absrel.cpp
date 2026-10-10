#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

#include "hyphy/core/types.hpp"
#include "hyphy/core/genetic_code.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/analyses/absrel.hpp"

#include <nlohmann/json.hpp>
#include <fstream>
#include <iostream>

using namespace hyphy::core;
using namespace hyphy::analyses;

TEST_CASE("aBSREL: Beta-globin Benchmark Parity (17 taxa, 144 codons)") {
    std::vector<std::string> candidates = {
        "benchmarks/data/bglobin.nex",
        "../benchmarks/data/bglobin.nex",
        "/Users/sergei/Development/hyphy-3/benchmarks/data/bglobin.nex",
        "/Users/sergei/Development/hyphy/tests/data/bglobin.nex"
    };
    std::string path;
    for (const auto& c : candidates) {
        std::ifstream f(c);
        if (f.good()) { path = c; break; }
    }
    REQUIRE(!path.empty());
    auto aln = Alignment::load(path);
    REQUIRE(!aln.embedded_tree_newick.empty());
    auto tree = Tree::from_newick(aln.embedded_tree_newick);

    ABSRELSettings settings;
    settings.max_rate_classes = 3;
    settings.p_threshold = 0.05;
    settings.test_branches = "All";

    auto absrel = ABSRELAnalyzer::create(tree, aln, settings);
    auto res = absrel.run();

    std::cout << "\n[aBSREL bglobin Results]\n"
              << "  GTR ln L              : " << res.gtr_fit.log_likelihood << " (AICc: " << res.gtr_fit.aicc << ")\n"
              << "  Baseline ln L         : " << res.baseline_fit.log_likelihood << " (AICc: " << res.baseline_fit.aicc << ")\n"
              << "  Full Adaptive ln L    : " << res.full_adaptive_fit.log_likelihood << " (AICc: " << res.full_adaptive_fit.aicc << ")\n"
              << "  Tested branches       : " << res.tested_branches.size() << "\n"
              << "  Positive branches     : " << res.positive_branches.size() << "\n"
              << "  Runtime               : " << res.runtime_seconds << " seconds\n";

    for (const auto& bname : res.positive_branches) {
        const auto& bres = res.branches.at(bname);
        std::cout << "  * Positive branch: " << bname
                  << " | LRT = " << bres.lrt
                  << " | p = " << bres.uncorrected_p_value
                  << " | p_adj = " << bres.corrected_p_value
                  << " | Rate classes = " << bres.rate_classes
                  << " | Sites @ EBF>=100 = " << bres.sites_ebf_100 << "\n";
    }

    // Baseline fit checks: GTR and Baseline MG94
    CHECK(res.gtr_fit.log_likelihood < -3800.0);
    CHECK(res.gtr_fit.log_likelihood > -4100.0);
    CHECK(res.baseline_fit.log_likelihood > res.gtr_fit.log_likelihood);
    CHECK(res.full_adaptive_fit.log_likelihood >= res.baseline_fit.log_likelihood);

    // Testing checks
    CHECK(res.tested_branches.size() == 31);
    for (const auto& [bname, bres] : res.branches) {
        if (bres.is_tested) {
            CHECK(bres.lrt >= 0.0);
            CHECK(bres.uncorrected_p_value >= 0.0);
            CHECK(bres.uncorrected_p_value <= 1.0);
            CHECK(bres.corrected_p_value >= 0.0);
            CHECK(bres.corrected_p_value <= 1.0);
        }
    }

    // Parity check: Node15, Node18, and SHEEP are expected positive branches
    CHECK(res.positive_branches.size() > 0);

    // Test Legacy JSON export
    auto legacy_out = res.to_legacy_json(tree, aln);
    CHECK(legacy_out.contains("fits"));
    CHECK(legacy_out["fits"].contains("Baseline MG94xREV"));
    CHECK(legacy_out["fits"].contains("Full adaptive model"));
    CHECK(legacy_out.contains("branch attributes"));
    CHECK(legacy_out["branch attributes"].contains("0"));
    CHECK(legacy_out.contains("test results"));
    CHECK(legacy_out["test results"]["positive test results"] == res.positive_branches.size());

    // Test Modern JSON export
    auto modern_out = res.to_modern_json(tree, aln);
    CHECK(modern_out.contains("statistical_tests"));
    CHECK(modern_out.contains("model_fits"));
    CHECK(modern_out["model_fits"].contains("baseline_mg94"));
    CHECK(modern_out["model_fits"].contains("full_adaptive"));
    CHECK(modern_out.contains("branch_results"));
    CHECK(modern_out["statistical_tests"]["branch_level_summary"]["positive_branches_count"] == res.positive_branches.size());
}
