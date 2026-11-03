#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "hyphy/core/types.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/core/rate_matrix.hpp"
#include "hyphy/analyses/contrast_fel.hpp"

#include <iostream>
#include <fstream>

using namespace hyphy::core;
using namespace hyphy::analyses;

inline std::string resolve_test_path(const std::string& path) {
    if (std::ifstream(path).good()) return path;
    if (std::ifstream("../" + path).good()) return "../" + path;
    if (std::ifstream("../../" + path).good()) return "../../" + path;
    return path;
}

TEST_CASE("Contrast-FEL: Multi-group differential selection on cd2.fna with tree tags") {
    std::string aln_path = resolve_test_path("benchmarks/data/cd2.fna");
    std::string tree_path = resolve_test_path("benchmarks/data/cd2.nwk");

    auto aln = Alignment::from_fasta(aln_path);
    auto tree = Tree::from_newick_file(tree_path);

    ContrastFELSettings settings;
    settings.p_value_threshold = 0.10;
    settings.q_value_threshold = 0.20;
    settings.shrinkage_lambda = 1.0;
    settings.full_model = false; // Fast fitting for unit test

    auto analyzer = ContrastFELAnalyzer::create_and_fit(tree, aln, settings);

    // Verify branch partitioning
    CHECK(analyzer.testable_groups.size() >= 2);
    std::cout << "Detected " << analyzer.testable_groups.size() << " testable branch sets:\n";
    for (const auto& grp : analyzer.testable_groups) {
        std::cout << "  - " << grp.name << ": " << grp.branch_count << " branches\n";
    }

    auto results = analyzer.run(false, false);
    CHECK(results.size() == aln.num_codons);

    // Check invariants
    size_t invariant_count = 0;
    for (const auto& r : results) {
        if (r.contrast_tag == "Invariant") {
            CHECK(r.alpha == 0.0);
            CHECK(r.p_value == 1.0);
            invariant_count++;
        }
    }
    CHECK(invariant_count > 0);
    std::cout << "Invariant sites correctly identified: " << invariant_count << "\n";

    // Verify Ridge Shrinkage effect
    size_t penalized_sites = 0;
    for (const auto& r : results) {
        if (r.shrinkage_penalty > 0.0) {
            penalized_sites++;
            // Penalized LRT must be <= raw LRT
            CHECK(r.penalized_lrt <= r.raw_lrt + 1e-4);
        }
    }
    std::cout << "Sites regularized by Ridge Shrinkage: " << penalized_sites << "\n";
    CHECK(penalized_sites > 0);

    // Verify summary and JSON serialization
    auto summary = analyzer.get_summary();
    CHECK(summary.total_sites == aln.num_codons);

    auto json = analyzer.to_modern_json(aln_path, tree_path);
    CHECK(json["$schema"] == "https://raw.githubusercontent.com/veg/hyphy-3/main/schemas/v3/contrast_fel.v3.schema.json");
    CHECK(json["schema_version"] == "3.0.0");
    CHECK(json["analysis"]["id"] == "contrast-fel");
    CHECK(json["analysis"]["shrinkage_lambda"] == 1.0);
    CHECK(json["software"]["name"] == "hyphy");
    CHECK(json.contains("provenance"));
    CHECK(json["site_results"].contains("columns"));
    CHECK(json["site_results"]["row_count"] == aln.num_codons);
    CHECK(json["site_results"]["data"]["site"].size() == aln.num_codons);
    // Dynamic per-group beta columns must exist for every testable group
    for (const auto& grp : analyzer.testable_groups) {
        CHECK(json["site_results"]["data"].contains("beta_" + grp.name));
    }
}

TEST_CASE("Contrast-FEL: Internal vs Leaves 2-group contrast") {
    std::string aln_path = resolve_test_path("benchmarks/data/cd2.fna");
    std::string tree_path = resolve_test_path("benchmarks/data/cd2.nwk");

    auto aln = Alignment::from_fasta(aln_path);
    auto tree = Tree::from_newick_file(tree_path);

    ContrastFELSettings settings;
    settings.p_value_threshold = 0.05;
    settings.q_value_threshold = 0.20;
    settings.shrinkage_lambda = 0.5;
    settings.branch_specs = "Internal,Leaves";
    settings.full_model = false;

    auto analyzer = ContrastFELAnalyzer::create_and_fit(tree, aln, settings);

    CHECK(analyzer.testable_groups.size() == 2);
    CHECK(analyzer.testable_groups[0].name == "Internal");
    CHECK(analyzer.testable_groups[1].name == "Leaves");

    auto results = analyzer.run(false, false);
    CHECK(results.size() == aln.num_codons);

    // Verify that every site has exactly 2 group betas
    for (const auto& r : results) {
        CHECK(r.beta.size() == 2);
        CHECK(r.q_value >= 0.0);
        CHECK(r.q_value <= 1.0);
    }
}

TEST_CASE("Issue #6: Contrast-FEL --branches subset selection, pairwise tests, and legacy JSON") {
    std::string aln_path = HYPHY3_SOURCE_DIR "/benchmarks/data/cd2_reduced.fna";
    std::string nwk = "((((Pig{G1}:0.147969,Cow{G1}:0.213430):0.085099,Horse{G2}:0.165787,Cat{G2}:0.264806):0.058611,((RhMonkey{G3}:0.002015,Baboon{G3}:0.003108):0.022733,(Human{G3}:0.004349,Chimp{G3}:0.000799):0.011873):0.101856):0.340802,Rat:0.050958,Mouse:0.097950);";

    auto aln = Alignment::from_fasta(aln_path);
    auto tree = Tree::from_newick(nwk);

    // 1. Selecting 2 of the 3 groups puts G3 + unlabeled into Background
    ContrastFELSettings s2;
    s2.full_model = false;
    s2.branch_specs = "G1,G2";
    ContrastFELAnalyzer a2(tree, aln, MG94Parameters{}, s2);
    REQUIRE(a2.testable_groups.size() == 2);
    CHECK(a2.testable_groups[0].name == "G1");
    CHECK(a2.testable_groups[1].name == "G2");
    CHECK(a2.has_background == true);

    // 2. Selecting all 3 groups produces C(3,2)=3 pairwise tests per site
    ContrastFELSettings s3;
    s3.full_model = false;
    s3.branch_specs = "G1,G2,G3";
    ContrastFELAnalyzer a3(tree, aln, MG94Parameters{}, s3);
    REQUIRE(a3.testable_groups.size() == 3);
    auto res3 = a3.run(false, false);
    REQUIRE(!res3.empty());
    CHECK(res3[0].pairwise.size() == 3);

    auto legacy = a3.to_legacy_json(aln_path, nwk);
    CHECK(legacy["MLE"]["headers"].is_array());
    CHECK(legacy["tested"]["0"].size() > 0);

    // 3. Unknown group throws std::runtime_error
    ContrastFELSettings s_bad;
    s_bad.branch_specs = "G1,NonExistentGroup";
    CHECK_THROWS_AS(ContrastFELAnalyzer(tree, aln, MG94Parameters{}, s_bad), std::runtime_error);
}

