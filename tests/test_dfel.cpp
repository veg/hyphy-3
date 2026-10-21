#include "hyphy/analyses/directional_fel.hpp"
#undef NDEBUG
#include <cassert>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace hyphy::core;
using namespace hyphy::analyses;

inline std::string resolve_test_path(const std::string& path) {
    if (std::ifstream(path).good()) return path;
    if (std::ifstream("../" + path).good()) return "../" + path;
    if (std::ifstream("../../" + path).good()) return "../../" + path;
    return path;
}

void test_chi2_p_value() {
    std::cout << "Testing Chi-Square p-value calculations...\n";

    // df = 1: p-value of 3.84146 is ~0.05
    assert(std::abs(chi2_p_value(3.84146, 1) - 0.05) < 0.001);

    // df = 2: p-value of 5.991 is ~0.05
    assert(std::abs(chi2_p_value(5.99146, 2) - 0.05) < 0.001);

    // df = 10: p-value of 18.307 is ~0.05
    assert(std::abs(chi2_p_value(18.307, 10) - 0.05) < 0.002);

    // Incomplete gamma boundary conditions
    assert(chi2_p_value(0.0, 5) == 1.0);
    assert(chi2_p_value(100.0, 2) < 1e-15);

    std::cout << "  Passed Chi-Square p-value tests.\n";
}

void test_compositional_heterogeneity() {
    std::cout << "Testing Compositional Heterogeneity Screening...\n";

    Alignment aln = Alignment::from_nexus(resolve_test_path("tests/data/adh.nex"));
    auto comp = DirectionalFELAnalyzer::run_composition_test(aln, 0.05);

    assert(comp.taxa.size() == aln.num_taxa);
    assert(comp.df_overall == static_cast<int>(3 * (aln.num_taxa - 1)));
    assert(comp.df_cp1 == static_cast<int>(3 * (aln.num_taxa - 1)));
    assert(comp.df_cp2 == static_cast<int>(3 * (aln.num_taxa - 1)));
    assert(comp.df_cp3 == static_cast<int>(3 * (aln.num_taxa - 1)));

    assert(comp.chi2_overall >= 0.0);
    assert(comp.p_value_overall >= 0.0 && comp.p_value_overall <= 1.0);
    assert(comp.p_value_cp3 >= 0.0 && comp.p_value_cp3 <= 1.0);

    for (const auto& tc : comp.taxa) {
        (void)tc;
        assert(tc.total_codons == aln.num_codons);
        assert(tc.gc_overall >= 0.0 && tc.gc_overall <= 1.0);
        assert(tc.gc_cp3 >= 0.0 && tc.gc_cp3 <= 1.0);
    }

    std::cout << "  Overall Chi2: " << comp.chi2_overall << ", df: " << comp.df_overall 
              << ", p: " << comp.p_value_overall << "\n";
    std::cout << "  CP3 Chi2: " << comp.chi2_cp3 << ", df: " << comp.df_cp3 
              << ", p: " << comp.p_value_cp3 << "\n";
    std::cout << "  Passed Compositional Heterogeneity tests.\n";
}

void test_tree_rerooting() {
    std::cout << "Testing Tree Re-rooting & Invariant Geometry...\n";

    Alignment aln = Alignment::from_nexus(resolve_test_path("tests/data/adh.nex"));
    Tree tree = Tree::from_newick(aln.embedded_tree_newick);

    Scalar orig_tree_len = 0.0;
    for (const auto& node : tree.nodes) {
        if (node.id != tree.root_id) orig_tree_len += node.branch_length;
    }

    // Pick an edge connecting child 0 of root to its neighbor
    int32_t c0 = tree.nodes[tree.root_id].children[0];
    int32_t c1 = tree.nodes[tree.root_id].children[1];

    Tree rerooted = DirectionalFELAnalyzer::reroot_on_edge(tree, c0, c1, 0.5);

    Scalar new_tree_len = 0.0;
    for (const auto& node : rerooted.nodes) {
        if (node.id != rerooted.root_id) new_tree_len += node.branch_length;
    }

    // Tree length must be conserved
    assert(std::abs(orig_tree_len - new_tree_len) < 1e-6);

    // Number of leaves must be conserved
    assert(rerooted.num_leaves() == tree.num_leaves());
    for (const auto& [name, id] : tree.leaf_name_to_id) {
        (void)id;
        assert(rerooted.leaf_name_to_id.count(name) > 0);
    }

    std::cout << "  Original tree length: " << orig_tree_len << ", Rerooted: " << new_tree_len << "\n";
    std::cout << "  Passed Tree Re-rooting tests.\n";
}

void test_auto_rooting_and_dfel_fit() {
    std::cout << "Testing Auto-Rooting Selection and dFEL Analyzer...\n";

    Alignment aln = Alignment::from_nexus(resolve_test_path("tests/data/adh.nex"));
    Tree tree = Tree::from_newick(aln.embedded_tree_newick);

    DirectionalFELSettings settings;
    settings.root_strategy = "auto";
    settings.full_model = false; // quick mode for test
    settings.num_regimes = 0; // auto

    auto dfel = DirectionalFELAnalyzer::create_and_fit(
        tree,
        aln,
        settings,
        [](const std::string& msg, double pct) {
            std::cout << "    [" << static_cast<int>(pct * 100) << "%] " << msg << "\n";
        }
    );

    // Verify auto-rooting results
    assert(!dfel.root_result.candidates.empty());
    assert(dfel.root_result.chosen_weight > 0.0);
    Scalar total_w = 0.0;
    for (const auto& c : dfel.root_result.candidates) total_w += c.weight;
    assert(std::abs(total_w - 1.0) < 1e-4);
    (void)total_w;

    std::cout << "  Auto-rooting selected branch: " << dfel.root_result.chosen_branch_name 
              << " with weight " << dfel.root_result.chosen_weight << "\n";

    // Verify regimes
    assert(dfel.regimes.size() >= 1);
    std::cout << "  Activated regimes count: " << dfel.regimes.size() << "\n";
    for (const auto& reg : dfel.regimes) {
        std::cout << "    Regime " << reg.id << " (" << reg.name << "): psi = " 
                  << reg.psi << ", branches = " << reg.branch_count << "\n";
    }

    // Run site testing
    auto results = dfel.run(false, false);
    assert(results.size() == aln.num_codons);

    size_t sig_count = 0;
    for (const auto& r : results) {
        assert(r.alpha >= 0.0);
        assert(r.beta >= 0.0);
        assert(r.lrt >= 0.0);
        assert(r.p_value >= 0.0 && r.p_value <= 1.0);
        if (r.p_value <= settings.p_value_threshold) sig_count++;
    }

    std::cout << "  Analyzed " << results.size() << " codons. Significant sites: " << sig_count << "\n";

    // Verify JSON export
    auto modern_j = dfel.to_modern_json(resolve_test_path("tests/data/adh.nex"), dfel.tree.to_newick(true));
    assert(modern_j["schema_version"] == "3.0.0");
    assert(modern_j["analysis"]["id"] == "dfel");
    assert(modern_j["compositional_heterogeneity"]["taxa"].size() == aln.num_taxa);
    assert(!modern_j["root_selection"]["candidates"].empty());
    assert(modern_j["site_results"]["data"]["site"].size() == aln.num_codons);

    auto legacy_j = dfel.to_legacy_json(resolve_test_path("tests/data/adh.nex"), dfel.tree.to_newick(true));
    assert(legacy_j["MLE"]["content"]["0"].size() == aln.num_codons);

    std::cout << "  Passed Auto-Rooting and dFEL fit tests.\n";
}

int main() {
    try {
        test_chi2_p_value();
        test_compositional_heterogeneity();
        test_tree_rerooting();
        test_auto_rooting_and_dfel_fit();
        std::cout << "\nALL DIRECTIONAL FEL TESTS PASSED!\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Test failed with exception: " << e.what() << "\n";
        return 1;
    }
}
