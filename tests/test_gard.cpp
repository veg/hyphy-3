#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "hyphy/core/types.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/core/neighbor_joining.hpp"
#include "hyphy/analyses/gard.hpp"

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

TEST_CASE("GARD: Nucleotide Alignment Loading & Slicing") {
    std::string adh_path = resolve_test_path("tests/data/adh.nex");
    auto aln = Alignment::load_nucleotide(adh_path);

    CHECK(aln.num_taxa == 23);
    CHECK(aln.num_nucleotides == 762);
    CHECK(aln.num_codons == 0);
    CHECK(aln.patterns.size() > 0);

    // Verify slicing
    auto slice1 = aln.slice_nucleotide(0, 460);
    CHECK(slice1.num_taxa == 23);
    CHECK(slice1.num_nucleotides == 460);

    auto slice2 = aln.slice_nucleotide(460, 762);
    CHECK(slice2.num_taxa == 23);
    CHECK(slice2.num_nucleotides == 302);
}

TEST_CASE("GARD: Neighbor-Joining Tree Inference") {
    std::string adh_path = resolve_test_path("tests/data/adh.nex");
    auto aln = Alignment::load_nucleotide(adh_path);

    Matrix D = NeighborJoining::compute_tn93_distances(aln);
    CHECK(D.rows() == 23);
    CHECK(D.cols() == 23);

    // Verify distance matrix symmetry and positive off-diagonals
    for (size_t i = 0; i < 23; ++i) {
        CHECK(D(i, i) == doctest::Approx(0.0));
        for (size_t j = i + 1; j < 23; ++j) {
            CHECK(D(i, j) == doctest::Approx(D(j, i)));
            CHECK(D(i, j) > 0.0);
        }
    }

    Tree tree = NeighborJoining::infer(aln);
    CHECK(tree.num_leaves() == 23);
    CHECK(tree.num_nodes() == 44); // 23 leaves + 21 internal nodes (including trifurcating root)
    CHECK(tree.root_id != INVALID_INDEX);

    // Verify all branch lengths are strictly positive
    for (const auto& node : tree.nodes) {
        if (node.id != tree.root_id) {
            CHECK(node.branch_length > 0.0);
        }
    }
}

TEST_CASE("GARD: Negative Control on HRVI.nex (0 Breakpoints Expected)") {
    std::string hrvi_path = resolve_test_path("tests/data/HRVI.nex");
    auto aln = Alignment::load_nucleotide(hrvi_path);

    CHECK(aln.num_taxa == 12);
    CHECK(aln.num_nucleotides == 428);

    GARDSettings settings;
    settings.input_file = hrvi_path;
    settings.run_mode = "Faster";
    settings.quiet = true;

    GARDAnalyzer analyzer(aln, settings);
    GARDResult result = analyzer.run();

    // Baseline c-AIC should be close to ~1450.97
    CHECK(result.baseline_aicc > 1400.0);
    CHECK(result.baseline_aicc < 1500.0);

    // No breakpoint should improve upon baseline
    CHECK(result.best_model.breakpoints.empty());
    CHECK(result.best_model.aicc == doctest::Approx(result.baseline_aicc));
    CHECK(result.topological_incongruence == false);
    CHECK(result.topology_test_result.run == false);
    CHECK(result.topology_test_result.breakpoint_tests.empty());
    CHECK(result.topology_test_result.any_incongruence_significant == false);

    // Check JSON serialization
    auto modern_json = result.to_json();
    CHECK(modern_json["analysis"]["id"] == "gard");
    CHECK(modern_json["dataset"]["taxa_count"] == 12);
    CHECK(modern_json["dataset"]["sites_count"] == 428);
    CHECK(modern_json["model_selection"]["best_model_breakpoints"].empty());

    auto legacy_json = result.to_legacy_json();
    CHECK(legacy_json["breakpointData"]["0"]["bps"][0][0] == "1");
    CHECK(legacy_json["breakpointData"]["0"]["bps"][0][1] == "428");
}

TEST_CASE("GARD: Positive Control on adh.nex (Breakpoint at Site 460)") {
    std::string adh_path = resolve_test_path("tests/data/adh.nex");
    auto aln = Alignment::load_nucleotide(adh_path);

    GARDSettings settings;
    settings.input_file = adh_path;
    settings.run_mode = "Faster";
    settings.max_breakpoints = 2;
    settings.quiet = true;

    GARDAnalyzer analyzer(aln, settings);
    GARDResult result = analyzer.run();

    CHECK(result.num_taxa == 23);
    CHECK(result.num_sites == 762);
    CHECK(result.baseline_aicc > 0.0);

    // Recombination should be detected
    CHECK(!result.best_model.breakpoints.empty());

    // The primary breakpoint on adh is known to be site 460 (0-based 459, or within small neighborhood)
    size_t primary_bp = result.best_model.breakpoints[0];
    CHECK(primary_bp >= 450);
    CHECK(primary_bp <= 470);

    // Significant improvement in c-AIC
    CHECK(result.best_model.aicc < result.baseline_aicc);
    CHECK(result.improvements.size() >= 1);
    CHECK(result.improvements[0].second > 20.0); // Substantial delta-AICc (empirically ~22.92)

    // Topological incongruence should be supported
    CHECK(result.topological_incongruence == true);
    CHECK(result.best_model.aicc < result.single_tree_aicc);

    // Topological incongruence tests (KH / SH with RELL bootstrap)
    CHECK(result.topology_test_result.run == true);
    CHECK(!result.topology_test_result.breakpoint_tests.empty());
    CHECK(result.topology_test_result.breakpoint_tests.size() == result.best_model.breakpoints.size());
    const auto& bp_test0 = result.topology_test_result.breakpoint_tests[0];
    CHECK(bp_test0.delta_lnl_left > 0.0);
    CHECK(bp_test0.delta_lnl_right > 0.0);
    CHECK(bp_test0.kh_p_value_left >= 0.0);
    CHECK(bp_test0.kh_p_value_left <= 1.0);
    CHECK(bp_test0.sh_p_value_left >= 0.0);
    CHECK(bp_test0.sh_p_value_left <= 1.0);
    CHECK(bp_test0.statistically_significant == true);
    CHECK(result.topology_test_result.any_incongruence_significant == true);
    CHECK(!result.topology_test_result.partition_matrix.empty());

    // Verify partitions
    CHECK(result.best_model.partitions.size() >= 2);
    CHECK(result.best_model.partitions[0].start_site == 1);
    CHECK(result.best_model.partitions[0].end_site == primary_bp + 1);
    CHECK(result.best_model.partitions[1].start_site == primary_bp + 2);

    // Modern JSON verification
    auto adh_modern_json = result.to_json();
    CHECK(adh_modern_json.contains("topology_test"));
    CHECK(adh_modern_json["topology_test"]["topological_incongruence_supported"] == true);
    CHECK(adh_modern_json["topology_test"]["any_incongruence_significant"] == true);
    CHECK(adh_modern_json["topology_test"]["breakpoints"].size() == result.best_model.breakpoints.size());

    // Issue #5: Verify legacy JSON has all partitions, trees, improvements, and topologyTests
    auto legacy_json = result.to_legacy_json();
    CHECK(legacy_json["breakpointData"].size() == result.best_model.partitions.size());
    CHECK(legacy_json["trees"].size() == result.best_model.partitions.size());
    CHECK(legacy_json["improvements"].size() == result.improvements.size());
    CHECK(legacy_json.contains("topologyTests"));
    CHECK(!legacy_json["topologyTests"].empty());
    for (auto& [k, imp] : legacy_json["improvements"].items()) {
        CHECK(imp.contains("breakpoints"));
        CHECK(imp["breakpoints"].is_array());
        CHECK(!imp["breakpoints"].empty());
    }
}

TEST_CASE("GARD: Dynamic Breakpoint Chart and Visualizer") {
    // 1. Test single-breakpoint dashboard generation
    GARDVisualizer viz(762, 23, false, true, true);
    viz.chart_width = 65;
    viz.start_single_bp(315, 10378.22);

    viz.record_candidate({459}, 10353.93);
    viz.update_best({459}, 10353.93);
    viz.tick(50);

    auto lines = viz.format_dashboard_lines();
    CHECK(lines.size() == 5);

    // Line 0: Status
    CHECK(lines[0].find("Status:") != std::string::npos);
    CHECK(lines[0].find("10353.93") != std::string::npos);
    CHECK(lines[0].find("460") != std::string::npos);

    // Line 1: Ruler
    CHECK(lines[1].find("Ruler:") != std::string::npos);
    CHECK(lines[1].find("1") != std::string::npos);
    CHECK(lines[1].find("460") != std::string::npos);
    CHECK(lines[1].find("762 bp") != std::string::npos);

    // Line 2: Track
    CHECK(lines[2].find("Track:") != std::string::npos);
    CHECK(lines[2].find("├") != std::string::npos);
    CHECK(lines[2].find("┤") != std::string::npos);
    CHECK(lines[2].find("▲") != std::string::npos);

    // Line 3: Density
    CHECK(lines[3].find("Density:") != std::string::npos);

    // Line 4: Partitions
    CHECK(lines[4].find("Parts:") != std::string::npos);
    CHECK(lines[4].find("[1..460] (460bp)") != std::string::npos);
    CHECK(lines[4].find("[461..762] (302bp)") != std::string::npos);

    // 2. Test multi-breakpoint dashboard generation
    viz.start_ga(2, 50, 10353.93, {459});
    viz.record_candidate({239, 519}, 10340.12);
    viz.update_best({239, 519}, 10340.12);
    viz.tick(100);

    auto ga_lines = viz.format_dashboard_lines();
    CHECK(ga_lines.size() == 5);
    CHECK(ga_lines[1].find("240") != std::string::npos);
    CHECK(ga_lines[1].find("520") != std::string::npos);
    CHECK(ga_lines[4].find("[1..240] (240bp)") != std::string::npos);
    CHECK(ga_lines[4].find("[241..520] (280bp)") != std::string::npos);
    CHECK(ga_lines[4].find("[521..762] (242bp)") != std::string::npos);
}

TEST_CASE("GARD: High Breakpoint Count GUI Stability & Adaptive Resolution Detail (B=8, 9, 12)") {
    GARDVisualizer viz(762, 23, false, true, true);
    viz.chart_width = 65;

    // Test ANSI truncation helper
    std::string styled = Console::bold(Console::brand("1234567890"));
    CHECK(Console::printable_len(styled) == 10);
    std::string truncated = GARDVisualizer::truncate_ansi_line(styled, 5);
    CHECK(Console::printable_len(truncated) == 5);

    // 1. Test B = 8 (9 partitions) - must not wrap in standard 80-col terminal
    std::vector<size_t> bps_8 = {85, 158, 239, 319, 399, 479, 559, 639};
    viz.start_ga(8, 50, 10300.0, bps_8);
    viz.update_ga_generation(25, 30, 8, 0.35, bps_8, 10250.0);

    auto lines_8 = viz.format_dashboard_lines();
    CHECK(lines_8.size() == 5);
    for (size_t i = 0; i < lines_8.size(); ++i) {
        CHECK(Console::printable_len(lines_8[i]) <= 80);
    }
    // Partition line should use adaptive compact/stats tier without overflow
    CHECK(lines_8[4].find("parts") != std::string::npos);

    // 2. Test B = 9 (10 partitions)
    std::vector<size_t> bps_9 = {75, 150, 225, 300, 375, 450, 525, 600, 675};
    viz.start_ga(9, 50, 10250.0, bps_9);
    viz.update_ga_generation(30, 30, 12, 0.28, bps_9, 10220.0);

    auto lines_9 = viz.format_dashboard_lines();
    CHECK(lines_9.size() == 5);
    for (size_t i = 0; i < lines_9.size(); ++i) {
        CHECK(Console::printable_len(lines_9[i]) <= 80);
    }
    CHECK(lines_9[4].find("parts") != std::string::npos);

    // 3. Test B = 12 (13 partitions)
    std::vector<size_t> bps_12 = {50, 110, 170, 230, 290, 350, 410, 470, 530, 590, 650, 710};
    viz.start_ga(12, 50, 10200.0, bps_12);
    viz.update_ga_generation(40, 30, 20, 0.20, bps_12, 10180.0);

    auto lines_12 = viz.format_dashboard_lines();
    CHECK(lines_12.size() == 5);
    for (size_t i = 0; i < lines_12.size(); ++i) {
        CHECK(Console::printable_len(lines_12[i]) <= 80);
    }
    CHECK(lines_12[4].find("parts") != std::string::npos);

    // 4. Test ruler collision avoidance with tightly clustered breakpoints
    std::vector<size_t> clustered_bps = {450, 452, 455};
    viz.start_ga(3, 50, 10200.0, clustered_bps);
    viz.update_best(clustered_bps, 10190.0);
    auto clustered_lines = viz.format_dashboard_lines();
    CHECK(clustered_lines.size() == 5);
    for (size_t i = 0; i < clustered_lines.size(); ++i) {
        CHECK(Console::printable_len(clustered_lines[i]) <= 80);
    }
}

TEST_CASE("GARD: Ohmic Spectral Trees, Rate Variation & Wheatstone Bridges on HRVI.nex") {
    std::string hrvi_path = resolve_test_path("tests/data/HRVI.nex");
    auto aln = Alignment::load_nucleotide(hrvi_path);

    GARDSettings settings;
    settings.input_file = hrvi_path;
    settings.run_mode = "Faster";
    settings.tree_method = GARDTreeMethod::Ohmic;
    settings.srv_mode = GARDSRVMode::DiscreteGamma;
    settings.allow_bridges = true;
    settings.quiet = true;

    GARDAnalyzer analyzer(aln, settings);
    GARDResult result = analyzer.run();

    CHECK(result.tree_method == GARDTreeMethod::Ohmic);
    CHECK(result.srv_mode == GARDSRVMode::DiscreteGamma);
    CHECK(result.gamma_rates.size() == 4);
    CHECK(result.gamma_alpha > 0.01);

    // Wheatstone bridge screening should detect active reticulation
    CHECK(result.ohmic_bridges_detected == true);
    CHECK(result.num_bridges >= 1);
    CHECK(result.ohmic_aicc < result.baseline_aicc);
    CHECK(result.ohmic_delta_aicc > 50.0);

    // Bridge trajectory should seed candidate breakpoint
    CHECK(!result.bridge_candidate_breakpoints.empty());

    // JSON export
    auto j = result.to_json();
    CHECK(j["rate_variation"]["mode"] == "discrete_gamma");
    CHECK(j["rate_variation"]["classes"] == 4);
    CHECK(j["tree_engine"]["method"] == "ohmic_spectral");
    CHECK(j["ohmic_network"]["bridges_detected"] == true);
    CHECK(j["ohmic_network"]["num_bridges"] >= 1);
}

TEST_CASE("GARD: Kishino-Hasegawa (KH) & Shimodaira-Hasegawa (SH) Topology Testing Engine") {
    // 1. Identical models
    std::vector<Scalar> lnl1 = {-1.5, -2.0, -0.8, -3.2, -1.1};
    std::vector<Scalar> lnl2 = {-1.5, -2.0, -0.8, -3.2, -1.1};
    auto kh_ident = GARDTopologyTester::run_kh_test(lnl1, lnl2, 500, 42);
    CHECK(kh_ident.delta_lnl == doctest::Approx(0.0));
    CHECK(kh_ident.p_asymptotic == doctest::Approx(1.0));
    CHECK(kh_ident.p_rell == doctest::Approx(1.0));

    auto sh_ident = GARDTopologyTester::run_sh_test({lnl1, lnl2}, 500, 42);
    CHECK(sh_ident.best_index == 0);
    CHECK(sh_ident.p_values[0] == doctest::Approx(1.0));
    CHECK(sh_ident.p_values[1] == doctest::Approx(1.0));

    // 2. Decisively superior tree across sites
    std::vector<Scalar> lnl_good(100, -1.0);
    std::vector<Scalar> lnl_bad(100, -3.0);
    auto kh_sup = GARDTopologyTester::run_kh_test(lnl_good, lnl_bad, 500, 42);
    CHECK(kh_sup.delta_lnl == doctest::Approx(200.0));
    CHECK(kh_sup.p_asymptotic < 1e-6);
    CHECK(kh_sup.p_rell < 0.01);

    auto sh_sup = GARDTopologyTester::run_sh_test({lnl_good, lnl_bad}, 500, 42);
    CHECK(sh_sup.best_index == 0);
    CHECK(sh_sup.p_values[0] == doctest::Approx(1.0));
    CHECK(sh_sup.p_values[1] < 0.01);

    // 3. Multi-tree candidate set (3 trees) with site variance
    std::vector<Scalar> t1(50, -1.0);
    std::vector<Scalar> t2 = t1;
    for (size_t i = 0; i < 50; ++i) {
        if (i % 2 == 0) t2[i] += 0.05; // slightly better at even sites
        else t2[i] -= 0.06;           // slightly worse at odd sites
    }
    std::vector<Scalar> t3(50, -3.0);  // decisively worse across all sites
    auto sh_multi = GARDTopologyTester::run_sh_test({t1, t2, t3}, 500, 42);
    CHECK(sh_multi.best_index == 0);
    CHECK(sh_multi.p_values[0] == doctest::Approx(1.0));
    CHECK(sh_multi.p_values[1] > 0.05); // Not rejected due to variance
    CHECK(sh_multi.p_values[2] < 0.01); // Decisively rejected
}

