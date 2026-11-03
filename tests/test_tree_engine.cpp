#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "hyphy/core/types.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/core/rate_matrix.hpp"
#include "hyphy/analyses/tree_engine.hpp"

using namespace hyphy::core;
using namespace hyphy::analyses;

TEST_CASE("TreeEngine: Spectral Scaffolding & Native Polytomy Detection") {
    std::string aln_path = HYPHY3_SOURCE_DIR "/benchmarks/data/cd2.fna";
    auto aln = Alignment::from_fasta(aln_path);

    CHECK(aln.num_taxa == 10);
    CHECK(aln.num_nucleotides == 561);

    // Build spectral scaffold
    Tree tree = TreeEngine::build_spectral_scaffold(aln, true, 0.05, 0.15);

    CHECK(tree.num_leaves() == 10);
    CHECK(tree.root_id != INVALID_INDEX);
    CHECK(tree.nodes.size() >= 11);

    // Verify all leaves are present
    for (const auto& seq : aln.sequences) {
        CHECK(tree.leaf_name_to_id.find(seq.name) != tree.leaf_name_to_id.end());
    }

    std::string nwk = tree.to_newick(false);
    CHECK(!nwk.empty());
    CHECK(nwk.back() == ';');
}

TEST_CASE("TreeEngine: Simultaneous All-Branch Quasi-Newton (L-BFGS) Optimization") {
    std::string aln_path = HYPHY3_SOURCE_DIR "/benchmarks/data/cd2.fna";
    auto aln = Alignment::from_fasta(aln_path).to_nucleotide_alignment();

    Tree tree = TreeEngine::build_spectral_scaffold(aln, true);

    GTRMatrix gtr;
    GTRParameters gtr_p;
    gtr.update(gtr_p, aln.nuc_frequencies);

    Scalar init_ll = LikelihoodEngine::compute_gtr_log_likelihood(tree, aln, gtr);
    CHECK(std::isfinite(init_ll));

    // Optimize all branches simultaneously
    Scalar opt_ll = TreeEngine::optimize_branches_gtr(tree, aln, gtr, 35, 1e-4);

    CHECK(std::isfinite(opt_ll));
    CHECK(opt_ll > init_ll); // Log-likelihood strictly improves!
    CHECK(opt_ll > -3850.0);  // Initial unoptimized GTR rates on CD2 gives ~ -3785
}

TEST_CASE("TreeEngine: End-to-End Pipeline Execution on CD2 and ADH") {
    TreeEngine engine;
    engine.settings.quiet = false;
    engine.settings.lbfgs_max_iterations = 30;

    SUBCASE("CD2 10 taxa") {
        std::string aln_path = HYPHY3_SOURCE_DIR "/benchmarks/data/cd2.fna";
        auto aln = Alignment::from_fasta(aln_path);

        auto res = engine.run(aln);

        CHECK(res.num_taxa == 10);
        CHECK(res.optimized_gtr_log_l > res.initial_gtr_log_l);
        CHECK(res.final_log_l > -4500.0);
        CHECK(res.num_final_nodes > 0);
        CHECK(!res.optimized_tree.to_newick().empty());

        auto json_out = engine.to_json(res);
        CHECK(json_out["analysis"]["name"] == "HyPhy-3 ML/Ohmic Tree Engine");
        CHECK(json_out["likelihoods"]["final"].get<double>() == doctest::Approx(res.final_log_l));
    }

    SUBCASE("ADH 23 taxa") {
        std::string aln_path = HYPHY3_SOURCE_DIR "/benchmarks/data/adh.fna";
        auto aln = Alignment::from_fasta(aln_path);

        auto res = engine.run(aln);

        CHECK(res.num_taxa == 23);
        CHECK(res.optimized_gtr_log_l > res.initial_gtr_log_l);
        CHECK(res.final_log_l > -7000.0);
        CHECK(res.runtime_seconds > 0.0);
    }
}

TEST_CASE("TreeEngine: Unified CLI Tree Resolution Helper (--tree auto)") {
    std::string aln_path = HYPHY3_SOURCE_DIR "/benchmarks/data/cd2.fna";
    auto aln = Alignment::from_fasta(aln_path);

    SUBCASE("Explicit --tree auto") {
        std::string source_desc;
        Tree tree = TreeEngine::resolve_tree(aln, "auto", false, false, "", &source_desc, true);
        CHECK(tree.num_leaves() == 10);
        CHECK(source_desc.find("TreeEngine") != std::string::npos);
        CHECK(source_desc.find("auto") != std::string::npos);
    }

    SUBCASE("Explicit --tree AUTO case-insensitivity") {
        std::string source_desc;
        Tree tree = TreeEngine::resolve_tree(aln, "AUTO", false, false, "", &source_desc, true);
        CHECK(tree.num_leaves() == 10);
        CHECK(source_desc.find("TreeEngine") != std::string::npos);
    }

    SUBCASE("Explicit tree file") {
        std::string tree_path = HYPHY3_SOURCE_DIR "/benchmarks/data/cd2.nwk";
        std::string source_desc;
        Tree tree = TreeEngine::resolve_tree(aln, tree_path, false, false, "", &source_desc, true);
        CHECK(tree.num_leaves() == 10);
        CHECK(source_desc == tree_path);
    }

    SUBCASE("Tree-free default fallback (no tree specified, runs TreeEngine)") {
        std::string source_desc;
        Tree tree = TreeEngine::resolve_tree(aln, "", false, false, "", &source_desc, true);
        CHECK(tree.num_leaves() == 10);
        CHECK(source_desc.find("TreeEngine") != std::string::npos);
    }

    SUBCASE("Tree-free with force_ohmic runs OhmicNetwork") {
        std::string source_desc;
        Tree tree = TreeEngine::resolve_tree(aln, "", true, false, "", &source_desc, true);
        CHECK(tree.num_leaves() == 10);
        CHECK(source_desc.find("OhmicNetwork") != std::string::npos);
    }
}

TEST_CASE("TreeEngine: Ohmic Divide-and-Conquer Subtree Optimization") {
    std::string aln_path = HYPHY3_SOURCE_DIR "/benchmarks/data/adh.fna";
    auto aln = Alignment::from_fasta(aln_path);

    TreeEngineSettings settings;
    settings.quiet = true;
    settings.enable_dnc = true;
    settings.dnc_min_taxa = 16;
    settings.dnc_target_island_size = 14;

    TreeEngine engine(settings);
    auto res = engine.run(aln);

    CHECK(res.dnc_used == true);
    CHECK(res.num_islands >= 2);
    CHECK(res.island_sizes.size() == res.num_islands);

    size_t total_island_taxa = 0;
    for (size_t sz : res.island_sizes) {
        CHECK(sz >= 4);
        CHECK(sz <= 14);
        total_island_taxa += sz;
    }
    CHECK(total_island_taxa <= 23);
    CHECK(total_island_taxa >= 18); // Bulk of tree covered by islands

    auto json_out = engine.to_json(res);
    CHECK(json_out.contains("divide_and_conquer"));
    CHECK(json_out["divide_and_conquer"]["enabled"] == true);
    CHECK(json_out["divide_and_conquer"]["num_islands"].get<size_t>() == res.num_islands);
    CHECK(!res.optimized_tree.to_newick().empty());
    CHECK(res.optimized_tree.num_leaves() == 23);
}

TEST_CASE("TreeEngine: Option C - Default GTR vs Explicit MG94 Codon Optimization") {
    std::string aln_path = HYPHY3_SOURCE_DIR "/benchmarks/data/cd2.fna";
    auto aln = Alignment::load(aln_path);

    SUBCASE("Default settings (Option C: pure GTR, no escalation)") {
        TreeEngineSettings settings;
        settings.quiet = true;
        TreeEngine engine(settings);
        auto res = engine.run(aln);

        CHECK(res.active_model == "GTR");
        CHECK(res.final_log_l == res.optimized_gtr_log_l);
    }

    SUBCASE("Explicit MG94 model (Option C: full MG94 branch length & omega optimization)") {
        TreeEngineSettings settings;
        settings.quiet = true;
        settings.model = TreeSearchModel::MG94;
        settings.lbfgs_max_iterations = 20;

        TreeEngine engine(settings);
        auto res = engine.run(aln);

        CHECK(res.active_model == "MG94");
        CHECK(res.mg94_omega > 0.0);
        CHECK(res.mg94_omega < 10.0);
        CHECK(res.final_log_l > -4000.0);
    }
}
