#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "hyphy/core/types.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/core/rate_matrix.hpp"
#include "hyphy/core/likelihood.hpp"

using namespace hyphy::core;

TEST_CASE("Alignment: Reading and frequency calculation on cd2_reduced") {
    std::string aln_path = HYPHY3_SOURCE_DIR "/benchmarks/data/cd2_reduced.fna";
    auto aln = Alignment::from_fasta(aln_path);

    CHECK(aln.num_taxa == 10);
    CHECK(aln.num_codons == 6);
    CHECK(aln.num_nucleotides == 18);

    // Non-gap empirical frequencies:
    // A: 0.456140, C: 0.114035, G: 0.254386, T: 0.175439
    CHECK(doctest::Approx(aln.nuc_frequencies(0)).epsilon(1e-4) == 0.456140);
    CHECK(doctest::Approx(aln.nuc_frequencies(1)).epsilon(1e-4) == 0.114035);
    CHECK(doctest::Approx(aln.nuc_frequencies(2)).epsilon(1e-4) == 0.254386);
    CHECK(doctest::Approx(aln.nuc_frequencies(3)).epsilon(1e-4) == 0.175439);

    // Verify pattern compression
    CHECK(aln.patterns.size() > 0);
    CHECK(aln.patterns.size() <= 6);
    size_t total_weight = 0;
    for (const auto& p : aln.patterns) {
        total_weight += p.weight;
    }
    CHECK(total_weight == 6);
}

TEST_CASE("Tree: Newick parsing and post-order traversal on cd2_reduced") {
    std::string nwk = "((((Pig:0.147969,Cow:0.213430):0.085099,Horse:0.165787,Cat:0.264806):0.058611,((RhMonkey:0.002015,Baboon:0.003108):0.022733,(Human:0.004349,Chimp:0.000799):0.011873):0.101856):0.340802,Rat:0.050958,Mouse:0.097950);";
    auto tree = Tree::from_newick(nwk);

    CHECK(tree.num_leaves() == 10);
    CHECK(tree.root_id != INVALID_INDEX);
    CHECK(tree.post_order.size() == tree.num_nodes());

    // Root should be the last element in post-order traversal
    CHECK(tree.post_order.back() == tree.root_id);

    // Leaves should be recognized
    CHECK(tree.leaf_name_to_id.count("Human") == 1);
    CHECK(tree.leaf_name_to_id.count("Chimp") == 1);
    CHECK(tree.leaf_name_to_id.count("Pig") == 1);
}

TEST_CASE("RateMatrix: GTR and MG94 generator properties") {
    Vector4 nuc_freqs;
    nuc_freqs << 0.3, 0.2, 0.2, 0.3;

    GTRParameters gtr_p;
    gtr_p.theta_CT = 2.0;

    GTRMatrix gtr;
    gtr.update(gtr_p, nuc_freqs);

    // Check row sum to 0
    for (int i = 0; i < 4; ++i) {
        CHECK(std::abs(gtr.Q.row(i).sum()) < 1e-10);
    }

    // Check transition matrix row sum to 1
    Matrix4 P = gtr.transition_matrix(0.1);
    for (int i = 0; i < 4; ++i) {
        CHECK(doctest::Approx(P.row(i).sum()).epsilon(1e-6) == 1.0);
    }

    // Check MG94
    Vector codon_freqs = Vector::Constant(61, 1.0 / 61.0);
    Eigen::Matrix<Scalar, 3, 4> pos_nuc = Eigen::Matrix<Scalar, 3, 4>::Constant(0.25);
    MG94Parameters mg_p;
    mg_p.alpha = 1.0;
    mg_p.beta = 0.5;

    MG94Matrix mg;
    mg.update(mg_p, pos_nuc, codon_freqs);

    for (size_t i = 0; i < NUM_SENSE_CODONS; ++i) {
        CHECK(std::abs(mg.Q.row(i).sum()) < 1e-9);
    }

    Matrix61 P_codon = mg.transition_matrix(0.05);
    for (size_t i = 0; i < NUM_SENSE_CODONS; ++i) {
        CHECK(doctest::Approx(P_codon.row(i).sum()).epsilon(1e-6) == 1.0);
    }
}

TEST_CASE("Likelihood: Point likelihood comparison against HyPhy 2.5 ground truth") {
    std::string aln_path = HYPHY3_SOURCE_DIR "/benchmarks/data/cd2_reduced.fna";
    std::string nwk = "((((Pig:0.147969,Cow:0.213430):0.085099,Horse:0.165787,Cat:0.264806):0.058611,((RhMonkey:0.002015,Baboon:0.003108):0.022733,(Human:0.004349,Chimp:0.000799):0.011873):0.101856):0.340802,Rat:0.050958,Mouse:0.097950);";

    auto aln = Alignment::from_fasta(aln_path);
    auto tree = Tree::from_newick(nwk);

    // MG94 Model with F3x4 frequencies
    MG94Parameters mg_p;
    mg_p.alpha = 1.0;
    mg_p.beta = 0.5;
    mg_p.theta_AC = 1.0;
    mg_p.theta_AT = 1.0;
    mg_p.theta_CG = 1.0;
    mg_p.theta_CT = 2.0;
    mg_p.theta_GT = 1.0;

    MG94Matrix mg94;
    mg94.update(mg_p, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4);

    Scalar mg94_log_l = LikelihoodEngine::compute_mg94_log_likelihood(tree, aln, mg94);
    std::cout << "[hyphy-next] Computed MG94 LogL: " << mg94_log_l << "\n";
    std::cout << "[HyPhy 2.5 ] Ground Truth MG94: -92.7774489\n";

    // Exact match assertion against HyPhy 2.5 ground truth:
    CHECK(doctest::Approx(mg94_log_l).epsilon(1e-4) == -92.77745);
}

#include "hyphy/opt/optimizer.hpp"
using namespace hyphy::opt;

TEST_CASE("Optimizer: L-BFGS-B optimization of omega on cd2_reduced") {
    std::string aln_path = HYPHY3_SOURCE_DIR "/benchmarks/data/cd2_reduced.fna";
    std::string nwk = "((((Pig:0.147969,Cow:0.213430):0.085099,Horse:0.165787,Cat:0.264806):0.058611,((RhMonkey:0.002015,Baboon:0.003108):0.022733,(Human:0.004349,Chimp:0.000799):0.011873):0.101856):0.340802,Rat:0.050958,Mouse:0.097950);";

    auto aln = Alignment::from_fasta(aln_path);
    auto tree = Tree::from_newick(nwk);

    MG94Parameters base_p;
    base_p.theta_CT = 2.0;

    MG94Fitter fitter(tree, aln);
    auto res = fitter.fit_omega(1.0, base_p);

    std::cout << "[hyphy-next] Optimized omega: " << res.x_opt(0) 
              << ", LogL: " << res.log_likelihood 
              << " in " << res.iterations << " iterations\n";
    std::cout << "[HyPhy 2.5 ] Ground Truth omega: 5.74349\n";

    CHECK(res.converged == true);
    CHECK(res.x_opt(0) > 0.1);
    CHECK(res.log_likelihood > -92.77745); // Strictly improved from initial -92.777
}

TEST_CASE("LikelihoodEngine: Exact Spectral Fréchet Joint Gradients vs Finite Differences") {
    std::string aln_path = HYPHY3_SOURCE_DIR "/benchmarks/data/cd2_reduced.fna";
    std::string nwk = "((((Pig:0.147969,Cow:0.213430):0.085099,Horse:0.165787,Cat:0.264806):0.058611,((RhMonkey:0.002015,Baboon:0.003108):0.022733,(Human:0.004349,Chimp:0.000799):0.011873):0.101856):0.340802,Rat:0.050958,Mouse:0.097950);";

    auto aln = Alignment::from_fasta(aln_path);
    auto tree = Tree::from_newick(nwk);

    MG94Parameters params;
    params.alpha = 1.0;
    params.beta = 0.45;
    params.theta_AC = 0.8;
    params.theta_AT = 1.2;
    params.theta_CG = 0.6;
    params.theta_CT = 1.8;
    params.theta_GT = 0.9;

    size_t num_nodes = tree.num_nodes();
    std::vector<Scalar> bl(num_nodes, 0.0);
    for (const auto& node : tree.nodes) {
        if (node.id != tree.root_id) bl[node.id] = std::max(1e-4, node.branch_length);
    }
    auto leaf_map = LikelihoodEngine::build_leaf_to_taxon_map(tree, aln);

    auto ana = LikelihoodEngine::compute_joint_frechet_gradients(tree, aln, params, bl, &leaf_map);

    auto eval_lnl = [&](const MG94Parameters& p, const std::vector<Scalar>& b) -> Scalar {
        MG94Matrix mg;
        mg.update(p, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *(aln.code ? aln.code : GeneticCode::universal()));
        Scalar conv = (mg.scale_factor > 1e-12) ? (3.0 / mg.scale_factor) : 1.0;
        Tree t_eval = tree;
        for (const auto& node : t_eval.nodes) {
            if (node.id != t_eval.root_id) {
                t_eval.nodes[node.id].branch_length = b[node.id] * conv;
            }
        }
        return LikelihoodEngine::compute_mg94_log_likelihood(t_eval, aln, mg, &leaf_map);
    };

    Scalar eps = 1e-6;

    // Test Global Parameters Gradients
    auto test_param = [&](Scalar ana_grad, auto mutate) {
        MG94Parameters p_plus = params, p_minus = params;
        mutate(p_plus, eps);
        mutate(p_minus, -eps);
        Scalar fd = (eval_lnl(p_plus, bl) - eval_lnl(p_minus, bl)) / (2.0 * eps);
        CHECK(std::abs(ana_grad - fd) / (std::abs(ana_grad) + 1.0) < 1e-4);
    };

    test_param(ana.grad_omega, [](MG94Parameters& p, Scalar d) { p.beta += d; });
    test_param(ana.grad_theta_AC, [](MG94Parameters& p, Scalar d) { p.theta_AC += d; });
    test_param(ana.grad_theta_AT, [](MG94Parameters& p, Scalar d) { p.theta_AT += d; });
    test_param(ana.grad_theta_CG, [](MG94Parameters& p, Scalar d) { p.theta_CG += d; });
    test_param(ana.grad_theta_CT, [](MG94Parameters& p, Scalar d) { p.theta_CT += d; });
    test_param(ana.grad_theta_GT, [](MG94Parameters& p, Scalar d) { p.theta_GT += d; });

    // Test Branch Length Gradients
    for (const auto& node : tree.nodes) {
        if (node.id != tree.root_id) {
            int32_t bid = node.id;
            auto bl_plus = bl, bl_minus = bl;
            bl_plus[bid] += eps;
            bl_minus[bid] -= eps;
            Scalar fd = (eval_lnl(params, bl_plus) - eval_lnl(params, bl_minus)) / (2.0 * eps);
            CHECK(std::abs(ana.grad_branches[bid] - fd) / (std::abs(ana.grad_branches[bid]) + 1.0) < 1e-4);
        }
    }

    // Verify Unified Joint Fit full model converges and improves likelihood
    MG94Fitter fitter(tree, aln);
    auto fit_res = fitter.fit_full_model(1.0, params);
    CHECK(fit_res.converged);
    CHECK(fit_res.log_likelihood > ana.log_likelihood);
    CHECK(fit_res.params.beta > 0.0);
}

TEST_CASE("Issue #2: NEXUS FORMAT NOLABELS with TAXLABELS block") {
    std::string tmp_nex = "test_nolabels_issue2.nex";
    {
        std::ofstream out(tmp_nex);
        out << "#NEXUS\n"
            << "BEGIN TAXA;\n"
            << "  DIMENSIONS NTAX=3;\n"
            << "  TAXLABELS 'Taxon A' Taxon_B \"Taxon C\";\n"
            << "END;\n"
            << "BEGIN CHARACTERS;\n"
            << "  DIMENSIONS NCHAR=9;\n"
            << "  FORMAT DATATYPE=DNA NOLABELS;\n"
            << "  MATRIX\n"
            << "    ATGAAACCC\n"
            << "    ATGAAACCA\n"
            << "    ATGAAGCCC\n"
            << "  ;\n"
            << "END;\n";
    }
    auto aln = Alignment::load(tmp_nex);
    std::remove(tmp_nex.c_str());

    REQUIRE(aln.num_taxa == 3);
    CHECK(aln.num_codons == 3);
    CHECK(aln.taxon_names[0] == "Taxon A");
    CHECK(aln.taxon_names[1] == "Taxon_B");
    CHECK(aln.taxon_names[2] == "Taxon C");
}

TEST_CASE("Issue #3: Tree::select_branches supports model tags, keywords, branch lists, and throws on empty match") {
    std::string nwk = "((A{FG}:0.1,B{FG}:0.2){FG}:0.3,(C{BG}:0.1,D:0.2):0.3);";
    auto tree = Tree::from_newick(nwk);

    // Tag selection (case-insensitive)
    auto fg_sel = tree.select_branches("fg");
    CHECK(fg_sel.num_selected == 3);

    // Unlabeled selection
    auto unl_sel = tree.select_branches("Unlabeled");
    CHECK(unl_sel.num_selected == 2); // D and parent of (C,D)

    // Comma-separated branch names + tags
    auto combo_sel = tree.select_branches("BG,D");
    CHECK(combo_sel.num_selected == 2); // C{BG} and D

    // Zero matches should throw std::runtime_error
    CHECK_THROWS_AS(tree.select_branches("NonExistentGroup"), std::runtime_error);
}
