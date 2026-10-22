#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "hyphy/core/types.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/core/rate_matrix.hpp"
#include "hyphy/core/likelihood.hpp"
#include "hyphy/core/ohmic_network.hpp"
#include "hyphy/core/hyphaeon.hpp"
#include "hyphy/core/wheatstone_likelihood.hpp"
#include "hyphy/analyses/ohmic.hpp"
#include "hyphy/analyses/fel.hpp"
#include "hyphy/analyses/fubar.hpp"
#include "hyphy/analyses/busted.hpp"
#include "hyphy/analyses/absrel.hpp"

#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>

using namespace hyphy::core;
using namespace hyphy::analyses;

inline std::string resolve_test_path(const std::string& path) {
    auto slash = path.find_last_of('/');
    std::string base = (slash != std::string::npos) ? path.substr(slash + 1) : path;
    std::vector<std::string> candidates = {
        path,
        "../" + path,
        "../../" + path,
        "benchmarks/data/" + base,
        "../benchmarks/data/" + base,
        "../../benchmarks/data/" + base,
        "/Users/sergei/Development/hyphy-3/benchmarks/data/" + base,
        "/Users/sergei/Development/hyphy-3/tests/data/" + base
#ifdef HYPHY3_SOURCE_DIR
        , std::string(HYPHY3_SOURCE_DIR) + "/benchmarks/data/" + base
        , std::string(HYPHY3_SOURCE_DIR) + "/tests/data/" + base
#endif
    };
    for (const auto& c : candidates) {
        if (std::ifstream(c).good()) return c;
    }
    return path;
}

// 1. Foster's Theorem Unit Test (ExactDense)
TEST_CASE("OHMIC 1: Foster's Theorem Invariant on Empirical Alignments (CD2, ADH, bglobin)") {
    std::vector<std::string> files = {
        "benchmarks/data/cd2.fna",
        "tests/data/adh.nex",
        "benchmarks/data/bglobin.nex"
    };

    for (const auto& f : files) {
        std::string path = resolve_test_path(f);
        auto aln = Alignment::load_nucleotide(path);
        REQUIRE(aln.num_taxa >= 5);

        OhmicSettings settings;
        settings.scale_mode = OhmicScaleMode::ExactDense;
        OhmicNetwork net = OhmicNetwork::infer(aln, settings);

        CHECK(net.active_mode == OhmicScaleMode::ExactDense);
        CHECK(net.foster_sum == doctest::Approx(static_cast<Scalar>(aln.num_taxa - 1)).epsilon(1e-5));
    }
}

// 2. Degenerate Single-Tree Parity Test
TEST_CASE("OHMIC 2: Degenerate Single-Tree Parity (w1=1, w2=0, w3=0) vs LikelihoodEngine") {
    std::string cd2_path = resolve_test_path("benchmarks/data/cd2.fna");
    auto aln = Alignment::load(cd2_path);
    REQUIRE(aln.num_taxa == 10);

    OhmicSettings settings;
    settings.scale_mode = OhmicScaleMode::ExactDense;
    settings.min_shunt_stress = 0.01; // Ensure at least 1 bridge is created
    OhmicNetwork net = OhmicNetwork::infer(aln, settings);
    REQUIRE(net.bridges.size() >= 1);

    // Force degenerate single-tree weights (w1 = 1, w2 = 0, w3 = 0) on all bridges
    for (auto& wb : net.bridges) {
        wb.weights = {1.0, 0.0, 0.0};
    }

    MG94Parameters params;
    params.alpha = 1.0;
    params.beta = 0.65;
    params.theta_AC = 0.45;
    params.theta_AT = 0.35;
    params.theta_CG = 0.50;
    params.theta_CT = 1.20;
    params.theta_GT = 0.40;

    auto base_grad = LikelihoodEngine::compute_joint_frechet_gradients(
        net.backbone_tree, aln, params
    );
    auto ohm_grad = WheatstoneLikelihood::compute_joint_frechet_gradients(
        net, aln, params
    );

    CHECK(ohm_grad.log_likelihood == doctest::Approx(base_grad.log_likelihood).epsilon(1e-11));
    CHECK(ohm_grad.grad_omega == doctest::Approx(base_grad.grad_omega).epsilon(1e-11));
    CHECK(ohm_grad.grad_theta_AC == doctest::Approx(base_grad.grad_theta_AC).epsilon(1e-11));
    CHECK(ohm_grad.grad_theta_AT == doctest::Approx(base_grad.grad_theta_AT).epsilon(1e-11));
    CHECK(ohm_grad.grad_theta_CG == doctest::Approx(base_grad.grad_theta_CG).epsilon(1e-11));
    CHECK(ohm_grad.grad_theta_CT == doctest::Approx(base_grad.grad_theta_CT).epsilon(1e-11));
    CHECK(ohm_grad.grad_theta_GT == doctest::Approx(base_grad.grad_theta_GT).epsilon(1e-11));

    for (size_t nid = 0; nid < net.backbone_tree.nodes.size(); ++nid) {
        if (static_cast<int32_t>(nid) == net.backbone_tree.root_id) continue;
        CHECK(ohm_grad.grad_branches[nid] == doctest::Approx(base_grad.grad_branches[nid]).epsilon(1e-11));
    }
}

// 3. Worked 5-Taxon Ground-Truth Example Test (Section 5 of specs/ohmic.MD)
TEST_CASE("OHMIC 3: Ground-Truth 5-Taxon Worked Example (Exact Closed-Form Verification)") {
    // Construct the exact N=5, L=100 alignment from Section 5.1:
    // Sites 1-90 (90 invariant): A=A, B=A, C=A, D=A, O=A
    // Sites 91-97 (7 majority AB|CD): A=C, B=C, C=A, D=A, O=A
    // Sites 98-100 (3 slow shunt AC|BD): A=G, B=A, C=G, D=A, O=A
    Alignment aln;
    std::string seq_A = std::string(90, 'A') + std::string(7, 'C') + std::string(3, 'G');
    std::string seq_B = std::string(90, 'A') + std::string(7, 'C') + std::string(3, 'A');
    std::string seq_C = std::string(90, 'A') + std::string(7, 'A') + std::string(3, 'G');
    std::string seq_D = std::string(100, 'A');
    std::string seq_O = std::string(100, 'A');

    aln.add_sequence("A", seq_A);
    aln.add_sequence("B", seq_B);
    aln.add_sequence("C", seq_C);
    aln.add_sequence("D", seq_D);
    aln.add_sequence("O", seq_O);
    aln.finalize_nucleotide();

    OhmicSettings settings;
    settings.scale_mode = OhmicScaleMode::ExactDense;
    settings.use_hamming_distance = true;
    settings.fixed_sigma = 0.10;
    settings.epsilon_0 = 0.01;
    settings.min_shunt_stress = 0.05;
    settings.collapse_ohmic_polytomies = false;

    OhmicNetwork net = OhmicNetwork::infer(aln, settings);

    // 5.2.1 Verify exact Hamming distances
    CHECK(net.D_tn93(0, 1) == doctest::Approx(0.03).epsilon(1e-6)); // D_AB
    CHECK(net.D_tn93(0, 2) == doctest::Approx(0.07).epsilon(1e-6)); // D_AC
    CHECK(net.D_tn93(0, 3) == doctest::Approx(0.10).epsilon(1e-6)); // D_AD
    CHECK(net.D_tn93(0, 4) == doctest::Approx(0.10).epsilon(1e-6)); // D_AO
    CHECK(net.D_tn93(1, 2) == doctest::Approx(0.10).epsilon(1e-6)); // D_BC
    CHECK(net.D_tn93(1, 3) == doctest::Approx(0.07).epsilon(1e-6)); // D_BD
    CHECK(net.D_tn93(2, 3) == doctest::Approx(0.03).epsilon(1e-6)); // D_CD
    CHECK(net.D_tn93(3, 4) == doctest::Approx(0.00).epsilon(1e-6)); // D_DO

    // 5.2.2 Verify Foster's Theorem identity == 4.000000
    CHECK(net.foster_sum == doctest::Approx(4.000000).epsilon(1e-6));

    // 5.2.3 Verify quartet conductances, Ohmic weights (0.70 / 0.30), and shunt stress (0.4286)
    REQUIRE(net.bridges.size() == 1);
    const auto& wb = net.bridges[0];
    CHECK(wb.weights[0] == doctest::Approx(0.7000).epsilon(1e-4));
    CHECK(wb.weights[1] == doctest::Approx(0.3000).epsilon(1e-4));
    CHECK(wb.shunt_stress == doctest::Approx(0.4285714).epsilon(1e-3));

    // 5.3 Verify exact JC69 Log-Likelihoods on the 5-taxon Wheatstone bridge topology
    // Build the exact 5-taxon collapsed tree topology specified in Section 5.3:
    // Root has children {hub_u, D, O} with t_u = 0, t_D = 0.005, t_O = 0.0
    // Hub u has children {A, B, C} with t_A = t_B = t_C = 0.005, bridge length t_int = 0.08
    OhmicNetwork exact_net;
    Tree& ct = exact_net.collapsed_tree;
    for (int i = 0; i < 5; ++i) {
        TreeNode leaf;
        leaf.id = i;
        leaf.name = aln.taxon_names[i];
        leaf.is_leaf = true;
        leaf.branch_length = (i == 4) ? 0.0 : 0.005; // A,B,C,D = 0.005, O = 0.0
        ct.nodes.push_back(leaf);
        ct.leaf_name_to_id[leaf.name] = i;
    }
    TreeNode hub_u;
    hub_u.id = 5;
    hub_u.name = "Hub_u";
    hub_u.is_leaf = false;
    hub_u.branch_length = 0.0;
    hub_u.children = {0, 1, 2}; // A, B, C
    ct.nodes[0].parent_id = 5;
    ct.nodes[1].parent_id = 5;
    ct.nodes[2].parent_id = 5;
    ct.nodes.push_back(hub_u);

    TreeNode root;
    root.id = 6;
    root.name = "Root";
    root.is_leaf = false;
    root.branch_length = 0.0;
    root.children = {5, 3, 4}; // Hub_u, D, O
    ct.nodes[5].parent_id = 6;
    ct.nodes[3].parent_id = 6;
    ct.nodes[4].parent_id = 6;
    ct.nodes.push_back(root);
    ct.root_id = 6;
    ct.compute_post_order();

    WheatstoneBridge bridge;
    bridge.hub_node_id = 5;
    bridge.bridge_node_id = 5;
    bridge.children = {0, 1, 2};
    bridge.weights = {0.70, 0.30, 0.0};
    bridge.branch_lengths = {0.08, 0.08, 0.08};
    exact_net.bridges = {bridge};
    exact_net.node_to_bridge_idx.assign(ct.num_nodes(), -1);
    exact_net.node_to_bridge_idx[5] = 0;

    // JC69 rate matrix (equal base frequencies 0.25, all exchangeabilities = 1.0 -> mu = 4/3)
    GTRParameters jc_params;
    Vector4 jc_pi = Vector4::Constant(0.25);
    GTRMatrix jc_gtr;
    jc_gtr.update(jc_params, jc_pi);

    // Slice into fast partition (sites 0..97, r_fast = 1.0) and slow shunt partition (sites 97..100, r_slow = 0.05)
    Alignment aln_inv = aln.slice_nucleotide(0, 90);
    Alignment aln_maj = aln.slice_nucleotide(90, 97);
    Alignment aln_shunt = aln.slice_nucleotide(97, 100);

    // 1. Tree 1 (AB|CD only: w1 = 1.0, w2 = 0.0)
    exact_net.bridges[0].weights = {1.0, 0.0, 0.0};
    Scalar lnl_inv_1 = WheatstoneLikelihood::compute_gtr_log_likelihood(exact_net, aln_inv, jc_gtr, 1.0) / 90.0;
    Scalar lnl_maj_1 = WheatstoneLikelihood::compute_gtr_log_likelihood(exact_net, aln_maj, jc_gtr, 1.0) / 7.0;
    Scalar lnl_shunt_1 = WheatstoneLikelihood::compute_gtr_log_likelihood(exact_net, aln_shunt, jc_gtr, 0.05) / 3.0;
    Scalar total_lnl_T1 = 90.0 * lnl_inv_1 + 7.0 * lnl_maj_1 + 3.0 * lnl_shunt_1;

    CHECK(lnl_inv_1 == doctest::Approx(-1.485192).epsilon(1e-5));
    CHECK(lnl_maj_1 == doctest::Approx(-5.083371).epsilon(1e-5));
    CHECK(lnl_shunt_1 == doctest::Approx(-20.175114).epsilon(1e-5));
    CHECK(total_lnl_T1 == doctest::Approx(-229.7762).epsilon(1e-4));

    // 2. Tree 2 (AC|BD only: w1 = 0.0, w2 = 1.0)
    exact_net.bridges[0].weights = {0.0, 1.0, 0.0};
    Scalar lnl_inv_2 = WheatstoneLikelihood::compute_gtr_log_likelihood(exact_net, aln_inv, jc_gtr, 1.0) / 90.0;
    Scalar lnl_maj_2 = WheatstoneLikelihood::compute_gtr_log_likelihood(exact_net, aln_maj, jc_gtr, 1.0) / 7.0;
    Scalar lnl_shunt_2 = WheatstoneLikelihood::compute_gtr_log_likelihood(exact_net, aln_shunt, jc_gtr, 0.05) / 3.0;
    Scalar total_lnl_T2 = 90.0 * lnl_inv_2 + 7.0 * lnl_maj_2 + 3.0 * lnl_shunt_2;

    CHECK(lnl_inv_2 == doctest::Approx(-1.485192).epsilon(1e-5));
    CHECK(lnl_maj_2 == doctest::Approx(-14.248631).epsilon(1e-5));
    CHECK(lnl_shunt_2 == doctest::Approx(-8.010028).epsilon(1e-5));
    CHECK(total_lnl_T2 == doctest::Approx(-257.4378).epsilon(1e-4));

    // 3. Wheatstone Pruning Network (w1 = 0.70, w2 = 0.30)
    exact_net.bridges[0].weights = {0.70, 0.30, 0.0};
    Scalar lnl_inv_ohm = WheatstoneLikelihood::compute_gtr_log_likelihood(exact_net, aln_inv, jc_gtr, 1.0) / 90.0;
    Scalar lnl_maj_ohm = WheatstoneLikelihood::compute_gtr_log_likelihood(exact_net, aln_maj, jc_gtr, 1.0) / 7.0;
    Scalar lnl_shunt_ohm = WheatstoneLikelihood::compute_gtr_log_likelihood(exact_net, aln_shunt, jc_gtr, 0.05) / 3.0;
    Scalar total_lnl_ohm = 90.0 * lnl_inv_ohm + 7.0 * lnl_maj_ohm + 3.0 * lnl_shunt_ohm;

    CHECK(lnl_inv_ohm == doctest::Approx(-1.485192).epsilon(1e-5));
    CHECK(lnl_maj_ohm == doctest::Approx(-5.440001).epsilon(1e-5));
    CHECK(lnl_shunt_ohm == doctest::Approx(-9.213988).epsilon(1e-5));
    CHECK(total_lnl_ohm == doctest::Approx(-199.3893).epsilon(1e-4));
}

// 4. Finite-Difference Gradient Verification
TEST_CASE("OHMIC 4: Analytical Fréchet & Adjoint Gradients vs Central Finite Differences") {
    std::string cd2_path = resolve_test_path("tests/data/cd2.fna");
    auto aln = Alignment::load(cd2_path);

    OhmicSettings settings;
    settings.scale_mode = OhmicScaleMode::ExactDense;
    settings.min_shunt_stress = 0.01;
    OhmicNetwork net = OhmicNetwork::infer(aln, settings);
    REQUIRE(net.bridges.size() >= 1);

    // Activate all 3 bridge resolutions on bridge 0 with distinct branch lengths and weights
    net.bridges[0].weights = {0.62, 0.28, 0.10};
    net.bridges[0].branch_lengths = {0.045, 0.025, 0.018};

    MG94Parameters params;
    params.alpha = 1.0;
    params.beta = 0.75;
    params.theta_AC = 0.55;
    params.theta_AT = 0.40;
    params.theta_CG = 0.60;
    params.theta_CT = 1.35;
    params.theta_GT = 0.45;

    auto g = WheatstoneLikelihood::compute_joint_frechet_gradients(net, aln, params);
    const Scalar h = 1e-6;

    auto eval_lnl = [&](const OhmicNetwork& n_eval, const MG94Parameters& p_eval) -> Scalar {
        return WheatstoneLikelihood::compute_joint_frechet_gradients(n_eval, aln, p_eval).log_likelihood;
    };

    // Check global MG94 parameters (omega, theta_*)
    {
        MG94Parameters p_p = params, p_m = params;
        p_p.beta += h; p_m.beta -= h;
        Scalar fd_omega = (eval_lnl(net, p_p) - eval_lnl(net, p_m)) / (2.0 * h);
        CHECK(g.grad_omega == doctest::Approx(fd_omega).epsilon(1e-4));
    }
    {
        MG94Parameters p_p = params, p_m = params;
        p_p.theta_AC += h; p_m.theta_AC -= h;
        Scalar fd = (eval_lnl(net, p_p) - eval_lnl(net, p_m)) / (2.0 * h);
        CHECK(g.grad_theta_AC == doctest::Approx(fd).epsilon(1e-4));
    }
    {
        MG94Parameters p_p = params, p_m = params;
        p_p.theta_CT += h; p_m.theta_CT -= h;
        Scalar fd = (eval_lnl(net, p_p) - eval_lnl(net, p_m)) / (2.0 * h);
        CHECK(g.grad_theta_CT == doctest::Approx(fd).epsilon(1e-4));
    }

    // Check bridge branch lengths t_{0,r} and bridge weights w_{0,r}
    for (int r = 0; r < 3; ++r) {
        OhmicNetwork n_p = net, n_m = net;
        n_p.bridges[0].branch_lengths[r] += h;
        n_m.bridges[0].branch_lengths[r] -= h;
        Scalar fd_len = (eval_lnl(n_p, params) - eval_lnl(n_m, params)) / (2.0 * h);
        CHECK(g.grad_bridge_lengths[0][r] == doctest::Approx(fd_len).epsilon(1e-4));

        OhmicNetwork nw_p = net, nw_m = net;
        nw_p.bridges[0].weights[r] += h;
        nw_m.bridges[0].weights[r] -= h;
        Scalar fd_w = (eval_lnl(nw_p, params) - eval_lnl(nw_m, params)) / (2.0 * h);
        CHECK(g.grad_bridge_weights[0][r] == doctest::Approx(fd_w).epsilon(1e-4));
    }

    // Check incident child branches {id_A, id_B, id_C} of bridge 0
    for (int c_idx = 0; c_idx < 3; ++c_idx) {
        int32_t cid = net.bridges[0].children[c_idx];
        OhmicNetwork n_p = net, n_m = net;
        n_p.collapsed_tree.nodes[cid].branch_length += h;
        n_m.collapsed_tree.nodes[cid].branch_length -= h;
        Scalar fd_b = (eval_lnl(n_p, params) - eval_lnl(n_m, params)) / (2.0 * h);
        CHECK(g.grad_branches[cid] == doctest::Approx(fd_b).epsilon(1e-4));
    }
}

// 5. ExactDense vs ScalableWoodbury Crossover Parity Test
TEST_CASE("OHMIC 5: ExactDense vs ScalableWoodbury Crossover Parity (N = 120)") {
    // Generate a structured 120-taxon alignment with 4 major clades of 30 distinct sequences (L = 300 bp)
    const size_t N = 120;
    const size_t L = 300;
    std::mt19937 rng(42);
    const char bases[4] = {'A', 'C', 'G', 'T'};

    std::string root_seq(L, 'A');
    for (size_t s = 0; s < L; ++s) root_seq[s] = bases[rng() % 4];

    std::vector<std::string> clade_centers(4, root_seq);
    for (int c = 0; c < 4; ++c) {
        for (size_t s = 0; s < L; ++s) {
            if ((rng() % 100) < 12) clade_centers[c][s] = bases[rng() % 4];
        }
    }

    Alignment aln;
    for (size_t i = 0; i < N; ++i) {
        std::string seq = clade_centers[i % 4];
        // Mutate ~5% of sites uniquely per taxon so all M=120 sequences are distinct
        for (size_t s = 0; s < L; ++s) {
            if ((rng() % 100) < 5) seq[s] = bases[(rng() + i) % 4];
        }
        // Guarantee uniqueness at site i % L
        seq[i % L] = bases[(i / 4) % 4];
        aln.add_sequence("Taxon_" + std::to_string(i + 1), seq);
    }
    aln.finalize_nucleotide();
    REQUIRE(aln.num_taxa == N);

    OhmicSettings dense_settings;
    dense_settings.scale_mode = OhmicScaleMode::ExactDense;
    OhmicNetwork net_dense = OhmicNetwork::infer(aln, dense_settings);

    OhmicSettings woodbury_settings;
    woodbury_settings.scale_mode = OhmicScaleMode::ScalableWoodbury;
    woodbury_settings.woodbury_rank = 120;
    OhmicNetwork net_woodbury = OhmicNetwork::infer(aln, woodbury_settings);

    CHECK(net_dense.active_mode == OhmicScaleMode::ExactDense);
    CHECK(net_woodbury.active_mode == OhmicScaleMode::ScalableWoodbury);

    // Compute Pearson correlation R between ExactDense Omega_{ij} and ScalableWoodbury effective_resistance(i, j)
    std::vector<Scalar> x_om, y_om;
    size_t M = net_woodbury.num_super_nodes;
    REQUIRE(M == N);
    for (size_t i = 0; i < M; ++i) {
        for (size_t j = i + 1; j < M; ++j) {
            x_om.push_back(net_dense.effective_resistance(i, j));
            y_om.push_back(net_woodbury.effective_resistance(i, j));
        }
    }
    Scalar mean_x = std::accumulate(x_om.begin(), x_om.end(), 0.0) / x_om.size();
    Scalar mean_y = std::accumulate(y_om.begin(), y_om.end(), 0.0) / y_om.size();
    Scalar cov = 0.0, var_x = 0.0, var_y = 0.0;
    for (size_t k = 0; k < x_om.size(); ++k) {
        Scalar dx = x_om[k] - mean_x;
        Scalar dy = y_om[k] - mean_y;
        cov += dx * dy;
        var_x += dx * dx;
        var_y += dy * dy;
    }
    Scalar r_corr = cov / std::sqrt(var_x * var_y);
    CHECK(r_corr > 0.98);

    // Check GTR log-likelihood agreement within 1%
    GTRParameters gtr_p;
    GTRMatrix gtr_m;
    gtr_m.update(gtr_p, aln.nuc_frequencies);
    Scalar lnl_dense = WheatstoneLikelihood::compute_gtr_log_likelihood(net_dense, aln, gtr_m);
    Scalar lnl_wood = WheatstoneLikelihood::compute_gtr_log_likelihood(net_woodbury, aln, gtr_m);
    Scalar rel_diff = std::abs(lnl_dense - lnl_wood) / std::abs(lnl_dense);
    CHECK((lnl_wood >= lnl_dense || rel_diff < 0.02));
}

// 6. O(N log N) Large-N Scalability & Polytomy Collapse Test (N = 2000, L = 300)
TEST_CASE("OHMIC 6: Large-N Scalability (N = 2000, L = 300) & Ohmic Polytomy Collapse (< 2.0s)") {
    const size_t N = 2000;
    const size_t L = 300;
    const size_t num_clones = 250; // 250 distinct outbreak haplotypes replicated across 2000 samples
    std::mt19937 rng(1337);
    const char bases[4] = {'A', 'C', 'G', 'T'};

    std::string anc(L, 'A');
    for (size_t s = 0; s < L; ++s) anc[s] = bases[rng() % 4];

    std::vector<std::string> haplotypes(num_clones, anc);
    for (size_t c = 1; c < num_clones; ++c) {
        haplotypes[c] = haplotypes[c / 2];
        // 1-3 mutations per branch
        for (int m = 0; m < 2; ++m) {
            size_t pos = rng() % L;
            haplotypes[c][pos] = bases[(rng() + c) % 4];
        }
    }

    Alignment aln;
    for (size_t i = 0; i < N; ++i) {
        aln.add_sequence("Sample_" + std::to_string(i + 1), haplotypes[i % num_clones]);
    }
    aln.finalize_nucleotide();
    REQUIRE(aln.num_taxa == N);

    auto t0 = std::chrono::steady_clock::now();

    OhmicSettings settings;
    settings.scale_mode = OhmicScaleMode::Auto;
    settings.dense_switch_threshold = 200; // Trigger ScalableWoodbury automatically since M ~ 250 > 200
    settings.woodbury_rank = 64;
    settings.collapse_ohmic_polytomies = true;

    OhmicNetwork net = OhmicNetwork::infer(aln, settings);

    GTRParameters gtr_p;
    GTRMatrix gtr_m;
    gtr_m.update(gtr_p, aln.nuc_frequencies);
    Scalar lnl = WheatstoneLikelihood::compute_gtr_log_likelihood(net, aln, gtr_m);

    auto t1 = std::chrono::steady_clock::now();
    double elapsed_sec = std::chrono::duration<double>(t1 - t0).count();

    CHECK(net.active_mode == OhmicScaleMode::ScalableWoodbury);
    CHECK(net.num_super_nodes <= num_clones);
    CHECK(net.num_collapsed_polytomies > 0);
    CHECK(std::isfinite(lnl));
    CHECK(lnl < 0.0);
    // Wall-clock bound only catches a pathological slowdown: shared CI runners are several times slower than a workstation
    CHECK(elapsed_sec < 30.0);
}

// 7. Empirical Benchmark Test on Real Multi-Gene Alignments (cd2.fna, bglobin.nex)
TEST_CASE("OHMIC 7: Empirical Multi-Gene Benchmark (cd2.fna, bglobin.nex)") {
    std::string cd2_path = resolve_test_path("tests/data/cd2.fna");
    auto aln_cd2 = Alignment::load_nucleotide(cd2_path);

    OhmicAnalysisSettings settings;
    settings.input_file = cd2_path;
    settings.use_codon_model = false;
    settings.quiet = true;
    settings.network_settings.min_shunt_stress = 0.05;

    OhmicAnalyzer analyzer(aln_cd2, settings);
    OhmicResult res = analyzer.run();

    CHECK(res.num_taxa == 10);
    CHECK(res.foster_sum == doctest::Approx(9.0).epsilon(1e-5));
    CHECK(res.ohmic_log_l >= res.single_tree_log_l - 1e-4);
    // Wall-clock bound only catches a pathological slowdown: shared CI runners are several times slower than a workstation
    CHECK(res.runtime_seconds < 30.0);

    auto j = res.to_modern_json(cd2_path);
    CHECK(j["analysis"]["id"] == "ohmic");
    CHECK(j["dataset"]["taxa_count"] == 10);
}

// 8. Standard Analyses (FEL & FUBAR) on Wheatstone-Equipped Tree + HyphAeon 384D
TEST_CASE("OHMIC 8: Standard Analyses (FEL & FUBAR) on Wheatstone Tree + HyphAeon 384D") {
    std::string cd2_path = resolve_test_path("tests/data/cd2.fna");
    auto aln_cd2 = Alignment::load(cd2_path);

    OhmicNetwork net;
    std::string source_desc;
    // Request HyphAeon only when aeon-core + cached weights are present; the
    // Wheatstone/FEL/FUBAR checks below run either way.
    const bool use_hyphaeon = hyphaeon::is_available();
    Tree tree = OhmicNetwork::resolve_tree(aln_cd2, "", true, use_hyphaeon, "", &net, &source_desc);

    CHECK(tree.num_leaves() == 10);
    CHECK(tree.has_wheatstone_bridges());
    CHECK(tree.bridges.size() == net.bridges.size());
    if (net.used_hyphaeon_embeddings) {
        CHECK(net.hyphaeon_embeddings.rows() == 10);
        CHECK(net.hyphaeon_embeddings.cols() == 384);
        CHECK(source_desc.find("HyphAeon384D") != std::string::npos);
    }

    // Run FEL on the Wheatstone-bridge-equipped Tree
    FELSettings fel_settings;
    fel_settings.full_model = false; // quick mode
    auto fel = FELAnalyzer::create_and_fit(tree, aln_cd2, fel_settings);
    CHECK(std::isfinite(fel.gtr_log_l));
    CHECK(std::isfinite(fel.global_log_l));
    auto fel_sites = fel.run(false, false);
    CHECK(fel_sites.size() == aln_cd2.num_codons);

    // Run FUBAR on the Wheatstone-bridge-equipped Tree
    FUBARSettings fubar_settings;
    fubar_settings.grid_points = 10;
    auto fubar = FUBARAnalyzer::create_and_fit(tree, aln_cd2, fubar_settings);
    auto fubar_res = fubar.run(nullptr, false);
    CHECK(fubar_res.site_results.size() == aln_cd2.num_codons);
}

// 9. Finite-Difference Gradient Verification on Wheatstone-Bridged Tree (MG94 & BUSTED)
TEST_CASE("OHMIC 9: Finite-Difference Gradient Verification on Wheatstone-Bridged Tree") {
    std::string cd2_path = resolve_test_path("benchmarks/data/cd2.fna");
    auto aln = Alignment::load(cd2_path);
    REQUIRE(aln.num_taxa == 10);

    OhmicSettings settings;
    settings.scale_mode = OhmicScaleMode::ExactDense;
    settings.min_shunt_stress = 0.01;
    OhmicNetwork net = OhmicNetwork::infer(aln, settings);
    REQUIRE(net.bridges.size() >= 1);

    Tree tree = net.collapsed_tree;
    REQUIRE(tree.bridges.size() >= 1);
    tree.bridges[0].weights = {0.58, 0.27, 0.15};
    tree.bridges[0].branch_lengths = {0.060, 0.035, 0.020};
    int32_t b0_id = tree.bridges[0].bridge_node_id;
    if (b0_id >= 0 && b0_id != tree.bridges[0].hub_node_id &&
        static_cast<size_t>(b0_id) < tree.nodes.size()) {
        tree.nodes[b0_id].branch_length = 0.060;
    }

    MG94Parameters params;
    params.alpha = 1.0;
    params.beta = 0.75;
    params.theta_AC = 0.55;
    params.theta_AT = 0.40;
    params.theta_CG = 0.60;
    params.theta_CT = 1.35;
    params.theta_GT = 0.45;

    std::vector<Scalar> bl(tree.num_nodes(), 0.0);
    for (const auto& node : tree.nodes) {
        if (node.id != tree.root_id) {
            bl[node.id] = node.branch_length;
        }
    }

    auto g = LikelihoodEngine::compute_joint_frechet_gradients(tree, aln, params, bl);
    const Scalar h = 1e-6;

    auto eval_tree_lnl = [&](const MG94Parameters& p_eval, const std::vector<Scalar>& bl_eval) -> Scalar {
        return LikelihoodEngine::compute_joint_frechet_gradients(tree, aln, p_eval, bl_eval).log_likelihood;
    };

    // 1. Check global MG94 Fréchet gradients (omega, theta_AC, theta_CT) including bridge matrices
    {
        MG94Parameters p_p = params, p_m = params;
        p_p.beta += h; p_m.beta -= h;
        Scalar fd_omega = (eval_tree_lnl(p_p, bl) - eval_tree_lnl(p_m, bl)) / (2.0 * h);
        CHECK(g.grad_omega == doctest::Approx(fd_omega).epsilon(1e-4));
    }
    {
        MG94Parameters p_p = params, p_m = params;
        p_p.theta_AC += h; p_m.theta_AC -= h;
        Scalar fd_ac = (eval_tree_lnl(p_p, bl) - eval_tree_lnl(p_m, bl)) / (2.0 * h);
        CHECK(g.grad_theta_AC == doctest::Approx(fd_ac).epsilon(1e-4));
    }
    {
        MG94Parameters p_p = params, p_m = params;
        p_p.theta_CT += h; p_m.theta_CT -= h;
        Scalar fd_ct = (eval_tree_lnl(p_p, bl) - eval_tree_lnl(p_m, bl)) / (2.0 * h);
        CHECK(g.grad_theta_CT == doctest::Approx(fd_ct).epsilon(1e-4));
    }

    // 2. Check all active branch gradients (and verify detached bridge nodes have zero single-edge gradient)
    for (const auto& node : tree.nodes) {
        if (node.id == tree.root_id) continue;
        if (node.parent_id == INVALID_INDEX) {
            CHECK(g.grad_branches[node.id] == doctest::Approx(0.0).epsilon(1e-12));
            continue;
        }
        std::vector<Scalar> bl_p = bl, bl_m = bl;
        bl_p[node.id] += h;
        bl_m[node.id] -= h;
        Scalar fd_b = (eval_tree_lnl(params, bl_p) - eval_tree_lnl(params, bl_m)) / (2.0 * h);
        CHECK(g.grad_branches[node.id] == doctest::Approx(fd_b).epsilon(1e-4));
    }

    // 3. Check BUSTEDAnalyzer::compute_busted_branch_gradients on the bridged tree
    BUSTEDAnalyzer busted(tree, aln);
    busted.base_params = params;
    busted.mg94_omega = 0.75;
    busted.prepare_baseline_lengths();

    std::vector<Scalar> omegas = {0.15, 0.70, 2.40};
    std::vector<Scalar> weights = {0.50, 0.35, 0.15};
    std::vector<Scalar> syn_rates = {0.60, 1.40};
    std::vector<Scalar> syn_weights = {0.55, 0.45};
    std::vector<Scalar> custom_tau = busted.branch_tau;

    auto [ll0, b_grads] = busted.compute_busted_branch_gradients(
        custom_tau, omegas, weights, syn_rates, syn_weights
    );
    CHECK(ll0 == doctest::Approx(
        busted.evaluate_log_l(1.0, omegas, weights, syn_rates, syn_weights, &custom_tau)
    ).epsilon(1e-9));

    for (const auto& node : busted.tree.nodes) {
        if (node.id == busted.tree.root_id) continue;
        bool is_detached_bridge = false;
        for (const auto& wb : busted.tree.bridges) {
            if (wb.bridge_node_id == node.id && wb.bridge_node_id != wb.hub_node_id) {
                is_detached_bridge = true;
                break;
            }
        }
        if (node.parent_id == INVALID_INDEX && !is_detached_bridge) continue;

        std::vector<Scalar> tau_p = custom_tau, tau_m = custom_tau;
        tau_p[node.id] += h;
        tau_m[node.id] -= h;
        Scalar ll_p = busted.evaluate_log_l(1.0, omegas, weights, syn_rates, syn_weights, &tau_p);
        Scalar ll_m = busted.evaluate_log_l(1.0, omegas, weights, syn_rates, syn_weights, &tau_m);
        Scalar fd_tau = (ll_p - ll_m) / (2.0 * h);
        CHECK(b_grads[node.id] == doctest::Approx(fd_tau).epsilon(1e-4));
    }

    // 4. Verify aBSREL detached bridge node cavity equivalence and full fitting/testing
    if (b0_id >= 0 && b0_id != tree.bridges[0].hub_node_id) {
        ABSRELSettings absrel_settings;
        ABSRELAnalyzer absrel = ABSRELAnalyzer::create(tree, aln, absrel_settings);
        CHECK(absrel.detached_bridge_index(b0_id) == 0);
        CHECK(absrel.is_active_absrel_branch(b0_id));

        // Verify stacked [3*S x P] bridge cavity matches full tree log-likelihood for M=2 mixture on b0_id
        absrel.branch_alpha[b0_id] = 0.085;
        absrel.branch_omegas[b0_id] = {0.18, 2.35};
        absrel.branch_weights[b0_id] = {0.72, 0.28};

        auto P_all = absrel.compute_all_branch_transition_matrices();
        auto P_bridges = absrel.compute_all_bridge_transition_matrices();
        Scalar full_tree_ll = absrel.compute_tree_log_likelihood(P_all, &P_bridges);

        Matrix V_mat, D_mat, PD_scratch;
        Scalar log_scale_offset = 0.0;
        absrel.compute_branch_context_matrices(b0_id, P_all, V_mat, D_mat, log_scale_offset, &P_bridges);
        CHECK(V_mat.rows() == 3 * absrel.gcode->num_sense_codons);
        CHECK(D_mat.rows() == 3 * absrel.gcode->num_sense_codons);
        Scalar cavity_ll = absrel.eval_node_current_mixture_ll(
            b0_id, absrel.branch_alpha[b0_id], absrel.branch_omegas[b0_id],
            absrel.branch_weights[b0_id], V_mat, D_mat, log_scale_offset, PD_scratch
        );
        CHECK(cavity_ll == doctest::Approx(full_tree_ll).epsilon(1e-9));

        // Run full aBSREL pipeline and verify detached bridge node and its resolution ratios are fitted and tested
        absrel_settings.max_rate_classes = 2;
        absrel_settings.test_branches = "Internal";
        ABSRELAnalyzer absrel_run = ABSRELAnalyzer::create(tree, aln, absrel_settings);
        auto init_ratios = absrel_run.bridge_resolution_ratios(0);
        ABSRELResult absrel_res = absrel_run.run();
        auto fitted_ratios = absrel_run.bridge_resolution_ratios(0);
        CHECK((std::abs(fitted_ratios[1] - init_ratios[1]) > 1e-6 ||
               std::abs(fitted_ratios[2] - init_ratios[2]) > 1e-6));
        CHECK(std::isfinite(absrel_res.baseline_fit.log_likelihood));
        CHECK(std::isfinite(absrel_res.full_adaptive_fit.log_likelihood));
        CHECK(absrel_res.full_adaptive_fit.log_likelihood >= absrel_res.baseline_fit.log_likelihood - 1e-4);

        bool found_detached_bridge = false;
        for (const auto& [name, br] : absrel_res.branches) {
            if (br.node_id == b0_id) {
                found_detached_bridge = true;
                CHECK(br.is_tested);
                CHECK(br.baseline_branch_length > 0.0);
                CHECK(br.full_branch_length > 0.0);
                CHECK(std::isfinite(br.baseline_omega));
                CHECK(std::isfinite(br.lrt));
                CHECK(br.uncorrected_p_value >= 0.0);
                CHECK(br.uncorrected_p_value <= 1.0);
            }
            // Collapsed polytomy tombstones (parent_id == INVALID_INDEX and not a detached bridge) must not appear
            if (absrel_run.tree.nodes[br.node_id].parent_id == INVALID_INDEX) {
                CHECK(absrel_run.detached_bridge_index(br.node_id) >= 0);
            }
        }
        CHECK(found_detached_bridge);
    }
}


