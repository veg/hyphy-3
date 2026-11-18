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

TEST_CASE("NEXUS: Embedded tree TRANSLATE resolution on CD2.nex and edge cases") {
    // 1. Alignment::load resolves TRANSLATE IDs (1->Pig, 2->Cow, 3->Horse, 8->Cat) inside embedded_tree_newick
    auto aln = Alignment::load(HYPHY3_SOURCE_DIR "/tests/data/CD2.nex");
    REQUIRE(!aln.embedded_tree_newick.empty());

    auto emb_tree = Tree::from_newick(aln.embedded_tree_newick);
    REQUIRE(emb_tree.num_leaves() == 10);
    for (const char* taxon : {"Pig", "Cow", "Horse", "Cat", "RhMonkey", "Baboon", "Human", "Chimp", "Rat", "Mouse"}) {
        CHECK(emb_tree.leaf_name_to_id.count(taxon) == 1);
    }
    for (const char* numeric_id : {"1", "2", "3", "8"}) {
        CHECK(emb_tree.leaf_name_to_id.count(numeric_id) == 0);
    }
    auto pr_sel = emb_tree.select_branches("PR");
    CHECK(pr_sel.num_selected == 8); // 4 primate tips + 4 annotated internal PR branches

    auto leaf_map = LikelihoodEngine::build_leaf_to_taxon_map(emb_tree, aln);
    CHECK(leaf_map.size() == emb_tree.num_nodes());

    // 2. Tree::from_newick_file directly parses a #NEXUS file with BEGIN TREES + TRANSLATE
    auto direct_tree = Tree::from_newick_file(HYPHY3_SOURCE_DIR "/tests/data/CD2.nex");
    REQUIRE(direct_tree.num_leaves() == 10);
    for (const char* taxon : {"Pig", "Cow", "Horse", "Cat", "RhMonkey", "Baboon", "Human", "Chimp", "Rat", "Mouse"}) {
        CHECK(direct_tree.leaf_name_to_id.count(taxon) == 1);
    }

    // 3. Quoted multi-word names, node model tags, internal bootstrap numbers, and trailing BEGIN HYPHY block
    std::string edge_nex = R"(#NEXUS
BEGIN TREES;
  TRANSLATE
    1 'Canis lupus',
    2 "Felis catus";
  UTREE * t1 = [&U] ((1{FG}:0.15,2:0.25)1:0.05,Outgroup:0.30);
END;
BEGIN HYPHY;
  Tree ignored_tree = (A:0.1,B:0.1);
END;
)";
    std::string nwk = extract_nexus_tree_newick(strip_nexus_comments(edge_nex), 0);
    REQUIRE(!nwk.empty());
    auto edge_tree = Tree::from_newick(nwk);
    REQUIRE(edge_tree.num_leaves() == 3);
    CHECK(edge_tree.leaf_name_to_id.count("Canis lupus") == 1);
    CHECK(edge_tree.leaf_name_to_id.count("Felis catus") == 1);
    CHECK(edge_tree.leaf_name_to_id.count("Outgroup") == 1);
    CHECK(edge_tree.leaf_name_to_id.count("1") == 0);
    CHECK(edge_tree.select_branches("FG").num_selected == 1);

    bool found_internal_bootstrap = false;
    for (const auto& nd : edge_tree.nodes) {
        if (!nd.is_leaf && nd.name == "1") found_internal_bootstrap = true;
    }
    CHECK(found_internal_bootstrap);
}

TEST_CASE("Issue #25: Tree::parse_label_map_file, apply_label_map, and apply_branch_labels") {
    std::string tmp_csv = "test_issue25_labels.csv";
    {
        std::ofstream out(tmp_csv);
        out << "# Taxon-to-group CSV with header and extra absent species\n"
            << "Species,Group\n"
            << "Pig,group_A\n"
            << "Cow,group_A\n"
            << "Horse,group_A\n"
            << "Cat,group_A\n"
            << "Human,group_B\n"
            << "Chimp,group_B\n"
            << "RhMonkey,group_B\n"
            << "Baboon,group_B\n"
            << "Rat,reference\n"
            << "Mouse,reference\n"
            << "Missing_Sp1,group_A\n";
    }
    auto groups = Tree::parse_label_map_file(tmp_csv);
    std::remove(tmp_csv.c_str());

    REQUIRE(groups.size() == 3);
    CHECK(groups[0].first == "group_A");
    CHECK(groups[1].first == "group_B");
    CHECK(groups[2].first == "reference");

    const std::string nwk =
        "((((Pig:0.14,Cow:0.21):0.08,Horse:0.16,Cat:0.26):0.05,"
        "((RhMonkey:0.002,Baboon:0.003):0.02,(Human:0.004,Chimp:0.001):0.01):0.10):0.34,"
        "Rat:0.05,Mouse:0.09);";

    // 1. All (monophyletic internal propagation)
    {
        Tree t = Tree::from_newick(nwk);
        size_t n_lbl = t.apply_label_map(groups, "All", tmp_csv);
        CHECK(n_lbl == 15); // 10 tips + 2 group_A internals + 3 group_B internals
        CHECK(t.select_branches("group_A").num_selected == 6); // 4 tips + (Pig,Cow) + ((Pig,Cow),Horse,Cat)
        CHECK(t.select_branches("group_B").num_selected == 7); // 4 tips + (Rh,Bab) + (Hum,Chi) + Primate ancestor
        CHECK(t.select_branches("reference").num_selected == 2); // Rat, Mouse (sibling root children)
    }

    // 2. Parsimony (Generalized Fitch over G+1 states)
    {
        Tree t = Tree::from_newick(nwk);
        t.apply_label_map(groups, "Parsimony", tmp_csv);
        CHECK(t.select_branches("group_A").num_selected >= 6);
        CHECK(t.select_branches("group_B").num_selected >= 7);
        CHECK(t.select_branches("reference").num_selected == 2);
    }

    // 3. None (tips only)
    {
        Tree t = Tree::from_newick(nwk);
        size_t n_lbl = t.apply_label_map(groups, "None", tmp_csv);
        CHECK(n_lbl == 10);
        CHECK(t.select_branches("group_A").num_selected == 4);
        CHECK(t.select_branches("group_B").num_selected == 4);
        CHECK(t.select_branches("reference").num_selected == 2);
    }

    // 4. apply_branch_labels (comma tip list with All vs None vs Parsimony)
    {
        Tree t_all = Tree::from_newick(nwk);
        t_all.apply_branch_labels("Foreground", "Human,Chimp,Baboon,RhMonkey", "All");
        CHECK(t_all.select_branches("Foreground").num_selected == 7);

        Tree t_none = Tree::from_newick(nwk);
        t_none.apply_branch_labels("Foreground", "Human,Chimp,Baboon,RhMonkey", "None");
        CHECK(t_none.select_branches("Foreground").num_selected == 4);

        Tree t_pars = Tree::from_newick(nwk);
        t_pars.apply_branch_labels("Foreground", "Human,Chimp,Baboon", "Parsimony");
        CHECK(t_pars.select_branches("Foreground").num_selected >= 4);

        CHECK_THROWS_AS(t_all.apply_branch_labels("Foreground", "Human", "InvalidMode"), std::runtime_error);
    }

    // 5. Error cases: conflicting duplicate taxon, and 0 matching taxa against tree
    {
        std::string tmp_conflict = "test_issue25_conflict.csv";
        {
            std::ofstream out(tmp_conflict);
            out << "Pig,group_A\nPig,group_B\n";
        }
        CHECK_THROWS_AS(Tree::parse_label_map_file(tmp_conflict), std::runtime_error);
        std::remove(tmp_conflict.c_str());

        Tree t = Tree::from_newick(nwk);
        std::vector<std::pair<std::string, std::vector<std::string>>> no_match = {
            {"group_A", {"NoSuchSpecies_1", "NoSuchSpecies_2"}}
        };
        CHECK_THROWS_AS(t.apply_label_map(no_match, "All", "nomatch.csv"), std::runtime_error);
    }
}
