#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/string_view.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/unordered_map.h>
#include <nanobind/eigen/dense.h>

#include "hyphy/core/types.hpp"
#include "hyphy/core/genetic_code.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/core/rate_matrix.hpp"
#include "hyphy/core/likelihood.hpp"
#include "hyphy/analyses/fel.hpp"
#include "hyphy/analyses/meme.hpp"
#include "hyphy/analyses/busted.hpp"
#include "hyphy/analyses/absrel.hpp"

namespace nb = nanobind;
using namespace nb::literals;
using namespace hyphy::core;
using namespace hyphy::analyses;

NB_MODULE(_hyphy3, m) {
    m.doc() = "HyPhy 3: Next-Generation High-Performance Molecular Evolution & Phylogenetics";

    // GeneticCode
    nb::class_<GeneticCode>(m, "GeneticCode")
        .def_static("universal", &GeneticCode::universal)
        .def_static("from_name", &GeneticCode::from_name, "name"_a)
        .def_ro("name", &GeneticCode::name)
        .def_ro("num_sense_codons", &GeneticCode::num_sense_codons);

    // TreeNode
    nb::class_<TreeNode>(m, "TreeNode")
        .def_ro("id", &TreeNode::id)
        .def_ro("name", &TreeNode::name)
        .def_ro("parent_id", &TreeNode::parent_id)
        .def_ro("children", &TreeNode::children)
        .def_rw("branch_length", &TreeNode::branch_length)
        .def_ro("is_leaf", &TreeNode::is_leaf);

    // Tree
    nb::class_<Tree>(m, "Tree")
        .def(nb::init<>())
        .def_static("from_newick", &Tree::from_newick, "newick"_a)
        .def_static("from_newick_file", &Tree::from_newick_file, "filepath"_a)
        .def("num_nodes", &Tree::num_nodes)
        .def("num_leaves", &Tree::num_leaves)
        .def_ro("root_id", &Tree::root_id)
        .def_ro("post_order", &Tree::post_order)
        .def_ro("nodes", &Tree::nodes)
        .def_ro("leaf_name_to_id", &Tree::leaf_name_to_id)
        .def("get_branch_lengths", [](const Tree& tree) {
            std::vector<Scalar> bl(tree.num_nodes(), 0.0);
            for (size_t i = 0; i < tree.nodes.size(); ++i) {
                bl[i] = tree.nodes[i].branch_length;
            }
            return bl;
        })
        .def("set_branch_lengths", [](Tree& tree, const std::vector<Scalar>& bl) {
            for (size_t i = 0; i < tree.nodes.size() && i < bl.size(); ++i) {
                tree.nodes[i].branch_length = bl[i];
            }
        })
        .def("get_branch_length", [](const Tree& tree, int32_t id) {
            return tree.nodes.at(id).branch_length;
        }, "id"_a)
        .def("set_branch_length", [](Tree& tree, int32_t id, Scalar length) {
            tree.nodes.at(id).branch_length = length;
        }, "id"_a, "length"_a);

    // Alignment
    nb::class_<Alignment>(m, "Alignment")
        .def(nb::init<>())
        .def_static("load", [](const std::string& filepath) { return Alignment::load(filepath); }, "filepath"_a)
        .def_static("from_fasta", [](const std::string& filepath) { return Alignment::from_fasta(filepath); }, "filepath"_a)
        .def_static("from_nexus", [](const std::string& filepath) { return Alignment::from_nexus(filepath); }, "filepath"_a)
        .def_ro("num_taxa", &Alignment::num_taxa)
        .def_prop_ro("num_sequences", [](const Alignment& a) { return a.num_taxa; })
        .def_ro("num_codons", &Alignment::num_codons)
        .def_ro("num_nucleotides", &Alignment::num_nucleotides)
        .def_prop_ro("num_sites", [](const Alignment& a) { return a.num_codons; })
        .def_prop_ro("num_patterns", [](const Alignment& a) { return a.patterns.size(); })
        .def_ro("embedded_tree_newick", &Alignment::embedded_tree_newick)
        .def_ro("taxon_names", &Alignment::taxon_names);

    // MG94Parameters
    nb::class_<MG94Parameters>(m, "MG94Parameters")
        .def(nb::init<>())
        .def_rw("alpha", &MG94Parameters::alpha)
        .def_rw("beta", &MG94Parameters::beta)
        .def_rw("theta_AC", &MG94Parameters::theta_AC)
        .def_rw("theta_AT", &MG94Parameters::theta_AT)
        .def_rw("theta_CG", &MG94Parameters::theta_CG)
        .def_rw("theta_CT", &MG94Parameters::theta_CT)
        .def_rw("theta_GT", &MG94Parameters::theta_GT);

    // FEL Results & Analyzer
    nb::class_<SiteResult>(m, "FELSiteResult")
        .def_ro("site_index", &SiteResult::site_index)
        .def_ro("alpha", &SiteResult::alpha)
        .def_ro("beta", &SiteResult::beta)
        .def_ro("alpha_null", &SiteResult::alpha_null)
        .def_ro("lrt", &SiteResult::lrt)
        .def_ro("p_value", &SiteResult::p_value)
        .def_ro("total_branch_length", &SiteResult::total_branch_length)
        .def_ro("log_l_alt", &SiteResult::log_l_alt)
        .def_ro("log_l_null", &SiteResult::log_l_null);

    nb::class_<FELAnalyzer>(m, "FELAnalyzer")
        .def("__init__", [](FELAnalyzer* self, Tree tree, Alignment aln, const MG94Parameters& params) {
            new (self) FELAnalyzer(std::move(tree), std::move(aln), params);
        }, "tree"_a, "alignment"_a, "params"_a = MG94Parameters{})
        .def_static("create_and_fit", [](Tree tree, Alignment aln, Scalar pvalue_threshold) {
            return FELAnalyzer::create_and_fit(std::move(tree), std::move(aln), pvalue_threshold);
        }, "tree"_a, "alignment"_a, "pvalue_threshold"_a = 0.1)
        .def("run", &FELAnalyzer::run)
        .def("to_json", [](const FELAnalyzer& fel) {
            return fel.to_json().dump();
        })
        .def_ro("global_log_l", &FELAnalyzer::global_log_l)
        .def_ro("global_aicc", &FELAnalyzer::global_aicc)
        .def_ro("site_results", &FELAnalyzer::site_results);

    // MEME Results & Analyzer
    nb::class_<MEMESiteResult>(m, "MEMESiteResult")
        .def_ro("site_index", &MEMESiteResult::site_index)
        .def_ro("alpha", &MEMESiteResult::alpha)
        .def_ro("beta1", &MEMESiteResult::beta1)
        .def_ro("p1", &MEMESiteResult::p1)
        .def_ro("beta_plus", &MEMESiteResult::beta_plus)
        .def_ro("p_plus", &MEMESiteResult::p_plus)
        .def_ro("lrt", &MEMESiteResult::lrt)
        .def_ro("p_value", &MEMESiteResult::p_value)
        .def_ro("branches_under_selection", &MEMESiteResult::branches_under_selection)
        .def_ro("total_branch_length", &MEMESiteResult::total_branch_length)
        .def_ro("log_l_meme", &MEMESiteResult::log_l_meme)
        .def_ro("log_l_fel", &MEMESiteResult::log_l_fel)
        .def_ro("p_value_meme_vs_fel", &MEMESiteResult::p_value_meme_vs_fel)
        .def_ro("branch_ebf", &MEMESiteResult::branch_ebf);

    nb::class_<MEMEAnalyzer>(m, "MEMEAnalyzer")
        .def("__init__", [](MEMEAnalyzer* self, Tree tree, Alignment aln, const MG94Parameters& params) {
            new (self) MEMEAnalyzer(std::move(tree), std::move(aln), params);
        }, "tree"_a, "alignment"_a, "params"_a = MG94Parameters{})
        .def_static("create_and_fit", [](Tree tree, Alignment aln, Scalar pvalue_threshold) {
            return MEMEAnalyzer::create_and_fit(std::move(tree), std::move(aln), pvalue_threshold);
        }, "tree"_a, "alignment"_a, "pvalue_threshold"_a = 0.1)
        .def("run", &MEMEAnalyzer::run)
        .def("to_json", [](const MEMEAnalyzer& meme) {
            return meme.to_json().dump();
        })
        .def_ro("global_log_l", &MEMEAnalyzer::global_log_l)
        .def_ro("global_aicc", &MEMEAnalyzer::global_aicc)
        .def_ro("site_results", &MEMEAnalyzer::site_results);

    // BUSTED Results & Analyzer
    nb::class_<BUSTEDSettings>(m, "BUSTEDSettings")
        .def(nb::init<>())
        .def_rw("num_rate_classes", &BUSTEDSettings::num_rate_classes)
        .def_rw("srv", &BUSTEDSettings::srv)
        .def_rw("num_syn_rate_classes", &BUSTEDSettings::num_syn_rate_classes)
        .def_rw("auto_select_k", &BUSTEDSettings::auto_select_k)
        .def_rw("max_k", &BUSTEDSettings::max_k)
        .def_rw("refine_branch_lengths", &BUSTEDSettings::refine_branch_lengths)
        .def_rw("p_value_threshold", &BUSTEDSettings::p_value_threshold);

    nb::class_<BUSTEDRateDistribution>(m, "BUSTEDRateDistribution")
        .def_ro("omegas", &BUSTEDRateDistribution::omegas)
        .def_ro("weights", &BUSTEDRateDistribution::weights)
        .def_ro("syn_rates", &BUSTEDRateDistribution::syn_rates)
        .def_ro("syn_weights", &BUSTEDRateDistribution::syn_weights)
        .def_ro("annotations", &BUSTEDRateDistribution::annotations);

    nb::class_<BUSTEDFit>(m, "BUSTEDFit")
        .def_ro("log_likelihood", &BUSTEDFit::log_likelihood)
        .def_ro("aicc", &BUSTEDFit::aicc)
        .def_ro("tree_scale", &BUSTEDFit::tree_scale)
        .def_ro("num_rate_classes", &BUSTEDFit::num_rate_classes)
        .def_ro("test_distribution", &BUSTEDFit::test_distribution)
        .def_ro("branch_lengths", &BUSTEDFit::branch_lengths);

    nb::class_<BUSTEDResult>(m, "BUSTEDResult")
        .def_ro("optimal_k", &BUSTEDResult::optimal_k)
        .def_ro("lrt", &BUSTEDResult::lrt)
        .def_ro("p_value", &BUSTEDResult::p_value)
        .def_ro("unconstrained", &BUSTEDResult::unconstrained)
        .def_ro("constrained", &BUSTEDResult::constrained)
        .def_ro("evidence_ratios", &BUSTEDResult::evidence_ratios)
        .def_ro("runtime_seconds", &BUSTEDResult::runtime_seconds);

    nb::class_<BUSTEDAnalyzer>(m, "BUSTEDAnalyzer")
        .def("__init__", [](BUSTEDAnalyzer* self, Tree tree, Alignment aln, const MG94Parameters& params) {
            new (self) BUSTEDAnalyzer(std::move(tree), std::move(aln), params);
        }, "tree"_a, "alignment"_a, "params"_a = MG94Parameters{})
        .def_static("create_and_fit", &BUSTEDAnalyzer::create_and_fit, "tree"_a, "alignment"_a)
        .def("run", &BUSTEDAnalyzer::run, "settings"_a = BUSTEDSettings{})
        .def("to_json", [](const BUSTEDAnalyzer& busted, const BUSTEDResult& res) {
            return busted.to_json(res).dump();
        });

    // aBSREL Results & Analyzer
    nb::class_<ABSRELSettings>(m, "ABSRELSettings")
        .def(nb::init<>())
        .def_rw("max_rate_classes", &ABSRELSettings::max_rate_classes)
        .def_rw("p_threshold", &ABSRELSettings::p_threshold)
        .def_rw("test_branches", &ABSRELSettings::test_branches)
        .def_rw("do_srv", &ABSRELSettings::do_srv)
        .def_rw("syn_rate_classes", &ABSRELSettings::syn_rate_classes)
        .def_rw("verbose", &ABSRELSettings::verbose);

    nb::class_<ABSRELRateDistribution>(m, "ABSRELRateDistribution")
        .def(nb::init<>())
        .def_ro("rates", &ABSRELRateDistribution::rates)
        .def_ro("weights", &ABSRELRateDistribution::weights);

    nb::class_<ABSRELBranchResult>(m, "ABSRELBranchResult")
        .def(nb::init<>())
        .def_ro("branch_name", &ABSRELBranchResult::branch_name)
        .def_ro("node_id", &ABSRELBranchResult::node_id)
        .def_ro("baseline_omega", &ABSRELBranchResult::baseline_omega)
        .def_ro("baseline_branch_length", &ABSRELBranchResult::baseline_branch_length)
        .def_ro("rate_classes", &ABSRELBranchResult::rate_classes)
        .def_ro("rate_distribution", &ABSRELBranchResult::rate_distribution)
        .def_ro("full_branch_length", &ABSRELBranchResult::full_branch_length)
        .def_ro("full_es", &ABSRELBranchResult::full_es)
        .def_ro("full_en", &ABSRELBranchResult::full_en)
        .def_ro("lrt", &ABSRELBranchResult::lrt)
        .def_ro("uncorrected_p_value", &ABSRELBranchResult::uncorrected_p_value)
        .def_ro("corrected_p_value", &ABSRELBranchResult::corrected_p_value)
        .def_ro("is_tested", &ABSRELBranchResult::is_tested)
        .def_ro("is_positive", &ABSRELBranchResult::is_positive)
        .def_ro("sites_ebf_100", &ABSRELBranchResult::sites_ebf_100)
        .def_ro("site_ebf", &ABSRELBranchResult::site_ebf);

    nb::class_<ABSRELFitSummary>(m, "ABSRELFitSummary")
        .def(nb::init<>())
        .def_ro("log_likelihood", &ABSRELFitSummary::log_likelihood)
        .def_ro("parameters", &ABSRELFitSummary::parameters)
        .def_ro("aicc", &ABSRELFitSummary::aicc);

    nb::class_<ABSRELResult>(m, "ABSRELResult")
        .def(nb::init<>())
        .def_ro("gtr_fit", &ABSRELResult::gtr_fit)
        .def_ro("baseline_fit", &ABSRELResult::baseline_fit)
        .def_ro("full_adaptive_fit", &ABSRELResult::full_adaptive_fit)
        .def_ro("branches", &ABSRELResult::branches)
        .def_ro("tested_branches", &ABSRELResult::tested_branches)
        .def_ro("positive_branches", &ABSRELResult::positive_branches)
        .def_ro("p_threshold", &ABSRELResult::p_threshold)
        .def_ro("runtime_seconds", &ABSRELResult::runtime_seconds)
        .def("to_json", [](const ABSRELResult& res, const Tree& tree, const Alignment& aln) {
            return res.to_json(tree, aln).dump();
        }, "tree"_a, "alignment"_a);

    nb::class_<ABSRELAnalyzer>(m, "ABSRELAnalyzer")
        .def("__init__", [](ABSRELAnalyzer* self, Tree tree, Alignment aln, const ABSRELSettings& settings) {
            new (self) ABSRELAnalyzer(std::move(tree), std::move(aln), settings);
        }, "tree"_a, "alignment"_a, "settings"_a = ABSRELSettings{})
        .def_static("create", &ABSRELAnalyzer::create, "tree"_a, "alignment"_a, "settings"_a = ABSRELSettings{})
        .def("run", [](ABSRELAnalyzer& absrel) {
            return absrel.run();
        });

    // Differentiable Engine: Inside-Outside Likelihood and Analytical Adjoint Gradients!
    m.def("compute_branch_length_gradients", [](
        const Tree& tree,
        const Alignment& aln,
        const MG94Parameters& params,
        const std::vector<Scalar>& branch_lengths
    ) -> std::pair<Scalar, std::vector<Scalar>> {
        return LikelihoodEngine::compute_branch_length_gradients(tree, aln, params, branch_lengths);
    }, "tree"_a, "aln"_a, "params"_a, "branch_lengths"_a = std::vector<Scalar>{}, "Compute log-likelihood and exact analytical branch gradients d ln L / d t_b via Inside-Outside");
}
