#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/string_view.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/unordered_map.h>
#include <nanobind/stl/map.h>
#include <nanobind/stl/function.h>
#include <nanobind/stl/array.h>
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
#include "hyphy/analyses/relax.hpp"
#include "hyphy/analyses/fubar.hpp"
#include "hyphy/analyses/prime.hpp"
#include "hyphy/core/neighbor_joining.hpp"
#include "hyphy/analyses/gard.hpp"
#include "hyphy/analyses/slac.hpp"
#include "hyphy/analyses/directional_fel.hpp"
#include "hyphy/core/ohmic_network.hpp"
#include "hyphy/core/wheatstone_likelihood.hpp"
#include "hyphy/analyses/ohmic.hpp"

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
        .def_static("from_newick", &Tree::from_newick, "newick"_a, "keep_rooted"_a = true)
        .def_static("from_newick_file", &Tree::from_newick_file, "filepath"_a, "keep_rooted"_a = true)
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
        .def_static("load_nucleotide", [](const std::string& filepath) { return Alignment::load_nucleotide(filepath); }, "filepath"_a)
        .def_static("from_fasta", [](const std::string& filepath) { return Alignment::from_fasta(filepath); }, "filepath"_a)
        .def_static("from_fasta_nucleotide", [](const std::string& filepath) { return Alignment::from_fasta_nucleotide(filepath); }, "filepath"_a)
        .def_static("from_nexus", [](const std::string& filepath) { return Alignment::from_nexus(filepath); }, "filepath"_a)
        .def_static("from_nexus_nucleotide", [](const std::string& filepath) { return Alignment::from_nexus_nucleotide(filepath); }, "filepath"_a)
        .def("slice_nucleotide", &Alignment::slice_nucleotide, "start"_a, "end"_a)
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
        .def_rw("theta_GT", &MG94Parameters::theta_GT)
        .def_rw("delta", &MG94Parameters::delta)
        .def_rw("psi", &MG94Parameters::psi);

    // FEL Settings, Results & Analyzer
    nb::class_<FELSettings>(m, "FELSettings")
        .def(nb::init<>())
        .def_rw("p_value_threshold", &FELSettings::p_value_threshold)
        .def_rw("full_model", &FELSettings::full_model)
        .def_rw("srv", &FELSettings::srv)
        .def_rw("test_branches", &FELSettings::test_branches)
        .def_rw("test_branch_regex", &FELSettings::test_branch_regex)
        .def_rw("test_branch_names", &FELSettings::test_branch_names);

    nb::class_<SubstitutionEvent>(m, "SubstitutionEvent")
        .def_ro("node_id", &SubstitutionEvent::node_id)
        .def_ro("branch_name", &SubstitutionEvent::branch_name)
        .def_ro("is_tested", &SubstitutionEvent::is_tested)
        .def_ro("branch_length", &SubstitutionEvent::branch_length)
        .def_ro("parent_codon", &SubstitutionEvent::parent_codon)
        .def_ro("child_codon", &SubstitutionEvent::child_codon)
        .def_prop_ro("parent_aa", [](const SubstitutionEvent& ev) { return std::string(1, ev.parent_aa); })
        .def_prop_ro("child_aa", [](const SubstitutionEvent& ev) { return std::string(1, ev.child_aa); })
        .def_ro("nt_distance", &SubstitutionEvent::nt_distance)
        .def_ro("is_synonymous", &SubstitutionEvent::is_synonymous)
        .def_ro("transition_prob", &SubstitutionEvent::transition_prob)
        .def_ro("surprise_bits", &SubstitutionEvent::surprise_bits);

    nb::class_<SiteSubstitutions>(m, "SiteSubstitutions")
        .def_ro("site_index", &SiteSubstitutions::site_index)
        .def_ro("root_codon", &SiteSubstitutions::root_codon)
        .def_prop_ro("root_aa", [](const SiteSubstitutions& s) { return std::string(1, s.root_aa); })
        .def_ro("events", &SiteSubstitutions::events)
        .def_ro("syn_test", &SiteSubstitutions::syn_test)
        .def_ro("nonsyn_test", &SiteSubstitutions::nonsyn_test)
        .def_ro("syn_bg", &SiteSubstitutions::syn_bg)
        .def_ro("nonsyn_bg", &SiteSubstitutions::nonsyn_bg)
        .def_ro("subs_1nt", &SiteSubstitutions::subs_1nt)
        .def_ro("subs_2nt", &SiteSubstitutions::subs_2nt)
        .def_ro("subs_3nt", &SiteSubstitutions::subs_3nt)
        .def_ro("mean_surprise_bits", &SiteSubstitutions::mean_surprise_bits)
        .def_ro("total_surprise_bits", &SiteSubstitutions::total_surprise_bits)
        .def_ro("mean_branch_length", &SiteSubstitutions::mean_branch_length)
        .def_ro("site_perplexity", &SiteSubstitutions::site_perplexity)
        .def_ro("summary_string", &SiteSubstitutions::summary_string)
        .def_ro("summary_string_partitioned", &SiteSubstitutions::summary_string_partitioned);

    nb::class_<SiteResult>(m, "FELSiteResult")
        .def_ro("site_index", &SiteResult::site_index)
        .def_ro("alpha", &SiteResult::alpha)
        .def_ro("beta", &SiteResult::beta)
        .def_ro("beta_nuisance", &SiteResult::beta_nuisance)
        .def_ro("beta_background", &SiteResult::beta_nuisance)
        .def_ro("alpha_null", &SiteResult::alpha_null)
        .def_ro("lrt", &SiteResult::lrt)
        .def_ro("p_value", &SiteResult::p_value)
        .def_ro("total_branch_length", &SiteResult::total_branch_length)
        .def_ro("log_l_alt", &SiteResult::log_l_alt)
        .def_ro("log_l_null", &SiteResult::log_l_null)
        .def_ro("subs", &SiteResult::subs);

    nb::class_<FELAnalyzer>(m, "FELAnalyzer")
        .def("__init__", [](FELAnalyzer* self, Tree tree, Alignment aln, const MG94Parameters& params, const FELSettings& settings) {
            new (self) FELAnalyzer(std::move(tree), std::move(aln), params, settings);
        }, "tree"_a, "alignment"_a, "params"_a = MG94Parameters{}, "settings"_a = FELSettings{})
        .def_static("create_and_fit", [](Tree tree, Alignment aln, const FELSettings& settings) {
            return FELAnalyzer::create_and_fit(std::move(tree), std::move(aln), settings);
        }, "tree"_a, "alignment"_a, "settings"_a = FELSettings{})
        .def_static("create_and_fit", [](Tree tree, Alignment aln, Scalar pvalue_threshold, bool full_model) {
            return FELAnalyzer::create_and_fit(std::move(tree), std::move(aln), pvalue_threshold, nullptr, full_model);
        }, "tree"_a, "alignment"_a, "pvalue_threshold"_a = 0.1, "full_model"_a = true)
        .def("run", &FELAnalyzer::run, "show_progress"_a = false, "force_progress"_a = false)
        .def("to_json", [](const FELAnalyzer& fel, const std::string& fmt) {
            auto format = (fmt == "legacy") ? JSONFormat::Legacy : JSONFormat::ModernV3;
            return fel.to_json("", "", format).dump();
        }, "format"_a = "modern")
        .def("to_legacy_json", [](const FELAnalyzer& fel) {
            return fel.to_legacy_json().dump();
        })
        .def("to_modern_json", [](const FELAnalyzer& fel) {
            return fel.to_modern_json().dump();
        })
        .def_ro("global_log_l", &FELAnalyzer::global_log_l)
        .def_ro("global_aicc", &FELAnalyzer::global_aicc)
        .def_ro("site_results", &FELAnalyzer::site_results)
        .def_ro("num_tested", &FELAnalyzer::num_tested)
        .def_ro("num_background", &FELAnalyzer::num_background)
        .def_ro("test_branch_names", &FELAnalyzer::test_branch_names)
        .def_ro("bg_branch_names", &FELAnalyzer::bg_branch_names)
        .def_ro("is_test_branch", &FELAnalyzer::is_test_branch)
        .def_ro("settings", &FELAnalyzer::settings);

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
        .def_static("create_and_fit", [](Tree tree, Alignment aln, Scalar pvalue_threshold, bool full_model) {
            return MEMEAnalyzer::create_and_fit(std::move(tree), std::move(aln), pvalue_threshold, nullptr, full_model);
        }, "tree"_a, "alignment"_a, "pvalue_threshold"_a = 0.1, "full_model"_a = true)
        .def("run", &MEMEAnalyzer::run, "show_progress"_a = false, "force_progress"_a = false)
        .def("to_json", [](const MEMEAnalyzer& meme, const std::string& fmt) {
            auto format = (fmt == "legacy") ? JSONFormat::Legacy : JSONFormat::ModernV3;
            return meme.to_json(format).dump();
        }, "format"_a = "modern")
        .def("to_legacy_json", [](const MEMEAnalyzer& meme) {
            return meme.to_legacy_json().dump();
        })
        .def("to_modern_json", [](const MEMEAnalyzer& meme) {
            return meme.to_modern_json().dump();
        })
        .def_ro("global_log_l", &MEMEAnalyzer::global_log_l)
        .def_ro("global_aicc", &MEMEAnalyzer::global_aicc)
        .def_ro("site_results", &MEMEAnalyzer::site_results);

    // BUSTED Results & Analyzer
    nb::class_<BUSTEDSettings>(m, "BUSTEDSettings")
        .def(nb::init<>())
        .def_rw("full_auto", &BUSTEDSettings::full_auto)
        .def_rw("num_rate_classes", &BUSTEDSettings::num_rate_classes)
        .def_rw("srv", &BUSTEDSettings::srv)
        .def_rw("num_syn_rate_classes", &BUSTEDSettings::num_syn_rate_classes)
        .def_rw("auto_select_k", &BUSTEDSettings::auto_select_k)
        .def_rw("max_k", &BUSTEDSettings::max_k)
        .def_rw("refine_branch_lengths", &BUSTEDSettings::refine_branch_lengths)
        .def_rw("multiple_hits", &BUSTEDSettings::multiple_hits)
        .def_rw("error_sink", &BUSTEDSettings::error_sink)
        .def_rw("error_sink_max_weight", &BUSTEDSettings::error_sink_max_weight)
        .def_rw("error_sink_bf_threshold", &BUSTEDSettings::error_sink_bf_threshold)
        .def_rw("test_branches", &BUSTEDSettings::test_branches)
        .def_rw("phenotype_mode", &BUSTEDSettings::phenotype_mode)
        .def_rw("relax_mode", &BUSTEDSettings::relax_mode)
        .def_rw("include_relax_in_ph", &BUSTEDSettings::include_relax_in_ph)
        .def_rw("comparison_branches", &BUSTEDSettings::comparison_branches)
        .def_rw("p_value_threshold", &BUSTEDSettings::p_value_threshold);

    nb::class_<BUSTEDRateDistribution>(m, "BUSTEDRateDistribution")
        .def_ro("omegas", &BUSTEDRateDistribution::omegas)
        .def_ro("weights", &BUSTEDRateDistribution::weights)
        .def_ro("syn_rates", &BUSTEDRateDistribution::syn_rates)
        .def_ro("syn_weights", &BUSTEDRateDistribution::syn_weights)
        .def_ro("error_sink_weight", &BUSTEDRateDistribution::error_sink_weight)
        .def_ro("annotations", &BUSTEDRateDistribution::annotations);

    nb::class_<BUSTEDFit>(m, "BUSTEDFit")
        .def_ro("log_likelihood", &BUSTEDFit::log_likelihood)
        .def_ro("aicc", &BUSTEDFit::aicc)
        .def_ro("tree_scale", &BUSTEDFit::tree_scale)
        .def_ro("num_rate_classes", &BUSTEDFit::num_rate_classes)
        .def_ro("test_distribution", &BUSTEDFit::test_distribution)
        .def_ro("has_background_distribution", &BUSTEDFit::has_background_distribution)
        .def_ro("background_distribution", &BUSTEDFit::background_distribution)
        .def_ro("branch_lengths", &BUSTEDFit::branch_lengths)
        .def_ro("delta", &BUSTEDFit::delta)
        .def_ro("psi", &BUSTEDFit::psi)
        .def_ro("frac_delta", &BUSTEDFit::frac_delta)
        .def_ro("frac_psi", &BUSTEDFit::frac_psi)
        .def_ro("error_sink_weight", &BUSTEDFit::error_sink_weight);

    nb::class_<BUSTEDResult>(m, "BUSTEDResult")
        .def_ro("optimal_k", &BUSTEDResult::optimal_k)
        .def_ro("lrt", &BUSTEDResult::lrt)
        .def_ro("p_value", &BUSTEDResult::p_value)
        .def_ro("unconstrained", &BUSTEDResult::unconstrained)
        .def_ro("constrained", &BUSTEDResult::constrained)
        .def_ro("has_background_test", &BUSTEDResult::has_background_test)
        .def_ro("background_lrt", &BUSTEDResult::background_lrt)
        .def_ro("background_p_value", &BUSTEDResult::background_p_value)
        .def_ro("shared_lrt", &BUSTEDResult::shared_lrt)
        .def_ro("shared_p_value", &BUSTEDResult::shared_p_value)
        .def_ro("shared_df", &BUSTEDResult::shared_df)
        .def_ro("union_intersection_p_value", &BUSTEDResult::union_intersection_p_value)
        .def_ro("background_delta_aicc", &BUSTEDResult::background_delta_aicc)
        .def_ro("background_neutral_supported", &BUSTEDResult::background_neutral_supported)
        .def_ro("ranking_score", &BUSTEDResult::ranking_score)
        .def_ro("phenotype_associated", &BUSTEDResult::phenotype_associated)
        .def_ro("phenotype_verdict", &BUSTEDResult::phenotype_verdict)
        .def_ro("constrained_background", &BUSTEDResult::constrained_background)
        .def_ro("shared_distributions", &BUSTEDResult::shared_distributions)
        .def_ro("has_relax_test", &BUSTEDResult::has_relax_test)
        .def_ro("relax_k", &BUSTEDResult::relax_k)
        .def_ro("relax_lrt", &BUSTEDResult::relax_lrt)
        .def_ro("relax_p_value", &BUSTEDResult::relax_p_value)
        .def_ro("relax_direction", &BUSTEDResult::relax_direction)
        .def_ro("relax_alternative", &BUSTEDResult::relax_alternative)
        .def_ro("settings", &BUSTEDResult::settings)
        .def_ro("evidence_ratios", &BUSTEDResult::evidence_ratios)
        .def_ro("runtime_seconds", &BUSTEDResult::runtime_seconds);

    nb::class_<BUSTEDAnalyzer>(m, "BUSTEDAnalyzer")
        .def("__init__", [](BUSTEDAnalyzer* self, Tree tree, Alignment aln, const MG94Parameters& params) {
            new (self) BUSTEDAnalyzer(std::move(tree), std::move(aln), params);
        }, "tree"_a, "alignment"_a, "params"_a = MG94Parameters{})
        .def_static("create_and_fit", [](Tree tree, Alignment aln) {
            return BUSTEDAnalyzer::create_and_fit(std::move(tree), std::move(aln));
        }, "tree"_a, "alignment"_a)
        .def("run", &BUSTEDAnalyzer::run, "settings"_a = BUSTEDSettings{}, "show_progress"_a = false, "force_progress"_a = false)
        .def("to_json", [](const BUSTEDAnalyzer& busted, const BUSTEDResult& res, const std::string& fmt) {
            auto format = (fmt == "legacy") ? JSONFormat::Legacy : JSONFormat::ModernV3;
            return busted.to_json(res, format).dump();
        }, "res"_a, "format"_a = "modern")
        .def("to_legacy_json", [](const BUSTEDAnalyzer& busted, const BUSTEDResult& res) {
            return busted.to_legacy_json(res).dump();
        })
        .def("to_modern_json", [](const BUSTEDAnalyzer& busted, const BUSTEDResult& res) {
            return busted.to_modern_json(res).dump();
        });

    // aBSREL Results & Analyzer
    nb::class_<ABSRELSettings>(m, "ABSRELSettings")
        .def(nb::init<>())
        .def_rw("max_rate_classes", &ABSRELSettings::max_rate_classes)
        .def_rw("p_threshold", &ABSRELSettings::p_threshold)
        .def_rw("test_branches", &ABSRELSettings::test_branches)
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
        .def("to_json", [](const ABSRELResult& res, const Tree& tree, const Alignment& aln, const std::string& fmt) {
            auto format = (fmt == "legacy") ? JSONFormat::Legacy : JSONFormat::ModernV3;
            return res.to_json(tree, aln, format).dump();
        }, "tree"_a, "alignment"_a, "format"_a = "modern")
        .def("to_legacy_json", [](const ABSRELResult& res, const Tree& tree, const Alignment& aln) {
            return res.to_legacy_json(tree, aln).dump();
        }, "tree"_a, "alignment"_a)
        .def("to_modern_json", [](const ABSRELResult& res, const Tree& tree, const Alignment& aln) {
            return res.to_modern_json(tree, aln).dump();
        }, "tree"_a, "alignment"_a);

    nb::class_<ABSRELAnalyzer>(m, "ABSRELAnalyzer")
        .def("__init__", [](ABSRELAnalyzer* self, Tree tree, Alignment aln, const ABSRELSettings& settings) {
            new (self) ABSRELAnalyzer(std::move(tree), std::move(aln), settings);
        }, "tree"_a, "alignment"_a, "settings"_a = ABSRELSettings{})
        .def_static("create", &ABSRELAnalyzer::create, "tree"_a, "alignment"_a, "settings"_a = ABSRELSettings{})
        .def("run", [](ABSRELAnalyzer& absrel, bool show_progress, bool force_progress) {
            return absrel.run(nullptr, show_progress, force_progress);
        }, "show_progress"_a = false, "force_progress"_a = false);

    // RELAX Results & Analyzer
    nb::class_<RELAXSettings>(m, "RELAXSettings")
        .def(nb::init<>())
        .def_rw("p_value_threshold", &RELAXSettings::p_value_threshold)
        .def_rw("test_branch_regex", &RELAXSettings::test_branch_regex)
        .def_rw("test_branch_names", &RELAXSettings::test_branch_names)
        .def_rw("reference_branch_regex", &RELAXSettings::reference_branch_regex)
        .def_rw("reference_branch_names", &RELAXSettings::reference_branch_names)
        .def_rw("refine_branch_lengths", &RELAXSettings::refine_branch_lengths)
        .def_rw("max_joint_cycles", &RELAXSettings::max_joint_cycles)
        .def_rw("joint_tolerance", &RELAXSettings::joint_tolerance)
        .def_rw("verbose", &RELAXSettings::verbose);

    nb::class_<RELAXRateDistribution>(m, "RELAXRateDistribution")
        .def(nb::init<>())
        .def_ro("omegas", &RELAXRateDistribution::omegas)
        .def_ro("weights", &RELAXRateDistribution::weights);

    nb::class_<RELAXModelFit>(m, "RELAXModelFit")
        .def(nb::init<>())
        .def_ro("log_likelihood", &RELAXModelFit::log_likelihood)
        .def_ro("aicc", &RELAXModelFit::aicc)
        .def_ro("parameters", &RELAXModelFit::parameters)
        .def_ro("tree_scale", &RELAXModelFit::tree_scale)
        .def_ro("k", &RELAXModelFit::k)
        .def_ro("realized_cycles", &RELAXModelFit::realized_cycles)
        .def_ro("converged", &RELAXModelFit::converged)
        .def_ro("reference_distribution", &RELAXModelFit::reference_distribution)
        .def_ro("test_distribution", &RELAXModelFit::test_distribution)
        .def_ro("unclassified_distribution", &RELAXModelFit::unclassified_distribution)
        .def_ro("branch_lengths", &RELAXModelFit::branch_lengths)
        .def_ro("site_log_likelihoods", &RELAXModelFit::site_log_likelihoods);

    nb::class_<RELAXResult>(m, "RELAXResult")
        .def(nb::init<>())
        .def_ro("gtr_log_likelihood", &RELAXResult::gtr_log_likelihood)
        .def_ro("gtr_aicc", &RELAXResult::gtr_aicc)
        .def_ro("gtr_parameters", &RELAXResult::gtr_parameters)
        .def_ro("mg94_log_likelihood", &RELAXResult::mg94_log_likelihood)
        .def_ro("mg94_aicc", &RELAXResult::mg94_aicc)
        .def_ro("mg94_parameters", &RELAXResult::mg94_parameters)
        .def_ro("mg94_omega_R", &RELAXResult::mg94_omega_R)
        .def_ro("mg94_omega_T", &RELAXResult::mg94_omega_T)
        .def_ro("alternative_fit", &RELAXResult::alternative_fit)
        .def_ro("null_fit", &RELAXResult::null_fit)
        .def_ro("branch_switched_evidence_ratios", &RELAXResult::branch_switched_evidence_ratios)
        .def_ro("lrt", &RELAXResult::lrt)
        .def_ro("p_value", &RELAXResult::p_value)
        .def_ro("k", &RELAXResult::k)
        .def_ro("is_relaxed", &RELAXResult::is_relaxed)
        .def_ro("is_intensified", &RELAXResult::is_intensified)
        .def_ro("is_significant", &RELAXResult::is_significant)
        .def_ro("test_branches", &RELAXResult::test_branches)
        .def_ro("reference_branches", &RELAXResult::reference_branches)
        .def_ro("branch_class", &RELAXResult::branch_class)
        .def_ro("runtime_seconds", &RELAXResult::runtime_seconds)
        .def("to_json", [](const RELAXResult& res, const Tree& tree, const Alignment& aln, const std::string& fmt) {
            auto format = (fmt == "legacy") ? JSONFormat::Legacy : JSONFormat::ModernV3;
            return res.to_json(tree, aln, format).dump();
        }, "tree"_a, "alignment"_a, "format"_a = "modern")
        .def("to_legacy_json", [](const RELAXResult& res, const Tree& tree, const Alignment& aln) {
            return res.to_legacy_json(tree, aln).dump();
        }, "tree"_a, "alignment"_a)
        .def("to_modern_json", [](const RELAXResult& res, const Tree& tree, const Alignment& aln) {
            return res.to_modern_json(tree, aln).dump();
        }, "tree"_a, "alignment"_a);

    nb::class_<RELAXAnalyzer>(m, "RELAXAnalyzer")
        .def("__init__", [](RELAXAnalyzer* self, Tree tree, Alignment aln, const RELAXSettings& settings) {
            new (self) RELAXAnalyzer(std::move(tree), std::move(aln), settings);
        }, "tree"_a, "alignment"_a, "settings"_a = RELAXSettings{})
        .def_static("create", &RELAXAnalyzer::create, "tree"_a, "alignment"_a, "settings"_a = RELAXSettings{})
        .def("run", [](RELAXAnalyzer& relax, bool show_progress, bool force_progress) {
            return relax.run(nullptr, show_progress, force_progress);
        }, "show_progress"_a = false, "force_progress"_a = false);

    // FUBAR
    nb::enum_<FUBARMethod>(m, "FUBARMethod")
        .value("VariationalBayes", FUBARMethod::VariationalBayes)
        .value("CollapsedGibbs", FUBARMethod::CollapsedGibbs)
        .value("MetropolisHastings", FUBARMethod::MetropolisHastings);

    nb::class_<FUBARSettings>(m, "FUBARSettings")
        .def(nb::init<>())
        .def_rw("grid_points", &FUBARSettings::grid_points)
        .def_rw("method", &FUBARSettings::method)
        .def_rw("concentration", &FUBARSettings::concentration)
        .def_rw("posterior_threshold", &FUBARSettings::posterior_threshold)
        .def_rw("non_zero", &FUBARSettings::non_zero)
        .def_rw("tree_scaling", &FUBARSettings::tree_scaling)
        .def_rw("mcmc_chains", &FUBARSettings::mcmc_chains)
        .def_rw("mcmc_chain_length", &FUBARSettings::mcmc_chain_length)
        .def_rw("mcmc_burn_in", &FUBARSettings::mcmc_burn_in)
        .def_rw("mcmc_samples", &FUBARSettings::mcmc_samples);

    nb::class_<FUBARSiteResult>(m, "FUBARSiteResult")
        .def_ro("site", &FUBARSiteResult::site)
        .def_ro("alpha", &FUBARSiteResult::alpha)
        .def_ro("beta", &FUBARSiteResult::beta)
        .def_ro("beta_minus_alpha", &FUBARSiteResult::beta_minus_alpha)
        .def_ro("prob_neg", &FUBARSiteResult::prob_neg)
        .def_ro("prob_pos", &FUBARSiteResult::prob_pos)
        .def_ro("bayes_factor", &FUBARSiteResult::bayes_factor)
        .def_ro("classification", &FUBARSiteResult::classification)
        .def_ro("posterior_weights", &FUBARSiteResult::posterior_weights);

    nb::class_<FUBARResult>(m, "FUBARResult")
        .def_ro("site_results", &FUBARResult::site_results)
        .def_ro("grid", &FUBARResult::grid)
        .def_ro("posterior_grid_weights", &FUBARResult::posterior_grid_weights)
        .def_ro("positive_sites", &FUBARResult::positive_sites)
        .def_ro("purifying_sites", &FUBARResult::purifying_sites)
        .def_ro("expected_false_positives", &FUBARResult::expected_false_positives)
        .def_ro("false_positives_ci_95", &FUBARResult::false_positives_ci_95)
        .def_ro("weight_non_positive", &FUBARResult::weight_non_positive)
        .def_ro("time_total_sec", &FUBARResult::time_total_sec);

    nb::class_<FUBARAnalyzer>(m, "FUBARAnalyzer")
        .def(nb::init<Tree, Alignment, MG94Parameters, FUBARSettings>(),
             "tree"_a, "alignment"_a, "params"_a = MG94Parameters{}, "settings"_a = FUBARSettings{})
        .def_static("create_and_fit", [](Tree tree, Alignment aln, const FUBARSettings& settings) {
            return FUBARAnalyzer::create_and_fit(std::move(tree), std::move(aln), settings, nullptr);
        }, "tree"_a, "alignment"_a, "settings"_a = FUBARSettings{})
        .def("run", [](FUBARAnalyzer& fubar, bool verbose) {
            return fubar.run(nullptr, verbose);
        }, "verbose"_a = false)
        .def("to_modern_json", &FUBARAnalyzer::to_modern_json, "result"_a, "cli_cmd"_a = "hyphy3 fubar", "input_filepath"_a = "")
        .def("to_legacy_json", &FUBARAnalyzer::to_legacy_json, "result"_a, "input_filepath"_a = "");
    // PRIME
    nb::class_<AminoAcidProperty>(m, "AminoAcidProperty")
        .def_ro("name", &AminoAcidProperty::name)
        .def_ro("description", &AminoAcidProperty::description)
        .def_ro("mean", &AminoAcidProperty::mean)
        .def_ro("std_dev", &AminoAcidProperty::std_dev);

    nb::class_<PropertySiteInference>(m, "PropertySiteInference")
        .def_ro("property_name", &PropertySiteInference::property_name)
        .def_ro("lambda_", &PropertySiteInference::lambda)
        .def_ro("lrt", &PropertySiteInference::lrt)
        .def_ro("p_value", &PropertySiteInference::p_value)
        .def_ro("adjusted_p_value", &PropertySiteInference::adjusted_p_value)
        .def_ro("log_l_null", &PropertySiteInference::log_l_null)
        .def_ro("classification", &PropertySiteInference::classification);

    nb::class_<PRIMESettings>(m, "PRIMESettings")
        .def(nb::init<>())
        .def_rw("p_value_threshold", &PRIMESettings::p_value_threshold)
        .def_rw("stage2_p_value", &PRIMESettings::stage2_p_value)
        .def_rw("property_set", &PRIMESettings::property_set)
        .def_rw("custom_property_file", &PRIMESettings::custom_property_file)
        .def_rw("shrinkage_lambda", &PRIMESettings::shrinkage_lambda)
        .def_rw("max_bound", &PRIMESettings::max_bound)
        .def_rw("full_model", &PRIMESettings::full_model)
        .def_rw("fit_global", &PRIMESettings::fit_global)
        .def_rw("impute_states", &PRIMESettings::impute_states)
        .def_rw("test_branches", &PRIMESettings::test_branches)
        .def_rw("test_branch_regex", &PRIMESettings::test_branch_regex)
        .def_rw("test_branch_names", &PRIMESettings::test_branch_names);

    nb::class_<PRIMESiteResult>(m, "PRIMESiteResult")
        .def_ro("site_index", &PRIMESiteResult::site_index)
        .def_ro("alpha", &PRIMESiteResult::alpha)
        .def_ro("beta", &PRIMESiteResult::beta)
        .def_ro("beta_nuisance", &PRIMESiteResult::beta_nuisance)
        .def_ro("fel_alpha", &PRIMESiteResult::fel_alpha)
        .def_ro("fel_beta", &PRIMESiteResult::fel_beta)
        .def_ro("log_l_alt", &PRIMESiteResult::log_l_alt)
        .def_ro("log_l_fel", &PRIMESiteResult::log_l_fel)
        .def_ro("omnibus_lrt", &PRIMESiteResult::omnibus_lrt)
        .def_ro("omnibus_p_value", &PRIMESiteResult::omnibus_p_value)
        .def_ro("omnibus_q_value", &PRIMESiteResult::omnibus_q_value)
        .def_ro("total_branch_length", &PRIMESiteResult::total_branch_length)
        .def_ro("properties", &PRIMESiteResult::properties)
        .def_ro("nonsyn_subs", &PRIMESiteResult::nonsyn_subs)
        .def_ro("unique_aa", &PRIMESiteResult::unique_aa)
        .def_ro("redundancy_ratio", &PRIMESiteResult::redundancy_ratio)
        .def("is_significant", &PRIMESiteResult::is_significant, "threshold"_a = 0.10);

    nb::class_<PRIMEGlobalResult>(m, "PRIMEGlobalResult")
        .def_ro("log_likelihood", &PRIMEGlobalResult::log_likelihood)
        .def_ro("aicc", &PRIMEGlobalResult::aicc)
        .def_ro("alpha", &PRIMEGlobalResult::alpha)
        .def_ro("beta", &PRIMEGlobalResult::beta)
        .def_ro("lambdas", &PRIMEGlobalResult::lambdas)
        .def_ro("property_names", &PRIMEGlobalResult::property_names)
        .def_ro("baseline_tree_length", &PRIMEGlobalResult::baseline_tree_length);

    nb::class_<PRIMEResult>(m, "PRIMEResult")
        .def_ro("settings", &PRIMEResult::settings)
        .def_ro("properties", &PRIMEResult::properties)
        .def_ro("global_fit", &PRIMEResult::global_fit)
        .def_ro("site_results", &PRIMEResult::site_results)
        .def_ro("total_sites", &PRIMEResult::total_sites)
        .def_ro("invariant_sites", &PRIMEResult::invariant_sites)
        .def_ro("significant_sites", &PRIMEResult::significant_sites)
        .def_ro("conserved_counts", &PRIMEResult::conserved_counts)
        .def_ro("changing_counts", &PRIMEResult::changing_counts)
        .def_ro("runtime_seconds", &PRIMEResult::runtime_seconds)
        .def("to_modern_json", [](const PRIMEResult& res, const Tree& tree, const Alignment& aln, const std::string& cli_cmd, const std::string& input_filepath) {
            return res.to_modern_json(tree, aln, cli_cmd, input_filepath).dump();
        }, "tree"_a, "alignment"_a, "cli_cmd"_a = "hyphy3 prime", "input_filepath"_a = "")
        .def("to_legacy_json", [](const PRIMEResult& res, const Tree& tree, const Alignment& aln, const std::string& input_filepath) {
            return res.to_legacy_json(tree, aln, input_filepath).dump();
        }, "tree"_a, "alignment"_a, "input_filepath"_a = "");

    nb::class_<PRIMEAnalyzer>(m, "PRIMEAnalyzer")
        .def(nb::init<Tree, Alignment, MG94Parameters, PRIMESettings>(),
             "tree"_a, "alignment"_a, "params"_a = MG94Parameters{}, "settings"_a = PRIMESettings{})
        .def_static("create_and_fit", [](Tree tree, Alignment aln, const PRIMESettings& settings) {
            return PRIMEAnalyzer::create_and_fit(std::move(tree), std::move(aln), settings);
        }, "tree"_a, "alignment"_a, "settings"_a = PRIMESettings{})
        .def("run", [](
            PRIMEAnalyzer& prime,
            bool show_progress,
            bool force_progress,
            std::function<void(size_t, const PRIMESiteResult&)> site_cb,
            std::function<void(size_t, size_t, size_t, const PRIMESiteResult&)> live_site_cb
        ) {
            return prime.run(nullptr, show_progress, force_progress, site_cb, live_site_cb);
        }, "show_progress"_a = false, "force_progress"_a = false,
           "site_callback"_a = nb::none(), "live_site_callback"_a = nb::none());

    // Differentiable Engine: Inside-Outside Likelihood and Analytical Adjoint Gradients!
    m.def("compute_branch_length_gradients", [](
        const Tree& tree,
        const Alignment& aln,
        const MG94Parameters& params,
        const std::vector<Scalar>& branch_lengths
    ) -> std::pair<Scalar, std::vector<Scalar>> {
        return LikelihoodEngine::compute_branch_length_gradients(tree, aln, params, branch_lengths);
    }, "tree"_a, "aln"_a, "params"_a, "branch_lengths"_a = std::vector<Scalar>{}, "Compute log-likelihood and exact analytical branch gradients d ln L / d t_b via Inside-Outside");

    // Neighbor Joining
    m.def("compute_tn93_distances", &NeighborJoining::compute_tn93_distances, "aln"_a);
    m.def("infer_nj_tree", &NeighborJoining::infer, "aln"_a);

    // GARD Settings, Partitions, Models, Results & Analyzer
    nb::class_<GARDSettings>(m, "GARDSettings")
        .def(nb::init<>())
        .def_rw("input_file", &GARDSettings::input_file)
        .def_rw("output_file", &GARDSettings::output_file)
        .def_rw("run_mode", &GARDSettings::run_mode)
        .def_rw("max_breakpoints", &GARDSettings::max_breakpoints)
        .def_rw("population_size", &GARDSettings::population_size)
        .def_rw("mutation_rate", &GARDSettings::mutation_rate)
        .def_rw("rate_small_shifts", &GARDSettings::rate_small_shifts)
        .def_rw("caic_diversity_threshold", &GARDSettings::caic_diversity_threshold)
        .def_rw("max_stagnant_generations", &GARDSettings::max_stagnant_generations)
        .def_rw("caic_improvement_threshold", &GARDSettings::caic_improvement_threshold)
        .def_rw("quiet", &GARDSettings::quiet)
        .def_rw("dynamic_chart", &GARDSettings::dynamic_chart)
        .def_rw("force_tty", &GARDSettings::force_tty)
        .def_rw("legacy_json", &GARDSettings::legacy_json);

    nb::class_<GARDPartition>(m, "GARDPartition")
        .def_ro("start_site", &GARDPartition::start_site)
        .def_ro("end_site", &GARDPartition::end_site)
        .def_ro("tree", &GARDPartition::tree)
        .def_ro("log_likelihood", &GARDPartition::log_likelihood);

    nb::class_<GARDModel>(m, "GARDModel")
        .def_ro("breakpoints", &GARDModel::breakpoints)
        .def_ro("log_likelihood", &GARDModel::log_likelihood)
        .def_ro("aicc", &GARDModel::aicc)
        .def_ro("partitions", &GARDModel::partitions);

    nb::class_<GARDResult>(m, "GARDResult")
        .def_ro("baseline_aicc", &GARDResult::baseline_aicc)
        .def_ro("baseline_log_l", &GARDResult::baseline_log_l)
        .def_ro("baseline_tree", &GARDResult::baseline_tree)
        .def_ro("num_sites", &GARDResult::num_sites)
        .def_ro("num_taxa", &GARDResult::num_taxa)
        .def_ro("potential_breakpoints", &GARDResult::potential_breakpoints)
        .def_ro("min_partition_size", &GARDResult::min_partition_size)
        .def_ro("variable_sites", &GARDResult::variable_sites)
        .def_ro("best_model", &GARDResult::best_model)
        .def_ro("single_tree_aicc", &GARDResult::single_tree_aicc)
        .def_ro("topological_incongruence", &GARDResult::topological_incongruence)
        .def_ro("site_breakpoint_support", &GARDResult::site_breakpoint_support)
        .def_ro("improvements", &GARDResult::improvements)
        .def_ro("total_models_evaluated", &GARDResult::total_models_evaluated)
        .def_ro("time_elapsed", &GARDResult::time_elapsed)
        .def("to_modern_json", [](const GARDResult& res) { return res.to_json().dump(); })
        .def("to_legacy_json", [](const GARDResult& res) { return res.to_legacy_json().dump(); });

    nb::class_<GARDAnalyzer>(m, "GARDAnalyzer")
        .def(nb::init<const Alignment&, GARDSettings>(), "alignment"_a, "settings"_a = GARDSettings{})
        .def("run", &GARDAnalyzer::run);

    // SLAC Settings, Results & Analyzer
    nb::class_<SLACSettings>(m, "SLACSettings")
        .def(nb::init<>())
        .def_rw("p_value_threshold", &SLACSettings::p_value_threshold)
        .def_rw("samples", &SLACSettings::samples)
        .def_rw("test_branches", &SLACSettings::test_branches)
        .def_rw("test_branch_regex", &SLACSettings::test_branch_regex)
        .def_rw("test_branch_names", &SLACSettings::test_branch_names)
        .def_rw("full_model", &SLACSettings::full_model)
        .def_rw("quiet", &SLACSettings::quiet)
        .def_rw("all_sites", &SLACSettings::all_sites)
        .def_rw("random_seed", &SLACSettings::random_seed);

    nb::class_<SLACSiteResult>(m, "SLACSiteResult")
        .def_ro("site_index", &SLACSiteResult::site_index)
        .def_ro("ES", &SLACSiteResult::ES)
        .def_ro("EN", &SLACSiteResult::EN)
        .def_ro("S", &SLACSiteResult::S)
        .def_ro("N", &SLACSiteResult::N)
        .def_ro("P_S", &SLACSiteResult::P_S)
        .def_ro("dS", &SLACSiteResult::dS)
        .def_ro("dN", &SLACSiteResult::dN)
        .def_ro("dN_minus_dS", &SLACSiteResult::dN_minus_dS)
        .def_ro("p_positive", &SLACSiteResult::p_positive)
        .def_ro("p_negative", &SLACSiteResult::p_negative)
        .def_ro("total_branch_length", &SLACSiteResult::total_branch_length)
        .def_ro("has_dS", &SLACSiteResult::has_dS)
        .def_ro("has_dN", &SLACSiteResult::has_dN)
        .def_ro("has_dN_minus_dS", &SLACSiteResult::has_dN_minus_dS)
        .def_ro("quantiles_2_5", &SLACSiteResult::quantiles_2_5)
        .def_ro("quantiles_median", &SLACSiteResult::quantiles_median)
        .def_ro("quantiles_97_5", &SLACSiteResult::quantiles_97_5);

    nb::class_<SLACBranchResult>(m, "SLACBranchResult")
        .def_ro("branch_name", &SLACBranchResult::branch_name)
        .def_ro("branch_length", &SLACBranchResult::branch_length)
        .def_ro("ES", &SLACBranchResult::ES)
        .def_ro("EN", &SLACBranchResult::EN)
        .def_ro("S", &SLACBranchResult::S)
        .def_ro("N", &SLACBranchResult::N)
        .def_ro("P_S", &SLACBranchResult::P_S)
        .def_ro("dS", &SLACBranchResult::dS)
        .def_ro("dN", &SLACBranchResult::dN)
        .def_ro("dN_minus_dS", &SLACBranchResult::dN_minus_dS)
        .def_ro("p_positive", &SLACBranchResult::p_positive)
        .def_ro("p_negative", &SLACBranchResult::p_negative)
        .def_ro("has_dS", &SLACBranchResult::has_dS)
        .def_ro("has_dN", &SLACBranchResult::has_dN)
        .def_ro("has_dN_minus_dS", &SLACBranchResult::has_dN_minus_dS);

    nb::class_<SLACResult>(m, "SLACResult")
        .def_ro("settings", &SLACResult::settings)
        .def_ro("baseline_log_l", &SLACResult::baseline_log_l)
        .def_ro("baseline_aicc", &SLACResult::baseline_aicc)
        .def_ro("gtr_log_l", &SLACResult::gtr_log_l)
        .def_ro("gtr_aicc", &SLACResult::gtr_aicc)
        .def_ro("global_omega", &SLACResult::global_omega)
        .def_ro("total_tree_length", &SLACResult::total_tree_length)
        .def_ro("total_tested_tree_length", &SLACResult::total_tested_tree_length)
        .def_ro("num_tested_branches", &SLACResult::num_tested_branches)
        .def_ro("num_background_branches", &SLACResult::num_background_branches)
        .def_ro("tested_branch_names", &SLACResult::tested_branch_names)
        .def_ro("background_branch_names", &SLACResult::background_branch_names)
        .def_ro("resolved_sites", &SLACResult::resolved_sites)
        .def_ro("averaged_sites", &SLACResult::averaged_sites)
        .def_ro("resolved_branches", &SLACResult::resolved_branches)
        .def_ro("averaged_branches", &SLACResult::averaged_branches)
        .def_ro("total_sites", &SLACResult::total_sites)
        .def_ro("diversifying_sites", &SLACResult::diversifying_sites)
        .def_ro("purifying_sites", &SLACResult::purifying_sites)
        .def_ro("neutral_sites", &SLACResult::neutral_sites)
        .def_ro("runtime_seconds", &SLACResult::runtime_seconds)
        .def("to_modern_json", [](const SLACResult& res, const std::string& input_filepath, const std::string& tree_string) {
            return res.to_modern_json(input_filepath, tree_string).dump();
        }, "input_filepath"_a = "", "tree_string"_a = "")
        .def("to_legacy_json", [](const SLACResult& res, const std::string& input_filepath, const std::string& tree_string) {
            return res.to_legacy_json(input_filepath, tree_string).dump();
        }, "input_filepath"_a = "", "tree_string"_a = "");

    nb::class_<SLACAnalyzer>(m, "SLACAnalyzer")
        .def(nb::init<Tree, Alignment, MG94Parameters, SLACSettings>(),
             "tree"_a, "alignment"_a, "params"_a = MG94Parameters{}, "settings"_a = SLACSettings{})
        .def_static("create_and_fit", [](Tree tree, Alignment aln, const SLACSettings& settings) {
            return SLACAnalyzer::create_and_fit(std::move(tree), std::move(aln), settings);
        }, "tree"_a, "alignment"_a, "settings"_a = SLACSettings{})
        .def("run", [](SLACAnalyzer& slac, bool show_progress, bool force_progress) {
            return slac.run(nullptr, show_progress, force_progress);
        }, "show_progress"_a = false, "force_progress"_a = false);

    // Directional FEL Results & Analyzer
    nb::class_<DirectionalFELSettings>(m, "DirectionalFELSettings")
        .def(nb::init<>())
        .def_rw("p_value_threshold", &DirectionalFELSettings::p_value_threshold)
        .def_rw("heterogeneity_threshold", &DirectionalFELSettings::heterogeneity_threshold)
        .def_rw("root_strategy", &DirectionalFELSettings::root_strategy)
        .def_rw("num_regimes", &DirectionalFELSettings::num_regimes)
        .def_rw("full_model", &DirectionalFELSettings::full_model);

    nb::class_<TaxonComposition>(m, "TaxonComposition")
        .def_ro("name", &TaxonComposition::name)
        .def_ro("total_codons", &TaxonComposition::total_codons)
        .def_ro("gc_overall", &TaxonComposition::gc_overall)
        .def_ro("gc_cp3", &TaxonComposition::gc_cp3);

    nb::class_<CompositionTestResult>(m, "CompositionTestResult")
        .def_ro("chi2_overall", &CompositionTestResult::chi2_overall)
        .def_ro("df_overall", &CompositionTestResult::df_overall)
        .def_ro("p_value_overall", &CompositionTestResult::p_value_overall)
        .def_ro("chi2_cp3", &CompositionTestResult::chi2_cp3)
        .def_ro("df_cp3", &CompositionTestResult::df_cp3)
        .def_ro("p_value_cp3", &CompositionTestResult::p_value_cp3)
        .def_ro("is_heterogeneous", &CompositionTestResult::is_heterogeneous)
        .def_ro("taxa", &CompositionTestResult::taxa);

    nb::class_<RootCandidate>(m, "RootCandidate")
        .def_ro("rank", &RootCandidate::rank)
        .def_ro("branch_name", &RootCandidate::branch_name)
        .def_ro("is_terminal", &RootCandidate::is_terminal)
        .def_ro("log_l", &RootCandidate::log_l)
        .def_ro("delta_log_l", &RootCandidate::delta_log_l)
        .def_ro("weight", &RootCandidate::weight);

    nb::class_<AutoRootResult>(m, "AutoRootResult")
        .def_ro("chosen_branch_name", &AutoRootResult::chosen_branch_name)
        .def_ro("chosen_log_l", &AutoRootResult::chosen_log_l)
        .def_ro("chosen_weight", &AutoRootResult::chosen_weight)
        .def_ro("candidates", &AutoRootResult::candidates);

    nb::class_<DirectionalRegime>(m, "DirectionalRegime")
        .def_ro("id", &DirectionalRegime::id)
        .def_ro("name", &DirectionalRegime::name)
        .def_ro("psi", &DirectionalRegime::psi)
        .def_ro("mean_gc", &DirectionalRegime::mean_gc)
        .def_ro("branch_count", &DirectionalRegime::branch_count);

    nb::class_<DirectionalSiteResult>(m, "DirectionalSiteResult")
        .def_ro("site_index", &DirectionalSiteResult::site_index)
        .def_ro("alpha", &DirectionalSiteResult::alpha)
        .def_ro("beta", &DirectionalSiteResult::beta)
        .def_ro("alpha_null", &DirectionalSiteResult::alpha_null)
        .def_ro("lrt", &DirectionalSiteResult::lrt)
        .def_ro("p_value", &DirectionalSiteResult::p_value)
        .def_ro("total_branch_length", &DirectionalSiteResult::total_branch_length);

    nb::class_<DirectionalFELAnalyzer>(m, "DirectionalFELAnalyzer")
        .def(nb::init<Tree, Alignment, MG94Parameters, DirectionalFELSettings>(),
             "tree"_a, "alignment"_a, "params"_a = MG94Parameters{}, "settings"_a = DirectionalFELSettings{})
        .def_static("run_composition_test", &DirectionalFELAnalyzer::run_composition_test, "alignment"_a, "threshold"_a = 0.05)
        .def_static("create_and_fit", [](Tree tree, Alignment aln, const DirectionalFELSettings& settings) {
            return DirectionalFELAnalyzer::create_and_fit(std::move(tree), std::move(aln), settings);
        }, "tree"_a, "alignment"_a, "settings"_a = DirectionalFELSettings{})
        .def("run", &DirectionalFELAnalyzer::run, "show_progress"_a = false, "force_progress"_a = false)
        .def("to_modern_json", [](const DirectionalFELAnalyzer& dfel, const std::string& input_filepath, const std::string& tree_string) {
            return dfel.to_modern_json(input_filepath, tree_string).dump();
        }, "input_filepath"_a = "", "tree_string"_a = "")
        .def("to_legacy_json", [](const DirectionalFELAnalyzer& dfel, const std::string& input_filepath, const std::string& tree_string) {
            return dfel.to_legacy_json(input_filepath, tree_string).dump();
        }, "input_filepath"_a = "", "tree_string"_a = "")
        .def_ro("comp_test", &DirectionalFELAnalyzer::comp_test)
        .def_ro("root_result", &DirectionalFELAnalyzer::root_result)
        .def_ro("regimes", &DirectionalFELAnalyzer::regimes)
        .def_ro("branch_regime_ids", &DirectionalFELAnalyzer::branch_regime_ids)
        .def_ro("branch_regime_weights", &DirectionalFELAnalyzer::branch_regime_weights)
        .def_ro("site_results", &DirectionalFELAnalyzer::site_results);

    // Ohmic Phylogenetic Networks & Wheatstone Pruning
    nb::enum_<OhmicScaleMode>(m, "OhmicScaleMode")
        .value("Auto", OhmicScaleMode::Auto)
        .value("ExactDense", OhmicScaleMode::ExactDense)
        .value("ScalableWoodbury", OhmicScaleMode::ScalableWoodbury);

    nb::class_<OhmicSettings>(m, "OhmicSettings")
        .def(nb::init<>())
        .def_rw("scale_mode", &OhmicSettings::scale_mode)
        .def_rw("dense_switch_threshold", &OhmicSettings::dense_switch_threshold)
        .def_rw("rank_ratio_threshold", &OhmicSettings::rank_ratio_threshold)
        .def_rw("woodbury_rank", &OhmicSettings::woodbury_rank)
        .def_rw("min_shunt_stress", &OhmicSettings::min_shunt_stress)
        .def_rw("max_bridges", &OhmicSettings::max_bridges)
        .def_rw("collapse_ohmic_polytomies", &OhmicSettings::collapse_ohmic_polytomies)
        .def_rw("include_tertiary", &OhmicSettings::include_tertiary)
        .def_rw("use_hamming_distance", &OhmicSettings::use_hamming_distance)
        .def_rw("fixed_sigma", &OhmicSettings::fixed_sigma)
        .def_rw("epsilon_0", &OhmicSettings::epsilon_0);

    nb::class_<WheatstoneBridge>(m, "WheatstoneBridge")
        .def(nb::init<>())
        .def_rw("hub_node_id", &WheatstoneBridge::hub_node_id)
        .def_rw("bridge_node_id", &WheatstoneBridge::bridge_node_id)
        .def_rw("children", &WheatstoneBridge::children)
        .def_rw("weights", &WheatstoneBridge::weights)
        .def_rw("branch_lengths", &WheatstoneBridge::branch_lengths)
        .def_rw("ohmic_support", &WheatstoneBridge::ohmic_support)
        .def_rw("shunt_stress", &WheatstoneBridge::shunt_stress);

    nb::class_<OhmicNetwork>(m, "OhmicNetwork")
        .def(nb::init<>())
        .def_ro("backbone_tree", &OhmicNetwork::backbone_tree)
        .def_ro("collapsed_tree", &OhmicNetwork::collapsed_tree)
        .def_rw("bridges", &OhmicNetwork::bridges)
        .def_ro("node_to_bridge_idx", &OhmicNetwork::node_to_bridge_idx)
        .def_ro("active_mode", &OhmicNetwork::active_mode)
        .def_ro("num_super_nodes", &OhmicNetwork::num_super_nodes)
        .def_ro("num_collapsed_polytomies", &OhmicNetwork::num_collapsed_polytomies)
        .def_ro("D_tn93", &OhmicNetwork::D_tn93)
        .def_ro("Omega", &OhmicNetwork::Omega)
        .def_ro("P_edge", &OhmicNetwork::P_edge)
        .def_ro("foster_sum", &OhmicNetwork::foster_sum)
        .def_ro("inv_aug_degrees", &OhmicNetwork::inv_aug_degrees)
        .def_ro("Z_resistance_coords", &OhmicNetwork::Z_resistance_coords)
        .def("effective_resistance", &OhmicNetwork::effective_resistance, "i"_a, "j"_a)
        .def_static("infer", [](const Alignment& aln, const OhmicSettings& settings) {
            return OhmicNetwork::infer(aln, settings, nullptr);
        }, "alignment"_a, "settings"_a = OhmicSettings{});

    nb::class_<OhmicAnalysisSettings>(m, "OhmicAnalysisSettings")
        .def(nb::init<>())
        .def_rw("network_settings", &OhmicAnalysisSettings::network_settings)
        .def_rw("use_codon_model", &OhmicAnalysisSettings::use_codon_model)
        .def_rw("optimize_bridge_weights", &OhmicAnalysisSettings::optimize_bridge_weights)
        .def_rw("optimize_branch_lengths", &OhmicAnalysisSettings::optimize_branch_lengths)
        .def_rw("quiet", &OhmicAnalysisSettings::quiet)
        .def_rw("input_file", &OhmicAnalysisSettings::input_file)
        .def_rw("output_file", &OhmicAnalysisSettings::output_file);

    nb::class_<OhmicBridgeSummary>(m, "OhmicBridgeSummary")
        .def_ro("bridge_index", &OhmicBridgeSummary::bridge_index)
        .def_ro("hub_name", &OhmicBridgeSummary::hub_name)
        .def_ro("bridge_name", &OhmicBridgeSummary::bridge_name)
        .def_ro("child_names", &OhmicBridgeSummary::child_names)
        .def_ro("initial_weights", &OhmicBridgeSummary::initial_weights)
        .def_ro("fitted_weights", &OhmicBridgeSummary::fitted_weights)
        .def_ro("fitted_branch_lengths", &OhmicBridgeSummary::fitted_branch_lengths)
        .def_ro("ohmic_support", &OhmicBridgeSummary::ohmic_support)
        .def_ro("shunt_stress", &OhmicBridgeSummary::shunt_stress);

    nb::class_<OhmicResult>(m, "OhmicResult")
        .def_ro("network", &OhmicResult::network)
        .def_ro("num_taxa", &OhmicResult::num_taxa)
        .def_ro("num_sites", &OhmicResult::num_sites)
        .def_ro("num_patterns", &OhmicResult::num_patterns)
        .def_ro("num_super_nodes", &OhmicResult::num_super_nodes)
        .def_ro("num_bridges", &OhmicResult::num_bridges)
        .def_ro("num_collapsed_polytomies", &OhmicResult::num_collapsed_polytomies)
        .def_ro("foster_sum", &OhmicResult::foster_sum)
        .def_ro("single_tree_log_l", &OhmicResult::single_tree_log_l)
        .def_ro("single_tree_aicc", &OhmicResult::single_tree_aicc)
        .def_ro("ohmic_log_l", &OhmicResult::ohmic_log_l)
        .def_ro("ohmic_aicc", &OhmicResult::ohmic_aicc)
        .def_ro("delta_log_l", &OhmicResult::delta_log_l)
        .def_ro("delta_aicc", &OhmicResult::delta_aicc)
        .def_ro("bridge_summaries", &OhmicResult::bridge_summaries)
        .def_ro("runtime_seconds", &OhmicResult::runtime_seconds)
        .def("to_modern_json", [](const OhmicResult& res, const std::string& filepath) {
            return res.to_modern_json(filepath).dump();
        }, "filepath"_a = "")
        .def("to_legacy_json", [](const OhmicResult& res, const std::string& filepath) {
            return res.to_legacy_json(filepath).dump();
        }, "filepath"_a = "");

    nb::class_<OhmicAnalyzer>(m, "OhmicAnalyzer")
        .def(nb::init<Alignment, OhmicAnalysisSettings>(), "alignment"_a, "settings"_a = OhmicAnalysisSettings{})
        .def("run", [](OhmicAnalyzer& analyzer) {
            return analyzer.run(nullptr);
        });
}
