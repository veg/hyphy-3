#pragma once

#include "hyphy/analyses/joint_align.hpp"
#include "hyphy/analyses/directional_fel.hpp"
#include "hyphy/analyses/prime.hpp"
#include "hyphy/core/console.hpp"
#include "hyphy/core/stats.hpp"
#include "nlohmann/json.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <optional>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace hyphy::analyses {

using namespace hyphy::core;

// ============================================================================
// One-Button Evolutionary Characterization ("Characterize / Auto / Omni")
//   From raw unaligned coding sequences (or an existing MSA) to complete
//   alignment-marginalized, reticulation-aware, compositionally audited
//   gene-, branch-, clade-, site-, and biochemical property-level selection
//   inference via a single shared Inside-Outside sufficiency tensor,
//   two-stage active-set continuous polish, and hierarchical Simes/BH FDR.
// ============================================================================

struct CharacterizeSettings {
    Scalar p_value_threshold = 0.10;
    Scalar fdr_threshold = 0.10;
    Scalar active_set_screen_p = 0.25; // Decision 1: Two-stage active-set continuous polish gate
    bool from_unaligned = false;
    bool enable_coupled_spr = true;
    bool enable_segment_routing = true;
    bool enable_coevolution = false;   // Tier 3 (Opt-in): Epistatic co-substitution pairs (--coevolution)
    bool emit_legacy_jsons = false;    // Decision 7: Optional Vision/Datamonkey legacy JSON bundle
    bool quiet = false;
    std::string genetic_code_name = "Universal";
    std::string foreground_tag = "";   // Decision 5: Optional user foreground tag (auto-discovered if empty)
};

struct SpatialDomainSegment {
    size_t start_codon = 0; // 0-indexed
    size_t end_codon = 0;   // 0-indexed inclusive
    std::string domain_type = "Purifying Core"; // "Purifying Core", "Permissive / Surface", "Hypervariable / IDR"
    Scalar mean_alpha = 1.0;
    Scalar mean_beta = 0.1;
    Scalar mean_omega = 0.1;
    Scalar mean_gap_fraction = 0.0;
    Scalar mean_n_eff = 1.0;
};

struct SynonymousConstraintIsland {
    size_t start_codon = 0; // 0-indexed
    size_t end_codon = 0;   // 0-indexed inclusive
    Scalar mean_alpha = 0.1;
    Scalar gene_mean_alpha = 1.0;
    Scalar suppression_ratio = 0.1;
    Scalar z_score = 0.0;
    Scalar p_value = 1.0;
};

struct MutationalContextAlert {
    bool triggered = false;
    size_t num_selected_sites_evaluated = 0;
    size_t cpg_driven_selected_sites = 0;
    size_t apobec_driven_selected_sites = 0;
    Scalar cpg_enrichment_fold = 1.0;
    Scalar cpg_fisher_p = 1.0;
    std::vector<size_t> flagged_cpg_sites;    // 1-indexed site numbers
    std::vector<size_t> flagged_apobec_sites; // 1-indexed site numbers
    std::string alert_summary;
};

struct CoevolutionSitePair {
    size_t site_i = 0; // 0-indexed
    size_t site_j = 0; // 0-indexed
    Scalar shared_branches = 0.0;
    Scalar expected_shared = 0.0;
    Scalar cosine_similarity = 0.0;
    Scalar p_value = 1.0;
    Scalar bh_q_value = 1.0;
    std::string top_co_mutating_branches;
};

struct RecombinationBridgeSegment {
    int32_t node_u = -1;
    int32_t node_v = -1;
    std::string clade_u_name;
    std::string target_v_name;
    Scalar bridge_conductance = 0.0;
    size_t peak_start_codon = 0;
    size_t peak_end_codon = 0;
    Scalar segment_delta_log_l = 0.0;
    Scalar p_value = 1.0;
    bool routed_to_spr_tree = false; // Decision 3: True if segment p <= 0.01 routed to local SPR tree T_{u->v}
};

struct CladeMacroRegimeSummary {
    // Tier 2 Triggered: Dual-level autonomous macro-regime partition (or user {Foreground} contrast)
    bool triggered = false; // True if user supplied {Foreground} tags OR autonomous split achieves Delta BIC >= 10
    std::string partition_source = "None (Homogeneous Tree-Wide Regime)";
    std::string test_regime_label = "Focal Clade";
    std::string reference_regime_label = "Background";
    size_t num_test_branches = 0;
    size_t num_reference_branches = 0;
    Scalar clade_relax_k = 1.0;
    Scalar clade_relax_lrt = 0.0;
    Scalar clade_relax_delta_bic = 0.0;
    Scalar clade_relax_p_value = 1.0;
    std::string clade_relax_call = "Baseline"; // "Intensified", "Relaxed", or "Baseline"
    size_t contrast_fel_differing_sites = 0;
};

struct BranchEvolutionaryProfile {
    int32_t node_id = -1;
    std::string branch_name;
    std::string model_tag;
    bool is_leaf = false;
    Scalar branch_length_subs = 0.0;
    std::string compositional_regime = "Neutral";
    std::array<Scalar, 3> directional_weights = {0.0, 1.0, 0.0}; // Decision 4: Softmax weights [w_AT, w_Neu, w_GC]
    Scalar gc_cp3 = 0.50;
    std::string macro_regime = "Reference"; // Decision 5: "Test" vs "Reference"

    // aBSREL-style branch episodic diversifying selection (with Decision 1 active-set continuous polish)
    Scalar omega_bg = 0.5;
    Scalar omega_plus = 1.0;
    Scalar weight_plus = 0.0;
    Scalar absrel_lrt = 0.0;
    Scalar absrel_p_value = 1.0;
    Scalar absrel_holm_p = 1.0;
    Scalar absrel_bh_q = 1.0;

    // Autonomous RELAX selection intensity parameter K_b (omega_b = omega_gene^{K_b})
    Scalar relax_k = 1.0;
    Scalar relax_lrt = 0.0;
    Scalar relax_p_value = 1.0;
    Scalar relax_bh_q = 1.0;
    std::string relax_call = "Baseline"; // "Intensified", "Relaxed", or "Baseline"

    // Tier 1 (#2): Lineage-Specific Loss-of-Constraint / Pseudogenization (Core-Invasion) Detector
    Scalar core_nonsyn_Expected = 0.0;
    Scalar core_nonsyn_Observed = 0.0;
    Scalar lof_bayes_factor = 0.0;
    Scalar lof_p_value = 1.0;
    bool lof_flagged = false;
    std::vector<size_t> core_invaded_sites; // 1-indexed core sites mutated on this branch

    size_t episodic_site_count = 0; // Sites with Empirical Bayes EBF >= 10 on this branch
};

struct UnifiedSiteAtlasEntry {
    size_t site_index = 0; // 0-indexed (1-indexed in display)
    Scalar column_homology_posterior = 1.0;
    Scalar gap_fraction = 0.0;
    bool in_indel_window = false;
    size_t ensemble_paths_evaluated = 1; // Decision 2: K=3 POA path realizations if PP_c < 0.99
    std::string active_topology_label = "Backbone"; // Decision 3: "Backbone" or "SPR(u->v)"
    std::string spatial_domain = "Purifying Core";  // Tier 1 (#3): 1D spatial segment label
    bool in_synonymous_island = false;              // Tier 1 (#3): Inside a synonymous constraint island

    // Pervasive selection (Alignment-Marginalized FEL & Compositional-Drift-Adjusted dFEL)
    Scalar alpha = 1.0;
    Scalar beta_pervasive = 1.0;
    Scalar fel_lrt = 0.0;
    Scalar fel_p_value = 1.0;
    Scalar fel_bh_q = 1.0;
    Scalar dfel_lrt = 0.0;
    Scalar dfel_p_value = 1.0;
    Scalar dfel_bh_q = 1.0;

    // Episodic selection (Exact Inside-Outside MEME + Active-Set Continuous Polish + EBF)
    Scalar beta_minus = 0.5;
    Scalar beta_plus = 1.0;
    Scalar weight_plus = 0.0;
    Scalar meme_lrt = 0.0;
    Scalar meme_p_value = 1.0;
    Scalar meme_bh_q = 1.0;
    bool meme_polished = false;
    size_t num_episodic_branches = 0;
    std::vector<std::string> top_episodic_branches;
    Scalar max_branch_ebf = 0.0;

    // Decision 5 (Tier 2 Triggered): Clade / Macro-Regime Differential Selection (Contrast-FEL)
    Scalar beta_reference = 1.0;
    Scalar beta_test = 1.0;
    Scalar contrast_fel_lrt = 0.0;
    Scalar contrast_fel_p_value = 1.0;
    Scalar contrast_fel_bh_q = 1.0;

    // Tier 1 (#1): Model-Free Amino-Acid Exchangeability, Toggle & Convergence Classifier
    // (Replaces locked-in 5-property scales with empirical state-space concentration & homoplasy)
    Scalar expected_nonsyn_subs = 0.0;       // M_s = sum_b E[N_{s,b} | D]
    Scalar effective_aa_alphabet = 1.0;      // N_eff(s) = exp(H(pi_s^AA))
    Scalar expected_neutral_n_eff = 1.0;     // E[N_eff | M_s] under neutral codon branching
    Scalar pair_concentration = 0.0;         // Top unordered AA pair fraction of nonsyn transitions
    size_t max_convergent_transitions = 0;   // Max independent branches transitioning into the same target AA
    char dominant_target_aa = '-';
    Scalar exchangeability_stat = 0.0;       // Concentration / convergence surprise statistic
    Scalar exchangeability_p_value = 1.0;
    Scalar exchangeability_bh_q = 1.0;
    std::string exchangeability_mode = "Baseline"; // "Toggle (A<->B)", "Convergent (->A)", "Focused Exchange", "Broadly Diversifying"
    std::string key_property_driver;         // Summary string (e.g., "V<->I Toggle (6.2 subs, p=0.002)" or "Convergent ->K (4 br, p=0.004)")

    // Tier 2 (#4): Mutational context annotations on selected sites
    bool cpg_hypermutable_context = false;
    bool apobec_motif_context = false;

    // Decision 6: Hierarchical Simes-adjusted omnibus p-value & BH q-value + conflict-free label
    Scalar simes_omnibus_p = 1.0;
    Scalar simes_omnibus_q = 1.0;
    std::string classification = "Neutral";
    bool rigid_alignment_artifact_eliminated = false;
};

struct GeneWideSelectionSummary {
    // Decision 2: BUSTED-E (Alignment-Marginalized 3-class ω mixture + Dedicated Error-Sink Class ω_err)
    std::array<Scalar, 3> omega_classes = {0.1, 0.8, 1.5};
    std::array<Scalar, 3> omega_weights = {0.7, 0.25, 0.05};
    Scalar error_sink_omega = 50.0;
    Scalar error_sink_weight = 0.0;
    Scalar unconstrained_log_l = 0.0;
    Scalar null_log_l = 0.0;
    Scalar sink_null_log_l_gain = 0.0; // Null lnL with the error sink minus without it
    Scalar busted_lrt = 0.0;
    Scalar busted_p_value = 1.0;
    bool has_gene_wide_positive_selection = false;
};

struct CharacterizeResult {
    // Metadata
    std::string genetic_code_name = "Universal";
    std::string input_mode = "aligned"; // "unaligned" or "aligned"
    std::string tree_source = "de_novo_inferred"; // "de_novo_inferred", "user_supplied", or "embedded_in_alignment"
    Scalar p_value_threshold = 0.10;
    Scalar fdr_threshold = 0.10;

    // Stage I: Co-Estimated Alignment, Tree, Network & Objective Quality
    JointAlignResult joint_align;
    JointAlignEngine::AlignmentObjectiveEvaluation alignment_objective;
    std::vector<RecombinationBridgeSegment> recombination_segments;

    // Stage II: Compositional Drift & Auto-Rooting Audit
    CompositionTestResult composition_test;
    bool directional_drift_active = false;
    std::string chosen_root_branch;
    Scalar chosen_root_weight = 0.0;

    // Stage III: Unified Gene, Clade, Branch, Domain & Site Selection Atlas
    GeneWideSelectionSummary gene_wide;
    CladeMacroRegimeSummary clade_macro_regime;
    std::vector<SpatialDomainSegment> spatial_domains;
    std::vector<SynonymousConstraintIsland> synonymous_islands;
    MutationalContextAlert mutational_context_alert;
    std::vector<CoevolutionSitePair> coevolution_pairs;
    std::vector<BranchEvolutionaryProfile> branch_profiles;
    std::vector<UnifiedSiteAtlasEntry> site_atlas;

    // Counts
    size_t sites_pervasive_diversifying = 0;
    size_t sites_pervasive_purifying = 0;
    size_t sites_episodic_diversifying = 0;
    size_t sites_regime_differential = 0;
    size_t sites_property_constrained_or_shifting = 0; // Toggle / Convergent / Focused Exchange sites
    size_t sites_simes_fdr_significant = 0;
    size_t sites_simes_diversifying = 0;
    size_t sites_simes_purifying = 0;
    size_t branches_episodic_diversifying = 0;
    size_t branches_intensified = 0;
    size_t branches_relaxed = 0;
    size_t branches_loss_of_constraint = 0;

    double runtime_seconds = 0.0;
};

class CharacterizeEngine {
public:
    CharacterizeSettings settings;

    explicit CharacterizeEngine(CharacterizeSettings s = CharacterizeSettings{})
        : settings(std::move(s)) {}

    // Benjamini-Hochberg step-up FDR q-values from a vector of raw p-values
    static std::vector<Scalar> compute_bh_q_values(const std::vector<Scalar>& p_vals) {
        size_t M = p_vals.size();
        std::vector<Scalar> q_vals(M, 1.0);
        if (M == 0) return q_vals;

        std::vector<size_t> order(M);
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
            return p_vals[a] < p_vals[b];
        });

        Scalar running_min = 1.0;
        for (size_t idx = M; idx-- > 0;) {
            size_t orig = order[idx];
            Scalar p = std::clamp<Scalar>(p_vals[orig], 0.0, 1.0);
            Scalar q = (p * static_cast<Scalar>(M)) / static_cast<Scalar>(idx + 1);
            running_min = std::min(running_min, q);
            q_vals[orig] = std::min<Scalar>(1.0, running_min);
        }
        return q_vals;
    }

    // Simes (1986) combined p-value across a family of K test p-values at a single site
    static Scalar compute_simes_p_value(std::vector<Scalar> p_family) {
        if (p_family.empty()) return 1.0;
        std::sort(p_family.begin(), p_family.end());
        size_t K = p_family.size();
        Scalar best = 1.0;
        for (size_t j = 0; j < K; ++j) {
            Scalar cand = (static_cast<Scalar>(K) * std::clamp<Scalar>(p_family[j], 0.0, 1.0))
                        / static_cast<Scalar>(j + 1);
            best = std::min(best, cand);
        }
        return std::min<Scalar>(1.0, best);
    }

    CharacterizeResult run_from_unaligned_file(
        const std::string& filepath,
        std::shared_ptr<const GeneticCode> gcode = nullptr
    ) {
        if (!gcode) gcode = GeneticCode::universal();
        auto seqs = JointAlignEngine::load_unaligned_fasta(filepath, gcode);
        return run_from_unaligned(seqs, gcode);
    }

    CharacterizeResult run_from_unaligned(
        const std::vector<Sequence>& unaligned_seqs,
        std::shared_ptr<const GeneticCode> gcode = nullptr,
        std::optional<Tree> input_tree = std::nullopt
    ) {
        auto t0 = std::chrono::steady_clock::now();
        if (!gcode) gcode = GeneticCode::universal();

        JointAlignSettings jas;
        jas.p_value_threshold = settings.p_value_threshold;
        jas.enable_coupled_topology_search = settings.enable_coupled_spr;
        jas.run_marginalized_fel = true;
        jas.quiet = settings.quiet;

        JointAlignEngine jae(jas);
        auto ja_res = jae.run_from_unaligned(unaligned_seqs, gcode, input_tree);

        auto res = run_unified_tensor_characterization(std::move(ja_res));
        res.input_mode = "unaligned";
        res.tree_source = input_tree.has_value() ? "user_supplied" : "de_novo_inferred";
        res.genetic_code_name = gcode->name.empty() ? settings.genetic_code_name : gcode->name;
        res.p_value_threshold = settings.p_value_threshold;
        res.fdr_threshold = settings.fdr_threshold;
        auto t1 = std::chrono::steady_clock::now();
        res.runtime_seconds = std::chrono::duration<double>(t1 - t0).count();
        return res;
    }

    CharacterizeResult run_from_alignment(
        const Alignment& input_aln,
        std::optional<Tree> input_tree = std::nullopt
    ) {
        auto t0 = std::chrono::steady_clock::now();

        JointAlignSettings jas;
        jas.p_value_threshold = settings.p_value_threshold;
        jas.enable_coupled_topology_search = settings.enable_coupled_spr;
        jas.run_marginalized_fel = true;
        jas.quiet = settings.quiet;

        JointAlignEngine jae(jas);
        auto ja_res = jae.run(input_aln, input_tree);

        auto res = run_unified_tensor_characterization(std::move(ja_res));
        res.input_mode = "aligned";
        if (input_tree.has_value()) {
            res.tree_source = "user_supplied";
        } else if (!input_aln.embedded_tree_newick.empty()) {
            res.tree_source = "embedded_in_alignment";
        } else {
            res.tree_source = "de_novo_inferred";
        }
        res.genetic_code_name = (input_aln.code && !input_aln.code->name.empty())
                              ? input_aln.code->name : settings.genetic_code_name;
        res.p_value_threshold = settings.p_value_threshold;
        res.fdr_threshold = settings.fdr_threshold;
        auto t1 = std::chrono::steady_clock::now();
        res.runtime_seconds = std::chrono::duration<double>(t1 - t0).count();
        return res;
    }

    static std::string to_json(const CharacterizeResult& res) {
        using json = nlohmann::json;
        json j;
        j["$schema"] = "https://hyphy.org/schemas/v3/characterize.v3.schema.json";
        j["analysis"] = {
            {"id", "characterize"},
            {"name", "HyPhy 3 One-Button Evolutionary Characterization"},
            {"version", "3.0.0"},
            {"schema_version", "3.0.0"},
            {"genetic_code", res.genetic_code_name},
            {"input_mode", res.input_mode},
            {"tree_source", res.tree_source},
            {"p_value_threshold", res.p_value_threshold},
            {"fdr_threshold", res.fdr_threshold},
            {"runtime_seconds", res.runtime_seconds}
        };

        const auto& adm = res.joint_align.admission_report;
        json adm_entries = json::array();
        for (const auto& e : adm.entries) {
            adm_entries.push_back({
                {"name", e.name},
                {"status", e.status},
                {"reason", e.reason},
                {"frame_offset", e.frame_offset},
                {"reverse_complemented", e.reverse_complemented},
                {"raw_nt_length", e.raw_nt_length},
                {"codon_length", e.codon_length},
                {"internal_stops_frame0", e.internal_stops_frame0},
                {"max_homology_z", e.max_homology_z},
                {"max_aa_identity", e.max_aa_identity}
            });
        }
        j["input_admission_audit"] = {
            {"total_input_sequences", adm.total_input_sequences},
            {"admitted_sequences", adm.admitted_sequences},
            {"rescued_frame_or_strand", adm.rescued_frame_or_strand},
            {"quarantined_sequences", adm.quarantined_sequences},
            {"median_pairwise_z", adm.median_pairwise_z},
            {"median_aa_identity", adm.median_aa_identity},
            {"taxa", adm_entries}
        };

        j["alignment_and_network"] = {
            {"num_taxa", res.joint_align.map_alignment.num_taxa},
            {"num_codons", res.joint_align.map_alignment.num_codons},
            {"total_gap_codons", res.alignment_objective.total_gap_codons},
            {"mg94_sub_log_l", res.alignment_objective.mg94_sub_log_l},
            {"branch_indel_events", res.alignment_objective.branch_indel_events},
            {"branch_indel_codons", res.alignment_objective.branch_indel_codons},
            {"branch_indel_log_l", res.alignment_objective.branch_indel_log_l},
            {"option1_joint_log_l", res.alignment_objective.joint_log_l},
            {"option3_loto_cv_log_l", res.alignment_objective.loto_cv_log_l},
            {"option5_mean_column_posterior", res.alignment_objective.mean_column_posterior},
            {"option5_min_column_posterior", res.alignment_objective.min_column_posterior},
            {"coupled_spr_moves_accepted", res.joint_align.total_coupled_moves_accepted},
            {"wheatstone_bridges", res.joint_align.ohmic_network.bridges.size()},
            {"polytomy_busbars", res.joint_align.ohmic_network.num_collapsed_polytomies},
            {"tree_newick", res.joint_align.optimized_tree.to_newick()},
            {"extended_newick", res.joint_align.extended_newick}
        };

        json bridges_arr = json::array();
        for (const auto& br : res.recombination_segments) {
            bridges_arr.push_back({
                {"clade_u", br.clade_u_name},
                {"target_v", br.target_v_name},
                {"conductance", br.bridge_conductance},
                {"peak_start_codon", br.peak_start_codon + 1},
                {"peak_end_codon", br.peak_end_codon + 1},
                {"segment_delta_log_l", br.segment_delta_log_l},
                {"p_value", br.p_value},
                {"routed_to_spr_tree", br.routed_to_spr_tree}
            });
        }
        j["recombination_bridges"] = bridges_arr;

        j["compositional_audit"] = {
            {"is_heterogeneous", res.composition_test.is_heterogeneous},
            {"chi2_overall", res.composition_test.chi2_overall},
            {"p_value_overall", res.composition_test.p_value_overall},
            {"chi2_cp3", res.composition_test.chi2_cp3},
            {"p_value_cp3", res.composition_test.p_value_cp3},
            {"directional_drift_active", res.directional_drift_active},
            {"chosen_root_branch", res.chosen_root_branch},
            {"chosen_root_weight", res.chosen_root_weight}
        };

        j["gene_wide_selection"] = {
            {"omega_classes", res.gene_wide.omega_classes},
            {"omega_weights", res.gene_wide.omega_weights},
            {"error_sink_omega", res.gene_wide.error_sink_omega},
            {"error_sink_weight", res.gene_wide.error_sink_weight},
            {"unconstrained_log_l", res.gene_wide.unconstrained_log_l},
            {"null_log_l", res.gene_wide.null_log_l},
            {"lrt", res.gene_wide.busted_lrt},
            {"p_value", res.gene_wide.busted_p_value},
            {"has_positive_selection", res.gene_wide.has_gene_wide_positive_selection}
        };

        j["clade_macro_regime"] = {
            {"triggered", res.clade_macro_regime.triggered},
            {"partition_source", res.clade_macro_regime.partition_source},
            {"test_regime_label", res.clade_macro_regime.test_regime_label},
            {"reference_regime_label", res.clade_macro_regime.reference_regime_label},
            {"num_test_branches", res.clade_macro_regime.num_test_branches},
            {"num_reference_branches", res.clade_macro_regime.num_reference_branches},
            {"clade_relax_k", res.clade_macro_regime.clade_relax_k},
            {"clade_relax_lrt", res.clade_macro_regime.clade_relax_lrt},
            {"clade_relax_delta_bic", res.clade_macro_regime.clade_relax_delta_bic},
            {"clade_relax_p_value", res.clade_macro_regime.clade_relax_p_value},
            {"clade_relax_call", res.clade_macro_regime.clade_relax_call},
            {"contrast_fel_differing_sites", res.clade_macro_regime.contrast_fel_differing_sites}
        };

        json dom_arr = json::array();
        for (const auto& d : res.spatial_domains) {
            dom_arr.push_back({
                {"start_codon", d.start_codon + 1},
                {"end_codon", d.end_codon + 1},
                {"domain_type", d.domain_type},
                {"mean_alpha", d.mean_alpha},
                {"mean_beta", d.mean_beta},
                {"mean_omega", d.mean_omega},
                {"mean_gap_fraction", d.mean_gap_fraction},
                {"mean_n_eff", d.mean_n_eff}
            });
        }
        j["spatial_domains"] = dom_arr;

        json isl_arr = json::array();
        for (const auto& isl : res.synonymous_islands) {
            isl_arr.push_back({
                {"start_codon", isl.start_codon + 1},
                {"end_codon", isl.end_codon + 1},
                {"mean_alpha", isl.mean_alpha},
                {"gene_mean_alpha", isl.gene_mean_alpha},
                {"suppression_ratio", isl.suppression_ratio},
                {"z_score", isl.z_score},
                {"p_value", isl.p_value}
            });
        }
        j["synonymous_constraint_islands"] = isl_arr;

        j["mutational_context_alert"] = {
            {"triggered", res.mutational_context_alert.triggered},
            {"num_selected_sites_evaluated", res.mutational_context_alert.num_selected_sites_evaluated},
            {"cpg_driven_selected_sites", res.mutational_context_alert.cpg_driven_selected_sites},
            {"apobec_driven_selected_sites", res.mutational_context_alert.apobec_driven_selected_sites},
            {"cpg_enrichment_fold", res.mutational_context_alert.cpg_enrichment_fold},
            {"cpg_fisher_p", res.mutational_context_alert.cpg_fisher_p},
            {"flagged_cpg_sites", res.mutational_context_alert.flagged_cpg_sites},
            {"flagged_apobec_sites", res.mutational_context_alert.flagged_apobec_sites},
            {"alert_summary", res.mutational_context_alert.alert_summary}
        };

        json coev_arr = json::array();
        for (const auto& cp : res.coevolution_pairs) {
            coev_arr.push_back({
                {"site_i", cp.site_i + 1},
                {"site_j", cp.site_j + 1},
                {"shared_branches", cp.shared_branches},
                {"expected_shared", cp.expected_shared},
                {"cosine_similarity", cp.cosine_similarity},
                {"p_value", cp.p_value},
                {"bh_q_value", cp.bh_q_value},
                {"top_co_mutating_branches", cp.top_co_mutating_branches}
            });
        }
        j["coevolution_pairs"] = coev_arr;

        json br_arr = json::array();
        for (const auto& bp : res.branch_profiles) {
            br_arr.push_back({
                {"branch_name", bp.branch_name},
                {"model_tag", bp.model_tag},
                {"is_leaf", bp.is_leaf},
                {"branch_length_subs", bp.branch_length_subs},
                {"compositional_regime", bp.compositional_regime},
                {"directional_weights", bp.directional_weights},
                {"gc_cp3", bp.gc_cp3},
                {"macro_regime", bp.macro_regime},
                {"omega_bg", bp.omega_bg},
                {"omega_plus", bp.omega_plus},
                {"weight_plus", bp.weight_plus},
                {"absrel_lrt", bp.absrel_lrt},
                {"absrel_p_value", bp.absrel_p_value},
                {"absrel_holm_p", bp.absrel_holm_p},
                {"absrel_bh_q", bp.absrel_bh_q},
                {"relax_k", bp.relax_k},
                {"relax_lrt", bp.relax_lrt},
                {"relax_p_value", bp.relax_p_value},
                {"relax_bh_q", bp.relax_bh_q},
                {"relax_call", bp.relax_call},
                {"core_nonsyn_expected", bp.core_nonsyn_Expected},
                {"core_nonsyn_observed", bp.core_nonsyn_Observed},
                {"lof_bayes_factor", bp.lof_bayes_factor},
                {"lof_p_value", bp.lof_p_value},
                {"lof_flagged", bp.lof_flagged},
                {"core_invaded_sites", bp.core_invaded_sites},
                {"episodic_site_count", bp.episodic_site_count}
            });
        }
        j["branch_profiles"] = br_arr;

        json site_arr = json::array();
        for (const auto& s : res.site_atlas) {
            site_arr.push_back({
                {"site", s.site_index + 1},
                {"column_homology_posterior", s.column_homology_posterior},
                {"gap_fraction", s.gap_fraction},
                {"in_indel_window", s.in_indel_window},
                {"ensemble_paths_evaluated", s.ensemble_paths_evaluated},
                {"active_topology_label", s.active_topology_label},
                {"spatial_domain", s.spatial_domain},
                {"in_synonymous_island", s.in_synonymous_island},
                {"alpha", s.alpha},
                {"beta_pervasive", s.beta_pervasive},
                {"fel_lrt", s.fel_lrt},
                {"fel_p_value", s.fel_p_value},
                {"fel_bh_q", s.fel_bh_q},
                {"dfel_lrt", s.dfel_lrt},
                {"dfel_p_value", s.dfel_p_value},
                {"dfel_bh_q", s.dfel_bh_q},
                {"beta_minus", s.beta_minus},
                {"beta_plus", s.beta_plus},
                {"weight_plus", s.weight_plus},
                {"meme_lrt", s.meme_lrt},
                {"meme_p_value", s.meme_p_value},
                {"meme_bh_q", s.meme_bh_q},
                {"meme_polished", s.meme_polished},
                {"num_episodic_branches", s.num_episodic_branches},
                {"top_episodic_branches", s.top_episodic_branches},
                {"max_branch_ebf", s.max_branch_ebf},
                {"beta_reference", s.beta_reference},
                {"beta_test", s.beta_test},
                {"contrast_fel_lrt", s.contrast_fel_lrt},
                {"contrast_fel_p_value", s.contrast_fel_p_value},
                {"contrast_fel_bh_q", s.contrast_fel_bh_q},
                {"expected_nonsyn_subs", s.expected_nonsyn_subs},
                {"effective_aa_alphabet", s.effective_aa_alphabet},
                {"expected_neutral_n_eff", s.expected_neutral_n_eff},
                {"pair_concentration", s.pair_concentration},
                {"max_convergent_transitions", s.max_convergent_transitions},
                {"dominant_target_aa", std::string(1, s.dominant_target_aa)},
                {"exchangeability_stat", s.exchangeability_stat},
                {"exchangeability_p_value", s.exchangeability_p_value},
                {"exchangeability_bh_q", s.exchangeability_bh_q},
                {"exchangeability_mode", s.exchangeability_mode},
                {"key_property_driver", s.key_property_driver},
                {"cpg_hypermutable_context", s.cpg_hypermutable_context},
                {"apobec_motif_context", s.apobec_motif_context},
                {"simes_omnibus_p", s.simes_omnibus_p},
                {"simes_omnibus_q", s.simes_omnibus_q},
                {"classification", s.classification},
                {"artifact_eliminated", s.rigid_alignment_artifact_eliminated}
            });
        }
        j["site_atlas"] = site_arr;

        j["summary_counts"] = {
            {"sites_pervasive_diversifying", res.sites_pervasive_diversifying},
            {"sites_pervasive_purifying", res.sites_pervasive_purifying},
            {"sites_episodic_diversifying", res.sites_episodic_diversifying},
            {"sites_regime_differential", res.sites_regime_differential},
            {"sites_exchangeability_or_convergent", res.sites_property_constrained_or_shifting},
            {"sites_simes_fdr_significant", res.sites_simes_fdr_significant},
            {"sites_simes_diversifying", res.sites_simes_diversifying},
            {"sites_simes_purifying", res.sites_simes_purifying},
            {"spatial_domains_detected", res.spatial_domains.size()},
            {"synonymous_constraint_islands", res.synonymous_islands.size()},
            {"branches_episodic_diversifying", res.branches_episodic_diversifying},
            {"branches_intensified", res.branches_intensified},
            {"branches_relaxed", res.branches_relaxed},
            {"branches_loss_of_constraint", res.branches_loss_of_constraint},
            {"false_positive_sites_eliminated", res.joint_align.false_positive_sites_eliminated},
            {"true_positive_sites_rescued", res.joint_align.true_positive_sites_rescued}
        };

        return j.dump(2);
    }

    // Decision 7: Optional legacy JSON bundle emission (.FEL.json, .MEME.json, .BUSTED.json)
    static void write_legacy_json_bundle(const CharacterizeResult& res, const std::string& base_prefix) {
        using json = nlohmann::json;
        const size_t N = res.joint_align.map_alignment.num_taxa;
        const size_t C = res.joint_align.map_alignment.num_codons;
        const std::string nwk = res.joint_align.optimized_tree.to_newick();

        json common_input = {
            {"file name", base_prefix},
            {"number of sequences", N},
            {"number of sites", C},
            {"partition count", 1},
            {"trees", {{"0", nwk}}}
        };
        json tested_map = json::object();
        json branch_attrs = json::object();
        for (const auto& bp : res.branch_profiles) {
            tested_map[bp.branch_name] = "test";
            branch_attrs[bp.branch_name] = {
                {"Global MG94xREV", bp.branch_length_subs},
                {"Nucleotide GTR", bp.branch_length_subs}
            };
        }

        json common_fits = {
            {"Global MG94xREV", {
                {"Log Likelihood", res.gene_wide.null_log_l},
                {"estimated parameters", res.branch_profiles.size() + 8}
            }},
            {"Nucleotide GTR", {
                {"Log Likelihood", res.gene_wide.null_log_l},
                {"estimated parameters", res.branch_profiles.size() + 5}
            }}
        };

        // 1. Legacy FEL JSON
        {
            json jf;
            jf["analysis"] = {
                {"info", "HyPhy 3 Characterize Legacy-Compatible FEL Export"},
                {"version", "3.0.0"}
            };
            jf["input"] = common_input;
            jf["tested"] = {{"0", tested_map}};
            jf["branch attributes"] = {{"0", branch_attrs}};
            jf["fits"] = common_fits;
            jf["MLE"]["headers"] = json::array({
                json::array({"alpha", "Synonymous substitution rate at a site"}),
                json::array({"beta", "Non-synonymous substitution rate at a site"}),
                json::array({"alpha=beta", "Neutral rate estimate under H0"}),
                json::array({"LRT", "Likelihood ratio test statistic for beta = alpha"}),
                json::array({"p-value", "Asymptotic p-value for FEL/dFEL test"}),
                json::array({"Total branch length", "Column homology posterior confidence"})
            });
            json content = json::array();
            for (const auto& s : res.site_atlas) {
                Scalar p = res.directional_drift_active ? s.dfel_p_value : s.fel_p_value;
                Scalar lrt = res.directional_drift_active ? s.dfel_lrt : s.fel_lrt;
                content.push_back({s.alpha, s.beta_pervasive, s.alpha, lrt, p, s.column_homology_posterior});
            }
            jf["MLE"]["content"]["0"] = content;
            std::ofstream out(base_prefix + ".FEL.json");
            out << jf.dump(2);
        }
        // 2. Legacy MEME JSON
        {
            json jm;
            jm["analysis"] = {
                {"info", "HyPhy 3 Characterize Legacy-Compatible MEME Export"},
                {"version", "3.0.0"}
            };
            jm["input"] = common_input;
            jm["tested"] = {{"0", tested_map}};
            jm["branch attributes"] = {{"0", branch_attrs}};
            jm["fits"] = common_fits;
            jm["MLE"]["headers"] = json::array({
                json::array({"alpha", "Synonymous substitution rate at a site"}),
                json::array({"beta-", "Non-synonymous rate for the purifying/neutral branch class"}),
                json::array({"p-", "Mixture weight for beta-"}),
                json::array({"beta+", "Non-synonymous rate for the diversifying branch class"}),
                json::array({"p+", "Mixture weight for beta+"}),
                json::array({"LRT", "Likelihood ratio test statistic for episodic diversifying selection"}),
                json::array({"p-value", "Asymptotic mixture chi-squared p-value"}),
                json::array({"# branches under selection", "Branches with Empirical Bayes Factor >= 10"})
            });
            json content = json::array();
            for (const auto& s : res.site_atlas) {
                content.push_back({
                    s.alpha, s.beta_minus, 1.0 - s.weight_plus,
                    s.beta_plus, s.weight_plus, s.meme_lrt, s.meme_p_value,
                    s.num_episodic_branches
                });
            }
            jm["MLE"]["content"]["0"] = content;
            std::ofstream out(base_prefix + ".MEME.json");
            out << jm.dump(2);
        }
        // 3. Legacy BUSTED JSON
        {
            json jb;
            jb["analysis"] = {
                {"info", "HyPhy 3 Characterize Legacy-Compatible BUSTED-E Export"},
                {"version", "3.0.0"}
            };
            jb["input"] = common_input;
            jb["tested"] = {{"0", tested_map}};
            jb["branch attributes"] = {{"0", branch_attrs}};
            jb["test results"] = {
                {"LRT", res.gene_wide.busted_lrt},
                {"p-value", res.gene_wide.busted_p_value}
            };
            jb["fits"]["Unconstrained model"] = {
                {"Log Likelihood", res.gene_wide.unconstrained_log_l},
                {"Rate Distributions", {
                    {"Test", {
                        {"0", {{"omega", res.gene_wide.omega_classes[0]}, {"proportion", res.gene_wide.omega_weights[0]}}},
                        {"1", {{"omega", res.gene_wide.omega_classes[1]}, {"proportion", res.gene_wide.omega_weights[1]}}},
                        {"2", {{"omega", res.gene_wide.omega_classes[2]}, {"proportion", res.gene_wide.omega_weights[2]}}}
                    }}
                }}
            };
            jb["fits"]["Constrained model"] = {
                {"Log Likelihood", res.gene_wide.null_log_l}
            };
            std::ofstream out(base_prefix + ".BUSTED.json");
            out << jb.dump(2);
        }
    }

    // Decision 7: Unified, self-contained interactive 6-stage multi-tab SPA (.CHARACTERIZE.html)
    static std::string render_standalone_html_from_json_string(
        const std::string& raw_json_payload,
        const std::string& dataset_label
    ) {
        // Safely escape any closing </script sequence inside the JSON payload
        std::string safe_json;
        safe_json.reserve(raw_json_payload.size() + 32);
        for (size_t i = 0; i < raw_json_payload.size(); ++i) {
            if (raw_json_payload[i] == '<' && i + 1 < raw_json_payload.size() && raw_json_payload[i + 1] == '/') {
                safe_json += "<\\/";
                ++i;
            } else {
                safe_json += raw_json_payload[i];
            }
        }

        std::string safe_label;
        for (char c : dataset_label) {
            if (c == '<') safe_label += "&lt;";
            else if (c == '>') safe_label += "&gt;";
            else if (c == '&') safe_label += "&amp;";
            else if (c == '"') safe_label += "&quot;";
            else safe_label += c;
        }

        std::ostringstream html;
        html << R"WALRUS_HTML(<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>WALRUS • Evolutionary Characterization — )WALRUS_HTML" << safe_label << R"WALRUS_HTML(</title>
<style>
  :root {
    --bg: #f4f5f7;
    --card: #ffffff;
    --ink: #0f172a;
    --ink-2: #334155;
    --muted: #64748b;
    --border: #e2e8f0;
    --border-strong: #cbd5e1;
    --accent: #ea580c;
    --pos: #dc2626;
    --pos-bg: #fef2f2;
    --epi: #7c3aed;
    --epi-bg: #f5f3ff;
    --neg: #2563eb;
    --neg-bg: #eff6ff;
    --prop: #0d9488;
    --prop-bg: #f0fdfa;
    --warn: #d97706;
    --warn-bg: #fffbeb;
    --ok: #16a34a;
    --ok-bg: #f0fdf4;
  }
  * { box-sizing: border-box; margin: 0; padding: 0; }
  body {
    background: var(--bg);
    color: var(--ink);
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Inter, Roboto, Helvetica, Arial, sans-serif;
    line-height: 1.48;
    padding: 20px 24px 60px;
    font-size: 13px;
  }
  .container { max-width: 1520px; margin: 0 auto; }

  /* Header */
  header.top-bar {
    background: linear-gradient(135deg, #0f172a 0%, #1e293b 100%);
    color: #fff;
    padding: 18px 24px;
    border-radius: 10px;
    display: flex;
    justify-content: space-between;
    align-items: center;
    gap: 16px;
    margin-bottom: 14px;
    box-shadow: 0 4px 14px rgba(15, 23, 42, 0.10);
  }
  header.top-bar h1 {
    font-size: 18px;
    font-weight: 800;
    letter-spacing: 0.4px;
    display: flex;
    align-items: center;
    gap: 10px;
  }
  header.top-bar .sub {
    color: #cbd5e1;
    font-size: 12.5px;
    margin-top: 4px;
  }
  .header-actions {
    display: flex;
    gap: 8px;
    align-items: center;
    flex-wrap: wrap;
  }
  .btn {
    background: rgba(255,255,255,0.12);
    color: #fff;
    border: 1px solid rgba(255,255,255,0.22);
    padding: 7px 12px;
    border-radius: 6px;
    font-size: 12px;
    font-weight: 600;
    cursor: pointer;
    transition: all 0.15s;
  }
  .btn:hover { background: rgba(255,255,255,0.22); }
  .btn-light {
    background: #fff;
    color: var(--ink);
    border: 1px solid var(--border-strong);
  }
  .btn-light:hover { background: #f8fafc; border-color: var(--ink-2); }

  /* 6-Stage Stepper Ribbon */
  .stage-ribbon {
    display: grid;
    grid-template-columns: repeat(6, 1fr);
    gap: 10px;
    margin-bottom: 14px;
  }
  @media (max-width: 1200px) {
    .stage-ribbon { grid-template-columns: repeat(3, 1fr); }
  }
  .stage-card {
    background: var(--card);
    border: 1px solid var(--border);
    border-top: 3px solid var(--muted);
    border-radius: 8px;
    padding: 11px 13px;
    cursor: pointer;
    transition: transform 0.12s, box-shadow 0.12s, border-color 0.12s;
    position: relative;
  }
  .stage-card:hover {
    transform: translateY(-1px);
    box-shadow: 0 4px 12px rgba(15, 23, 42, 0.07);
    border-color: var(--ink-2);
  }
  .stage-card.st-ok { border-top-color: var(--ok); }
  .stage-card.st-active { border-top-color: var(--prop); }
  .stage-card.st-warn { border-top-color: var(--warn); }
  .stage-card.st-sig { border-top-color: var(--pos); }
  .stage-card .st-head {
    display: flex;
    justify-content: space-between;
    align-items: center;
    font-size: 10.5px;
    font-weight: 700;
    text-transform: uppercase;
    letter-spacing: 0.6px;
    color: var(--muted);
    margin-bottom: 4px;
  }
  .stage-card .st-title {
    font-size: 12.5px;
    font-weight: 700;
    color: var(--ink);
    margin-bottom: 3px;
  }
  .stage-card .st-val {
    font-size: 15px;
    font-weight: 800;
    color: var(--ink);
    font-variant-numeric: tabular-nums;
  }
  .stage-card .st-sub {
    font-size: 11px;
    color: var(--muted);
    margin-top: 3px;
    line-height: 1.35;
  }

  /* What WALRUS Changed Strip */
  .delta-banner {
    background: linear-gradient(90deg, #fff7ed 0%, #f0fdfa 100%);
    border: 1px solid #fed7aa;
    border-left: 4px solid var(--accent);
    border-radius: 8px;
    padding: 10px 16px;
    margin-bottom: 16px;
    display: flex;
    align-items: center;
    justify-content: space-between;
    flex-wrap: wrap;
    gap: 12px;
  }
  .delta-banner .delta-title {
    font-weight: 800;
    color: #9a3412;
    font-size: 12.5px;
    display: flex;
    align-items: center;
    gap: 8px;
  }
  .delta-chips {
    display: flex;
    flex-wrap: wrap;
    gap: 8px;
  }
  .delta-chip {
    background: #ffffff;
    border: 1px solid var(--border-strong);
    border-radius: 999px;
    padding: 3px 10px;
    font-size: 11.5px;
    font-weight: 600;
    color: var(--ink-2);
  }
  .delta-chip b { color: var(--ink); }

  /* Tabs */
  .tabs {
    display: flex;
    gap: 8px;
    margin-bottom: 16px;
    border-bottom: 2px solid var(--border);
    padding-bottom: 8px;
    flex-wrap: wrap;
  }
  .tab-btn {
    background: var(--card);
    border: 1px solid var(--border-strong);
    border-radius: 7px;
    padding: 9px 16px;
    font-size: 13px;
    font-weight: 700;
    color: var(--ink-2);
    cursor: pointer;
    transition: all 0.12s;
  }
  .tab-btn:hover { background: #f8fafc; color: var(--ink); }
  .tab-btn.active {
    background: var(--ink);
    color: #fff;
    border-color: var(--ink);
  }
  .tab-panel { display: none; }
  .tab-panel.active { display: block; }

  /* Cards & Grids */
  .card {
    background: var(--card);
    border: 1px solid var(--border);
    border-radius: 10px;
    padding: 16px 18px;
    margin-bottom: 16px;
    box-shadow: 0 1px 3px rgba(15, 23, 42, 0.03);
  }
  .card-header {
    display: flex;
    justify-content: space-between;
    align-items: baseline;
    flex-wrap: wrap;
    gap: 8px;
    margin-bottom: 10px;
  }
  .card-title {
    font-size: 14px;
    font-weight: 800;
    color: var(--ink);
  }
  .card-subtitle {
    font-size: 12px;
    color: var(--muted);
  }
  .grid-2 { display: grid; grid-template-columns: 1.25fr 1fr; gap: 16px; }
  .grid-2-eq { display: grid; grid-template-columns: 1fr 1fr; gap: 16px; }
  .grid-3 { display: grid; grid-template-columns: repeat(3, 1fr); gap: 14px; }
  @media (max-width: 1024px) {
    .grid-2, .grid-2-eq, .grid-3 { grid-template-columns: 1fr; }
  }

  /* Pills & Badges */
  .pill {
    display: inline-block;
    padding: 2px 8px;
    border-radius: 999px;
    font-size: 11px;
    font-weight: 700;
    white-space: nowrap;
  }
  .pill-pos { background: #fee2e2; color: #b91c1c; }
  .pill-epi { background: #ede9fe; color: #6d28d9; }
  .pill-diff { background: #ffedd5; color: #c2410c; }
  .pill-neg { background: #dbeafe; color: #1d4ed8; }
  .pill-prop { background: #ccfbf1; color: #0f766e; }
  .pill-warn { background: #fef3c7; color: #b45309; }
  .pill-ok { background: #dcfce7; color: #15803d; }
  .pill-neu { background: #f1f5f9; color: #475569; }
  .site-chip {
    display: inline-block;
    padding: 2px 7px;
    margin: 2px;
    border-radius: 5px;
    font-size: 11px;
    font-weight: 700;
    background: #f1f5f9;
    border: 1px solid var(--border-strong);
    color: var(--ink);
    cursor: pointer;
  }
  .site-chip:hover { background: var(--ink); color: #fff; }

  /* Filter Toolbar & Preset Pills */
  .toolbar {
    display: flex;
    gap: 10px;
    align-items: center;
    flex-wrap: wrap;
    margin-bottom: 12px;
  }
  .toolbar input[type="text"], .toolbar select {
    padding: 6px 11px;
    border: 1px solid var(--border-strong);
    border-radius: 6px;
    font-size: 12.5px;
    background: #fff;
    color: var(--ink);
  }
  .filter-pills {
    display: flex;
    gap: 6px;
    flex-wrap: wrap;
    margin-bottom: 12px;
  }
  .filter-pill {
    background: #f8fafc;
    border: 1px solid var(--border-strong);
    border-radius: 999px;
    padding: 4px 11px;
    font-size: 11.5px;
    font-weight: 600;
    color: var(--ink-2);
    cursor: pointer;
  }
  .filter-pill:hover { border-color: var(--ink); color: var(--ink); }
  .filter-pill.active {
    background: var(--ink);
    color: #fff;
    border-color: var(--ink);
  }

  /* Tables */
  .table-wrap {
    max-height: 580px;
    overflow: auto;
    border: 1px solid var(--border);
    border-radius: 8px;
  }
  table {
    width: 100%;
    border-collapse: collapse;
    font-size: 12px;
  }
  th, td {
    padding: 7px 10px;
    border-bottom: 1px solid var(--border);
    text-align: left;
  }
  th {
    background: #f8fafc;
    font-weight: 700;
    font-size: 11px;
    text-transform: uppercase;
    letter-spacing: 0.4px;
    color: var(--ink-2);
    position: sticky;
    top: 0;
    z-index: 2;
    cursor: pointer;
    user-select: none;
    white-space: nowrap;
  }
  th:hover { background: #e2e8f0; }
  td.num, th.num {
    text-align: right;
    font-variant-numeric: tabular-nums;
  }
  tbody tr { cursor: pointer; transition: background 0.08s; }
  tbody tr:hover { background: #f1f5f9; }
  tbody tr.row-selected { background: #fff7ed !important; outline: 1px solid var(--accent); }

  /* Floating Tooltip */
  #walrus-tooltip {
    position: fixed;
    pointer-events: none;
    background: rgba(15, 23, 42, 0.96);
    color: #fff;
    padding: 10px 13px;
    border-radius: 8px;
    font-size: 12px;
    line-height: 1.45;
    max-width: 340px;
    z-index: 9999;
    display: none;
    box-shadow: 0 8px 24px rgba(0, 0, 0, 0.25);
    border: 1px solid rgba(255,255,255,0.15);
  }

  /* Slide-Over Site Forensics Drawer */
  .drawer-backdrop {
    position: fixed;
    inset: 0;
    background: rgba(15, 23, 42, 0.35);
    z-index: 1000;
    display: none;
  }
  .drawer-backdrop.open { display: block; }
  .site-drawer {
    position: fixed;
    top: 0;
    right: -520px;
    width: 480px;
    max-width: 92vw;
    height: 100vh;
    background: #fff;
    z-index: 1001;
    box-shadow: -8px 0 30px rgba(15, 23, 42, 0.18);
    transition: right 0.22s ease;
    display: flex;
    flex-direction: column;
  }
  .site-drawer.open { right: 0; }
  .drawer-head {
    background: var(--ink);
    color: #fff;
    padding: 16px 20px;
    display: flex;
    justify-content: space-between;
    align-items: center;
  }
  .drawer-body {
    padding: 18px 20px;
    overflow-y: auto;
    flex: 1;
  }
  .kv-grid {
    display: grid;
    grid-template-columns: 1fr 1fr;
    gap: 8px 12px;
    margin-top: 8px;
  }
  .kv-box {
    background: #f8fafc;
    border: 1px solid var(--border);
    border-radius: 6px;
    padding: 8px 10px;
  }
  .kv-box .k { font-size: 10.5px; text-transform: uppercase; color: var(--muted); font-weight: 700; }
  .kv-box .v { font-size: 14px; font-weight: 800; color: var(--ink); margin-top: 2px; font-variant-numeric: tabular-nums; }
  .svg-container {
    width: 100%;
    overflow-x: auto;
  }
  svg text { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; }
</style>
</head>
<body>
<div id="walrus-tooltip"></div>

<!-- Slide-Over Site Forensics Drawer -->
<div id="drawerBackdrop" class="drawer-backdrop" onclick="closeSiteDrawer()"></div>
<aside id="siteDrawer" class="site-drawer">
  <div class="drawer-head">
    <div>
      <div style="font-size:11px;color:#94a3b8;text-transform:uppercase;font-weight:700;">Stage 5 • Site Forensics Inspector</div>
      <div id="drawerTitle" style="font-size:18px;font-weight:800;margin-top:2px;">Codon —</div>
    </div>
    <div style="display:flex;gap:6px;align-items:center;">
      <button class="btn" onclick="stepSiteDrawer(-1)" title="Previous Codon">&larr;</button>
      <button class="btn" onclick="stepSiteDrawer(1)" title="Next Codon">&rarr;</button>
      <button class="btn" onclick="closeSiteDrawer()">&#10005;</button>
    </div>
  </div>
  <div id="drawerContent" class="drawer-body"></div>
</aside>

<div class="container">
  <header class="top-bar">
    <div>
      <h1>
        <span>WALRUS &bull; UNIFIED EVOLUTIONARY CHARACTERIZATION</span>
        <span class="pill" style="background:#ea580c;color:#fff;font-size:11px;">HyPhy 3</span>
      </h1>
      <div class="sub" id="headerSubtitle">Loading dataset...</div>
    </div>
    <div class="header-actions">
      <button class="btn" onclick="exportSiteCSV()">Export Site Atlas CSV</button>
      <button class="btn" onclick="downloadEmbeddedJSON()">Download Full JSON</button>
    </div>
  </header>

  <!-- 6-Stage Pipeline Ribbon -->
  <div id="stageRibbon" class="stage-ribbon"></div>

  <!-- What WALRUS Changed vs Standard Pipeline Strip -->
  <div id="deltaBanner" class="delta-banner"></div>

  <!-- Top-Level Navigation Tabs -->
  <div class="tabs">
    <button class="tab-btn active" id="btn-tab-overview" onclick="switchTab('tab-overview')">1. Overview &amp; 5-Second CDS Map</button>
    <button class="tab-btn" id="btn-tab-integrity" onclick="switchTab('tab-integrity')">2. Stage 0–2 Integrity, Homology &amp; Reticulation</button>
    <button class="tab-btn" id="btn-tab-tree" onclick="switchTab('tab-tree')">3. Stage 3–4 Interactive Phylogram &amp; Lineages</button>
    <button class="tab-btn" id="btn-tab-sites" onclick="switchTab('tab-sites')">4. Stage 5 Site, dFEL &amp; Convergence Explorer</button>
  </div>

  <!-- ===================================================================== -->
  <!-- TAB 1: OVERVIEW & 5-SECOND INTERACTIVE CDS ARCHITECTURE MAP           -->
  <!-- ===================================================================== -->
  <div id="tab-overview" class="tab-panel active">
    <div class="card">
      <div class="card-header">
        <div>
          <div class="card-title">Interactive 5-Second CDS Architecture, Reticulation &amp; Bidirectional Selection Profile</div>
          <div class="card-subtitle">Upward bars: Diversifying (Red/Purple) &amp; Convergent/Toggle (Teal) selection (&minus;log<sub>10</sub> q). Downward bars: Purifying constraint (Blue). Hover or click any codon to inspect.</div>
        </div>
        <div style="display:flex;gap:10px;align-items:center;font-size:11.5px;font-weight:600;">
          <span style="display:inline-flex;align-items:center;gap:4px;"><span style="width:10px;height:10px;background:#dc2626;display:inline-block;border-radius:2px;"></span>Episodic/Pervasive Div.</span>
          <span style="display:inline-flex;align-items:center;gap:4px;"><span style="width:10px;height:10px;background:#7c3aed;display:inline-block;border-radius:2px;"></span>Convergent/Toggle</span>
          <span style="display:inline-flex;align-items:center;gap:4px;"><span style="width:10px;height:10px;background:#2563eb;display:inline-block;border-radius:2px;"></span>Purifying</span>
          <span style="display:inline-flex;align-items:center;gap:4px;"><span style="width:12px;height:10px;background:#fed7aa;border:1px solid #f97316;display:inline-block;border-radius:2px;"></span>Stage 2 Recombination Bridge</span>
        </div>
      </div>
      <div id="cdsPlotContainer" class="svg-container"></div>
    </div>

    <div class="grid-2">
      <div class="card">
        <div class="card-header">
          <div class="card-title">Automated 6-Stage Evolutionary Synthesis</div>
          <span class="pill pill-neu">Plain-English Forensic Summary</span>
        </div>
        <div id="narrativeContainer" style="font-size:13px;line-height:1.6;color:var(--ink-2);"></div>
      </div>
      <div class="card">
        <div class="card-header">
          <div class="card-title">Stage 4 Gene-Wide Selection (BUSTED-E) &amp; Error-Sink Decomposition</div>
          <span id="bustedVerdictPill" class="pill pill-neu"></span>
        </div>
        <div id="bustedBreakdownContainer"></div>
      </div>
    </div>
  </div>

  <!-- ===================================================================== -->
  <!-- TAB 2: STAGE 0-2 INTEGRITY, HOMOLOGY & RETICULATION FORENSICS         -->
  <!-- ===================================================================== -->
  <div id="tab-integrity" class="tab-panel">
    <div class="grid-2-eq">
      <div class="card">
        <div class="card-header">
          <div>
            <div class="card-title">Stage 2: Recombination &amp; Paralog Conversion Bridges</div>
            <div class="card-subtitle">Localized Ohmic Wheatstone bridges and alternate SPR topologies routed during cavity evaluation</div>
          </div>
          <span id="bridgeCountPill" class="pill pill-neu"></span>
        </div>
        <div id="bridgesContainer"></div>
      </div>

      <div class="card">
        <div class="card-header">
          <div>
            <div class="card-title">Stage 1: Alignment Objective &amp; Indel Error-Sink Suppressed Sites</div>
            <div class="card-subtitle">Sites where POA homology uncertainty or flanking indels triggered false-positive calls under rigid alignments</div>
          </div>
          <span id="artifactCountPill" class="pill pill-neu"></span>
        </div>
        <div id="artifactContainer"></div>
      </div>
    </div>

    <div class="grid-2-eq">
      <div class="card">
        <div class="card-header">
          <div>
            <div class="card-title">Stage 4C: 1D Spatial Selective Domains &amp; Synonymous Constraint Islands</div>
            <div class="card-subtitle">Penalized BIC segmentation of log(&omega;) and binomial scan for dual-coding / RNA structure islands</div>
          </div>
        </div>
        <div id="domainsAndIslandsContainer"></div>
      </div>

      <div class="card">
        <div class="card-header">
          <div>
            <div class="card-title">Stage 4D: Mutational Context (CpG / APOBEC) &amp; Coevolution</div>
            <div class="card-subtitle">Dinucleotide hypermutability enrichment among selected sites and epistatic co-substitution pairs</div>
          </div>
        </div>
        <div id="contextAndCoevContainer"></div>
      </div>
    </div>

    <div class="card">
      <div class="card-header">
        <div>
          <div class="card-title">Stage 0: Input Sequence Admission, Reading-Frame Triage &amp; Quarantine Audit</div>
          <div class="card-subtitle">Pre-alignment screening for frameshifts, premature stop codons, fragment coverage, and non-homologous outliers</div>
        </div>
        <div class="toolbar" style="margin-bottom:0;">
          <input type="text" id="admSearch" placeholder="Search sequence name or reason..." oninput="renderAdmissionTable()">
          <select id="admStatusFilter" onchange="renderAdmissionTable()">
            <option value="ALL">All Input Sequences</option>
            <option value="QUARANTINED">Quarantined Only</option>
            <option value="RESCUED">Frame/Strand Rescued</option>
            <option value="ADMITTED">Admitted</option>
          </select>
        </div>
      </div>
      <div id="admissionTableContainer" class="table-wrap"></div>
    </div>
  </div>

  <!-- ===================================================================== -->
  <!-- TAB 3: STAGE 3-4 INTERACTIVE PHYLOGRAM & LINEAGE REGIMES              -->
  <!-- ===================================================================== -->
  <div id="tab-tree" class="tab-panel">
    <div class="card">
      <div class="card-header">
        <div>
          <div class="card-title">Interactive Phylogenetic Tree &amp; Lineage Regime Explorer</div>
          <div class="card-subtitle">Hover or click any branch to inspect its GC3 drift regime, episodic diversifying selection (aBSREL), RELAX intensity parameter K<sub>b</sub>, and core-invasion LoF audit.</div>
        </div>
        <div class="toolbar" style="margin-bottom:0;">
          <label style="font-weight:700;font-size:12px;">Color Branches By:</label>
          <select id="treeColorMode" onchange="renderPhylogram()">
            <option value="gc3">Stage 3: 3rd-Position GC3% &amp; Drift Regime</option>
            <option value="absrel">Stage 4: Episodic Diversifying (aBSREL p &le; 0.05)</option>
            <option value="relax">Stage 4: Selection Intensity Shift (RELAX K_b)</option>
            <option value="lof">Stage 4B: Loss-of-Constraint / LoF Alert</option>
            <option value="macro">Stage 4: Macro-Regime Partition (Test vs Ref)</option>
          </select>
          <label style="display:inline-flex;align-items:center;gap:4px;font-weight:600;font-size:12px;cursor:pointer;">
            <input type="checkbox" id="treeShowBridges" checked onchange="renderPhylogram()"> Show Recombination Bridges
          </label>
          <label style="display:inline-flex;align-items:center;gap:4px;font-weight:600;font-size:12px;cursor:pointer;">
            <input type="checkbox" id="treeCladogram" onchange="renderPhylogram()"> Cladogram (Equal Branch Lengths)
          </label>
          <label style="font-weight:600;font-size:12px;">Row Height:</label>
          <input type="range" id="treeRowHeight" min="11" max="26" value="15" oninput="renderPhylogram()" style="width:90px;">
        </div>
      </div>
      <div id="treeLegendBar" style="margin-bottom:8px;font-size:12px;color:var(--ink-2);display:flex;gap:14px;flex-wrap:wrap;align-items:center;"></div>
      <div id="treeSvgContainer" class="svg-container" style="max-height:640px;overflow:auto;border:1px solid var(--border);border-radius:8px;background:#fafaf9;"></div>
    </div>

    <div class="card">
      <div class="card-header">
        <div>
          <div class="card-title">Lineage Evolutionary Profiles (aBSREL / RELAX / GC3 / Loss-of-Constraint)</div>
          <div class="card-subtitle">Click any column header to sort; click any invaded core site pill to open the Site Forensics Drawer.</div>
        </div>
        <div class="toolbar" style="margin-bottom:0;">
          <input type="text" id="branchSearch" placeholder="Search branch or taxon name..." oninput="renderBranchTable(); renderPhylogram();" style="width:240px;">
          <select id="branchFilter" onchange="renderBranchTable()">
            <option value="ALL">All Branches</option>
            <option value="ABSREL_SIG">aBSREL Significant (p &le; 0.05)</option>
            <option value="LOF_ALERT">LoF / Core-Invasion Flagged</option>
            <option value="RELAX_SHIFT">RELAX Intensified or Relaxed</option>
            <option value="LEAVES">Terminal Taxa (Leaves) Only</option>
            <option value="INTERNAL">Internal Branches Only</option>
          </select>
        </div>
      </div>
      <div class="table-wrap">
        <table id="branchTable">
          <thead>
            <tr>
              <th onclick="sortBranches('branch_name')">Branch</th>
              <th onclick="sortBranches('macro_regime')">Macro-Regime</th>
              <th class="num" onclick="sortBranches('branch_length_subs')">Length (subs/codon)</th>
              <th class="num" onclick="sortBranches('gc_cp3')">GC3% (Regime)</th>
              <th class="num" onclick="sortBranches('omega_plus')">&omega;<sup>+</sup> (w<sup>+</sup>)</th>
              <th class="num" onclick="sortBranches('absrel_p_value')">aBSREL p (Holm / q)</th>
              <th class="num" onclick="sortBranches('relax_k')">RELAX K<sub>b</sub></th>
              <th class="num" onclick="sortBranches('relax_p_value')">RELAX p (Call)</th>
              <th class="num" onclick="sortBranches('core_nonsyn_observed')">Core Subs Obs (Exp)</th>
              <th class="num" onclick="sortBranches('lof_bayes_factor')">LoF BF (p)</th>
              <th>Invaded Core Sites</th>
            </tr>
          </thead>
          <tbody id="branchTableBody"></tbody>
        </table>
      </div>
    </div>
  </div>

  <!-- ===================================================================== -->
  <!-- TAB 4: STAGE 5 SITE, dFEL & CONVERGENCE EXPLORER                      -->
  <!-- ===================================================================== -->
  <div id="tab-sites" class="tab-panel">
    <div class="grid-2-eq">
      <div class="card">
        <div class="card-header">
          <div>
            <div class="card-title">Stationary FEL vs. Non-Reversible dFEL (Threshold-Crossing Sites)</div>
            <div class="card-subtitle">Sites where directional GC3 drift correction unmasks hidden selection or suppresses stationary artifacts (dashed line: p = 0.05, LRT = 3.84)</div>
          </div>
        </div>
        <div id="dfelPlotContainer" class="svg-container"></div>
      </div>

      <div class="card">
        <div class="card-header">
          <div>
            <div class="card-title">Convergent Evolution &amp; Amino-Acid Toggle Lollipop Atlas</div>
            <div class="card-subtitle">Independent phylogenetic lineages converging on the same target amino acid (colored by &minus;log<sub>10</sub> p; click any lollipop to inspect)</div>
          </div>
        </div>
        <div id="convPlotContainer" class="svg-container"></div>
      </div>
    </div>

    <div class="card">
      <div class="card-header">
        <div>
          <div class="card-title">Unified Site Atlas (FEL / dFEL / MEME / Exchangeability / Homology Forensics)</div>
          <div class="card-subtitle">Click any row to open the Site Forensics Drawer; click any column header to sort.</div>
        </div>
        <div class="toolbar" style="margin-bottom:0;">
          <input type="text" id="siteSearchInput" placeholder="Search site #, classification, branch, AA, driver..." oninput="renderSiteTable()" style="width:280px;">
          <button class="btn btn-light" onclick="exportSiteCSV()">Export Filtered CSV</button>
        </div>
      </div>

      <div id="siteFilterPills" class="filter-pills"></div>

      <div class="table-wrap">
        <table id="siteAtlasTable">
          <thead>
            <tr>
              <th class="num" onclick="sortSites('site')">Codon</th>
              <th onclick="sortSites('classification')">Unified Classification</th>
              <th onclick="sortSites('spatial_domain')">1D Domain / Topo</th>
              <th class="num" onclick="sortSites('column_homology_posterior')">POA PP (Gap%)</th>
              <th class="num" onclick="sortSites('simes_omnibus_p')">Simes p (BH q)</th>
              <th class="num" onclick="sortSites('alpha')">&alpha;</th>
              <th class="num" onclick="sortSites('beta_pervasive')">&beta;<sub>perv</sub></th>
              <th class="num" onclick="sortSites('fel_p_value')">Stationary FEL p</th>
              <th class="num" onclick="sortSites('dfel_p_value')">Directional dFEL p</th>
              <th class="num" onclick="sortSites('meme_p_value')">MEME p (&beta;<sup>+</sup>)</th>
              <th class="num" onclick="sortSites('effective_aa_alphabet')">N<sub>eff</sub> / E[N<sub>eff</sub>]</th>
              <th class="num" onclick="sortSites('max_convergent_transitions')">Conv. Lineages (AA)</th>
              <th>Driving Branches / Context</th>
            </tr>
          </thead>
          <tbody id="siteAtlasTableBody"></tbody>
        </table>
      </div>
    </div>
  </div>
</div>

<script id="walrus-data" type="application/json">)WALRUS_HTML"
             << safe_json
             << R"WALRUS_HTML(</script>

<script>
const DATA = JSON.parse(document.getElementById('walrus-data').textContent);
const DATASET_LABEL = )WALRUS_HTML" << nlohmann::json(dataset_label).dump() << R"WALRUS_HTML(;

// Global UI state
let currentSiteFilter = 'ALL';
let siteSortKey = 'site';
let siteSortAsc = true;
let branchSortKey = 'absrel_p_value';
let branchSortAsc = true;
let selectedCodon = null;
let selectedBranchName = null;

const siteByNumber = {};
(DATA.site_atlas || []).forEach(s => { siteByNumber[s.site] = s; });
const branchByName = {};
(DATA.branch_profiles || []).forEach(b => { branchByName[b.branch_name] = b; });

function fmtP(p) {
  if (p === null || p === undefined || isNaN(p)) return '—';
  if (p < 1e-4) return p.toExponential(1);
  return p.toFixed(4);
}
function fmtNum(x, d = 2) {
  if (x === null || x === undefined || isNaN(x)) return '—';
  return Number(x).toFixed(d);
}
function escHtml(s) {
  return String(s ?? '').replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;').replace(/"/g,'&quot;');
}

function showTooltip(evt, htmlContent) {
  const tt = document.getElementById('walrus-tooltip');
  tt.innerHTML = htmlContent;
  tt.style.display = 'block';
  const pad = 14;
  let x = evt.clientX + pad;
  let y = evt.clientY + pad;
  if (x + 350 > window.innerWidth) x = evt.clientX - 350;
  if (y + 220 > window.innerHeight) y = Math.max(10, evt.clientY - 200);
  tt.style.left = x + 'px';
  tt.style.top = y + 'px';
}
function hideTooltip() {
  document.getElementById('walrus-tooltip').style.display = 'none';
}

function switchTab(tabId) {
  document.querySelectorAll('.tab-panel').forEach(el => el.classList.remove('active'));
  document.querySelectorAll('.tab-btn').forEach(el => el.classList.remove('active'));
  const panel = document.getElementById(tabId);
  if (panel) panel.classList.add('active');
  const btn = document.getElementById('btn-' + tabId);
  if (btn) btn.classList.add('active');
  if (tabId === 'tab-tree') renderPhylogram();
}

function getClassificationPillClass(cls) {
  if (!cls) return 'pill-neu';
  if (cls.includes('Pervasive') && cls.includes('Diversifying')) return 'pill-pos';
  if (cls.includes('Episodic')) return 'pill-epi';
  if (cls.includes('Regime')) return 'pill-diff';
  if (cls.includes('Purifying')) return 'pill-neg';
  if (cls.includes('Toggle') || cls.includes('Convergent')) return 'pill-prop';
  return 'pill-neu';
}

function isDfelDiscordant(s) {
  const fSig = s.fel_p_value <= 0.05;
  const dSig = s.dfel_p_value <= 0.05;
  return fSig !== dSig;
}

// ============================================================================
// 1. HEADER, 6-STAGE RIBBON & WHAT WALRUS CHANGED BANNER
// ============================================================================
function renderHeaderAndRibbon() {
  const aln = DATA.alignment_and_network || {};
  const adm = DATA.input_admission_audit || {};
  const comp = DATA.compositional_audit || {};
  const gw = DATA.gene_wide_selection || {};
  const counts = DATA.summary_counts || {};
  const bridges = DATA.recombination_bridges || [];
  const rt = (DATA.analysis && DATA.analysis.runtime_seconds) ? DATA.analysis.runtime_seconds : 0;

  document.getElementById('headerSubtitle').innerHTML =
    `Dataset: <b>${escHtml(DATASET_LABEL)}</b> &bull; <b>${aln.num_taxa || 0}</b> admitted taxa &times; <b>${aln.num_codons || 0}</b> codons &bull; Mean Homology PP = <b>${fmtNum(aln.option5_mean_column_posterior, 4)}</b> &bull; Runtime: <b>${fmtNum(rt, 2)} s</b>`;

  const nQuar = adm.quarantined_sequences || 0;
  const nResc = adm.rescued_frame_or_strand || 0;
  const nArt = (DATA.site_atlas || []).filter(s => s.artifact_eliminated).length;
  const nBridges = bridges.length;
  const nDisc = (DATA.site_atlas || []).filter(isDfelDiscordant).length;
  const nUnmasked = (DATA.site_atlas || []).filter(s => s.dfel_p_value <= 0.05 && s.fel_p_value > 0.05).length;

  const stages = [
    {
      num: 'Stage 0',
      badge: nQuar > 0 ? 'Quarantined' : 'Clean',
      cls: nQuar > 0 ? 'st-warn' : 'st-ok',
      title: 'Input Admission & Triage',
      val: `${adm.admitted_sequences || aln.num_taxa}/${adm.total_input_sequences || aln.num_taxa} Taxa`,
      sub: `${nQuar} quarantined &bull; ${nResc} frame-rescued`,
      tab: 'tab-integrity'
    },
    {
      num: 'Stage 1',
      badge: nArt > 0 ? `${nArt} Filtered` : 'High PP',
      cls: nArt > 0 ? 'st-active' : 'st-ok',
      title: 'POA Homology & Indels',
      val: `PP = ${fmtNum(aln.option5_mean_column_posterior, 3)}`,
      sub: `${nArt} indel-artifact site(s) suppressed`,
      tab: 'tab-integrity'
    },
    {
      num: 'Stage 2',
      badge: nBridges > 0 ? `${nBridges} Bridge(s)` : 'Tree-Like',
      cls: nBridges > 0 ? 'st-warn' : 'st-ok',
      title: 'Reticulation & SPR Routing',
      val: nBridges > 0 ? `${nBridges} Active Tract(s)` : 'Single Tree',
      sub: nBridges > 0 ? `Peak ΔlnL = +${fmtNum(bridges[0].segment_delta_log_l, 1)} (p=${fmtP(bridges[0].p_value)})` : 'No recombination bridges (p > 0.01)',
      tab: 'tab-integrity'
    },
    {
      num: 'Stage 3',
      badge: comp.directional_drift_active ? '3-Regime dFEL' : 'Stationary',
      cls: comp.directional_drift_active ? 'st-active' : 'st-ok',
      title: 'GC3 Compositional Drift',
      val: `CP3 p = ${fmtP(comp.p_value_cp3)}`,
      sub: comp.directional_drift_active ? `χ² = ${fmtNum(comp.chi2_cp3, 1)} &bull; ${nUnmasked} sites unmasked` : `χ² = ${fmtNum(comp.chi2_cp3, 1)} (Stationary FEL)`,
      tab: 'tab-tree'
    },
    {
      num: 'Stage 4',
      badge: gw.has_positive_selection ? 'POSITIVE' : 'Purifying/Neutral',
      cls: gw.has_positive_selection ? 'st-sig' : 'st-ok',
      title: 'Gene-Wide & 1D Domains',
      val: `p = ${fmtP(gw.p_value)}`,
      sub: `BUSTED-E LRT = ${fmtNum(gw.lrt, 2)} &bull; ${counts.spatial_domains_detected || 1} domain(s)`,
      tab: 'tab-overview'
    },
    {
      num: 'Stage 5',
      badge: `${counts.sites_simes_fdr_significant || 0} FDR Sites`,
      cls: (counts.sites_episodic_diversifying || 0) > 0 ? 'st-sig' : 'st-active',
      title: 'Site & Convergence Atlas',
      val: `${counts.sites_episodic_diversifying || 0} Epi / ${counts.sites_pervasive_diversifying || 0} Perv`,
      sub: `${counts.sites_pervasive_purifying || 0} purifying &bull; ${counts.sites_exchangeability_or_convergent || 0} toggle/conv`,
      tab: 'tab-sites'
    }
  ];

  document.getElementById('stageRibbon').innerHTML = stages.map(s => `
    <div class="stage-card ${s.cls}" onclick="switchTab('${s.tab}')">
      <div class="st-head"><span>${s.num}</span><span class="pill pill-neu" style="font-size:10px;padding:1px 6px;">${s.badge}</span></div>
      <div class="st-title">${s.title}</div>
      <div class="st-val">${s.val}</div>
      <div class="st-sub">${s.sub}</div>
    </div>
  `).join('');

  document.getElementById('deltaBanner').innerHTML = `
    <div class="delta-title">
      <span>&#9889; What WALRUS Resolved on ${escHtml(DATASET_LABEL)}:</span>
    </div>
    <div class="delta-chips">
      <span class="delta-chip">Quarantined Sequences: <b>${nQuar}</b></span>
      <span class="delta-chip">Indel-Artifact Sites Suppressed: <b>${nArt}</b></span>
      <span class="delta-chip">Recombination / Conversion Bridges: <b>${nBridges}</b></span>
      <span class="delta-chip">dFEL vs Stationary FEL Discordant Sites: <b>${nDisc}</b> (<b>${nUnmasked}</b> unmasked)</span>
      <span class="delta-chip">Episodic Selected Branches (aBSREL): <b>${counts.branches_episodic_diversifying || 0}</b></span>
      <span class="delta-chip">Loss-of-Constraint (LoF) Lineages: <b>${counts.branches_loss_of_constraint || 0}</b></span>
    </div>
  `;
}

// ============================================================================
// 2. TAB 1: INTERACTIVE 5-SECOND CDS ARCHITECTURE & SELECTION PLOT
// ============================================================================
function renderOverviewTab() {
  const sites = DATA.site_atlas || [];
  const domains = DATA.spatial_domains || [];
  const islands = DATA.synonymous_constraint_islands || [];
  const bridges = DATA.recombination_bridges || [];
  const C = sites.length || 1;

  const W = 1380, H = 370;
  const padL = 62, padR = 26, padT = 56, padB = 64;
  const plotW = W - padL - padR;
  const plotH = H - padT - padB;
  const zeroY = padT + plotH * 0.48;
  const maxUp = 4.2, maxDown = 4.2;

  const xOf = c => padL + ((c - 0.5) / C) * plotW;
  const wCod = Math.max(1.5, plotW / C);
  const yUp = val => zeroY - (Math.min(maxUp, Math.max(0, val)) / maxUp) * (plotH * 0.46);
  const yDown = val => zeroY + (Math.min(maxDown, Math.max(0, val)) / maxDown) * (plotH * 0.46);

  let svg = `<svg viewBox="0 0 ${W} ${H}" width="100%" height="auto" style="display:block;user-select:none;">`;

  // Horizontal grid lines
  [-4, -2, 0, 2, 4].forEach(gv => {
    const yy = gv >= 0 ? yUp(gv) : yDown(-gv);
    const lbl = gv === 4 ? '4 (Div.)' : (gv === -4 ? '4 (Pur.)' : String(Math.abs(gv)));
    svg += `<line x1="${padL}" y1="${yy}" x2="${W - padR}" y2="${yy}" stroke="${gv === 0 ? '#334155' : '#e2e8f0'}" stroke-width="${gv === 0 ? 1.3 : 1}" stroke-dasharray="${gv === 0 ? '' : '3,3'}"/>`;
    svg += `<text x="${padL - 8}" y="${yy + 4}" text-anchor="end" font-size="11" fill="#475569" font-weight="${gv === 0 ? '700' : '500'}">${lbl}</text>`;
  });
  svg += `<text x="16" y="${zeroY}" transform="rotate(-90 16 ${zeroY})" text-anchor="middle" font-size="11.5" font-weight="700" fill="#1e293b">Selection Signal (−log10 q)</text>`;

  // Stage 2 Recombination Bridge shaded regions
  bridges.forEach(b => {
    const x1 = xOf(b.peak_start_codon) - wCod * 0.5;
    const x2 = xOf(b.peak_end_codon) + wCod * 0.5;
    svg += `<rect x="${x1}" y="${padT}" width="${Math.max(4, x2 - x1)}" height="${plotH * 0.94}" fill="#fed7aa" fill-opacity="0.42" stroke="#f97316" stroke-width="1" stroke-dasharray="4,3" rx="4"/>`;
    svg += `<text x="${0.5 * (x1 + x2)}" y="${padT + 14}" text-anchor="middle" font-size="10.5" font-weight="700" fill="#9a3412">Bridge: ${escHtml(b.clade_u)} ↔ ${escHtml(b.target_v)} (Codons ${b.peak_start_codon}–${b.peak_end_codon}, ΔlnL=+${fmtNum(b.segment_delta_log_l, 1)})</text>`;
  });

  // Top 1D Spatial Domain Architecture Pills
  domains.forEach(d => {
    const x1 = xOf(d.start_codon) - wCod * 0.4;
    const x2 = xOf(d.end_codon) + wCod * 0.4;
    let col = '#0284c7';
    if ((d.domain_type || '').includes('Hypervariable')) col = '#ea580c';
    else if ((d.domain_type || '').includes('Permissive')) col = '#0d9488';
    const dw = Math.max(6, x2 - x1 - 2);
    svg += `<rect x="${x1 + 1}" y="10" width="${dw}" height="22" rx="6" fill="${col}"/>`;
    if (dw > 75) {
      svg += `<text x="${x1 + 1 + dw / 2}" y="25" text-anchor="middle" font-size="10.5" font-weight="700" fill="#fff">${escHtml(d.domain_type)} (${d.start_codon}–${d.end_codon}, ω=${fmtNum(d.mean_omega, 2)})</text>`;
    }
  });

  // Synonymous Constraint Island Pills (between domains and plot)
  islands.forEach(isl => {
    const x1 = xOf(isl.start_codon) - wCod * 0.4;
    const x2 = xOf(isl.end_codon) + wCod * 0.4;
    svg += `<rect x="${x1}" y="35" width="${Math.max(8, x2 - x1)}" height="14" rx="4" fill="#16a34a"/>`;
    svg += `<text x="${0.5 * (x1 + x2)}" y="45" text-anchor="middle" font-size="9.5" font-weight="700" fill="#fff">SynIsl ${isl.start_codon}–${isl.end_codon}</text>`;
  });

  // Site bars + bottom POA Homology Posterior strip
  const ppStripY = H - 36;
  svg += `<text x="${padL - 8}" y="${ppStripY + 8}" text-anchor="end" font-size="10" font-weight="700" fill="#475569">POA PP</text>`;

  sites.forEach(s => {
    const cx = xOf(s.site);
    const qLog = Math.min(maxUp, -Math.log10(Math.max(s.simes_omnibus_q || 1.0, 1e-5)));
    const cls = s.classification || 'Neutral';

    let upVal = 0, downVal = 0, barCol = '#94a3b8';
    if (cls.includes('Diversifying')) {
      upVal = Math.max(0.9, qLog);
      barCol = cls.includes('Pervasive') ? '#dc2626' : '#7c3aed';
    } else if (cls.includes('Toggle') || cls.includes('Convergent') || (s.max_convergent_transitions >= 3 && s.exchangeability_p_value <= 0.05)) {
      upVal = Math.max(0.7, Math.min(3.4, -Math.log10(Math.max(s.exchangeability_p_value || 1.0, 1e-4))));
      barCol = '#0d9488';
    } else if (cls.includes('Purifying')) {
      downVal = Math.max(0.5, qLog);
      barCol = '#2563eb';
    }

    if (upVal > 0) {
      const yy = yUp(upVal);
      svg += `<rect x="${cx - wCod * 0.42}" y="${yy}" width="${Math.max(1.6, wCod * 0.84)}" height="${zeroY - yy}" fill="${barCol}" opacity="0.92"/>`;
    }
    if (downVal > 0) {
      const yy = yDown(downVal);
      svg += `<rect x="${cx - wCod * 0.38}" y="${zeroY}" width="${Math.max(1.4, wCod * 0.76)}" height="${yy - zeroY}" fill="${barCol}" opacity="0.80"/>`;
    }

    // Bottom POA Homology Posterior cell
    const pp = s.column_homology_posterior ?? 1.0;
    const ppCol = s.artifact_eliminated ? '#dc2626' : (pp >= 0.99 ? '#16a34a' : (pp >= 0.90 ? '#f59e0b' : '#ef4444'));
    svg += `<rect x="${cx - wCod * 0.5}" y="${ppStripY}" width="${Math.max(1.2, wCod)}" height="10" fill="${ppCol}"/>`;

    // Invisible full-height interactive hitbox for smooth hovering & clicking
    svg += `<rect x="${cx - wCod * 0.5}" y="${padT}" width="${Math.max(2.0, wCod)}" height="${H - padT - 20}" fill="transparent" style="cursor:pointer;" onmousemove="hoverCodon(event, ${s.site})" onmouseout="hideTooltip()" onclick="openSiteDrawer(${s.site})"/>`;
  });

  // X-axis ticks
  const step = C <= 120 ? 10 : (C <= 300 ? 25 : (C <= 700 ? 50 : 100));
  for (let c = step; c <= C; c += step) {
    const xx = xOf(c);
    svg += `<line x1="${xx}" y1="${ppStripY + 10}" x2="${xx}" y2="${ppStripY + 15}" stroke="#64748b"/>`;
    svg += `<text x="${xx}" y="${H - 8}" text-anchor="middle" font-size="10.5" fill="#475569">${c}</text>`;
  }
  svg += `<text x="${padL + plotW / 2}" y="${H - 2}" text-anchor="middle" font-size="11.5" font-weight="700" fill="#1e293b">Codon Position (1 to ${C})</text>`;
  svg += `</svg>`;
  document.getElementById('cdsPlotContainer').innerHTML = svg;

  // Render Narrative & BUSTED-E Breakdown
  const aln = DATA.alignment_and_network || {};
  const comp = DATA.compositional_audit || {};
  const gw = DATA.gene_wide_selection || {};
  const counts = DATA.summary_counts || {};
  const mca = DATA.mutational_context_alert || {};

  document.getElementById('narrativeContainer').innerHTML = `
    <p style="margin-bottom:8px;"><b>1. Homology, Indel &amp; Reticulation Integrity (Stages 0–2):</b> Across <b>${aln.num_taxa}</b> admitted sequences and <b>${aln.num_codons}</b> codons, POA column homology achieved mean posterior <b>PP = ${fmtNum(aln.option5_mean_column_posterior, 4)}</b> (min <b>${fmtNum(aln.option5_min_column_posterior, 3)}</b>), with <b>${aln.branch_indel_events || 0}</b> branch indel events. ${bridges.length > 0 ? `Ohmic Wheatstone screening localized <b>${bridges.length} recombination / gene-conversion bridge(s)</b> and routed their site likelihood cavities through local SPR topologies.` : `No non-tree-like recombination bridges exceeded the routing threshold.`}</p>
    <p style="margin-bottom:8px;"><b>2. Compositional Drift &amp; Lineage Regimes (Stage 3):</b> 3rd-codon-position nucleotide heterogeneity reported <b>CP3 &chi;&sup2; = ${fmtNum(comp.chi2_cp3, 1)}</b> (p = <b>${fmtP(comp.p_value_cp3)}</b>), activating <b>${comp.directional_drift_active ? 'Non-Reversible 3-Regime Directional Softmax Kernels (dFEL)' : 'Stationary Reversible Codon Kernels (FEL)'}</b>. Lineage-specific core-invasion screening flagged <b>${counts.branches_loss_of_constraint || 0} branch(es)</b> for potential relaxation of purifying core constraint.</p>
    <p style="margin-bottom:8px;"><b>3. Gene-Wide Selection &amp; 1D Architecture (Stage 4):</b> Error-compensated <b>BUSTED-E</b> ${gw.has_positive_selection ? '<b style="color:#dc2626;">detected statistically significant gene-wide episodic diversifying selection</b>' : '<b>did not detect gene-wide diversifying selection</b>'} (LRT = <b>${fmtNum(gw.lrt, 2)}</b>, p = <b>${fmtP(gw.p_value)}</b>). 1D spatial segmentation partitioned the protein into <b>${domains.length} selective domain(s)</b> and <b>${islands.length} synonymous constraint island(s)</b>.</p>
    <p><b>4. Site &amp; Mechanism Resolution (Stage 5):</b> Active-set polish identified <b>${counts.sites_episodic_diversifying || 0} episodic diversifying</b> (MEME), <b>${counts.sites_pervasive_diversifying || 0} pervasive diversifying</b> (${comp.directional_drift_active ? 'dFEL' : 'FEL'}), <b>${counts.sites_pervasive_purifying || 0} pervasive purifying</b>, and <b>${counts.sites_exchangeability_or_convergent || 0} convergent / amino-acid toggle</b> sites (${counts.sites_simes_fdr_significant || 0} significant at Simes FDR q &le; 0.10).</p>
    ${mca.triggered ? `<div style="margin-top:10px;padding:8px 12px;background:#fef2f2;border:1px solid #fecaca;border-radius:6px;color:#991b1b;font-weight:600;font-size:12px;">&#9888; Mutational Context Alert: ${escHtml(mca.alert_summary)}</div>` : ''}
  `;

  const vp = document.getElementById('bustedVerdictPill');
  vp.className = 'pill ' + (gw.has_positive_selection ? 'pill-pos' : 'pill-neu');
  vp.textContent = gw.has_positive_selection ? `POSITIVE SELECTION (p = ${fmtP(gw.p_value)})` : `NO GENE-WIDE SIGNAL (p = ${fmtP(gw.p_value)})`;

  const omC = gw.omega_classes || [0, 1, 1];
  const omW = gw.omega_weights || [0, 0, 0];
  document.getElementById('bustedBreakdownContainer').innerHTML = `
    <div class="kv-grid" style="margin-bottom:12px;">
      <div class="kv-box"><div class="k">BUSTED-E LRT (p-value)</div><div class="v">${fmtNum(gw.lrt, 2)} (p = ${fmtP(gw.p_value)})</div></div>
      <div class="kv-box"><div class="k">Error-Sink &omega;<sub>err</sub> (Weight)</div><div class="v">${fmtNum(gw.error_sink_omega, 1)} (${fmtNum((gw.error_sink_weight || 0) * 100, 2)}%)</div></div>
    </div>
    <table style="font-size:12px;">
      <thead><tr><th>Rate Class</th><th class="num">&omega; = dN/dS</th><th class="num">Proportion</th><th>Visual Weight</th></tr></thead>
      <tbody>
        <tr><td><span class="pill pill-neg">Purifying (&omega;<sub>1</sub>)</span></td><td class="num">${fmtNum(omC[0], 3)}</td><td class="num">${fmtNum((omW[0] || 0) * 100, 1)}%</td><td><div style="background:#dbeafe;height:10px;border-radius:4px;width:${Math.min(100, (omW[0] || 0) * 100)}%;"></div></td></tr>
        <tr><td><span class="pill pill-neu">Neutral / Relaxed (&omega;<sub>2</sub>)</span></td><td class="num">${fmtNum(omC[1], 3)}</td><td class="num">${fmtNum((omW[1] || 0) * 100, 1)}%</td><td><div style="background:#cbd5e1;height:10px;border-radius:4px;width:${Math.min(100, (omW[1] || 0) * 100)}%;"></div></td></tr>
        <tr><td><span class="pill pill-pos">Diversifying (&omega;<sub>3</sub>)</span></td><td class="num"><b>${fmtNum(omC[2], 2)}</b></td><td class="num"><b>${fmtNum((omW[2] || 0) * 100, 2)}%</b></td><td><div style="background:#dc2626;height:10px;border-radius:4px;width:${Math.min(100, Math.max(2, (omW[2] || 0) * 100))}%;"></div></td></tr>
      </tbody>
    </table>
  `;
}

function hoverCodon(evt, siteNum) {
  const s = siteByNumber[siteNum];
  if (!s) return;
  const html = `
    <div style="font-weight:800;font-size:13px;margin-bottom:4px;">Codon ${s.site} &bull; ${escHtml(s.classification)}</div>
    <div>Domain: <b>${escHtml(s.spatial_domain)}</b> ${s.in_synonymous_island ? '<span style="color:#4ade80;">[Synonymous Island]</span>' : ''}</div>
    <div>Simes Omnibus: <b>p = ${fmtP(s.simes_omnibus_p)}</b> (BH q = ${fmtNum(s.simes_omnibus_q, 3)})</div>
    <div>Rates: <b>&alpha; = ${fmtNum(s.alpha, 2)}</b>, <b>&beta;<sub>perv</sub> = ${fmtNum(s.beta_pervasive, 2)}</b>, <b>&beta;<sup>+</sup> = ${fmtNum(s.beta_plus, 2)}</b></div>
    <div>FEL p = <b>${fmtP(s.fel_p_value)}</b> &bull; dFEL p = <b>${fmtP(s.dfel_p_value)}</b> &bull; MEME p = <b>${fmtP(s.meme_p_value)}</b></div>
    <div>Convergence: <b>${s.max_convergent_transitions || 0} lineages (&rarr; ${escHtml(s.dominant_target_aa)})</b>, N<sub>eff</sub> = ${fmtNum(s.effective_aa_alphabet, 2)}</div>
    <div>Homology PP = <b>${fmtNum(s.column_homology_posterior, 3)}</b> &bull; Gaps = <b>${fmtNum((s.gap_fraction || 0) * 100, 1)}%</b> ${s.artifact_eliminated ? '<br><b style="color:#f87171;">[Indel Artifact Suppressed]</b>' : ''}</div>
    <div style="margin-top:4px;color:#94a3b8;font-size:11px;">Click to open full Site Forensics Drawer</div>
  `;
  showTooltip(evt, html);
}

// ============================================================================
// 3. TAB 2: STAGE 0-2 INTEGRITY, HOMOLOGY & RETICULATION FORENSICS
// ============================================================================
function renderIntegrityTab() {
  const bridges = DATA.recombination_bridges || [];
  document.getElementById('bridgeCountPill').textContent = `${bridges.length} Bridge(s)`;
  document.getElementById('bridgeCountPill').className = 'pill ' + (bridges.length > 0 ? 'pill-warn' : 'pill-ok');

  if (bridges.length === 0) {
    document.getElementById('bridgesContainer').innerHTML = `<div style="padding:14px;color:var(--muted);">No localized recombination or paralog conversion bridges exceeded the routing threshold (all segments consistent with the single species/gene tree).</div>`;
  } else {
    document.getElementById('bridgesContainer').innerHTML = `
      <div class="table-wrap">
        <table>
          <thead><tr><th>Recipient Clade (u)</th><th>Donor Target (v)</th><th class="num">Codon Window</th><th class="num">Segment &Delta;lnL</th><th class="num">p-value</th><th>SPR Routed</th></tr></thead>
          <tbody>
            ${bridges.map(b => `
              <tr>
                <td><b>${escHtml(b.clade_u)}</b></td>
                <td><b>${escHtml(b.target_v)}</b></td>
                <td class="num">${b.peak_start_codon}–${b.peak_end_codon}</td>
                <td class="num" style="color:#c2410c;font-weight:700;">+${fmtNum(b.segment_delta_log_l, 2)}</td>
                <td class="num">${fmtP(b.p_value)}</td>
                <td><span class="pill ${b.routed_to_spr_tree ? 'pill-warn' : 'pill-neu'}">${b.routed_to_spr_tree ? 'Routed' : 'Audit Only'}</span></td>
              </tr>
            `).join('')}
          </tbody>
        </table>
      </div>
    `;
  }

  const artSites = (DATA.site_atlas || []).filter(s => s.artifact_eliminated);
  const aln = DATA.alignment_and_network || {};
  document.getElementById('artifactCountPill').textContent = `${artSites.length} Suppressed`;
  document.getElementById('artifactCountPill').className = 'pill ' + (artSites.length > 0 ? 'pill-prop' : 'pill-ok');

  document.getElementById('artifactContainer').innerHTML = `
    <div class="kv-grid" style="margin-bottom:10px;">
      <div class="kv-box"><div class="k">Option 1 Joint lnL (MG94 + Indel)</div><div class="v">${fmtNum(aln.option1_joint_log_l, 1)}</div></div>
      <div class="kv-box"><div class="k">Mean / Min Homology Posterior PP</div><div class="v">${fmtNum(aln.option5_mean_column_posterior, 4)} / ${fmtNum(aln.option5_min_column_posterior, 3)}</div></div>
    </div>
    <div style="font-size:12px;color:var(--ink-2);margin-bottom:6px;">
      <b>Indel-Window / Low-Homology False-Positive Sites Eliminated (${artSites.length}):</b>
    </div>
    <div>
      ${artSites.length === 0
        ? `<span style="color:var(--muted);font-size:12px;">Zero alignment-artifact sites detected.</span>`
        : artSites.map(s => `<span class="site-chip" onclick="openSiteDrawer(${s.site})" title="PP=${fmtNum(s.column_homology_posterior,3)}, Gap=${fmtNum(s.gap_fraction*100,1)}%">Codon ${s.site} (PP=${fmtNum(s.column_homology_posterior, 2)})</span>`).join('')}
    </div>
  `;

  // Domains & Synonymous Constraint Islands
  const domains = DATA.spatial_domains || [];
  const islands = DATA.synonymous_constraint_islands || [];
  document.getElementById('domainsAndIslandsContainer').innerHTML = `
    <div style="font-weight:700;font-size:12px;margin-bottom:6px;">1D Selective Domains (${domains.length}):</div>
    <div class="table-wrap" style="margin-bottom:12px;">
      <table>
        <thead><tr><th>Codons</th><th>Domain Type</th><th class="num">Mean &omega;</th><th class="num">Mean &alpha;</th><th class="num">Mean &beta;</th><th class="num">Gap%</th><th class="num">Mean N<sub>eff</sub></th></tr></thead>
        <tbody>
          ${domains.map(d => `
            <tr>
              <td class="num"><b>${d.start_codon}–${d.end_codon}</b></td>
              <td><span class="pill ${(d.domain_type||'').includes('Purifying') ? 'pill-neg' : ((d.domain_type||'').includes('Hyper') ? 'pill-diff' : 'pill-prop')}">${escHtml(d.domain_type)}</span></td>
              <td class="num">${fmtNum(d.mean_omega, 2)}</td>
              <td class="num">${fmtNum(d.mean_alpha, 2)}</td>
              <td class="num">${fmtNum(d.mean_beta, 2)}</td>
              <td class="num">${fmtNum((d.mean_gap_fraction || 0) * 100, 1)}%</td>
              <td class="num">${fmtNum(d.mean_n_eff, 2)}</td>
            </tr>
          `).join('')}
        </tbody>
      </table>
    </div>
    <div style="font-weight:700;font-size:12px;margin-bottom:6px;">Synonymous Constraint Islands (${islands.length}):</div>
    ${islands.length === 0
      ? `<div style="color:var(--muted);font-size:12px;">No statistically significant synonymous constraint islands detected.</div>`
      : `<div class="table-wrap"><table><thead><tr><th>Codons</th><th class="num">Island &alpha;</th><th class="num">Gene &alpha;</th><th class="num">Suppression Ratio</th><th class="num">Z-Score</th><th class="num">p-value</th></tr></thead><tbody>
          ${islands.map(isl => `<tr><td class="num"><b>${isl.start_codon}–${isl.end_codon}</b></td><td class="num">${fmtNum(isl.mean_alpha, 3)}</td><td class="num">${fmtNum(isl.gene_mean_alpha, 2)}</td><td class="num">${fmtNum(isl.suppression_ratio, 3)}</td><td class="num">${fmtNum(isl.z_score, 2)}</td><td class="num">${fmtP(isl.p_value)}</td></tr>`).join('')}
        </tbody></table></div>`}
  `;

  // Mutational Context & Coevolution
  const mca = DATA.mutational_context_alert || {};
  const coev = DATA.coevolution_pairs || [];
  document.getElementById('contextAndCoevContainer').innerHTML = `
    <div class="kv-grid" style="margin-bottom:10px;">
      <div class="kv-box"><div class="k">CpG Selected Sites (Enrichment)</div><div class="v">${mca.cpg_driven_selected_sites || 0} / ${mca.num_selected_sites_evaluated || 0} (${fmtNum(mca.cpg_enrichment_fold, 1)}x, p=${fmtP(mca.cpg_fisher_p)})</div></div>
      <div class="kv-box"><div class="k">APOBEC Motif Selected Sites</div><div class="v">${mca.apobec_driven_selected_sites || 0}</div></div>
    </div>
    <div style="margin-bottom:10px;font-size:12px;">
      <b>Flagged CpG Selected Codons:</b> ${(mca.flagged_cpg_sites || []).length === 0 ? 'None' : mca.flagged_cpg_sites.map(c => `<span class="site-chip" onclick="openSiteDrawer(${c})">Codon ${c}</span>`).join('')}
    </div>
    <div style="font-weight:700;font-size:12px;margin-bottom:6px;">Epistatic Co-Mutating Site Pairs (${coev.length}):</div>
    ${coev.length === 0
      ? `<div style="color:var(--muted);font-size:12px;">No epistatic coevolution pairs reported (enable with <code>--coevolution</code>).</div>`
      : `<div class="table-wrap"><table><thead><tr><th>Site Pair</th><th class="num">Shared Branches</th><th class="num">Cosine r</th><th class="num">p (BH q)</th><th>Top Co-Mutating Lineages</th></tr></thead><tbody>
          ${coev.map(cp => `<tr><td><span class="site-chip" onclick="openSiteDrawer(${cp.site_i})">${cp.site_i}</span> &harr; <span class="site-chip" onclick="openSiteDrawer(${cp.site_j})">${cp.site_j}</span></td><td class="num">${fmtNum(cp.shared_branches, 1)}</td><td class="num">${fmtNum(cp.cosine_similarity, 2)}</td><td class="num">${fmtP(cp.p_value)} (${fmtNum(cp.bh_q_value, 2)})</td><td>${escHtml(cp.top_co_mutating_branches)}</td></tr>`).join('')}
        </tbody></table></div>`}
  `;

  renderAdmissionTable();
}

function renderAdmissionTable() {
  const adm = DATA.input_admission_audit || {};
  const taxa = adm.taxa || [];
  const q = (document.getElementById('admSearch')?.value || '').toLowerCase();
  const st = document.getElementById('admStatusFilter')?.value || 'ALL';

  const filtered = taxa.filter(t => {
    const statusUp = (t.status || '').toUpperCase();
    if (st === 'QUARANTINED' && !statusUp.includes('QUARANTIN')) return false;
    if (st === 'RESCUED' && !t.reverse_complemented && (t.frame_offset || 0) === 0) return false;
    if (st === 'ADMITTED' && statusUp.includes('QUARANTIN')) return false;
    if (q && !`${t.name} ${t.status} ${t.reason}`.toLowerCase().includes(q)) return false;
    return true;
  });

  document.getElementById('admissionTableContainer').innerHTML = `
    <table>
      <thead>
        <tr>
          <th>Sequence / Assembly</th>
          <th>Admission Status</th>
          <th class="num">Codon Length</th>
          <th class="num">Frame Offset (RC)</th>
          <th class="num">Internal Stops</th>
          <th class="num">Max Homology Z</th>
          <th class="num">Max AA Identity</th>
          <th>Diagnostic Reason</th>
        </tr>
      </thead>
      <tbody>
        ${filtered.map(t => {
          const isQ = (t.status || '').toUpperCase().includes('QUARANTIN');
          return `
            <tr>
              <td><b>${escHtml(t.name)}</b></td>
              <td><span class="pill ${isQ ? 'pill-pos' : 'pill-ok'}">${escHtml(t.status)}</span></td>
              <td class="num">${t.codon_length || 0}</td>
              <td class="num">${t.frame_offset || 0}${t.reverse_complemented ? ' (RC)' : ''}</td>
              <td class="num">${t.internal_stops_frame0 || 0}</td>
              <td class="num">${fmtNum(t.max_homology_z, 1)}</td>
              <td class="num">${fmtNum((t.max_aa_identity || 0) * 100, 1)}%</td>
              <td>${escHtml(t.reason || 'Admitted in-frame')}</td>
            </tr>
          `;
        }).join('')}
      </tbody>
    </table>
  `;
}

// ============================================================================
// 4. TAB 3: INTERACTIVE ZERO-DEPENDENCY SVG PHYLOGRAM & LINEAGE TABLE
// ============================================================================
function parseNewick(nwk) {
  if (!nwk) return null;
  let i = 0;
  const n = nwk.length;
  let autoId = 1;
  function skipWs() { while (i < n && /\s/.test(nwk[i])) i++; }
  function parseSubtree() {
    skipWs();
    const node = { children: [], name: '', length: 0.0, tag: '' };
    if (i < n && nwk[i] === '(') {
      i++; // '('
      while (true) {
        node.children.push(parseSubtree());
        skipWs();
        if (i < n && nwk[i] === ',') { i++; continue; }
        if (i < n && nwk[i] === ')') { i++; break; }
        break;
      }
    }
    skipWs();
    // Parse name (and optional {TAG} or [&&NHX...])
    let rawName = '';
    while (i < n && !/[,():;]/.test(nwk[i])) {
      if (nwk[i] === '[') {
        while (i < n && nwk[i] !== ']') i++;
        if (i < n) i++;
      } else {
        rawName += nwk[i++];
      }
    }
    rawName = rawName.trim();
    const tagMatch = rawName.match(/\{([^}]+)\}/);
    if (tagMatch) {
      node.tag = tagMatch[1];
      rawName = rawName.replace(/\{[^}]+\}/g, '').trim();
    }
    node.name = rawName || `Node_${autoId++}`;
    skipWs();
    if (i < n && nwk[i] === ':') {
      i++;
      let lenStr = '';
      while (i < n && !/[,();]/.test(nwk[i])) {
        if (nwk[i] === '[') {
          while (i < n && nwk[i] !== ']') i++;
          if (i < n) i++;
        } else {
          lenStr += nwk[i++];
        }
      }
      const v = parseFloat(lenStr);
      node.length = isNaN(v) ? 0.0 : Math.max(0.0, v);
    }
    return node;
  }
  const root = parseSubtree();
  root.isRoot = true;
  return root;
}

const parsedTreeRoot = parseNewick((DATA.alignment_and_network || {}).tree_newick || (DATA.alignment_and_network || {}).extended_newick || '');

function gc3Color(gcPct) {
  // Smooth blue (30% GC3) -> teal (55%) -> amber (75%) -> crimson (95%)
  const t = Math.max(0, Math.min(1, (gcPct - 30) / 65));
  if (t < 0.4) {
    const u = t / 0.4;
    const r = Math.round(37 + u * (13 - 37));
    const g = Math.round(99 + u * (148 - 99));
    const b = Math.round(235 + u * (136 - 235));
    return `rgb(${r},${g},${b})`;
  } else {
    const u = (t - 0.4) / 0.6;
    const r = Math.round(13 + u * (220 - 13));
    const g = Math.round(148 + u * (38 - 148));
    const b = Math.round(136 + u * (38 - 136));
    return `rgb(${r},${g},${b})`;
  }
}

function renderPhylogram() {
  if (!parsedTreeRoot) return;
  const mode = document.getElementById('treeColorMode')?.value || 'gc3';
  const showBridges = document.getElementById('treeShowBridges')?.checked ?? true;
  const useCladogram = document.getElementById('treeCladogram')?.checked ?? false;
  const rowH = parseInt(document.getElementById('treeRowHeight')?.value || '15', 10);
  const searchQ = (document.getElementById('branchSearch')?.value || '').toLowerCase();

  // Legend update
  const leg = document.getElementById('treeLegendBar');
  if (mode === 'gc3') {
    leg.innerHTML = `<b>GC3% Color Scale:</b> <span style="color:#2563eb;font-weight:700;">&#9644; 30% (AT-Rich)</span> <span style="color:#0d9488;font-weight:700;">&#9644; 55% (Mid)</span> <span style="color:#ea580c;font-weight:700;">&#9644; 75% (GC-Rich)</span> <span style="color:#dc2626;font-weight:700;">&#9644; 95% (Ultra-GC3)</span>`;
  } else if (mode === 'absrel') {
    leg.innerHTML = `<b>aBSREL Episodic Selection:</b> <span style="color:#dc2626;font-weight:700;">&#9644; Holm-corrected q &le; 0.05</span> <span style="color:#ea580c;font-weight:700;">&#9644; Nominal p &le; 0.05</span> <span style="color:#94a3b8;font-weight:700;">&#9644; Background (p &gt; 0.05)</span>`;
  } else if (mode === 'relax') {
    leg.innerHTML = `<b>RELAX Selection Intensity K<sub>b</sub>:</b> <span style="color:#dc2626;font-weight:700;">&#9644; Intensified (K &gt; 1, p &le; 0.05)</span> <span style="color:#2563eb;font-weight:700;">&#9644; Relaxed (K &lt; 1, p &le; 0.05)</span> <span style="color:#94a3b8;font-weight:700;">&#9644; Neutral / Unchanged</span>`;
  } else if (mode === 'lof') {
    leg.innerHTML = `<b>Core-Invasion / Loss-of-Constraint Audit:</b> <span style="color:#dc2626;font-weight:700;">&#9644; LoF Alert (BF &ge; 25 &amp; Bonferroni p &le; 0.05)</span> <span style="color:#f59e0b;font-weight:700;">&#9644; Elevated Core Subs (BF &ge; 10)</span> <span style="color:#94a3b8;font-weight:700;">&#9644; Intact Core Constraint</span>`;
  } else {
    leg.innerHTML = `<b>Macro-Regime Partition:</b> <span style="color:#7c3aed;font-weight:700;">&#9644; Test / Foreground Clade</span> <span style="color:#0d9488;font-weight:700;">&#9644; Reference Clade</span>`;
  }

  // Assign leaf y-positions and x-depths
  let leafIdx = 0;
  let maxDepth = 0;
  const nodeByName = {};

  function layoutPass(nd, depthSubs, depthClado) {
    const bp = branchByName[nd.name];
    const bl = bp ? bp.branch_length_subs : nd.length;
    nd.effLen = useCladogram ? (nd.isRoot ? 0 : 1.0) : (nd.isRoot ? 0 : Math.max(0.0005, bl));
    nd.xStart = depthSubs;
    nd.xEnd = depthSubs + nd.effLen;
    if (nd.xEnd > maxDepth) maxDepth = nd.xEnd;
    nodeByName[nd.name] = nd;

    if (!nd.children || nd.children.length === 0) {
      nd.y = leafIdx * rowH + 24;
      leafIdx++;
    } else {
      nd.children.forEach(ch => layoutPass(ch, nd.xEnd, depthClado + 1));
      const ys = nd.children.map(c => c.y);
      nd.y = 0.5 * (Math.min(...ys) + Math.max(...ys));
    }
  }
  layoutPass(parsedTreeRoot, 0, 0);

  const numLeaves = Math.max(1, leafIdx);
  const W = 1320;
  const H = Math.max(280, numLeaves * rowH + 44);
  const padL = 24, padR = 260;
  const treeW = W - padL - padR;
  const xScale = maxDepth > 0 ? treeW / maxDepth : 1;

  let svg = `<svg viewBox="0 0 ${W} ${H}" width="100%" height="${H}" style="display:block;">`;

  function branchStroke(nd) {
    const bp = branchByName[nd.name];
    if (!bp) return { col: '#94a3b8', w: 1.6 };
    if (mode === 'gc3') {
      return { col: gc3Color((bp.gc_cp3 || 0.5) * 100), w: 2.3 };
    } else if (mode === 'absrel') {
      if (bp.absrel_holm_p <= 0.05) return { col: '#dc2626', w: 3.4 };
      if (bp.absrel_p_value <= 0.05) return { col: '#ea580c', w: 2.8 };
      return { col: '#cbd5e1', w: 1.5 };
    } else if (mode === 'relax') {
      if ((bp.relax_call || '').includes('Intensified')) return { col: '#dc2626', w: 3.0 };
      if ((bp.relax_call || '').includes('Relaxed')) return { col: '#2563eb', w: 3.0 };
      return { col: '#cbd5e1', w: 1.5 };
    } else if (mode === 'lof') {
      if (bp.lof_flagged) return { col: '#dc2626', w: 3.6 };
      if ((bp.lof_bayes_factor || 0) >= 10) return { col: '#f59e0b', w: 2.6 };
      return { col: '#cbd5e1', w: 1.5 };
    } else {
      const isTest = (bp.macro_regime || '').toLowerCase().includes('test') || (bp.macro_regime || '').toLowerCase().includes('fg');
      return { col: isTest ? '#7c3aed' : '#0d9488', w: isTest ? 2.8 : 1.8 };
    }
  }

  function drawSubtree(nd) {
    const x1 = padL + nd.xStart * xScale;
    const x2 = padL + nd.xEnd * xScale;
    const y = nd.y;
    const st = branchStroke(nd);
    const isSel = selectedBranchName === nd.name;
    const isMatch = searchQ && nd.name.toLowerCase().includes(searchQ);

    if (!nd.isRoot) {
      svg += `<line x1="${x1}" y1="${y}" x2="${x2}" y2="${y}" stroke="${isSel ? '#0f172a' : st.col}" stroke-width="${isSel ? st.w + 1.8 : st.w}" stroke-linecap="round" style="cursor:pointer;" onmousemove="hoverBranch(event, '${escHtml(nd.name)}')" onmouseout="hideTooltip()" onclick="selectBranch('${escHtml(nd.name)}')"/>`;
      svg += `<circle cx="${x2}" cy="${y}" r="${isSel || isMatch ? 4.2 : 2.4}" fill="${isMatch ? '#ea580c' : st.col}" style="cursor:pointer;" onmousemove="hoverBranch(event, '${escHtml(nd.name)}')" onmouseout="hideTooltip()" onclick="selectBranch('${escHtml(nd.name)}')"/>`;
    }

    if (nd.children && nd.children.length > 0) {
      const yMin = Math.min(...nd.children.map(c => c.y));
      const yMax = Math.max(...nd.children.map(c => c.y));
      svg += `<line x1="${x2}" y1="${yMin}" x2="${x2}" y2="${yMax}" stroke="#94a3b8" stroke-width="1.4"/>`;
      nd.children.forEach(drawSubtree);
    } else {
      const bp = branchByName[nd.name];
      let badge = '';
      if (bp) {
        if (mode === 'gc3') badge = ` (${fmtNum((bp.gc_cp3 || 0) * 100, 0)}% GC3)`;
        else if (mode === 'absrel' && bp.absrel_p_value <= 0.05) badge = ` [ω+=${fmtNum(bp.omega_plus, 1)}, p=${fmtP(bp.absrel_p_value)}]`;
        else if (mode === 'relax' && bp.relax_p_value <= 0.05) badge = ` [K=${fmtNum(bp.relax_k, 2)}]`;
        else if (mode === 'lof' && bp.lof_flagged) badge = ` [LoF Alert!]`;
      }
      const fontW = (isSel || isMatch) ? '800' : '500';
      const fillCol = isMatch ? '#c2410c' : '#0f172a';
      svg += `<text x="${x2 + 6}" y="${y + 3.5}" font-size="${Math.min(12, Math.max(9.5, rowH - 3))}" font-weight="${fontW}" fill="${fillCol}" style="cursor:pointer;" onmousemove="hoverBranch(event, '${escHtml(nd.name)}')" onmouseout="hideTooltip()" onclick="selectBranch('${escHtml(nd.name)}')">${escHtml(nd.name)}${escHtml(badge)}</text>`;
    }
  }
  drawSubtree(parsedTreeRoot);

  // Stage 2 Recombination Bridge Arcs overlay
  if (showBridges) {
    (DATA.recombination_bridges || []).forEach(b => {
      const uNd = nodeByName[b.clade_u];
      const vNd = nodeByName[b.target_v];
      if (uNd && vNd) {
        const xU = padL + uNd.xEnd * xScale, yU = uNd.y;
        const xV = padL + vNd.xEnd * xScale, yV = vNd.y;
        const cx = Math.max(xU, xV) + 45;
        const cy = 0.5 * (yU + yV);
        svg += `<path d="M ${xU} ${yU} Q ${cx} ${cy} ${xV} ${yV}" fill="none" stroke="#ea580c" stroke-width="2.2" stroke-dasharray="5,4"/>`;
      }
    });
  }

  svg += `</svg>`;
  document.getElementById('treeSvgContainer').innerHTML = svg;
}

function hoverBranch(evt, brName) {
  const bp = branchByName[brName];
  if (!bp) {
    showTooltip(evt, `<b>${escHtml(brName)}</b>`);
    return;
  }
  const html = `
    <div style="font-weight:800;font-size:13px;margin-bottom:4px;">Branch: ${escHtml(bp.branch_name)} ${bp.is_leaf ? '(Leaf)' : '(Internal)'}</div>
    <div>Length: <b>${fmtNum(bp.branch_length_subs, 4)} subs/codon</b> &bull; Macro-Regime: <b>${escHtml(bp.macro_regime)}</b></div>
    <div>GC3 Drift Regime: <b>${escHtml(bp.compositional_regime)}</b> (CP3 GC = <b>${fmtNum((bp.gc_cp3 || 0) * 100, 1)}%</b>)</div>
    <div>aBSREL: <b>&omega;<sup>+</sup> = ${fmtNum(bp.omega_plus, 2)}</b> (${fmtNum((bp.weight_plus || 0) * 100, 1)}%), p = <b>${fmtP(bp.absrel_p_value)}</b> (Holm p = ${fmtP(bp.absrel_holm_p)})</div>
    <div>RELAX: <b>K<sub>b</sub> = ${fmtNum(bp.relax_k, 2)}</b>, p = <b>${fmtP(bp.relax_p_value)}</b> (<b>${escHtml(bp.relax_call)}</b>)</div>
    <div>Core Subs: Obs = <b>${fmtNum(bp.core_nonsyn_observed, 1)}</b> vs Exp = <b>${fmtNum(bp.core_nonsyn_expected, 2)}</b>, LoF BF = <b>${fmtNum(bp.lof_bayes_factor, 1)}</b> ${bp.lof_flagged ? '<b style="color:#f87171;">[LoF ALERT]</b>' : ''}</div>
  `;
  showTooltip(evt, html);
}

function selectBranch(brName) {
  selectedBranchName = brName;
  renderPhylogram();
  renderBranchTable();
  const row = document.getElementById('br-row-' + brName);
  if (row) row.scrollIntoView({ behavior: 'smooth', block: 'center' });
}

function sortBranches(key) {
  if (branchSortKey === key) branchSortAsc = !branchSortAsc;
  else { branchSortKey = key; branchSortAsc = (key === 'branch_name' || key.includes('p_value')); }
  renderBranchTable();
}

function renderBranchTable() {
  const q = (document.getElementById('branchSearch')?.value || '').toLowerCase();
  const flt = document.getElementById('branchFilter')?.value || 'ALL';
  let list = (DATA.branch_profiles || []).filter(bp => {
    if (flt === 'ABSREL_SIG' && bp.absrel_p_value > 0.05) return false;
    if (flt === 'LOF_ALERT' && !bp.lof_flagged && (bp.lof_bayes_factor || 0) < 10) return false;
    if (flt === 'RELAX_SHIFT' && bp.relax_p_value > 0.05) return false;
    if (flt === 'LEAVES' && !bp.is_leaf) return false;
    if (flt === 'INTERNAL' && bp.is_leaf) return false;
    if (q && !`${bp.branch_name} ${bp.compositional_regime} ${bp.macro_regime} ${bp.relax_call}`.toLowerCase().includes(q)) return false;
    return true;
  });

  list.sort((a, b) => {
    const va = a[branchSortKey], vb = b[branchSortKey];
    if (typeof va === 'string') return branchSortAsc ? va.localeCompare(vb) : vb.localeCompare(va);
    return branchSortAsc ? (va - vb) : (vb - va);
  });

  document.getElementById('branchTableBody').innerHTML = list.map(bp => {
    const inv = (bp.core_invaded_sites || []).map(c => `<span class="site-chip" onclick="event.stopPropagation(); openSiteDrawer(${c})">${c}</span>`).join('');
    const isSel = selectedBranchName === bp.branch_name;
    return `
      <tr id="br-row-${escHtml(bp.branch_name)}" class="${isSel ? 'row-selected' : ''}" onclick="selectBranch('${escHtml(bp.branch_name)}')">
        <td><b>${escHtml(bp.branch_name)}</b> ${bp.lof_flagged ? '<span class="pill pill-pos">LoF Alert</span>' : ''}</td>
        <td>${escHtml(bp.macro_regime)}</td>
        <td class="num">${fmtNum(bp.branch_length_subs, 4)}</td>
        <td class="num">${fmtNum((bp.gc_cp3 || 0) * 100, 1)}% <span class="pill pill-neu">${escHtml(bp.compositional_regime)}</span></td>
        <td class="num">${fmtNum(bp.omega_plus, 2)} (${fmtNum((bp.weight_plus || 0) * 100, 0)}%)</td>
        <td class="num">${bp.absrel_p_value <= 0.05 ? `<b style="color:#dc2626;">${fmtP(bp.absrel_p_value)}</b>` : fmtP(bp.absrel_p_value)} (${fmtNum(bp.absrel_holm_p, 3)} / ${fmtNum(bp.absrel_bh_q, 3)})</td>
        <td class="num">${fmtNum(bp.relax_k, 2)}</td>
        <td class="num">${fmtP(bp.relax_p_value)} <span class="pill ${(bp.relax_call||'').includes('Intens') ? 'pill-pos' : ((bp.relax_call||'').includes('Relax') ? 'pill-neg' : 'pill-neu')}">${escHtml(bp.relax_call)}</span></td>
        <td class="num">${fmtNum(bp.core_nonsyn_observed, 1)} (${fmtNum(bp.core_nonsyn_expected, 2)})</td>
        <td class="num">${fmtNum(bp.lof_bayes_factor, 1)} (${fmtP(bp.lof_p_value)})</td>
        <td>${inv || '—'}</td>
      </tr>
    `;
  }).join('');
}

// ============================================================================
// 5. TAB 4: STAGE 5 SITE, dFEL & CONVERGENCE EXPLORER
// ============================================================================
function renderSitesTab() {
  const sites = DATA.site_atlas || [];
  const presets = [
    { id: 'ALL', label: `All Codons (${sites.length})` },
    { id: 'SIMES_FDR', label: `Simes FDR q ≤ 0.10 (${sites.filter(s => s.simes_omnibus_q <= 0.10).length})` },
    { id: 'EPISODIC', label: `Episodic Diversifying (${sites.filter(s => (s.classification || '').includes('Episodic')).length})` },
    { id: 'PERVASIVE_POS', label: `Pervasive Diversifying (${sites.filter(s => (s.classification || '').includes('Pervasive') && (s.classification || '').includes('Diversifying')).length})` },
    { id: 'DFEL_DISCORD', label: `dFEL vs FEL Discordant (${sites.filter(isDfelDiscordant).length})` },
    { id: 'ARTIFACT_ELIM', label: `Artifact Eliminated (${sites.filter(s => s.artifact_eliminated).length})` },
    { id: 'CONVERGENT', label: `Convergent / Toggle (${sites.filter(s => (s.max_convergent_transitions || 0) >= 3 || (s.classification || '').includes('Toggle')).length})` },
    { id: 'PURIFYING', label: `Purifying (${sites.filter(s => (s.classification || '').includes('Purifying')).length})` },
    { id: 'CPG_APOBEC', label: `CpG / APOBEC Context (${sites.filter(s => s.cpg_hypermutable_context || s.apobec_motif_context).length})` }
  ];

  document.getElementById('siteFilterPills').innerHTML = presets.map(p => `
    <button class="filter-pill ${currentSiteFilter === p.id ? 'active' : ''}" onclick="setSiteFilter('${p.id}')">${escHtml(p.label)}</button>
  `).join('');

  renderDfelComparisonPlot();
  renderConvergenceLollipopPlot();
  renderSiteTable();
}

function setSiteFilter(fid) {
  currentSiteFilter = fid;
  renderSitesTab();
}

function renderDfelComparisonPlot() {
  const sites = DATA.site_atlas || [];
  const W = 640, H = 260, padL = 54, padR = 24, padT = 20, padB = 42;
  const pw = W - padL - padR, ph = H - padT - padB;
  const maxLrt = Math.max(8.0, ...sites.map(s => Math.min(25, Math.max(s.fel_lrt || 0, s.dfel_lrt || 0))));
  const xOf = v => padL + (Math.min(maxLrt, v) / maxLrt) * pw;
  const yOf = v => padT + ph - (Math.min(maxLrt, v) / maxLrt) * ph;
  const thrX = xOf(3.84), thrY = yOf(3.84);

  let svg = `<svg viewBox="0 0 ${W} ${H}" width="100%" height="auto" style="display:block;">`;
  svg += `<line x1="${padL}" y1="${padT + ph}" x2="${padL + pw}" y2="${padT}" stroke="#cbd5e1" stroke-width="1.2" stroke-dasharray="4,4"/>`;
  svg += `<line x1="${thrX}" y1="${padT}" x2="${thrX}" y2="${padT + ph}" stroke="#dc2626" stroke-width="1" stroke-dasharray="3,3"/>`;
  svg += `<line x1="${padL}" y1="${thrY}" x2="${padL + pw}" y2="${thrY}" stroke="#dc2626" stroke-width="1" stroke-dasharray="3,3"/>`;
  svg += `<text x="${padL + 6}" y="${thrY - 6}" font-size="10" font-weight="700" fill="#15803d">Unmasked by dFEL (dFEL p ≤ 0.05, FEL > 0.05)</text>`;
  svg += `<text x="${thrX + 6}" y="${padT + ph - 8}" font-size="10" font-weight="700" fill="#b45309">Suppressed by dFEL</text>`;

  sites.forEach(s => {
    const fx = xOf(s.fel_lrt || 0);
    const dy = yOf(s.dfel_lrt || 0);
    const unmasked = s.dfel_p_value <= 0.05 && s.fel_p_value > 0.05;
    const suppressed = s.fel_p_value <= 0.05 && s.dfel_p_value > 0.05;
    const isDiv = (s.beta_pervasive || 0) > (s.alpha || 0);
    let col = '#94a3b8', r = 3.2;
    if (unmasked) { col = isDiv ? '#dc2626' : '#2563eb'; r = 5.2; }
    else if (suppressed) { col = '#f59e0b'; r = 4.8; }
    svg += `<circle cx="${fx}" cy="${dy}" r="${r}" fill="${col}" stroke="#0f172a" stroke-width="0.6" opacity="0.88" style="cursor:pointer;" onmousemove="hoverCodon(event, ${s.site})" onmouseout="hideTooltip()" onclick="openSiteDrawer(${s.site})"/>`;
  });

  svg += `<line x1="${padL}" y1="${padT + ph}" x2="${padL + pw}" y2="${padT + ph}" stroke="#334155"/>`;
  svg += `<line x1="${padL}" y1="${padT}" x2="${padL}" y2="${padT + ph}" stroke="#334155"/>`;
  svg += `<text x="${padL + pw / 2}" y="${H - 8}" text-anchor="middle" font-size="11" font-weight="700" fill="#1e293b">Stationary FEL LRT (capped at ${fmtNum(maxLrt, 0)})</text>`;
  svg += `<text x="15" y="${padT + ph / 2}" transform="rotate(-90 15 ${padT + ph / 2})" text-anchor="middle" font-size="11" font-weight="700" fill="#1e293b">Directional dFEL LRT</text>`;
  svg += `</svg>`;
  document.getElementById('dfelPlotContainer').innerHTML = svg;
}

function renderConvergenceLollipopPlot() {
  const sites = DATA.site_atlas || [];
  const C = sites.length || 1;
  const convSites = sites.filter(s => (s.max_convergent_transitions || 0) >= 2);
  const maxConv = Math.max(5, ...convSites.map(s => s.max_convergent_transitions || 0));

  const W = 640, H = 260, padL = 48, padR = 20, padT = 24, padB = 42;
  const pw = W - padL - padR, ph = H - padT - padB;
  const xOf = c => padL + ((c - 0.5) / C) * pw;
  const yOf = m => padT + ph - (m / (maxConv * 1.12)) * ph;

  let svg = `<svg viewBox="0 0 ${W} ${H}" width="100%" height="auto" style="display:block;">`;
  svg += `<line x1="${padL}" y1="${padT + ph}" x2="${padL + pw}" y2="${padT + ph}" stroke="#334155"/>`;
  svg += `<line x1="${padL}" y1="${padT}" x2="${padL}" y2="${padT + ph}" stroke="#334155"/>`;

  // Top 5 convergent sites get automatic callout labels
  const top5 = [...convSites].sort((a, b) => (b.max_convergent_transitions - a.max_convergent_transitions) || (a.exchangeability_p_value - b.exchangeability_p_value)).slice(0, 5);
  const topSet = new Set(top5.map(s => s.site));

  convSites.forEach(s => {
    const xx = xOf(s.site);
    const yy = yOf(s.max_convergent_transitions);
    const nlp = Math.min(10, -Math.log10(Math.max(s.exchangeability_p_value || 1.0, 1e-10)));
    const col = nlp >= 4 ? '#dc2626' : (nlp >= 2 ? '#7c3aed' : '#0d9488');
    svg += `<line x1="${xx}" y1="${padT + ph}" x2="${xx}" y2="${yy}" stroke="#94a3b8" stroke-width="1.2"/>`;
    svg += `<circle cx="${xx}" cy="${yy}" r="${topSet.has(s.site) ? 5.2 : 3.6}" fill="${col}" stroke="#0f172a" stroke-width="0.7" style="cursor:pointer;" onmousemove="hoverCodon(event, ${s.site})" onmouseout="hideTooltip()" onclick="openSiteDrawer(${s.site})"/>`;
    if (topSet.has(s.site)) {
      svg += `<text x="${xx}" y="${Math.max(12, yy - 7)}" text-anchor="middle" font-size="9.5" font-weight="800" fill="#0f172a">C${s.site}→${escHtml(s.dominant_target_aa)} (${s.max_convergent_transitions})</text>`;
    }
  });

  svg += `<text x="${padL + pw / 2}" y="${H - 8}" text-anchor="middle" font-size="11" font-weight="700" fill="#1e293b">Codon Position (1 to ${C})</text>`;
  svg += `<text x="15" y="${padT + ph / 2}" transform="rotate(-90 15 ${padT + ph / 2})" text-anchor="middle" font-size="11" font-weight="700" fill="#1e293b">Convergent Lineages</text>`;
  svg += `</svg>`;
  document.getElementById('convPlotContainer').innerHTML = svg;
}

function getFilteredSortedSites() {
  const q = (document.getElementById('siteSearchInput')?.value || '').toLowerCase();
  let list = (DATA.site_atlas || []).filter(s => {
    const cls = s.classification || '';
    if (currentSiteFilter === 'SIMES_FDR' && s.simes_omnibus_q > 0.10) return false;
    if (currentSiteFilter === 'EPISODIC' && !cls.includes('Episodic')) return false;
    if (currentSiteFilter === 'PERVASIVE_POS' && !(cls.includes('Pervasive') && cls.includes('Diversifying'))) return false;
    if (currentSiteFilter === 'DFEL_DISCORD' && !isDfelDiscordant(s)) return false;
    if (currentSiteFilter === 'ARTIFACT_ELIM' && !s.artifact_eliminated) return false;
    if (currentSiteFilter === 'CONVERGENT' && (s.max_convergent_transitions || 0) < 3 && !cls.includes('Toggle')) return false;
    if (currentSiteFilter === 'PURIFYING' && !cls.includes('Purifying')) return false;
    if (currentSiteFilter === 'CPG_APOBEC' && !s.cpg_hypermutable_context && !s.apobec_motif_context) return false;
    if (q) {
      const hay = `${s.site} ${s.classification} ${s.spatial_domain} ${s.key_property_driver} ${(s.top_episodic_branches || []).join(' ')} ${s.dominant_target_aa}`.toLowerCase();
      if (!hay.includes(q)) return false;
    }
    return true;
  });

  list.sort((a, b) => {
    const va = a[siteSortKey], vb = b[siteSortKey];
    if (typeof va === 'string') return siteSortAsc ? va.localeCompare(vb) : vb.localeCompare(va);
    return siteSortAsc ? (va - vb) : (vb - va);
  });
  return list;
}

function sortSites(key) {
  if (siteSortKey === key) siteSortAsc = !siteSortAsc;
  else { siteSortKey = key; siteSortAsc = (key === 'site' || key.includes('p_value') || key.includes('omnibus')); }
  renderSiteTable();
}

function renderSiteTable() {
  const list = getFilteredSortedSites();
  document.getElementById('siteAtlasTableBody').innerHTML = list.map(s => {
    const pillCls = getClassificationPillClass(s.classification);
    const unmasked = s.dfel_p_value <= 0.05 && s.fel_p_value > 0.05;
    const suppressed = s.fel_p_value <= 0.05 && s.dfel_p_value > 0.05;
    const brs = (s.top_episodic_branches || []).slice(0, 3).join(', ');
    let badges = '';
    if (s.artifact_eliminated) badges += `<span class="pill pill-prop" style="margin-right:4px;">Artifact Suppressed</span>`;
    if (unmasked) badges += `<span class="pill pill-ok" style="margin-right:4px;">dFEL Unmasked</span>`;
    if (suppressed) badges += `<span class="pill pill-warn" style="margin-right:4px;">dFEL Suppressed</span>`;
    if (s.cpg_hypermutable_context) badges += `<span class="pill pill-warn" style="margin-right:4px;">CpG</span>`;
    if (s.apobec_motif_context) badges += `<span class="pill pill-warn" style="margin-right:4px;">APOBEC</span>`;

    return `
      <tr class="${selectedCodon === s.site ? 'row-selected' : ''}" onclick="openSiteDrawer(${s.site})">
        <td class="num"><b>${s.site}</b></td>
        <td><span class="pill ${pillCls}">${escHtml(s.classification)}</span></td>
        <td>${escHtml(s.spatial_domain)}${s.in_synonymous_island ? ' <span class="pill pill-ok">SynIsl</span>' : ''}</td>
        <td class="num">${fmtNum(s.column_homology_posterior, 3)} (${fmtNum((s.gap_fraction || 0) * 100, 0)}%)</td>
        <td class="num">${fmtP(s.simes_omnibus_p)} (${fmtNum(s.simes_omnibus_q, 3)})</td>
        <td class="num">${fmtNum(s.alpha, 2)}</td>
        <td class="num">${fmtNum(s.beta_pervasive, 2)}</td>
        <td class="num">${s.fel_p_value <= 0.05 ? `<b>${fmtP(s.fel_p_value)}</b>` : fmtP(s.fel_p_value)}</td>
        <td class="num">${s.dfel_p_value <= 0.05 ? `<b style="color:#2563eb;">${fmtP(s.dfel_p_value)}</b>` : fmtP(s.dfel_p_value)}</td>
        <td class="num">${s.meme_p_value <= 0.05 ? `<b style="color:#dc2626;">${fmtP(s.meme_p_value)}</b>` : fmtP(s.meme_p_value)} (${fmtNum(s.beta_plus, 1)})</td>
        <td class="num">${fmtNum(s.effective_aa_alphabet, 2)} / ${fmtNum(s.expected_neutral_n_eff, 2)}</td>
        <td class="num">${s.max_convergent_transitions || 0} (${escHtml(s.dominant_target_aa)})</td>
        <td>${badges}${escHtml(brs || s.key_property_driver || '—')}</td>
      </tr>
    `;
  }).join('');
}

// ============================================================================
// 6. SLIDE-OVER SITE FORENSICS DRAWER & EXPORTS
// ============================================================================
function openSiteDrawer(siteNum) {
  const s = siteByNumber[siteNum];
  if (!s) return;
  selectedCodon = siteNum;
  document.getElementById('drawerTitle').innerHTML = `Codon ${s.site} &bull; <span class="pill ${getClassificationPillClass(s.classification)}">${escHtml(s.classification)}</span>`;

  const unmasked = s.dfel_p_value <= 0.05 && s.fel_p_value > 0.05;
  const suppressed = s.fel_p_value <= 0.05 && s.dfel_p_value > 0.05;

  document.getElementById('drawerContent').innerHTML = `
    <div style="margin-bottom:14px;">
      <div style="font-size:11.5px;font-weight:800;text-transform:uppercase;color:var(--muted);">1. Omnibus Verdict &amp; Spatial Context</div>
      <div class="kv-grid">
        <div class="kv-box"><div class="k">Simes Omnibus p (BH q)</div><div class="v">${fmtP(s.simes_omnibus_p)} (q = ${fmtNum(s.simes_omnibus_q, 3)})</div></div>
        <div class="kv-box"><div class="k">1D Selective Domain</div><div class="v">${escHtml(s.spatial_domain)}</div></div>
        <div class="kv-box"><div class="k">Active Topology</div><div class="v">${escHtml(s.active_topology_label || 'Primary Tree')}</div></div>
        <div class="kv-box"><div class="k">POA Homology PP (Gaps)</div><div class="v">${fmtNum(s.column_homology_posterior, 3)} (${fmtNum((s.gap_fraction || 0) * 100, 1)}%)</div></div>
      </div>
      ${s.artifact_eliminated ? `<div style="margin-top:8px;padding:8px 10px;background:#f0fdfa;border:1px solid #5eead4;border-radius:6px;color:#0f766e;font-weight:700;font-size:12px;">&#10003; Alignment / Indel Error-Sink Mitigation: Spurious rigid-alignment positive selection call was suppressed at this codon.</div>` : ''}
    </div>

    <div style="margin-bottom:14px;">
      <div style="font-size:11.5px;font-weight:800;text-transform:uppercase;color:var(--muted);">2. Pervasive Selection: Stationary FEL vs. Non-Reversible dFEL</div>
      <div class="kv-grid">
        <div class="kv-box"><div class="k">Synonymous Rate &alpha;</div><div class="v">${fmtNum(s.alpha, 3)}</div></div>
        <div class="kv-box"><div class="k">Pervasive Nonsyn &beta;<sub>perv</sub></div><div class="v">${fmtNum(s.beta_pervasive, 3)}</div></div>
        <div class="kv-box"><div class="k">Stationary FEL LRT (p / q)</div><div class="v">${fmtNum(s.fel_lrt, 2)} (${fmtP(s.fel_p_value)} / ${fmtNum(s.fel_bh_q, 2)})</div></div>
        <div class="kv-box"><div class="k">Directional dFEL LRT (p / q)</div><div class="v">${fmtNum(s.dfel_lrt, 2)} (${fmtP(s.dfel_p_value)} / ${fmtNum(s.dfel_bh_q, 2)})</div></div>
      </div>
      ${unmasked ? `<div style="margin-top:8px;padding:7px 10px;background:#eff6ff;border:1px solid #93c5fd;border-radius:6px;color:#1d4ed8;font-weight:700;font-size:12px;">&#9733; Unmasked by Non-Reversible dFEL: Missed by stationary FEL (p = ${fmtP(s.fel_p_value)}), significant after GC3 drift correction (p = ${fmtP(s.dfel_p_value)}).</div>` : ''}
      ${suppressed ? `<div style="margin-top:8px;padding:7px 10px;background:#fffbeb;border:1px solid #fcd34d;border-radius:6px;color:#b45309;font-weight:700;font-size:12px;">&#9888; Stationary FEL False Call Suppressed by dFEL: Stationary FEL p = ${fmtP(s.fel_p_value)} &rarr; Non-reversible dFEL p = ${fmtP(s.dfel_p_value)}.</div>` : ''}
    </div>

    <div style="margin-bottom:14px;">
      <div style="font-size:11.5px;font-weight:800;text-transform:uppercase;color:var(--muted);">3. Episodic Diversifying Selection (MEME)</div>
      <div class="kv-grid">
        <div class="kv-box"><div class="k">Episodic &beta;<sup>+</sup> (Weight w<sup>+</sup>)</div><div class="v">${fmtNum(s.beta_plus, 2)} (${fmtNum((s.weight_plus || 0) * 100, 1)}%)</div></div>
        <div class="kv-box"><div class="k">MEME LRT (p / BH q)</div><div class="v">${fmtNum(s.meme_lrt, 2)} (${fmtP(s.meme_p_value)} / ${fmtNum(s.meme_bh_q, 2)})</div></div>
        <div class="kv-box"><div class="k">Episodic Branches (EBF &ge; 10)</div><div class="v">${s.num_episodic_branches || 0} (Max EBF = ${fmtNum(s.max_branch_ebf, 1)})</div></div>
        <div class="kv-box"><div class="k">Background &beta;<sup>&minus;</sup></div><div class="v">${fmtNum(s.beta_minus, 3)}</div></div>
      </div>
      <div style="margin-top:8px;font-size:12px;">
        <b>Driving Branches (click to locate on Phylogram):</b><br>
        ${(s.top_episodic_branches || []).length === 0 ? '<span style="color:var(--muted);">None with EBF &ge; 10</span>' : s.top_episodic_branches.map(b => `<span class="site-chip" onclick="closeSiteDrawer(); switchTab('tab-tree'); selectBranch('${escHtml(b)}');">${escHtml(b)}</span>`).join('')}
      </div>
    </div>

    <div style="margin-bottom:14px;">
      <div style="font-size:11.5px;font-weight:800;text-transform:uppercase;color:var(--muted);">4. Exchangeability, Toggles &amp; Convergent Evolution</div>
      <div class="kv-grid">
        <div class="kv-box"><div class="k">Exchangeability Mode</div><div class="v">${escHtml(s.exchangeability_mode || '—')}</div></div>
        <div class="kv-box"><div class="k">Nonsyn Substitutions (M<sub>s</sub>)</div><div class="v">${fmtNum(s.expected_nonsyn_subs, 1)}</div></div>
        <div class="kv-box"><div class="k">Observed N<sub>eff</sub> vs Neutral E[N<sub>eff</sub>]</div><div class="v">${fmtNum(s.effective_aa_alphabet, 2)} vs ${fmtNum(s.expected_neutral_n_eff, 2)}</div></div>
        <div class="kv-box"><div class="k">Max Convergent Lineages (AA)</div><div class="v">${s.max_convergent_transitions || 0} (&rarr; ${escHtml(s.dominant_target_aa)})</div></div>
        <div class="kv-box"><div class="k">Pair Concentration</div><div class="v">${fmtNum((s.pair_concentration || 0) * 100, 1)}%</div></div>
        <div class="kv-box"><div class="k">Exchangeability p (BH q)</div><div class="v">${fmtP(s.exchangeability_p_value)} (${fmtNum(s.exchangeability_bh_q, 2)})</div></div>
      </div>
      <div style="margin-top:8px;font-size:12px;"><b>Mechanism Summary:</b> ${escHtml(s.key_property_driver || 'Consistent with neutral amino-acid breadth')}</div>
    </div>
  `;

  document.getElementById('drawerBackdrop').classList.add('open');
  document.getElementById('siteDrawer').classList.add('open');
}

function closeSiteDrawer() {
  document.getElementById('drawerBackdrop').classList.remove('open');
  document.getElementById('siteDrawer').classList.remove('open');
}

function stepSiteDrawer(delta) {
  if (!selectedCodon) return;
  const C = (DATA.site_atlas || []).length;
  let next = selectedCodon + delta;
  if (next < 1) next = C;
  if (next > C) next = 1;
  openSiteDrawer(next);
}

function exportSiteCSV() {
  const list = getFilteredSortedSites();
  const cols = ['site','classification','spatial_domain','column_homology_posterior','gap_fraction','simes_omnibus_p','simes_omnibus_q','alpha','beta_pervasive','fel_lrt','fel_p_value','dfel_lrt','dfel_p_value','beta_plus','weight_plus','meme_lrt','meme_p_value','effective_aa_alphabet','expected_neutral_n_eff','max_convergent_transitions','dominant_target_aa','exchangeability_p_value','artifact_eliminated'];
  const rows = [cols.join(',')];
  list.forEach(s => {
    rows.push(cols.map(k => JSON.stringify(s[k] ?? '')).join(','));
  });
  const blob = new Blob([rows.join('\n')], { type: 'text/csv;charset=utf-8;' });
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = `${DATASET_LABEL}.walrus_sites.csv`;
  a.click();
}

function downloadEmbeddedJSON() {
  const raw = document.getElementById('walrus-data').textContent;
  const blob = new Blob([raw], { type: 'application/json' });
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = `${DATASET_LABEL}.CHARACTERIZE.json`;
  a.click();
}

// Initialize all views on load
renderHeaderAndRibbon();
renderOverviewTab();
renderIntegrityTab();
renderBranchTable();
renderSitesTab();
</script>
</body>
</html>
)WALRUS_HTML";
        return html.str();
    }

    static std::string render_standalone_html(const CharacterizeResult& res, const std::string& dataset_label) {
        return render_standalone_html_from_json_string(to_json(res), dataset_label);
    }

private:
    CharacterizeResult run_unified_tensor_characterization(JointAlignResult&& ja_res) {
        CharacterizeResult res;
        res.joint_align = std::move(ja_res);

        const Alignment& aln = res.joint_align.map_alignment;
        const Tree& tree = res.joint_align.optimized_tree;
        size_t N = aln.num_taxa;
        size_t C = aln.num_codons;
        const auto& gcode = *(aln.code ? aln.code : GeneticCode::universal());
        int S = gcode.num_sense_codons;

        // ---------------------------------------------------------------------
        // Step 1: Objective Alignment Quality Evaluation (Options 1, 3 & 5)
        //         + Ohmic Wheatstone Recombination Segment Localization & Routing
        // ---------------------------------------------------------------------
        if (!settings.quiet) {
            Panel::print_step(
                1, 4,
                "Stage 1-3 (Phase I): Homology, Indel & Reticulation Audit",
                "Options 1/3/5 Alignment Objective • Ohmic Wheatstone Segment Routing"
            );
        }
        {
            std::vector<std::pair<std::string, Alignment>> cands = {{"MAP_Alignment", aln}};
            JointAlignSettings eval_s;
            eval_s.quiet = true;
            auto evals = JointAlignEngine::evaluate_candidate_alignments(cands, tree, eval_s);
            if (!evals.empty()) {
                res.alignment_objective = evals[0];
                if (!res.joint_align.site_homology_posteriors.empty()) {
                    Scalar sum_pp = 0.0, min_pp = 1.0;
                    for (Scalar p : res.joint_align.site_homology_posteriors) {
                        sum_pp += p;
                        min_pp = std::min(min_pp, p);
                    }
                    res.alignment_objective.mean_column_posterior =
                        sum_pp / static_cast<Scalar>(res.joint_align.site_homology_posteriors.size());
                    res.alignment_objective.min_column_posterior = min_pp;
                }
            }
        }

        // Fit shared baseline GTR + MG94 parameters once
        GTRFitter gtr_fitter(tree, aln);
        auto gtr_fit = gtr_fitter.fit();
        MG94Parameters mg_params;
        mg_params.theta_AC = gtr_fit.params.theta_AC;
        mg_params.theta_AT = gtr_fit.params.theta_AT;
        mg_params.theta_CG = gtr_fit.params.theta_CG;
        mg_params.theta_CT = gtr_fit.params.theta_CT;
        mg_params.theta_GT = gtr_fit.params.theta_GT;

        MG94Fitter mg_fitter(gtr_fit.tree, aln);
        auto mg_fit = mg_fitter.fit_omega_and_scale(1.0, mg_params);
        mg_params.alpha = 1.0;
        mg_params.beta = mg_fit.x_opt(0);
        Tree work_tree = mg_fit.tree;

        MG94Matrix mg94_base;
        mg94_base.update(mg_params, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, gcode);
        Scalar conv = (mg94_base.scale_factor > 1e-12) ? (3.0 / mg94_base.scale_factor) : 1.0;
        Vector pi_codon = aln.codon_frequencies_f3x4;
        std::vector<size_t> leaf_to_taxon = LikelihoodEngine::build_leaf_to_taxon_map(work_tree, aln);

        // Decision 3: Localize Ohmic Wheatstone bridges and store alternate SPR topologies T_{u->v}
        // for any segment with p <= 0.01 so Inside-Outside cavities route through T_{u->v}!
        struct RoutedSegmentTree {
            Tree spr_tree;
            std::vector<size_t> leaf_to_taxon;
            std::string label;
            size_t start_c = 0;
            size_t end_c = 0;
            Scalar segment_delta_log_l = 0.0;
        };
        std::vector<RoutedSegmentTree> routed_spr_trees;
        std::vector<int> site_routed_tree_idx(C, -1);

        if (N >= 4 && C >= 10) {
            std::vector<Matrix> P_base(work_tree.num_nodes());
            for (const auto& nd : work_tree.nodes) {
                if (nd.id != work_tree.root_id && nd.parent_id != INVALID_INDEX) {
                    P_base[nd.id] = mg94_base.transition_matrix(nd.branch_length * conv);
                }
            }
            std::vector<SitePattern> site_pats(C);
            for (size_t c = 0; c < C; ++c) {
                site_pats[c].states.resize(N);
                for (size_t t = 0; t < N; ++t) site_pats[c].states[t] = aln.sequences[t].codons[c];
                site_pats[c].weight = 1;
            }
            std::vector<Scalar> base_site_ll(C, 0.0);
            std::vector<Vector> node_L(work_tree.num_nodes(), Vector::Zero(S));
            std::vector<double> node_sc(work_tree.num_nodes(), 0.0);
            for (size_t c = 0; c < C; ++c) {
                base_site_ll[c] = FELAnalyzer::compute_pattern_log_likelihood(
                    work_tree, site_pats[c], leaf_to_taxon, mg94_base, P_base, node_L, node_sc
                );
            }

            // Build work_tree-native candidate SPR moves (u -> v):
            // Map each node on work_tree to its representative leaf taxon index and subtree leaf set.
            std::vector<std::vector<size_t>> wt_leaves(work_tree.num_nodes());
            for (int32_t nid : work_tree.post_order) {
                if (nid < 0 || static_cast<size_t>(nid) >= work_tree.num_nodes()) continue;
                if (work_tree.nodes[nid].is_leaf) {
                    if (nid < static_cast<int32_t>(leaf_to_taxon.size()) && leaf_to_taxon[nid] < N) {
                        wt_leaves[nid].push_back(leaf_to_taxon[nid]);
                    }
                } else {
                    for (int32_t ch : work_tree.nodes[nid].children) {
                        if (ch >= 0 && static_cast<size_t>(ch) < wt_leaves.size()) {
                            wt_leaves[nid].insert(wt_leaves[nid].end(), wt_leaves[ch].begin(), wt_leaves[ch].end());
                        }
                    }
                }
            }

            // Compute pairwise tree path distance between nodes on work_tree
            const Matrix& Omega_mat = res.joint_align.ohmic_network.Omega;
            bool has_omega = (Omega_mat.rows() == static_cast<Eigen::Index>(N) &&
                              Omega_mat.cols() == static_cast<Eigen::Index>(N));

            struct BridgeProposal {
                int32_t u_id;
                int32_t v_id;
                Scalar score;
            };
            std::vector<BridgeProposal> proposals;

            // Source A: Sliding-window codon affinity + Ohmic stretch chords between non-sister nodes on work_tree
            const size_t win_w = std::min<size_t>(30, std::max<size_t>(10, C / 6));
            const size_t win_step = std::max<size_t>(5, win_w / 2);

            for (const auto& u_nd : work_tree.nodes) {
                int32_t u = u_nd.id;
                if (u == work_tree.root_id || u_nd.parent_id == INVALID_INDEX) continue;
                if (wt_leaves[u].empty() || wt_leaves[u].size() > N / 2) continue;

                // Find sibling leaves of u on work_tree
                int32_t p_u = u_nd.parent_id;
                std::vector<size_t> sib_leaves;
                for (int32_t sib : work_tree.nodes[p_u].children) {
                    if (sib != u && sib >= 0 && static_cast<size_t>(sib) < wt_leaves.size()) {
                        sib_leaves.insert(sib_leaves.end(), wt_leaves[sib].begin(), wt_leaves[sib].end());
                    }
                }
                if (sib_leaves.empty()) continue;

                int32_t best_v = INVALID_INDEX;
                Scalar best_v_score = 0.0;

                for (const auto& v_nd : work_tree.nodes) {
                    int32_t v = v_nd.id;
                    if (v == work_tree.root_id || v == u || v == p_u || v_nd.parent_id == INVALID_INDEX) continue;
                    if (v_nd.parent_id == p_u) continue; // already sibling
                    if (TreeEngine::is_in_subtree(work_tree, v, u) || TreeEngine::is_in_subtree(work_tree, u, v)) continue;
                    if (wt_leaves[v].empty()) continue;

                    // Find sibling leaves of v on work_tree (to distinguish mobile recipient u from stationary donor v)
                    int32_t p_v = v_nd.parent_id;
                    std::vector<size_t> v_sib_leaves;
                    if (p_v >= 0 && static_cast<size_t>(p_v) < work_tree.num_nodes()) {
                        for (int32_t vsib : work_tree.nodes[p_v].children) {
                            if (vsib != v && vsib >= 0 && static_cast<size_t>(vsib) < wt_leaves.size()) {
                                v_sib_leaves.insert(v_sib_leaves.end(), wt_leaves[vsib].begin(), wt_leaves[vsib].end());
                            }
                        }
                    }

                    // Scan windows for clade-average excess codon identity between Cu and Cv vs Cu and Csib(u),
                    // requiring every leaf in Cu to participate so parent clades do not mask a single recombinant leaf.
                    Scalar max_win_gain = 0.0;
                    for (size_t ws = 0; ws + win_w <= C; ws += win_step) {
                        Scalar min_leaf_gain = 1e9;
                        Scalar sum_leaf_gain = 0.0;
                        bool all_leaves_valid = true;
                        for (size_t lu : wt_leaves[u]) {
                            Scalar sum_uv = 0.0, cnt_uv = 0.0;
                            for (size_t lv : wt_leaves[v]) {
                                for (size_t c = ws; c < ws + win_w; ++c) {
                                    int cu = aln.sequences[lu].codons[c];
                                    int cv = aln.sequences[lv].codons[c];
                                    if (cu >= 0 && cv >= 0) {
                                        cnt_uv += 1.0;
                                        if (cu == cv) sum_uv += 1.0;
                                    }
                                }
                            }
                            Scalar sum_us = 0.0, cnt_us = 0.0;
                            for (size_t ls : sib_leaves) {
                                for (size_t c = ws; c < ws + win_w; ++c) {
                                    int cu = aln.sequences[lu].codons[c];
                                    int cs = aln.sequences[ls].codons[c];
                                    if (cu >= 0 && cs >= 0) {
                                        cnt_us += 1.0;
                                        if (cu == cs) sum_us += 1.0;
                                    }
                                }
                            }
                            if (cnt_uv < 8.0 || cnt_us < 8.0) {
                                all_leaves_valid = false;
                                break;
                            }
                            Scalar id_uv = sum_uv / cnt_uv;
                            Scalar id_us = sum_us / cnt_us;
                            Scalar lg = id_uv - id_us;
                            min_leaf_gain = std::min(min_leaf_gain, lg);
                            sum_leaf_gain += lg;
                        }
                        if (!all_leaves_valid || min_leaf_gain <= 0.08) continue;

                        Scalar mean_gain = sum_leaf_gain / static_cast<Scalar>(wt_leaves[u].size());
                        // Stationary donor check: if v is still concordant with its own sibling in this window
                        // while u is discordant with its sibling, u -> v is the true recipient -> donor direction.
                        Scalar donor_concordance_bonus = 0.0;
                        if (!v_sib_leaves.empty()) {
                            size_t lv0 = wt_leaves[v][0];
                            size_t lvs0 = v_sib_leaves[0];
                            int m_vv = 0, v_cnt = 0;
                            for (size_t c = ws; c < ws + win_w; ++c) {
                                int cv = aln.sequences[lv0].codons[c];
                                int cvs = aln.sequences[lvs0].codons[c];
                                if (cv >= 0 && cvs >= 0) {
                                    v_cnt++;
                                    if (cv == cvs) m_vv++;
                                }
                            }
                            if (v_cnt >= 8) {
                                donor_concordance_bonus = 0.25 * (static_cast<Scalar>(m_vv) / static_cast<Scalar>(v_cnt));
                            }
                        }
                        Scalar win_score = mean_gain + donor_concordance_bonus;
                        if (win_score > max_win_gain) max_win_gain = win_score;
                    }

                    Scalar ohm_bonus = 0.0;
                    if (has_omega && max_win_gain > 0.0) {
                        size_t rep_u = wt_leaves[u][0];
                        size_t rep_v = wt_leaves[v][0];
                        size_t rep_sib = sib_leaves[0];
                        Scalar om_uv = Omega_mat(static_cast<Eigen::Index>(rep_u), static_cast<Eigen::Index>(rep_v));
                        Scalar om_us = Omega_mat(static_cast<Eigen::Index>(rep_u), static_cast<Eigen::Index>(rep_sib));
                        if (om_uv > 1e-6 && om_uv < 1.25 * om_us) {
                            ohm_bonus = 0.10 * (om_us / om_uv);
                        }
                    }

                    Scalar sc = max_win_gain + ohm_bonus;
                    if (sc > best_v_score) {
                        best_v_score = sc;
                        best_v = v;
                    }
                }
                if (best_v != INVALID_INDEX && best_v_score > 0.15) {
                    proposals.push_back({u, best_v, best_v_score});
                }
            }

            std::sort(proposals.begin(), proposals.end(), [](const auto& a, const auto& b) {
                return a.score > b.score;
            });
            if (proposals.size() > 4) proposals.resize(4);

            const Scalar bonf_factor = static_cast<Scalar>(std::max<size_t>(1, proposals.size()));

            for (const auto& prop : proposals) {
                int32_t u_id = prop.u_id;
                int32_t v_id = prop.v_id;
                auto spr_opt = TreeEngine::apply_spr_move(work_tree, u_id, v_id, /*reuse_collapsed_parent=*/false);
                if (!spr_opt.has_value()) continue;
                Tree br_tree = std::move(spr_opt->tree);
                int32_t w_id = spr_opt->w;

                // Build symmetric 2-point local branch-scale transition matrices (standard 50/50 vs near-tip 85/15 transfer)
                // on both work_tree and br_tree so u_id is not penalized by its backbone-inflated length
                // while rate heterogeneity treats both topologies identically.
                std::vector<Matrix> P_base_short = P_base;
                if (u_id >= 0 && static_cast<size_t>(u_id) < work_tree.num_nodes() &&
                    u_id != work_tree.root_id && work_tree.nodes[u_id].parent_id != INVALID_INDEX) {
                    P_base_short[u_id] = mg94_base.transition_matrix(0.20 * work_tree.nodes[u_id].branch_length * conv);
                    int32_t p_u = work_tree.nodes[u_id].parent_id;
                    if (p_u >= 0 && p_u != work_tree.root_id && work_tree.nodes[p_u].children.size() == 2) {
                        int32_t s_u = (work_tree.nodes[p_u].children[0] == u_id)
                                    ? work_tree.nodes[p_u].children[1] : work_tree.nodes[p_u].children[0];
                        Scalar old_ts = work_tree.nodes[s_u].branch_length;
                        Scalar old_tp = work_tree.nodes[p_u].branch_length;
                        P_base_short[s_u] = mg94_base.transition_matrix(0.15 * old_ts * conv);
                        P_base_short[p_u] = mg94_base.transition_matrix((old_tp + 0.85 * old_ts) * conv);
                    }
                }

                std::vector<Matrix> P_br(br_tree.num_nodes());
                std::vector<Matrix> P_br_short(br_tree.num_nodes());
                for (const auto& nd : br_tree.nodes) {
                    if (nd.id != br_tree.root_id && nd.parent_id != INVALID_INDEX) {
                        P_br[nd.id] = mg94_base.transition_matrix(nd.branch_length * conv);
                        P_br_short[nd.id] = P_br[nd.id];
                    }
                }
                if (u_id >= 0 && static_cast<size_t>(u_id) < br_tree.num_nodes() &&
                    u_id != br_tree.root_id && br_tree.nodes[u_id].parent_id != INVALID_INDEX) {
                    P_br_short[u_id] = mg94_base.transition_matrix(0.20 * br_tree.nodes[u_id].branch_length * conv);
                }
                if (v_id >= 0 && static_cast<size_t>(v_id) < work_tree.num_nodes() &&
                    w_id >= 0 && static_cast<size_t>(w_id) < br_tree.num_nodes()) {
                    Scalar old_tv = std::max<Scalar>(1e-4, work_tree.nodes[v_id].branch_length);
                    P_br_short[v_id] = mg94_base.transition_matrix(0.15 * old_tv * conv);
                    P_br_short[w_id] = mg94_base.transition_matrix(0.85 * old_tv * conv);
                }

                auto l2t_br = LikelihoodEngine::build_leaf_to_taxon_map(br_tree, aln);
                std::vector<Vector> nL_br(br_tree.num_nodes(), Vector::Zero(S));
                std::vector<double> nS_br(br_tree.num_nodes(), 0.0);

                std::vector<Scalar> site_diffs(C, 0.0);
                Scalar best_sum = 0.0, cur_sum = 0.0;
                size_t best_s = 0, best_e = 0, cur_s = 0;
                for (size_t c = 0; c < C; ++c) {
                    Scalar ll_b1 = base_site_ll[c];
                    Scalar ll_b2 = FELAnalyzer::compute_pattern_log_likelihood(
                        work_tree, site_pats[c], leaf_to_taxon, mg94_base, P_base_short, node_L, node_sc
                    );
                    Scalar ll_base_sym = std::max(ll_b1, ll_b2);

                    Scalar ll_r1 = FELAnalyzer::compute_pattern_log_likelihood(
                        br_tree, site_pats[c], l2t_br, mg94_base, P_br, nL_br, nS_br
                    );
                    Scalar ll_r2 = FELAnalyzer::compute_pattern_log_likelihood(
                        br_tree, site_pats[c], l2t_br, mg94_base, P_br_short, nL_br, nS_br
                    );
                    Scalar ll_br_sym = std::max(ll_r1, ll_r2);

                    Scalar diff = ll_br_sym - ll_base_sym;
                    site_diffs[c] = diff;
                    if (cur_sum + diff > 0.0) {
                        cur_sum += diff;
                        if (cur_sum > best_sum) {
                            best_sum = cur_sum;
                            best_s = cur_s;
                            best_e = c;
                        }
                    } else {
                        cur_sum = 0.0;
                        cur_s = c + 1;
                    }
                }

                size_t seg_len = (best_e >= best_s) ? (best_e - best_s + 1) : 0;
                size_t pos_sites = 0;
                for (size_t c = best_s; c <= best_e && c < C; ++c) {
                    if (site_diffs[c] > 0.05) pos_sites++;
                }

                if (best_sum > 3.0 && seg_len >= 8 && pos_sites >= 3) {
                    // Calibrate scan statistic S_max via 2,000 site-order permutations of site_diffs
                    // to control for look-elsewhere random-walk excursions across O(C^2) intervals.
                    constexpr int kNumPerms = 2000;
                    int exceed_count = 0;
                    std::vector<Scalar> perm_diffs = site_diffs;
                    uint64_t rng_state = 0x9E3779B97F4A7C15ULL ^
                        (static_cast<uint64_t>(std::max<int32_t>(0, u_id)) * 1315423911ULL +
                         static_cast<uint64_t>(std::max<int32_t>(0, v_id)) * 2654435761ULL);
                    auto next_u64 = [&]() -> uint64_t {
                        rng_state ^= rng_state >> 12;
                        rng_state ^= rng_state << 25;
                        rng_state ^= rng_state >> 27;
                        return rng_state * 2685821657736338717ULL;
                    };
                    for (int p_it = 0; p_it < kNumPerms; ++p_it) {
                        for (size_t i = C - 1; i > 0; --i) {
                            size_t j = static_cast<size_t>(next_u64() % (i + 1));
                            std::swap(perm_diffs[i], perm_diffs[j]);
                        }
                        Scalar p_best = 0.0, p_cur = 0.0;
                        for (size_t c = 0; c < C; ++c) {
                            p_cur = std::max<Scalar>(0.0, p_cur + perm_diffs[c]);
                            if (p_cur > p_best) p_best = p_cur;
                        }
                        if (p_best >= best_sum - 1e-12) {
                            exceed_count++;
                        }
                    }
                    Scalar p_perm = static_cast<Scalar>(exceed_count + 1) / static_cast<Scalar>(kNumPerms + 1);
                    Scalar p_corr = std::min<Scalar>(1.0, bonf_factor * p_perm);
                    Scalar bic_thresh = 3.0 * std::log(std::max<Scalar>(10.0, static_cast<Scalar>(C)));

                    // Polish local branch lengths (u_id, v_id, w_id) on [best_s, best_e] for downstream cavity routing
                    if (u_id >= 0 && static_cast<size_t>(u_id) < br_tree.num_nodes() &&
                        v_id >= 0 && static_cast<size_t>(v_id) < work_tree.num_nodes() &&
                        w_id >= 0 && static_cast<size_t>(w_id) < br_tree.num_nodes()) {
                        Scalar orig_tu = br_tree.nodes[u_id].branch_length;
                        Scalar orig_tv = std::max<Scalar>(1e-4, work_tree.nodes[v_id].branch_length);
                        Scalar best_seg_ll = -1e30;
                        Scalar best_tu = orig_tu, best_tv = 0.5 * orig_tv, best_tw = 0.5 * orig_tv;
                        for (Scalar v_split : {Scalar(0.15), Scalar(0.50)}) {
                            P_br[v_id] = mg94_base.transition_matrix(v_split * orig_tv * conv);
                            P_br[w_id] = mg94_base.transition_matrix((1.0 - v_split) * orig_tv * conv);
                            for (Scalar mult : {Scalar(0.10), Scalar(0.20), Scalar(0.50), Scalar(1.00)}) {
                                Scalar cand_tu = std::max<Scalar>(1e-4, orig_tu * mult);
                                P_br[u_id] = mg94_base.transition_matrix(cand_tu * conv);
                                Scalar seg_ll = 0.0;
                                for (size_t c = best_s; c <= best_e && c < C; ++c) {
                                    seg_ll += FELAnalyzer::compute_pattern_log_likelihood(
                                        br_tree, site_pats[c], l2t_br, mg94_base, P_br, nL_br, nS_br
                                    );
                                }
                                if (seg_ll > best_seg_ll) {
                                    best_seg_ll = seg_ll;
                                    best_tu = cand_tu;
                                    best_tv = v_split * orig_tv;
                                    best_tw = (1.0 - v_split) * orig_tv;
                                }
                            }
                        }
                        br_tree.nodes[u_id].branch_length = best_tu;
                        br_tree.nodes[v_id].branch_length = best_tv;
                        br_tree.nodes[w_id].branch_length = best_tw;
                    }

                    RecombinationBridgeSegment seg;
                    seg.node_u = u_id;
                    seg.node_v = v_id;
                    seg.clade_u_name = (u_id >= 0 && static_cast<size_t>(u_id) < work_tree.nodes.size())
                                     ? work_tree.nodes[u_id].name : ("Node_" + std::to_string(u_id));
                    seg.target_v_name = (v_id >= 0 && static_cast<size_t>(v_id) < work_tree.nodes.size())
                                      ? work_tree.nodes[v_id].name : ("Node_" + std::to_string(v_id));
                    seg.bridge_conductance = prop.score;
                    seg.peak_start_codon = best_s;
                    seg.peak_end_codon = best_e;
                    seg.segment_delta_log_l = best_sum;
                    seg.p_value = p_corr;
                    if (settings.enable_segment_routing && seg.p_value <= 0.01 && (2.0 * best_sum > bic_thresh)) {
                        // Check if an overlapping segment with higher delta_log_l already covers >50% of [best_s, best_e]
                        size_t overlap_higher = 0;
                        for (size_t c = best_s; c <= best_e && c < C; ++c) {
                            int prev_r = site_routed_tree_idx[c];
                            if (prev_r >= 0 && routed_spr_trees[prev_r].segment_delta_log_l >= best_sum) {
                                overlap_higher++;
                            }
                        }
                        if (overlap_higher * 2 <= seg_len) {
                            seg.routed_to_spr_tree = true;
                            int r_idx = static_cast<int>(routed_spr_trees.size());
                            routed_spr_trees.push_back({
                                br_tree, l2t_br,
                                "SPR(" + seg.clade_u_name + "->" + seg.target_v_name + ")",
                                best_s, best_e, best_sum
                            });
                            for (size_t c = best_s; c <= best_e && c < C; ++c) {
                                int prev_r = site_routed_tree_idx[c];
                                if (prev_r < 0 || best_sum > routed_spr_trees[prev_r].segment_delta_log_l) {
                                    site_routed_tree_idx[c] = r_idx;
                                }
                            }
                        }
                    }
                    res.recombination_segments.push_back(seg);
                }
            }
        }
        if (!settings.quiet) {
            std::ostringstream ss;
            ss << std::fixed << std::setprecision(2)
               << "      Option 1 Joint lnP = " << res.alignment_objective.joint_log_l
               << "  │  Option 3 LOTO CV = " << res.alignment_objective.loto_cv_log_l
               << "  │  Mean PP = " << std::setprecision(4) << res.alignment_objective.mean_column_posterior
               << "  │  Localized Bridges = " << res.recombination_segments.size()
               << " (" << routed_spr_trees.size() << " SPR-routed)\n";
            std::cout << Console::muted(ss.str()) << std::flush;
        }

        // ---------------------------------------------------------------------
        // Step 2: Compositional Heterogeneity & Directional Drift Audit (dFEL)
        //         (Decision 4: Unified 3-Regime Directional Softmax Kernel)
        // ---------------------------------------------------------------------
        if (!settings.quiet) {
            Panel::print_step(
                2, 4,
                "Stage 2 (Phase II): Compositional Drift & Root Polarity Audit",
                "Phylogenetic Cavity Flux & BIC Gate • Unified 3-Regime Softmax Transition Kernel"
            );
        }
        res.composition_test = DirectionalFELAnalyzer::run_composition_test(aln, 0.05);
        res.directional_drift_active = false;

        std::unique_ptr<DirectionalFELAnalyzer> dfel_analyzer;
        if (res.composition_test.is_heterogeneous) {
            DirectionalFELSettings dfs;
            dfs.p_value_threshold = settings.p_value_threshold;
            dfs.root_strategy = "auto";
            dfs.num_regimes = 0; // Auto: requires phylogenetic directional LRT/BIC + cavity flux on the tree
            AutoRootResult root_res = DirectionalFELAnalyzer::auto_select_root(work_tree, aln, gtr_fit.params, dfs);
            dfel_analyzer = std::make_unique<DirectionalFELAnalyzer>(
                root_res.best_tree, aln, mg_params, dfs, res.composition_test, root_res
            );
            res.composition_test = dfel_analyzer->comp_test;
            res.directional_drift_active = (dfel_analyzer->regimes.size() == 3);
        }

        if (res.directional_drift_active && dfel_analyzer) {
            res.chosen_root_branch = dfel_analyzer->root_result.chosen_branch_name;
            res.chosen_root_weight = dfel_analyzer->root_result.chosen_weight;
#ifndef __EMSCRIPTEN__
            dfel_analyzer->run();
#endif

            // Map dfel_analyzer's BFS-indexed branch_regime_weights and branch_regime_ids back onto work_tree's
            // node IDs via leaf-bipartition matching so work_tree node IDs (and Stage 2 routed_spr_trees)
            // are preserved and work_tree branch lengths are never double-scaled by conversion_factor!
            auto build_node_leaf_sets = [&](const Tree& t_in) -> std::vector<std::vector<size_t>> {
                auto l2t = LikelihoodEngine::build_leaf_to_taxon_map(t_in, aln);
                std::vector<std::vector<size_t>> lsets(t_in.num_nodes());
                for (int32_t nid : t_in.post_order) {
                    if (nid < 0 || static_cast<size_t>(nid) >= t_in.num_nodes()) continue;
                    if (t_in.nodes[nid].is_leaf) {
                        if (nid < static_cast<int32_t>(l2t.size()) && l2t[nid] < N) {
                            lsets[nid].push_back(l2t[nid]);
                        }
                    } else {
                        for (int32_t ch : t_in.nodes[nid].children) {
                            if (ch >= 0 && static_cast<size_t>(ch) < lsets.size()) {
                                lsets[nid].insert(lsets[nid].end(), lsets[ch].begin(), lsets[ch].end());
                            }
                        }
                        std::sort(lsets[nid].begin(), lsets[nid].end());
                    }
                }
                return lsets;
            };

            auto wt_lsets = build_node_leaf_sets(work_tree);
            auto df_lsets = build_node_leaf_sets(dfel_analyzer->tree);
            std::vector<int> mapped_ids(work_tree.num_nodes(), 1);
            std::vector<std::vector<Scalar>> mapped_weights(work_tree.num_nodes(), {0.0, 1.0, 0.0});

            for (size_t bid = 0; bid < work_tree.num_nodes(); ++bid) {
                if (static_cast<int32_t>(bid) == work_tree.root_id || wt_lsets[bid].empty()) continue;
                std::vector<size_t> comp_set;
                comp_set.reserve(N);
                for (size_t t = 0; t < N; ++t) {
                    if (!std::binary_search(wt_lsets[bid].begin(), wt_lsets[bid].end(), t)) {
                        comp_set.push_back(t);
                    }
                }
                std::array<Scalar, 3> w_acc = {0.0, 0.0, 0.0};
                Scalar match_cnt = 0.0;
                for (size_t did = 0; did < dfel_analyzer->tree.num_nodes(); ++did) {
                    if (static_cast<int32_t>(did) == dfel_analyzer->tree.root_id || df_lsets[did].empty()) continue;
                    if (df_lsets[did] == wt_lsets[bid] || df_lsets[did] == comp_set) {
                        const auto& dw = dfel_analyzer->branch_regime_weights[did];
                        if (dw.size() == 3) {
                            w_acc[0] += dw[0];
                            w_acc[1] += dw[1];
                            w_acc[2] += dw[2];
                            match_cnt += 1.0;
                        }
                    }
                }
                if (match_cnt > 0.0) {
                    w_acc[0] /= match_cnt;
                    w_acc[1] /= match_cnt;
                    w_acc[2] /= match_cnt;
                    mapped_weights[bid] = {w_acc[0], w_acc[1], w_acc[2]};
                    int r_id = 1;
                    if (w_acc[2] > w_acc[1] && w_acc[2] > w_acc[0]) r_id = 2;
                    else if (w_acc[0] > w_acc[1] && w_acc[0] > w_acc[2]) r_id = 0;
                    mapped_ids[bid] = r_id;
                }
            }
            dfel_analyzer->branch_regime_ids = std::move(mapped_ids);
            dfel_analyzer->branch_regime_weights = std::move(mapped_weights);
        } else {
            dfel_analyzer.reset();
            res.chosen_root_branch = "Stationary (Unrooted)";
            res.chosen_root_weight = 1.0;
        }
        if (!settings.quiet) {
            std::ostringstream ss;
            ss << std::fixed << std::setprecision(3)
               << "      Phylo Dir LRT = " << res.composition_test.phylo_dir_lrt
               << " (p = " << res.composition_test.phylo_dir_p_value
               << ", CP3 p = " << res.composition_test.p_value_cp3 << ")"
               << "  │  Kernel Mode: " << (res.directional_drift_active ? "3-Regime Directional Softmax (AT/Neutral/GC)" : "1-Regime Stationary")
               << "  │  Root: " << res.chosen_root_branch << "\n";
            std::cout << Console::muted(ss.str()) << std::flush;
        }

        // ---------------------------------------------------------------------
        // Step 3: Shared Spectral Rate Grid + Ensemble Inside-Outside Tensor
        //         + Two-Stage Active-Set Continuous Polish (Decisions 1, 2, 3, 4, 5)
        // ---------------------------------------------------------------------
        if (!settings.quiet) {
            Panel::print_step(
                3, 4,
                "Stage 4 (Phase III): Shared Inside-Outside Tensor & Active-Set Polish",
                "Ensemble K=3 POA Marginalization • BUSTED-E + aBSREL + RELAX + MEME + Contrast-FEL"
            );
        }

        const std::vector<Scalar> omega_grid = {
            0.0, 0.02, 0.05, 0.10, 0.20, 0.35, 0.50, 0.75,
            1.00, 1.35, 1.80, 2.50, 4.00, 7.00, 12.00, 25.00
        };
        const size_t G = omega_grid.size();

        struct SpectralDecomp {
            Vector evals;
            Matrix R; // D^{-1/2} V
            Matrix L; // V^T D^{1/2}
        };

        // Decision 4: Build spectral grids for K=1 (Stationary) or K=3 (AT/Neutral/GC Directional Regimes)
        size_t num_k_regimes = (res.directional_drift_active && dfel_analyzer && dfel_analyzer->regimes.size() == 3) ? 3 : 1;
        std::vector<std::vector<SpectralDecomp>> regime_spectral_grids(num_k_regimes, std::vector<SpectralDecomp>(G));
        std::vector<MG94Matrix> m_syn_reg(num_k_regimes), m_non_reg(num_k_regimes);
        std::vector<Vector> sqrt_pi_reg(num_k_regimes), inv_sqrt_pi_reg(num_k_regimes);

        MG94Parameters p_syn = mg_params; p_syn.alpha = 1.0; p_syn.beta = 0.0;
        MG94Parameters p_non = mg_params; p_non.alpha = 0.0; p_non.beta = 1.0;
        MG94Matrix m_syn_base, m_non_base;
        m_syn_base.update(p_syn, aln.pos_nuc_frequencies, pi_codon, gcode);
        m_non_base.update(p_non, aln.pos_nuc_frequencies, pi_codon, gcode);
        const Matrix& Q_non = m_non_base.Q;

        for (size_t k_reg = 0; k_reg < num_k_regimes; ++k_reg) {
            Eigen::Matrix<Scalar, 3, 4> pos_freqs = (num_k_regimes == 3)
                ? dfel_analyzer->regimes[k_reg].pos_nuc_freqs : aln.pos_nuc_frequencies;
            Vector pi_k = (num_k_regimes == 3)
                ? dfel_analyzer->regimes[k_reg].codon_freqs : pi_codon;

            m_syn_reg[k_reg].update(p_syn, pos_freqs, pi_k, gcode);
            m_non_reg[k_reg].update(p_non, pos_freqs, pi_k, gcode);

            sqrt_pi_reg[k_reg] = pi_k.array().max(1e-15).sqrt();
            inv_sqrt_pi_reg[k_reg] = 1.0 / sqrt_pi_reg[k_reg].array();

            for (size_t g = 0; g < G; ++g) {
                Matrix Q_g = m_syn_reg[k_reg].Q + omega_grid[g] * m_non_reg[k_reg].Q;
                for (int i = 0; i < S; ++i) {
                    Q_g(i, i) = 0.0;
                    Q_g(i, i) = -Q_g.row(i).sum();
                }
                Matrix sym = sqrt_pi_reg[k_reg].asDiagonal() * Q_g * inv_sqrt_pi_reg[k_reg].asDiagonal();
                sym = 0.5 * (sym + sym.transpose());
                Eigen::SelfAdjointEigenSolver<Matrix> es(sym);
                regime_spectral_grids[k_reg][g].evals = es.eigenvalues();
                regime_spectral_grids[k_reg][g].R = inv_sqrt_pi_reg[k_reg].asDiagonal() * es.eigenvectors();
                regime_spectral_grids[k_reg][g].L = es.eigenvectors().transpose() * sqrt_pi_reg[k_reg].asDiagonal();
            }
        }

        // Extract active branches (all non-root nodes with a parent)
        std::vector<int32_t> active_branches;
        for (const auto& nd : work_tree.nodes) {
            if (nd.id != work_tree.root_id && nd.parent_id != INVALID_INDEX) {
                active_branches.push_back(nd.id);
            }
        }
        size_t B = active_branches.size();
        size_t num_nodes = work_tree.num_nodes();
        size_t max_tree_nodes = num_nodes;
        for (const auto& rt : routed_spr_trees) {
            max_tree_nodes = std::max(max_tree_nodes, rt.spr_tree.num_nodes());
        }
        const Vector& root_pi = (num_k_regimes == 3 && dfel_analyzer && dfel_analyzer->root_frequencies.size() == S)
                              ? dfel_analyzer->root_frequencies : pi_codon;

        // Decision 4: Helper to evaluate branch transition matrix for grid index g (or continuous ω via log-linear spectral interpolation)
        auto branch_grid_transition_matrix = [&](int32_t bid, size_t g_idx, Scalar eff_t) -> Matrix {
            Scalar t_pos = std::max<Scalar>(0.0, eff_t);
            if (num_k_regimes == 1 || !dfel_analyzer || bid < 0 || bid >= static_cast<int32_t>(dfel_analyzer->branch_regime_weights.size())) {
                Vector exp_d = (regime_spectral_grids[0][g_idx].evals.array() * t_pos).exp();
                Matrix P = regime_spectral_grids[0][g_idx].R * exp_d.asDiagonal() * regime_spectral_grids[0][g_idx].L;
                return P.cwiseMax(1e-15);
            }
            const auto& w_vec = dfel_analyzer->branch_regime_weights[bid];
            Matrix P_mix = Matrix::Zero(S, S);
            for (size_t k_reg = 0; k_reg < 3; ++k_reg) {
                Scalar wk = (k_reg < w_vec.size()) ? w_vec[k_reg] : (k_reg == 1 ? 1.0 : 0.0);
                if (wk <= 1e-6) continue;
                Vector exp_d = (regime_spectral_grids[k_reg][g_idx].evals.array() * t_pos).exp();
                P_mix.noalias() += wk * (regime_spectral_grids[k_reg][g_idx].R * exp_d.asDiagonal() * regime_spectral_grids[k_reg][g_idx].L);
            }
            return P_mix.cwiseMax(1e-15);
        };
        auto branch_grid_apply_P = [&](int32_t bid, size_t g_idx, Scalar eff_t, const Vector& vec) -> Vector {
            Scalar t_pos = std::max<Scalar>(0.0, eff_t);
            if (num_k_regimes == 1 || !dfel_analyzer || bid < 0 || bid >= static_cast<int32_t>(dfel_analyzer->branch_regime_weights.size())) {
                const auto& sd = regime_spectral_grids[0][g_idx];
                Vector exp_d = (sd.evals.array() * t_pos).exp();
                Vector out = sd.R * (exp_d.cwiseProduct(sd.L * vec));
                return out.cwiseMax(1e-15);
            }
            const auto& w_vec = dfel_analyzer->branch_regime_weights[bid];
            Vector out = Vector::Zero(S);
            for (size_t k_reg = 0; k_reg < 3; ++k_reg) {
                Scalar wk = (k_reg < w_vec.size()) ? w_vec[k_reg] : (k_reg == 1 ? 1.0 : 0.0);
                if (wk <= 1e-6) continue;
                const auto& sd = regime_spectral_grids[k_reg][g_idx];
                Vector exp_d = (sd.evals.array() * t_pos).exp();
                out.noalias() += wk * (sd.R * (exp_d.cwiseProduct(sd.L * vec)));
            }
            return out.cwiseMax(1e-15);
        };
        auto branch_grid_apply_PT = [&](int32_t bid, size_t g_idx, Scalar eff_t, const Vector& vec) -> Vector {
            Scalar t_pos = std::max<Scalar>(0.0, eff_t);
            if (num_k_regimes == 1 || !dfel_analyzer || bid < 0 || bid >= static_cast<int32_t>(dfel_analyzer->branch_regime_weights.size())) {
                const auto& sd = regime_spectral_grids[0][g_idx];
                Vector exp_d = (sd.evals.array() * t_pos).exp();
                Vector out = sd.L.transpose() * (exp_d.cwiseProduct(sd.R.transpose() * vec));
                return out.cwiseMax(1e-15);
            }
            const auto& w_vec = dfel_analyzer->branch_regime_weights[bid];
            Vector out = Vector::Zero(S);
            for (size_t k_reg = 0; k_reg < 3; ++k_reg) {
                Scalar wk = (k_reg < w_vec.size()) ? w_vec[k_reg] : (k_reg == 1 ? 1.0 : 0.0);
                if (wk <= 1e-6) continue;
                const auto& sd = regime_spectral_grids[k_reg][g_idx];
                Vector exp_d = (sd.evals.array() * t_pos).exp();
                out.noalias() += wk * (sd.L.transpose() * (exp_d.cwiseProduct(sd.R.transpose() * vec)));
            }
            return out.cwiseMax(1e-15);
        };

        // Compute Inside L_{v, c} and Outside U_{u, c} for every site c
        // incorporating Decision 2 (K=3 POA path ensemble marginalization for PP_c < 0.99 with Bayes-updated posterior weights)
        // and Decision 3 (Segment-specific SPR tree routing for recombinant segments p <= 0.01)
        std::vector<std::vector<Vector>> site_inside(C, std::vector<Vector>(max_tree_nodes, Vector::Ones(S)));
        std::vector<std::vector<Vector>> site_outside(C, std::vector<Vector>(max_tree_nodes, Vector::Zero(S)));
        std::vector<Scalar> site_base_ll(C, 0.0);
        std::vector<Scalar> site_alpha(C, 1.0);
        std::vector<size_t> site_g_near(C, 6);
        std::vector<size_t> site_ensemble_k(C, 1);

        const Alignment& init_aln = res.joint_align.initial_alignment;
        bool has_init_aln = (init_aln.num_taxa == N && init_aln.num_codons == C);

        // Build lookup of local 3-way NNI WheatstoneBridge hubs on work_tree (both explicit 3-child
        // bridges on work_tree.bridges and binary quartet bridges from res.joint_align.ohmic_network.bridges)
        struct LocalWheatstoneHub {
            int32_t hub_u = INVALID_INDEX;
            int32_t bridge_v = INVALID_INDEX;
            std::array<int32_t, 3> children{INVALID_INDEX, INVALID_INDEX, INVALID_INDEX};
            std::array<Scalar, 3> weights{1.0, 0.0, 0.0};
            std::array<Scalar, 3> branch_lengths{0.0, 0.0, 0.0};
        };
        std::vector<LocalWheatstoneHub> wt_wheatstone_hubs(max_tree_nodes);
        auto register_wb = [&](const WheatstoneBridge& wb) {
            int32_t id_A = wb.children[0];
            int32_t id_B = wb.children[1];
            int32_t id_C = wb.children[2];
            if (id_A < 0 || id_B < 0 || id_C < 0 ||
                static_cast<size_t>(id_A) >= work_tree.num_nodes() ||
                static_cast<size_t>(id_B) >= work_tree.num_nodes() ||
                static_cast<size_t>(id_C) >= work_tree.num_nodes()) {
                return;
            }
            int32_t v_id = wb.bridge_node_id;
            int32_t u_id = wb.hub_node_id;
            if (v_id >= 0 && static_cast<size_t>(v_id) < work_tree.num_nodes() &&
                work_tree.nodes[v_id].parent_id != INVALID_INDEX &&
                work_tree.nodes[v_id].children.size() == 2) {
                // Binary work_tree node v_id with parent u_id = parent(v_id)
                u_id = work_tree.nodes[v_id].parent_id;
            }
            if (u_id < 0 || static_cast<size_t>(u_id) >= work_tree.num_nodes()) return;
            if (wt_wheatstone_hubs[u_id].hub_u != INVALID_INDEX) return;

            LocalWheatstoneHub hub;
            hub.hub_u = u_id;
            hub.bridge_v = (v_id >= 0 && static_cast<size_t>(v_id) < work_tree.num_nodes()) ? v_id : u_id;
            hub.children = {id_A, id_B, id_C};
            Scalar w_sum = wb.weights[0] + wb.weights[1] + wb.weights[2];
            if (w_sum > 1e-12) {
                hub.weights = {wb.weights[0] / w_sum, wb.weights[1] / w_sum, wb.weights[2] / w_sum};
            }
            Scalar scale_r = 1.0;
            if (hub.bridge_v != hub.hub_u && wb.branch_lengths[0] > 1e-6 &&
                work_tree.nodes[hub.bridge_v].branch_length > 1e-6) {
                scale_r = work_tree.nodes[hub.bridge_v].branch_length / wb.branch_lengths[0];
            }
            for (int r = 0; r < 3; ++r) {
                hub.branch_lengths[r] = std::max<Scalar>(1e-6, wb.branch_lengths[r] * scale_r);
            }
            wt_wheatstone_hubs[u_id] = hub;
        };
        for (const auto& wb : work_tree.bridges) register_wb(wb);
        for (const auto& wb : res.joint_align.ohmic_network.bridges) register_wb(wb);

        // Map each column c to its Stage 1 indel window and offset (if any) for exact Forward-Backward
        // transducer candidate placements and posterior weights
        std::vector<int> site_to_win(C, -1);
        std::vector<size_t> site_to_win_offset(C, 0);
        for (size_t wi = 0; wi < res.joint_align.windows.size(); ++wi) {
            const auto& win = res.joint_align.windows[wi].window;
            for (size_t c = win.start_codon; c < win.end_codon && c < C; ++c) {
                site_to_win[c] = static_cast<int>(wi);
                site_to_win_offset[c] = c - win.start_codon;
            }
        }

        Scalar total_tree_subs = 0.0;
        for (const auto& nd : work_tree.nodes) {
            if (nd.id != work_tree.root_id && nd.parent_id != INVALID_INDEX) {
                total_tree_subs += nd.branch_length;
            }
        }

        // Empirical Bayes prior distribution over the 16-point ω grid from the gene's own site MLEs
        // (with Laplace Dirichlet smoothing +1 per grid point; purifying/neutral classes g <= 8 have full weight 1.0)
        std::array<Scalar, 16> emp_grid_prior{};
        for (size_t g = 0; g < G; ++g) emp_grid_prior[g] = 1.0;

        for (size_t c = 0; c < C; ++c) {
            Scalar raw_alpha_c = (c < res.joint_align.site_results.size())
                               ? std::max<Scalar>(0.0, res.joint_align.site_results[c].marginalized_alpha) : 1.0;
            Scalar beta_c = (c < res.joint_align.site_results.size())
                          ? std::max<Scalar>(0.0, res.joint_align.site_results[c].marginalized_beta) : mg_params.beta;
            if (dfel_analyzer && c < dfel_analyzer->site_results.size() && site_routed_tree_idx[c] < 0) {
                raw_alpha_c = std::max<Scalar>(0.0, dfel_analyzer->site_results[c].alpha);
                beta_c = std::max<Scalar>(0.0, dfel_analyzer->site_results[c].beta);
            }

            size_t non_gap_taxa = 0;
            for (size_t t = 0; t < N; ++t) {
                if (aln.sequences[t].codons[c] >= 0) non_gap_taxa++;
            }
            Scalar occ_frac = (N > 0) ? (static_cast<Scalar>(non_gap_taxa) / static_cast<Scalar>(N)) : 1.0;
            Scalar t_syn_c = (total_tree_subs * occ_frac) / 3.0;

            // Conjugate Gamma(1,1)-Poisson Empirical Bayes posterior mean E[α_c | S_c = α̂_c * T_{syn,c}]
            Scalar eb_alpha_c = (1.0 + std::min<Scalar>(raw_alpha_c, 50.0) * t_syn_c) / (1.0 + t_syn_c);
            site_alpha[c] = std::clamp<Scalar>(eb_alpha_c, 0.05, 20.0);
            Scalar omega_c = beta_c / site_alpha[c];

            size_t g_mle = 0;
            Scalar best_d = 1e9;
            for (size_t g = 0; g < G; ++g) {
                Scalar d = std::abs(omega_grid[g] - omega_c);
                if (d < best_d) { best_d = d; g_mle = g; }
            }
            emp_grid_prior[g_mle] += 1.0;
            // Background tree cavity (U_bar, L_bar) and reference denominator use the purifying/neutral baseline (g <= 8, ω <= 1)
            site_g_near[c] = std::min<size_t>(8, g_mle);
        }

        Scalar max_emp_grid_count = *std::max_element(emp_grid_prior.begin(), emp_grid_prior.end());
        for (size_t g = 0; g < G; ++g) {
            emp_grid_prior[g] = (g <= 8) ? 1.0 : (emp_grid_prior[g] / std::max<Scalar>(1.0, max_emp_grid_count));
        }

        // Per-site Inside-Outside cavity builder. When node_mask_gamma != nullptr (Pass 1 recalculation on sites
        // where a branch error-sink tract was detected), any edge `ch` with mask gamma_ch > 0 has its upward
        // message blended toward the flat missing-data vector 1: (1 - gamma_ch)*(P_ch * L_ch) + gamma_ch * 1,
        // severing corrupted subtrees so sibling and ancestor Outside vectors are recomputed cleanly.
        auto compute_site_cavity = [&](size_t c, const std::vector<Scalar>* node_mask_gamma) {
            Scalar alpha_c = site_alpha[c];
            size_t g_near = site_g_near[c];

            // Decision 3: Choose active tree for site c (Backbone vs Routed Recombinant SPR Tree)
            int r_idx = site_routed_tree_idx[c];
            const Tree& active_tree = (r_idx >= 0) ? routed_spr_trees[r_idx].spr_tree : work_tree;
            const auto& active_l2t = (r_idx >= 0) ? routed_spr_trees[r_idx].leaf_to_taxon : leaf_to_taxon;

            size_t act_nodes = std::max(max_tree_nodes, active_tree.num_nodes());
            auto apply_site_P = [&](int32_t nid, const Vector& v) -> Vector {
                return branch_grid_apply_P(nid, g_near, alpha_c * active_tree.nodes[nid].branch_length * conv, v);
            };
            auto apply_site_PT = [&](int32_t nid, const Vector& u) -> Vector {
                return branch_grid_apply_PT(nid, g_near, alpha_c * active_tree.nodes[nid].branch_length * conv, u);
            };

            auto edge_mask_gamma = [&](int32_t nid) -> Scalar {
                if (!node_mask_gamma || nid < 0 || static_cast<size_t>(nid) >= node_mask_gamma->size()) return 0.0;
                return (*node_mask_gamma)[nid];
            };
            auto masked_upward_msg = [&](int32_t nid, const Vector& loc_in_nid) -> Vector {
                Vector raw_m = apply_site_P(nid, loc_in_nid);
                Scalar gm = edge_mask_gamma(nid);
                if (gm <= 1e-6) return raw_m;
                return (1.0 - gm) * raw_m + gm * Vector::Ones(S);
            };
            auto masked_upward_sc = [&](int32_t nid, Scalar raw_sc) -> Scalar {
                Scalar gm = edge_mask_gamma(nid);
                return (gm <= 1e-6) ? raw_sc : ((1.0 - gm) * raw_sc);
            };

            // Decision 2: Construct up to K=3 POA path realizations from Stage 1 Forward-Backward transducer
            // candidate placements (`window_taxon_results`) with exact transducer log-posteriors (zero hardcoded weights)
            Scalar pp_c = (c < res.joint_align.site_homology_posteriors.size())
                        ? res.joint_align.site_homology_posteriors[c] : 1.0;
            struct ColumnRealization {
                std::vector<int> codons;
                Scalar prior_weight;
                Scalar log_l = -1e30;
                Scalar post_weight = 0.0;
                std::vector<Vector> loc_in;
                std::vector<Vector> loc_out;
            };
            std::vector<ColumnRealization> realizations;
            std::vector<int> map_col(N, -1);
            for (size_t t = 0; t < N; ++t) map_col[t] = aln.sequences[t].codons[c];

            int wi = site_to_win[c];
            bool has_transducer_cands = (
                pp_c < 0.99 &&
                wi >= 0 &&
                static_cast<size_t>(wi) < res.joint_align.window_taxon_results.size() &&
                res.joint_align.window_taxon_results[wi].size() == N
            );

            if (pp_c >= 0.99) {
                realizations.push_back({map_col, 1.0, -1e30, 1.0, {}, {}});
            } else if (has_transducer_cands) {
                size_t w_off = site_to_win_offset[c];
                const auto& win_tr = res.joint_align.window_taxon_results[wi];
                std::vector<Scalar> cand_log_priors;
                for (size_t rank = 0; rank < 3; ++rank) {
                    std::vector<int> col_r(N, -1);
                    Scalar log_prior_r = 0.0;
                    for (size_t t = 0; t < N; ++t) {
                        const auto& cands = win_tr[t].top_candidates;
                        if (!cands.empty()) {
                            size_t idx = std::min(rank, cands.size() - 1);
                            if (w_off < cands[idx].codons.size()) {
                                col_r[t] = cands[idx].codons[w_off];
                            } else {
                                col_r[t] = map_col[t];
                            }
                            log_prior_r += (cands[idx].log_score - win_tr[t].log_partition_Z);
                        } else {
                            col_r[t] = map_col[t];
                        }
                    }
                    bool duplicate = false;
                    for (const auto& ex : realizations) {
                        if (ex.codons == col_r) { duplicate = true; break; }
                    }
                    if (!duplicate) {
                        realizations.push_back({std::move(col_r), 1.0, -1e30, 0.0, {}, {}});
                        cand_log_priors.push_back(log_prior_r);
                    }
                }
                if (!cand_log_priors.empty()) {
                    Scalar max_lp = *std::max_element(cand_log_priors.begin(), cand_log_priors.end());
                    Scalar sum_pw = 0.0;
                    for (size_t m = 0; m < realizations.size(); ++m) {
                        realizations[m].prior_weight = std::exp(std::clamp<Scalar>(cand_log_priors[m] - max_lp, -50.0, 0.0));
                        sum_pw += realizations[m].prior_weight;
                    }
                    if (sum_pw > 0.0) {
                        for (auto& r : realizations) r.prior_weight /= sum_pw;
                    }
                }
            } else {
                // Fallback when window_taxon_results was compacted: uniform prior (1.0) so posterior weights
                // are governed strictly by the exact phylogenetic column likelihood exp(ell_c^{(m)})
                std::vector<int> alt1_col = map_col;
                std::vector<int> alt2_col = map_col;
                for (size_t t = 0; t < N; ++t) {
                    if (has_init_aln && init_aln.sequences[t].codons[c] != map_col[t]) {
                        alt1_col[t] = init_aln.sequences[t].codons[c];
                    } else if (map_col[t] < 0 && c + 1 < C && aln.sequences[t].codons[c + 1] >= 0) {
                        alt1_col[t] = aln.sequences[t].codons[c + 1];
                    }
                    if (map_col[t] < 0 && c > 0 && aln.sequences[t].codons[c - 1] >= 0) {
                        alt2_col[t] = aln.sequences[t].codons[c - 1];
                    } else if (map_col[t] >= 0 && ((c + 1 < C && aln.sequences[t].codons[c + 1] < 0) || (c > 0 && aln.sequences[t].codons[c - 1] < 0))) {
                        alt2_col[t] = -1;
                    }
                }
                realizations.push_back({map_col, 1.0, -1e30, 0.0, {}, {}});
                if (alt1_col != map_col) {
                    realizations.push_back({alt1_col, 1.0, -1e30, 0.0, {}, {}});
                }
                if (alt2_col != map_col && alt2_col != alt1_col) {
                    realizations.push_back({alt2_col, 1.0, -1e30, 0.0, {}, {}});
                }
            }
            site_ensemble_k[c] = realizations.size();

            Scalar max_realiz_ll = -1e30;
            for (auto& realiz : realizations) {
                realiz.loc_in.assign(act_nodes, Vector::Ones(S));
                realiz.loc_out.assign(act_nodes, Vector::Zero(S));
                std::vector<Scalar> log_sc(act_nodes, 0.0);

                for (int32_t nid : active_tree.post_order) {
                    const auto& nd = active_tree.nodes[nid];
                    if (nd.is_leaf) {
                        size_t tid = (nid < static_cast<int32_t>(active_l2t.size())) ? active_l2t[nid] : INVALID_INDEX;
                        if (tid < N && realiz.codons[tid] >= 0) {
                            realiz.loc_in[nid].setZero();
                            realiz.loc_in[nid](realiz.codons[tid]) = 1.0;
                        } else {
                            realiz.loc_in[nid].setOnes();
                        }
                    } else if (r_idx < 0 && nid >= 0 && static_cast<size_t>(nid) < wt_wheatstone_hubs.size() &&
                               wt_wheatstone_hubs[nid].hub_u == nid) {
                        // Local 3-way NNI WheatstoneBridge integration at hub nid
                        const auto& wh = wt_wheatstone_hubs[nid];
                        int32_t id_A = wh.children[0];
                        int32_t id_B = wh.children[1];
                        int32_t id_C = wh.children[2];
                        int32_t v_id = wh.bridge_v;

                        Vector M_A = masked_upward_msg(id_A, realiz.loc_in[id_A]);
                        Vector M_B = masked_upward_msg(id_B, realiz.loc_in[id_B]);
                        Vector M_C = masked_upward_msg(id_C, realiz.loc_in[id_C]);
                        Scalar sc = masked_upward_sc(id_A, log_sc[id_A])
                                  + masked_upward_sc(id_B, log_sc[id_B])
                                  + masked_upward_sc(id_C, log_sc[id_C]);

                        Vector W_abc = Vector::Zero(S);
                        if (wh.weights[0] > 0.0) {
                            Matrix P0 = branch_grid_transition_matrix(v_id, g_near, alpha_c * wh.branch_lengths[0] * conv);
                            W_abc.noalias() += wh.weights[0] * (P0 * M_A.cwiseProduct(M_B)).cwiseProduct(M_C);
                        }
                        if (wh.weights[1] > 0.0) {
                            Matrix P1 = branch_grid_transition_matrix(v_id, g_near, alpha_c * wh.branch_lengths[1] * conv);
                            W_abc.noalias() += wh.weights[1] * (P1 * M_A.cwiseProduct(M_C)).cwiseProduct(M_B);
                        }
                        if (wh.weights[2] > 0.0) {
                            Matrix P2 = branch_grid_transition_matrix(v_id, g_near, alpha_c * wh.branch_lengths[2] * conv);
                            W_abc.noalias() += wh.weights[2] * (P2 * M_B.cwiseProduct(M_C)).cwiseProduct(M_A);
                        }

                        Vector prod = W_abc;
                        for (int32_t ch : nd.children) {
                            if (ch != v_id && ch != id_A && ch != id_B && ch != id_C) {
                                prod.array() *= masked_upward_msg(ch, realiz.loc_in[ch]).array();
                                sc += masked_upward_sc(ch, log_sc[ch]);
                            }
                        }
                        Scalar mx = prod.maxCoeff();
                        if (mx > 0.0) {
                            prod /= mx;
                            sc += std::log(mx);
                        }
                        realiz.loc_in[nid] = prod;
                        log_sc[nid] = sc;
                    } else {
                        Vector prod = Vector::Ones(S);
                        Scalar sc = 0.0;
                        for (int32_t ch : nd.children) {
                            prod.array() *= masked_upward_msg(ch, realiz.loc_in[ch]).array();
                            sc += masked_upward_sc(ch, log_sc[ch]);
                        }
                        Scalar mx = prod.maxCoeff();
                        if (mx > 0.0) {
                            prod /= mx;
                            sc += std::log(mx);
                        }
                        realiz.loc_in[nid] = prod;
                        log_sc[nid] = sc;
                    }
                }
                Scalar root_dot = std::max<Scalar>(1e-300, root_pi.dot(realiz.loc_in[active_tree.root_id]));
                realiz.log_l = log_sc[active_tree.root_id] + std::log(root_dot);
                if (realiz.log_l > max_realiz_ll) max_realiz_ll = realiz.log_l;

                std::vector<Vector> node_post(act_nodes, Vector::Zero(S));
                std::vector<bool> outside_done(act_nodes, false);
                node_post[active_tree.root_id] = root_pi;
                for (auto it = active_tree.post_order.rbegin(); it != active_tree.post_order.rend(); ++it) {
                    int32_t uid = *it;
                    const auto& u_nd = active_tree.nodes[uid];
                    if (u_nd.children.empty() || outside_done[uid]) continue;

                    if (r_idx < 0 && uid >= 0 && static_cast<size_t>(uid) < wt_wheatstone_hubs.size() &&
                        wt_wheatstone_hubs[uid].hub_u == uid) {
                        const auto& wh = wt_wheatstone_hubs[uid];
                        int32_t id_A = wh.children[0];
                        int32_t id_B = wh.children[1];
                        int32_t id_C = wh.children[2];
                        int32_t v_id = wh.bridge_v;

                        Vector M_A = masked_upward_msg(id_A, realiz.loc_in[id_A]);
                        Vector M_B = masked_upward_msg(id_B, realiz.loc_in[id_B]);
                        Vector M_C = masked_upward_msg(id_C, realiz.loc_in[id_C]);

                        Matrix P0 = branch_grid_transition_matrix(v_id, g_near, alpha_c * wh.branch_lengths[0] * conv);
                        Matrix P1 = (wh.weights[1] > 0.0)
                                  ? branch_grid_transition_matrix(v_id, g_near, alpha_c * wh.branch_lengths[1] * conv)
                                  : Matrix::Identity(S, S);
                        Matrix P2 = (wh.weights[2] > 0.0)
                                  ? branch_grid_transition_matrix(v_id, g_near, alpha_c * wh.branch_lengths[2] * conv)
                                  : Matrix::Identity(S, S);

                        Vector PB0 = P0 * M_A.cwiseProduct(M_B);
                        Vector PB1 = (wh.weights[1] > 0.0) ? Vector(P1 * M_A.cwiseProduct(M_C)) : Vector(Vector::Zero(S));
                        Vector PB2 = (wh.weights[2] > 0.0) ? Vector(P2 * M_B.cwiseProduct(M_C)) : Vector(Vector::Zero(S));
                        Vector W_abc = wh.weights[0] * PB0.cwiseProduct(M_C)
                                     + wh.weights[1] * PB1.cwiseProduct(M_B)
                                     + wh.weights[2] * PB2.cwiseProduct(M_A);

                        std::vector<int32_t> extra_ch;
                        for (int32_t ch : u_nd.children) {
                            if (ch != v_id && ch != id_A && ch != id_B && ch != id_C) {
                                extra_ch.push_back(ch);
                            }
                        }
                        Vector M_rest = Vector::Ones(S);
                        for (int32_t d : extra_ch) {
                            M_rest.array() *= masked_upward_msg(d, realiz.loc_in[d]).array();
                        }
                        for (size_t idx = 0; idx < extra_ch.size(); ++idx) {
                            int32_t d = extra_ch[idx];
                            Vector cav_d = node_post[uid].cwiseProduct(W_abc);
                            for (size_t j = 0; j < extra_ch.size(); ++j) {
                                if (j != idx) {
                                    int32_t dj = extra_ch[j];
                                    cav_d.array() *= masked_upward_msg(dj, realiz.loc_in[dj]).array();
                                }
                            }
                            Scalar s_sum = cav_d.sum();
                            if (s_sum > 0.0) cav_d /= s_sum;
                            else cav_d = root_pi;
                            realiz.loc_out[d] = cav_d;

                            Vector post_d = apply_site_PT(d, cav_d);
                            Scalar pd_sum = post_d.sum();
                            if (pd_sum > 0.0) post_d /= pd_sum;
                            else post_d = root_pi;
                            node_post[d] = post_d;
                        }

                        Vector U_hub = node_post[uid].cwiseProduct(M_rest);

                        if (v_id != uid && v_id >= 0 && static_cast<size_t>(v_id) < act_nodes) {
                            Vector cav_v = U_hub.cwiseProduct(M_C);
                            Scalar sv = cav_v.sum();
                            if (sv > 0.0) cav_v /= sv;
                            else cav_v = root_pi;
                            realiz.loc_out[v_id] = cav_v;

                            Vector post_v = P0.transpose() * cav_v;
                            Scalar pv_s = post_v.sum();
                            if (pv_s > 0.0) post_v /= pv_s;
                            else post_v = root_pi;
                            node_post[v_id] = post_v;
                            outside_done[v_id] = true;
                        }

                        Vector QT0 = P0.transpose() * U_hub.cwiseProduct(M_C);
                        Vector QT1 = (wh.weights[1] > 0.0) ? Vector(P1.transpose() * U_hub.cwiseProduct(M_B)) : Vector(Vector::Zero(S));
                        Vector QT2 = (wh.weights[2] > 0.0) ? Vector(P2.transpose() * U_hub.cwiseProduct(M_A)) : Vector(Vector::Zero(S));

                        Vector V_A = wh.weights[0] * QT0.cwiseProduct(M_B)
                                   + wh.weights[1] * QT1.cwiseProduct(M_C)
                                   + wh.weights[2] * U_hub.cwiseProduct(PB2);
                        Vector V_B = wh.weights[0] * QT0.cwiseProduct(M_A)
                                   + wh.weights[1] * U_hub.cwiseProduct(PB1)
                                   + wh.weights[2] * QT2.cwiseProduct(M_C);
                        Vector V_C = wh.weights[0] * U_hub.cwiseProduct(PB0)
                                   + wh.weights[1] * QT1.cwiseProduct(M_A)
                                   + wh.weights[2] * QT2.cwiseProduct(M_B);

                        std::array<int32_t, 3> q_ids = {id_A, id_B, id_C};
                        std::array<Vector, 3> q_Vs = {std::move(V_A), std::move(V_B), std::move(V_C)};
                        for (int q = 0; q < 3; ++q) {
                            int32_t qid = q_ids[q];
                            Vector cav_q = std::move(q_Vs[q]);
                            Scalar sq = cav_q.sum();
                            if (sq > 0.0) cav_q /= sq;
                            else cav_q = root_pi;
                            realiz.loc_out[qid] = cav_q;

                            Vector post_q = apply_site_PT(qid, cav_q);
                            Scalar pq_s = post_q.sum();
                            if (pq_s > 0.0) post_q /= pq_s;
                            else post_q = root_pi;
                            node_post[qid] = post_q;
                        }
                        continue;
                    }

                    std::vector<Vector> child_msgs(u_nd.children.size());
                    for (size_t k = 0; k < u_nd.children.size(); ++k) {
                        int32_t vid = u_nd.children[k];
                        child_msgs[k] = masked_upward_msg(vid, realiz.loc_in[vid]);
                    }
                    for (size_t k = 0; k < u_nd.children.size(); ++k) {
                        int32_t vid = u_nd.children[k];
                        Vector cav_u = node_post[uid];
                        for (size_t j = 0; j < u_nd.children.size(); ++j) {
                            if (j != k) cav_u.array() *= child_msgs[j].array();
                        }
                        Scalar s_sum = cav_u.sum();
                        if (s_sum > 0.0) cav_u /= s_sum;
                        else cav_u = root_pi;
                        realiz.loc_out[vid] = cav_u;

                        Vector post_v = apply_site_PT(vid, cav_u);
                        Scalar pv_sum = post_v.sum();
                        if (pv_sum > 0.0) post_v /= pv_sum;
                        else post_v = root_pi;
                        node_post[vid] = post_v;
                    }
                }
            }

            // Bayes-update the realization weights: w_m \propto w_m^{prior} * exp(log_l_m - max_ll)
            Scalar sum_post_w = 0.0;
            for (auto& realiz : realizations) {
                realiz.post_weight = realiz.prior_weight * std::exp(std::clamp<Scalar>(realiz.log_l - max_realiz_ll, -50.0, 0.0));
                sum_post_w += realiz.post_weight;
            }
            if (sum_post_w > 0.0) {
                for (auto& realiz : realizations) realiz.post_weight /= sum_post_w;
            } else {
                realizations[0].post_weight = 1.0;
            }

            std::vector<Vector> acc_inside(max_tree_nodes, Vector::Zero(S));
            std::vector<Vector> acc_outside(max_tree_nodes, Vector::Zero(S));
            Scalar acc_ll = 0.0;
            for (const auto& realiz : realizations) {
                if (realiz.post_weight <= 1e-9) continue;
                acc_ll += realiz.post_weight * realiz.log_l;
                for (size_t nid = 0; nid < max_tree_nodes; ++nid) {
                    if (nid >= active_tree.num_nodes() ||
                        (static_cast<int32_t>(nid) != active_tree.root_id &&
                         active_tree.nodes[nid].parent_id == INVALID_INDEX)) {
                        acc_inside[nid].noalias() += realiz.post_weight * (Vector::Ones(S) / static_cast<Scalar>(S));
                        acc_outside[nid].noalias() += realiz.post_weight * root_pi;
                    } else {
                        Scalar in_sum = realiz.loc_in[nid].sum();
                        Vector norm_in = (in_sum > 0.0) ? (realiz.loc_in[nid] / in_sum).eval() : realiz.loc_in[nid];
                        acc_inside[nid].noalias() += realiz.post_weight * norm_in;
                        acc_outside[nid].noalias() += realiz.post_weight * realiz.loc_out[nid];
                    }
                }
            }
            site_inside[c] = std::move(acc_inside);
            site_outside[c] = std::move(acc_outside);
            site_base_ll[c] = acc_ll;
        };

        for (size_t c = 0; c < C; ++c) {
            compute_site_cavity(c, nullptr);
        }

        // Phase 4.4: Erlang-3 Minimum-Duration Baum-Welch Markov Error-Sink HMM with Schwarz BIC Selection
        // and Pass 1 Masked Cavity Recalculation (zero hardcoded grid subsets, transition rates, or sibling hacks)
        std::vector<std::vector<std::array<Scalar, 16>>> W_tensor(
            C, std::vector<std::array<Scalar, 16>>(B)
        );
        std::vector<std::vector<Scalar>> W_sink(C, std::vector<Scalar>(B, 0.0));
        std::vector<std::vector<Scalar>> markov_sink_gamma(C, std::vector<Scalar>(B, 0.0));

        Vector inv_one_minus_pi = Vector::Zero(S);
        for (int i = 0; i < S; ++i) {
            inv_one_minus_pi(i) = 1.0 / std::max<Scalar>(1e-12, 1.0 - pi_codon(i));
        }
        // Off-diagonal random-jump kernel likelihood of site c on branch bid, from its current cavities
        auto site_sink_likelihood = [&](size_t c, int32_t bid) -> Scalar {
            const Vector& U = site_outside[c][bid];
            const Vector& L = site_inside[c][bid];
            Scalar pi_dot_L = pi_codon.dot(L);
            Scalar sink_val = 0.0;
            for (int i = 0; i < S; ++i) {
                Scalar p_sink_L_i = (pi_dot_L - pi_codon(i) * L(i)) * inv_one_minus_pi(i);
                sink_val += U(i) * std::max<Scalar>(0.0, p_sink_L_i);
            }
            return std::max<Scalar>(0.0, sink_val);
        };
        Matrix P_mh_jump = Matrix::Zero(S, S);
        std::vector<std::vector<uint8_t>> codon_dnuc(S, std::vector<uint8_t>(S, 0));
        for (int i = 0; i < S; ++i) {
            Scalar row_mh_sum = 0.0;
            for (int j = 0; j < S; ++j) {
                if (i == j) continue;
                uint8_t d = static_cast<uint8_t>(gcode.diff_matrix[static_cast<size_t>(i)][static_cast<size_t>(j)].num_diffs);
                codon_dnuc[i][j] = d;
                if (d >= 2) {
                    P_mh_jump(i, j) = pi_codon(j);
                    row_mh_sum += pi_codon(j);
                }
            }
            if (row_mh_sum > 1e-12) P_mh_jump.row(i) /= row_mh_sum;
        }

        const Scalar bic_branch_sink_threshold = 2.0 * std::log(static_cast<Scalar>(std::max<size_t>(2, B * C)));
        bool any_branch_sink_detected = false;

        auto build_branch_hadamard = [&](int32_t bid, std::array<Scalar, 3>& br_w) -> std::vector<std::vector<Matrix>> {
            Matrix U_mat(S, C), L_mat(S, C);
            for (size_t c = 0; c < C; ++c) {
                U_mat.col(static_cast<Eigen::Index>(c)) = site_outside[c][bid];
                L_mat.col(static_cast<Eigen::Index>(c)) = site_inside[c][bid];
            }
            br_w = {0.0, 1.0, 0.0};
            if (num_k_regimes == 3 && dfel_analyzer && bid >= 0 && bid < static_cast<int32_t>(dfel_analyzer->branch_regime_weights.size())) {
                const auto& wv = dfel_analyzer->branch_regime_weights[bid];
                if (wv.size() == 3) { br_w[0] = wv[0]; br_w[1] = wv[1]; br_w[2] = wv[2]; }
            }
            std::vector<std::vector<Matrix>> H_reg(num_k_regimes, std::vector<Matrix>(G));
            for (size_t k_reg = 0; k_reg < num_k_regimes; ++k_reg) {
                Scalar wk = (num_k_regimes == 1) ? 1.0 : br_w[k_reg];
                if (wk <= 1e-6) continue;
                for (size_t g = 0; g < G; ++g) {
                    Matrix RT_U = regime_spectral_grids[k_reg][g].R.transpose() * U_mat;
                    Matrix L_L  = regime_spectral_grids[k_reg][g].L * L_mat;
                    H_reg[k_reg][g] = RT_U.cwiseProduct(L_L);
                }
            }
            return H_reg;
        };

        auto eval_hadamard_at = [&](const std::vector<std::vector<Matrix>>& H_reg,
                                    const std::array<Scalar, 3>& br_w,
                                    size_t c, size_t g_idx, Scalar eff_t) -> Scalar {
            Scalar t_pos = std::max<Scalar>(0.0, eff_t);
            Scalar val = 0.0;
            for (size_t k_reg = 0; k_reg < num_k_regimes; ++k_reg) {
                Scalar wk = (num_k_regimes == 1) ? 1.0 : br_w[k_reg];
                if (wk <= 1e-6) continue;
                Vector exp_d = (regime_spectral_grids[k_reg][g_idx].evals.array() * t_pos).exp();
                val += wk * H_reg[k_reg][g_idx].col(static_cast<Eigen::Index>(c)).dot(exp_d);
            }
            return std::max<Scalar>(1e-30, val);
        };

        for (size_t bi = 0; bi < B; ++bi) {
            int32_t bid = active_branches[bi];
            Scalar raw_bl = work_tree.nodes[bid].branch_length * conv;
            std::array<Scalar, 3> br_w;
            auto H_reg = build_branch_hadamard(bid, br_w);

            // Compute off-diagonal random-jump kernel likelihood f_sink[c] for each site c on branch bi
            std::vector<Scalar> f_sink(C, 0.0);
            std::vector<bool> is_homologous_site(C, true);
            size_t n_homologous = 0;
            for (size_t c = 0; c < C; ++c) {
                f_sink[c] = site_sink_likelihood(c, bid);
                Scalar f_raw = eval_hadamard_at(H_reg, br_w, c, site_g_near[c], site_alpha[c] * raw_bl);
                if (f_raw >= f_sink[c]) {
                    is_homologous_site[c] = true;
                    n_homologous++;
                } else {
                    is_homologous_site[c] = false;
                }
            }

            // 1D golden-section cavity MLE of the uncontaminated branch length ref_bl in [1e-4, raw_bl]
            Scalar ref_bl = raw_bl;
            if (n_homologous > 0 && n_homologous < C && raw_bl > 1e-4) {
                auto eval_clean_bl_ll = [&](Scalar cand_bl) -> Scalar {
                    Scalar ll = 0.0;
                    for (size_t c = 0; c < C; ++c) {
                        if (!is_homologous_site[c]) continue;
                        ll += std::log(eval_hadamard_at(H_reg, br_w, c, site_g_near[c], site_alpha[c] * cand_bl));
                    }
                    return ll;
                };
                Scalar lo_b = 1e-4, hi_b = raw_bl;
                for (int it = 0; it < 12; ++it) {
                    Scalar m1 = lo_b + 0.381966 * (hi_b - lo_b);
                    Scalar m2 = lo_b + 0.618034 * (hi_b - lo_b);
                    if (eval_clean_bl_ll(m1) > eval_clean_bl_ll(m2)) hi_b = m2;
                    else lo_b = m1;
                }
                ref_bl = 0.5 * (lo_b + hi_b);
            }

            // Evaluate local emission odds W_sink[c][bi] against the full 16-point Empirical Bayes biological denominator
            for (size_t c = 0; c < C; ++c) {
                Scalar site_eff_bl = site_alpha[c] * ref_bl;
                Scalar ref_val = eval_hadamard_at(H_reg, br_w, c, site_g_near[c], site_eff_bl);
                for (size_t g = 0; g < G; ++g) {
                    Scalar bio_l = emp_grid_prior[g] * eval_hadamard_at(H_reg, br_w, c, g, site_eff_bl);
                    if (bio_l > ref_val) ref_val = bio_l;
                }
                W_sink[c][bi] = std::max<Scalar>(0.0, f_sink[c] / std::max<Scalar>(1e-30, ref_val));
            }

            // Erlang-3 Minimum-Duration (0 -> 1_1 -> 1_2 -> 1_3 self-loop rho_b) Baum-Welch EM + Schwarz BIC gate
            if (C >= 3) {
                auto lse = [](Scalar a, Scalar b) -> Scalar {
                    if (a <= -1e20) return b;
                    if (b <= -1e20) return a;
                    Scalar mx = std::max(a, b);
                    return mx + std::log(std::exp(a - mx) + std::exp(b - mx));
                };
                std::vector<Scalar> log_em1(C, -30.0);
                for (size_t c = 0; c < C; ++c) {
                    log_em1[c] = std::min<Scalar>(8.0, std::log(std::max<Scalar>(1e-15, W_sink[c][bi])));
                }

                Scalar eps_b = 1.0 / static_cast<Scalar>(std::max<size_t>(100, C));
                Scalar rho_b = 0.75;
                Scalar log_Z = 0.0;
                const Scalar NEG_INF = -1e25;
                std::vector<std::array<Scalar, 4>> f(C, {NEG_INF, NEG_INF, NEG_INF, NEG_INF});
                std::vector<std::array<Scalar, 4>> b(C, {NEG_INF, NEG_INF, NEG_INF, NEG_INF});

                for (int em_it = 0; em_it < 8; ++em_it) {
                    Scalar log_1m_eps = std::log(1.0 - eps_b);
                    Scalar log_eps = std::log(eps_b);
                    Scalar log_1m_rho = std::log(1.0 - rho_b);
                    Scalar log_rho = std::log(rho_b);

                    f[0] = {log_1m_eps, log_eps + log_em1[0], NEG_INF, NEG_INF};
                    for (size_t c = 1; c < C; ++c) {
                        f[c][0] = lse(f[c - 1][0] + log_1m_eps, f[c - 1][3] + log_1m_rho);
                        f[c][1] = f[c - 1][0] + log_eps + log_em1[c];
                        f[c][2] = f[c - 1][1] + log_em1[c];
                        f[c][3] = lse(f[c - 1][2], f[c - 1][3] + log_rho) + log_em1[c];
                    }
                    log_Z = lse(f[C - 1][0], f[C - 1][3]);

                    b[C - 1] = {0.0, NEG_INF, NEG_INF, 0.0};
                    for (int c = static_cast<int>(C) - 2; c >= 0; --c) {
                        b[c][0] = lse(log_1m_eps + b[c + 1][0], log_eps + log_em1[c + 1] + b[c + 1][1]);
                        b[c][1] = log_em1[c + 1] + b[c + 1][2];
                        b[c][2] = log_em1[c + 1] + b[c + 1][3];
                        b[c][3] = lse(log_1m_rho + b[c + 1][0], log_rho + log_em1[c + 1] + b[c + 1][3]);
                    }

                    Scalar xi_01 = 0.0, sum_g0 = 0.0, xi_33 = 0.0, sum_g3 = 0.0;
                    for (size_t c = 0; c + 1 < C; ++c) {
                        sum_g0 += std::exp(std::clamp<Scalar>(f[c][0] + b[c][0] - log_Z, -60.0, 0.0));
                        xi_01  += std::exp(std::clamp<Scalar>(f[c][0] + log_eps + log_em1[c + 1] + b[c + 1][1] - log_Z, -60.0, 0.0));
                        sum_g3 += std::exp(std::clamp<Scalar>(f[c][3] + b[c][3] - log_Z, -60.0, 0.0));
                        xi_33  += std::exp(std::clamp<Scalar>(f[c][3] + log_rho + log_em1[c + 1] + b[c + 1][3] - log_Z, -60.0, 0.0));
                    }
                    Scalar min_eps = 1.0 / static_cast<Scalar>(std::max<size_t>(100, B * C));
                    eps_b = std::clamp<Scalar>(xi_01 / std::max<Scalar>(1.0, sum_g0), min_eps, 0.05);
                    if (sum_g3 > 1e-4) {
                        rho_b = std::clamp<Scalar>(xi_33 / sum_g3, 0.25, 1.0 - 1.0 / static_cast<Scalar>(C));
                    }
                }

                // Schwarz BIC gate: 2 * ln Z_b(hat{eps}_b, hat{rho}_b) > 2 * ln(B * C)
                if (2.0 * log_Z > bic_branch_sink_threshold) {
                    Scalar sum_gamma = 0.0;
                    for (size_t c = 0; c < C; ++c) {
                        Scalar g1 = std::exp(std::clamp<Scalar>(f[c][1] + b[c][1] - log_Z, -60.0, 0.0))
                                  + std::exp(std::clamp<Scalar>(f[c][2] + b[c][2] - log_Z, -60.0, 0.0))
                                  + std::exp(std::clamp<Scalar>(f[c][3] + b[c][3] - log_Z, -60.0, 0.0));
                        markov_sink_gamma[c][bi] = std::clamp<Scalar>(g1, 0.0, 1.0);
                        sum_gamma += markov_sink_gamma[c][bi];
                    }
                    if (sum_gamma >= 0.50) {
                        any_branch_sink_detected = true;
                        // Update branch bi's length via 1D golden-section cavity MLE on the (1 - gamma_{c,b}) weighted sites
                        auto eval_unmasked_bl_ll = [&](Scalar cand_bl) -> Scalar {
                            Scalar ll = 0.0;
                            for (size_t c = 0; c < C; ++c) {
                                Scalar w_clean = 1.0 - markov_sink_gamma[c][bi];
                                if (w_clean <= 1e-4) continue;
                                ll += w_clean * std::log(eval_hadamard_at(H_reg, br_w, c, site_g_near[c], site_alpha[c] * cand_bl));
                            }
                            return ll;
                        };
                        Scalar lo_b = 1e-4, hi_b = raw_bl;
                        for (int it = 0; it < 12; ++it) {
                            Scalar m1 = lo_b + 0.381966 * (hi_b - lo_b);
                            Scalar m2 = lo_b + 0.618034 * (hi_b - lo_b);
                            if (eval_unmasked_bl_ll(m1) > eval_unmasked_bl_ll(m2)) hi_b = m2;
                            else lo_b = m1;
                        }
                        work_tree.nodes[bid].branch_length = (0.5 * (lo_b + hi_b)) / conv;
                    }
                }
            }
        }

        // Pass 1 Cavity Recalculation: If any branch had an error-sink tract passing the Schwarz BIC gate,
        // recompute the Inside-Outside cavities for affected sites with the corrupted edge messages masked to 1.
        // This directly removes error leakage into sibling, parent, and cousin Outside vectors without any sibling quarantine hack.
        if (any_branch_sink_detected) {
            for (size_t c = 0; c < C; ++c) {
                bool site_has_mask = false;
                std::vector<Scalar> node_mask(max_tree_nodes, 0.0);
                for (size_t bi = 0; bi < B; ++bi) {
                    if (markov_sink_gamma[c][bi] > 1e-3) {
                        int32_t bid = active_branches[bi];
                        if (bid >= 0 && static_cast<size_t>(bid) < max_tree_nodes) {
                            node_mask[bid] = markov_sink_gamma[c][bi];
                            site_has_mask = true;
                        }
                    }
                }
                if (site_has_mask) {
                    compute_site_cavity(c, &node_mask);
                }
            }
        }

        // Populate the final sufficiency tensor W_tensor[c][bi][g] from the clean cavities
        for (size_t bi = 0; bi < B; ++bi) {
            int32_t bid = active_branches[bi];
            std::array<Scalar, 3> br_w;
            auto H_reg = build_branch_hadamard(bid, br_w);
            Scalar bl = work_tree.nodes[bid].branch_length * conv;

            for (size_t c = 0; c < C; ++c) {
                Scalar c_bl = bl;
                int r_idx = site_routed_tree_idx[c];
                if (r_idx >= 0) {
                    const Tree& spr_t = routed_spr_trees[r_idx].spr_tree;
                    if (bid >= static_cast<int32_t>(spr_t.num_nodes()) ||
                        (bid != spr_t.root_id && spr_t.nodes[bid].parent_id == INVALID_INDEX)) {
                        for (size_t g = 0; g < G; ++g) W_tensor[c][bi][g] = 1.0;
                        W_sink[c][bi] = 1.0;
                        continue;
                    }
                    c_bl = spr_t.nodes[bid].branch_length * conv;
                }
                Scalar eff_site_bl = site_alpha[c] * c_bl;
                for (size_t g = 0; g < G; ++g) {
                    W_tensor[c][bi][g] = eval_hadamard_at(H_reg, br_w, c, g, eff_site_bl);
                }
                Scalar ref = std::max<Scalar>(1e-30, W_tensor[c][bi][site_g_near[c]]);
                for (size_t g = 0; g < G; ++g) {
                    W_tensor[c][bi][g] /= ref;
                }
                Scalar g_mkt = markov_sink_gamma[c][bi];
                if (g_mkt >= 0.50) {
                    for (size_t g = 0; g < G; ++g) {
                        W_tensor[c][bi][g] = (1.0 - g_mkt) * W_tensor[c][bi][g] + g_mkt * 1.0;
                    }
                    W_sink[c][bi] = 1.0;
                } else {
                    // The BUSTED-E mixture adds W_sink to W_tensor, so both need the same cavities and ref
                    W_sink[c][bi] = site_sink_likelihood(c, bid) / ref;
                }
            }
        }

        // Helper for continuous ω evaluation on site c, branch bi via log-linear interpolation of W_tensor[c][bi]
        auto eval_branch_site_omega = [&](size_t c, size_t bi, Scalar omega_val) -> Scalar {
            Scalar w = std::clamp<Scalar>(omega_val, omega_grid.front(), omega_grid.back());
            if (w <= omega_grid[0]) return W_tensor[c][bi][0];
            if (w >= omega_grid[G - 1]) return W_tensor[c][bi][G - 1];
            size_t hi = 1;
            while (hi < G && omega_grid[hi] < w) ++hi;
            size_t lo = hi - 1;
            Scalar span = std::max<Scalar>(1e-12, omega_grid[hi] - omega_grid[lo]);
            Scalar frac = (w - omega_grid[lo]) / span;
            Scalar v_lo = W_tensor[c][bi][lo];
            Scalar v_hi = W_tensor[c][bi][hi];
            return (1.0 - frac) * v_lo + frac * v_hi;
        };

        // Pre-diagonalize Q_syn and Q_non or use a fast spectral decomposition helper for arbitrary (α, β)
        // Notice: Q(α, β) = α (Q_syn + (β/α) Q_non). For any ratio ω = β/α, we can either diagonalize
        // symmetric(Q(1, ω)) once per Nelder-Mead step (only ONE 61x61 eigendecomposition for all B branches!)
        // instead of computing B separate matrix exponentials per evaluation!
        Vector sqrt_pi_base = pi_codon.array().max(1e-15).sqrt();
        Vector inv_sqrt_pi_base = 1.0 / sqrt_pi_base.array();

        auto spectral_decompose_regime_omega = [&](size_t k_reg, Scalar omega_val) -> SpectralDecomp {
            Matrix Q_w = m_syn_reg[k_reg].Q + omega_val * m_non_reg[k_reg].Q;
            for (int i = 0; i < S; ++i) {
                Q_w(i, i) = 0.0;
                Q_w(i, i) = -Q_w.row(i).sum();
            }
            Matrix sym = sqrt_pi_reg[k_reg].asDiagonal() * Q_w * inv_sqrt_pi_reg[k_reg].asDiagonal();
            sym = 0.5 * (sym + sym.transpose());
            Eigen::SelfAdjointEigenSolver<Matrix> es(sym);
            SpectralDecomp sd;
            sd.evals = es.eigenvalues();
            sd.R = inv_sqrt_pi_reg[k_reg].asDiagonal() * es.eigenvectors();
            sd.L = es.eigenvectors().transpose() * sqrt_pi_reg[k_reg].asDiagonal();
            return sd;
        };

        // ---------------------------------------------------------------------
        // 3A. MEME Episodic Site Scan + Decision 1 Active-Set Continuous Polish
        // ---------------------------------------------------------------------
        res.site_atlas.resize(C);
        std::vector<size_t> branch_epi_counts(B, 0);

        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic)
        #endif
        for (size_t c = 0; c < C; ++c) {
            auto& entry = res.site_atlas[c];
            entry.site_index = c;
            entry.column_homology_posterior = (c < res.joint_align.site_homology_posteriors.size())
                                            ? res.joint_align.site_homology_posteriors[c] : 1.0;
            entry.ensemble_paths_evaluated = site_ensemble_k[c];
            if (site_routed_tree_idx[c] >= 0) {
                entry.active_topology_label = routed_spr_trees[site_routed_tree_idx[c]].label;
            }

            if (c < res.joint_align.site_results.size()) {
                const auto& sr = res.joint_align.site_results[c];
                entry.gap_fraction = sr.gap_fraction;
                entry.in_indel_window = sr.in_indel_window;
                entry.alpha = sr.marginalized_alpha;
                entry.beta_pervasive = sr.marginalized_beta;
                entry.fel_lrt = sr.marginalized_lrt;
                entry.fel_p_value = sr.marginalized_p_value;
                entry.rigid_alignment_artifact_eliminated =
                    (sr.rigid_p_value <= settings.p_value_threshold && sr.marginalized_p_value > settings.p_value_threshold);
            }
            if (dfel_analyzer && c < dfel_analyzer->site_results.size()) {
                entry.dfel_lrt = dfel_analyzer->site_results[c].lrt;
                entry.dfel_p_value = dfel_analyzer->site_results[c].p_value;
                if (site_routed_tree_idx[c] < 0) {
                    entry.alpha = dfel_analyzer->site_results[c].alpha;
                    entry.beta_pervasive = dfel_analyzer->site_results[c].beta;
                }
            } else if (res.directional_drift_active && dfel_analyzer) {
                Scalar df_null_ll = 0.0;
                for (size_t bi = 0; bi < B; ++bi) df_null_ll += std::log(std::max<Scalar>(1e-30, W_tensor[c][bi][8]));
                Scalar df_alt_ll = df_null_ll;
                size_t best_df_g = 8;
                for (size_t g = 0; g < G; ++g) {
                    Scalar ll = 0.0;
                    for (size_t bi = 0; bi < B; ++bi) ll += std::log(std::max<Scalar>(1e-30, W_tensor[c][bi][g]));
                    if (ll > df_alt_ll) { df_alt_ll = ll; best_df_g = g; }
                }
                Scalar w_lo = (best_df_g > 0) ? omega_grid[best_df_g - 1] : omega_grid[0];
                Scalar w_hi = (best_df_g + 1 < G) ? omega_grid[best_df_g + 1] : omega_grid[G - 1];
                Scalar best_df_w = omega_grid[best_df_g];
                for (int it = 0; it < 10; ++it) {
                    Scalar m1 = w_lo + 0.381966 * (w_hi - w_lo);
                    Scalar m2 = w_lo + 0.618034 * (w_hi - w_lo);
                    Scalar l1 = 0.0, l2 = 0.0;
                    for (size_t bi = 0; bi < B; ++bi) {
                        l1 += std::log(std::max<Scalar>(1e-30, eval_branch_site_omega(c, bi, m1)));
                        l2 += std::log(std::max<Scalar>(1e-30, eval_branch_site_omega(c, bi, m2)));
                    }
                    if (l1 > l2) { w_hi = m2; if (l1 > df_alt_ll) { df_alt_ll = l1; best_df_w = m1; } }
                    else { w_lo = m1; if (l2 > df_alt_ll) { df_alt_ll = l2; best_df_w = m2; } }
                }
                entry.dfel_lrt = std::max<Scalar>(0.0, 2.0 * (df_alt_ll - df_null_ll)) * entry.column_homology_posterior;
                entry.dfel_p_value = std::erfc(std::sqrt(entry.dfel_lrt / 2.0));
                if (site_routed_tree_idx[c] < 0) {
                    entry.beta_pervasive = best_df_w * std::max<Scalar>(0.05, entry.alpha);
                }
            } else {
                entry.dfel_lrt = entry.fel_lrt;
                entry.dfel_p_value = entry.fel_p_value;
            }

            // Count non-gap taxa and distinct amino acids in the MAP column (E6, E11):
            // 1. Columns with < 3 non-gap sequences lack sufficient evolutionary branch support for site testing.
            // 2. Columns with <= 1 distinct amino acid across all non-gap sequences in the MAP alignment have zero
            //    observed amino-acid variation and cannot support positive diversifying selection (beta > alpha or beta+ > alpha).
            size_t non_gap_taxa = 0;
            uint32_t aa_bitmask = 0;
            size_t distinct_map_aa = 0;
            for (size_t t = 0; t < N; ++t) {
                int8_t cod = aln.sequences[t].codons[c];
                if (cod >= 0 && cod < S) {
                    non_gap_taxa++;
                    int a_idx = aa_char_to_index(gcode.sense_codon_aa[cod]);
                    if (a_idx >= 0 && a_idx < 20) {
                        uint32_t bit = (1u << static_cast<uint32_t>(a_idx));
                        if ((aa_bitmask & bit) == 0) {
                            aa_bitmask |= bit;
                            distinct_map_aa++;
                        }
                    }
                }
            }
            if (non_gap_taxa < 3) {
                entry.alpha = 0.0;
                entry.beta_pervasive = 0.0;
                entry.fel_lrt = 0.0;
                entry.fel_p_value = 1.0;
                entry.dfel_lrt = 0.0;
                entry.dfel_p_value = 1.0;
                entry.beta_minus = 0.0;
                entry.beta_plus = 0.0;
                entry.weight_plus = 0.0;
                entry.meme_lrt = 0.0;
                entry.meme_p_value = 1.0;
                entry.max_branch_ebf = 0.0;
                entry.num_episodic_branches = 0;
                entry.top_episodic_branches.clear();
                continue;
            }
            if (distinct_map_aa <= 1) {
                bool was_false_pos = (entry.beta_pervasive > entry.alpha);
                entry.beta_pervasive = 0.0;
                if (was_false_pos) {
                    entry.fel_lrt = 0.0;
                    entry.fel_p_value = 1.0;
                    entry.dfel_lrt = 0.0;
                    entry.dfel_p_value = 1.0;
                }
                entry.beta_minus = 0.0;
                entry.beta_plus = 0.0;
                entry.weight_plus = 0.0;
                entry.meme_lrt = 0.0;
                entry.meme_p_value = 1.0;
                entry.max_branch_ebf = 0.0;
                entry.num_episodic_branches = 0;
                entry.top_episodic_branches.clear();
                continue;
            }

            // Stage 1 coarse grid screen on W_tensor[c]
            Scalar best_null_ll = -1e18;
            size_t best_g_null = 6;
            for (size_t g0 = 0; g0 <= 8; ++g0) {
                Scalar ll = 0.0;
                for (size_t bi = 0; bi < B; ++bi) ll += std::log(W_tensor[c][bi][g0]);
                if (ll > best_null_ll) {
                    best_null_ll = ll;
                    best_g_null = g0;
                }
            }

            Scalar best_alt_ll = best_null_ll;
            size_t best_g_minus = best_g_null;
            size_t best_g_plus = 8;
            Scalar best_w_plus = 0.0;

            const std::array<Scalar, 7> w_grid = {0.02, 0.05, 0.10, 0.20, 0.35, 0.55, 0.80};
            for (size_t gm = 0; gm <= 8; ++gm) {
                for (size_t gp = 9; gp < G; ++gp) {
                    for (Scalar wp : w_grid) {
                        Scalar ll = 0.0;
                        for (size_t bi = 0; bi < B; ++bi) {
                            Scalar mix = (1.0 - wp) * W_tensor[c][bi][gm] + wp * W_tensor[c][bi][gp];
                            ll += std::log(std::max<Scalar>(1e-30, mix));
                        }
                        if (ll > best_alt_ll) {
                            best_alt_ll = ll;
                            best_g_minus = gm;
                            best_g_plus = gp;
                            best_w_plus = wp;
                        }
                    }
                }
            }

            Scalar w_minus_opt = omega_grid[best_g_minus];
            Scalar w_plus_opt = omega_grid[best_g_plus];
            Scalar wp_opt = best_w_plus;

            Scalar screen_lrt = std::max<Scalar>(0.0, 2.0 * (best_alt_ll - best_null_ll));
            Scalar screen_p = (screen_lrt <= 1e-6) ? 1.0
                : (0.30 * stats::chi2_p_value(screen_lrt, 1.0) + 0.37 * stats::chi2_p_value(screen_lrt, 2.0));

            Scalar eff_alpha = std::max<Scalar>(0.10, entry.alpha);
            Scalar raw_meme_lrt = screen_lrt;

            // Decision 1 (Upgraded + Spectral-Accelerated): Exact 4D/3D Nelder-Mead Felsenstein Pruning Polish
            // Diagonalizes Q^{(k)}(1, w1) and Q^{(k)}(1, w+) ONCE per Nelder-Mead evaluation and computes all B branch
            // matrices via Softmax-weighted R^{(k)} * exp(alpha * t_b * evals^{(k)}) * L^{(k)}!
            if (screen_p < settings.active_set_screen_p || entry.beta_pervasive > entry.alpha) {
                entry.meme_polished = true;
#ifdef __EMSCRIPTEN__
                // In single-threaded WebAssembly, polish (w_minus, w_plus, wp) directly on the continuous
                // Inside-Outside sufficiency surface eval_branch_site_omega(c, bi, w)
                auto eval_io_meme_ll = [&](Scalar wm, Scalar wp_val, Scalar wpl) -> Scalar {
                    Scalar ll = 0.0;
                    for (size_t bi = 0; bi < B; ++bi) {
                        Scalar m1 = eval_branch_site_omega(c, bi, wm);
                        Scalar mp = eval_branch_site_omega(c, bi, wpl);
                        ll += std::log(std::max<Scalar>(1e-30, (1.0 - wp_val) * m1 + wp_val * mp));
                    }
                    return ll;
                };
                std::array<Scalar, 3> init_io = {w_minus_opt, wp_opt, std::max<Scalar>(1.05, w_plus_opt)};
                std::array<Scalar, 3> lb_io = {0.0, 0.0, 1.0};
                std::array<Scalar, 3> ub_io = {1.0, 0.99, 25.0};
                auto opt_io = opt::NelderMeadND<3>::minimize(
                    [&](const std::array<Scalar, 3>& x) { return -eval_io_meme_ll(x[0], x[1], x[2]); },
                    init_io, lb_io, ub_io, 1e-4, 60
                );
                Scalar polished_alt_ll = std::max(best_alt_ll, -opt_io.f);
                if (-opt_io.f >= best_alt_ll) {
                    w_minus_opt = opt_io.x[0];
                    wp_opt = opt_io.x[1];
                    w_plus_opt = opt_io.x[2];
                }
                Scalar lo_n = 0.0, hi_n = 1.0;
                Scalar polished_null_ll = best_null_ll;
                for (int it = 0; it < 12; ++it) {
                    Scalar m1 = lo_n + 0.381966 * (hi_n - lo_n);
                    Scalar m2 = lo_n + 0.618034 * (hi_n - lo_n);
                    Scalar l1 = eval_io_meme_ll(m1, 0.0, m1);
                    Scalar l2 = eval_io_meme_ll(m2, 0.0, m2);
                    if (l1 > l2) { hi_n = m2; polished_null_ll = std::max(polished_null_ll, l1); }
                    else { lo_n = m1; polished_null_ll = std::max(polished_null_ll, l2); }
                }
                raw_meme_lrt = (w_plus_opt > 1.0 && wp_opt > 1e-5)
                    ? std::max<Scalar>(0.0, 2.0 * (polished_alt_ll - polished_null_ll))
                    : 0.0;
                entry.beta_minus = w_minus_opt * eff_alpha;
                entry.beta_plus = w_plus_opt * eff_alpha;
                entry.weight_plus = wp_opt;
#else
                int r_idx = site_routed_tree_idx[c];
                const Tree& active_tree = (r_idx >= 0) ? routed_spr_trees[r_idx].spr_tree : work_tree;
                const auto& active_l2t = (r_idx >= 0) ? routed_spr_trees[r_idx].leaf_to_taxon : leaf_to_taxon;

                SitePattern pat;
                pat.states.resize(N);
                std::vector<bool> taxon_masked(N, false);
                for (size_t bi = 0; bi < B; ++bi) {
                    if (markov_sink_gamma[c][bi] >= 0.50) {
                        int32_t bid = active_branches[bi];
                        if (bid >= 0 && static_cast<size_t>(bid) < active_l2t.size()) {
                            size_t tid = active_l2t[bid];
                            if (tid < N) taxon_masked[tid] = true;
                        }
                    }
                }
                bool is_invar = true;
                int8_t first_st = -1;
                for (size_t t = 0; t < N; ++t) {
                    int8_t st = taxon_masked[t] ? static_cast<int8_t>(-1) : aln.sequences[t].codons[c];
                    pat.states[t] = st;
                    if (st >= 0) {
                        if (first_st == -1) first_st = st;
                        else if (st != first_st) is_invar = false;
                    }
                }
                pat.weight = 1;

                if (is_invar) {
                    raw_meme_lrt = 0.0;
                    wp_opt = 0.0;
                    entry.beta_minus = w_minus_opt * eff_alpha;
                    entry.beta_plus = w_minus_opt * eff_alpha;
                    entry.weight_plus = 0.0;
                } else {
                    std::vector<Vector> t_node_L(active_tree.num_nodes(), Vector::Zero(S));
                    std::vector<double> t_node_sc(active_tree.num_nodes(), 0.0);
                    std::vector<Matrix> P_br(active_tree.num_nodes());
                    std::vector<SpectralDecomp> sd1_reg(num_k_regimes), sdp_reg(num_k_regimes);
                    Scalar cached_w1 = -1.0, cached_wp = -1.0;

                    auto eval_exact_meme_ll = [&](Scalar a_val, Scalar w1_val, Scalar p1_val, Scalar b_plus_val) -> Scalar {
                        Scalar a_pos = std::max<Scalar>(1e-5, a_val);
                        Scalar w_plus_val = b_plus_val / a_pos;
                        if (std::abs(w1_val - cached_w1) > 1e-12) {
                            for (size_t k_reg = 0; k_reg < num_k_regimes; ++k_reg) {
                                sd1_reg[k_reg] = spectral_decompose_regime_omega(k_reg, w1_val);
                            }
                            cached_w1 = w1_val;
                        }
                        if (std::abs(w_plus_val - cached_wp) > 1e-12) {
                            for (size_t k_reg = 0; k_reg < num_k_regimes; ++k_reg) {
                                sdp_reg[k_reg] = spectral_decompose_regime_omega(k_reg, w_plus_val);
                            }
                            cached_wp = w_plus_val;
                        }
                        Scalar p_plus_wt = 1.0 - p1_val;
                        for (const auto& nd : active_tree.nodes) {
                            if (nd.id != active_tree.root_id && nd.parent_id != INVALID_INDEX) {
                                Scalar eff_t = a_pos * nd.branch_length * conv;
                                if (num_k_regimes == 1 || !dfel_analyzer || nd.id < 0 ||
                                    nd.id >= static_cast<int32_t>(dfel_analyzer->branch_regime_weights.size())) {
                                    Vector e1 = (sd1_reg[0].evals.array() * eff_t).exp();
                                    Vector ep = (sdp_reg[0].evals.array() * eff_t).exp();
                                    Matrix m1 = sd1_reg[0].R * e1.asDiagonal() * sd1_reg[0].L;
                                    Matrix mp = sdp_reg[0].R * ep.asDiagonal() * sdp_reg[0].L;
                                    P_br[nd.id] = (p1_val * m1 + p_plus_wt * mp).cwiseMax(1e-15);
                                } else {
                                    const auto& w_vec = dfel_analyzer->branch_regime_weights[nd.id];
                                    Matrix P_mix = Matrix::Zero(S, S);
                                    for (size_t k_reg = 0; k_reg < 3; ++k_reg) {
                                        Scalar wk = (k_reg < w_vec.size()) ? w_vec[k_reg] : (k_reg == 1 ? 1.0 : 0.0);
                                        if (wk <= 1e-6) continue;
                                        Vector e1 = (sd1_reg[k_reg].evals.array() * eff_t).exp();
                                        Vector ep = (sdp_reg[k_reg].evals.array() * eff_t).exp();
                                        Matrix m1 = sd1_reg[k_reg].R * e1.asDiagonal() * sd1_reg[k_reg].L;
                                        Matrix mp = sdp_reg[k_reg].R * ep.asDiagonal() * sdp_reg[k_reg].L;
                                        P_mix.noalias() += wk * (p1_val * m1 + p_plus_wt * mp);
                                    }
                                    P_br[nd.id] = P_mix.cwiseMax(1e-15);
                                }
                            }
                        }
                        const MG94Matrix& root_mg = (num_k_regimes == 3 && dfel_analyzer)
                            ? dfel_analyzer->regimes[1].base_matrix : mg94_base;
                        return FELAnalyzer::compute_pattern_log_likelihood(
                            active_tree, pat, active_l2t, root_mg, P_br, t_node_L, t_node_sc
                        );
                    };

                    Scalar a_init = std::clamp<Scalar>(entry.alpha, 1e-3, 50.0);
                    Scalar b_init = std::clamp<Scalar>(entry.beta_pervasive, 1e-3, 200.0);
                    Scalar w1_init = std::clamp<Scalar>(b_init / std::max<Scalar>(a_init, 1e-3), 0.0, 1.0);
                    Scalar p1_init = std::clamp<Scalar>(1.0 - wp_opt, 0.05, 0.98);
                    Scalar bp_init = std::max<Scalar>(b_init * 1.5, w_plus_opt * a_init);

                    std::vector<std::array<Scalar, 4>> start_grid = {
                        {a_init, w1_init, p1_init, bp_init},
                        {a_init, w1_init, 0.75, b_init},
                        {a_init, 0.0, 0.85, std::max<Scalar>(b_init * 3.0, a_init * 2.0)},
                        {a_init, 0.25, 0.75, std::max<Scalar>(b_init * 5.0, a_init * 3.0)}
                    };
                    std::array<Scalar, 4> lb_alt = {1e-4, 0.0, 1e-5, 1e-4};
                    std::array<Scalar, 4> ub_alt = {100.0, 1.0, 1.0 - 1e-5, 1000.0};
                    auto obj_alt = [&](const std::array<Scalar, 4>& x) -> Scalar {
                        return -eval_exact_meme_ll(x[0], x[1], x[2], x[3]);
                    };
                    auto opt_alt = opt::NelderMeadND<4>::minimize(obj_alt, start_grid[0], lb_alt, ub_alt, 1e-4, 150, start_grid);

                    eff_alpha = opt_alt.x[0];
                    w_minus_opt = opt_alt.x[1];
                    Scalar p1_opt = opt_alt.x[2];
                    Scalar bp_abs = opt_alt.x[3];
                    wp_opt = 1.0 - p1_opt;
                    w_plus_opt = bp_abs / std::max<Scalar>(1e-4, eff_alpha);
                    Scalar exact_alt_ll = -opt_alt.f;

                    if (bp_abs > eff_alpha && wp_opt > 1e-5) {
                        std::vector<std::array<Scalar, 3>> null_grid = {
                            {eff_alpha, w_minus_opt, p1_opt},
                            {a_init, 1.0, 0.5},
                            {eff_alpha, 1.0, 0.5}
                        };
                        std::array<Scalar, 3> lb_null = {1e-4, 0.0, 1e-5};
                        std::array<Scalar, 3> ub_null = {100.0, 1.0, 1.0 - 1e-5};
                        auto obj_null = [&](const std::array<Scalar, 3>& x) -> Scalar {
                            return -eval_exact_meme_ll(x[0], x[1], x[2], x[0]);
                        };
                        auto opt_null = opt::NelderMeadND<3>::minimize(obj_null, null_grid[0], lb_null, ub_null, 1e-4, 120, null_grid);
                        Scalar exact_null_ll = -opt_null.f;
                        raw_meme_lrt = std::max<Scalar>(0.0, 2.0 * (exact_alt_ll - exact_null_ll));
                    } else {
                        raw_meme_lrt = 0.0;
                    }
                    entry.alpha = eff_alpha;
                    entry.beta_minus = eff_alpha * w_minus_opt;
                    entry.beta_plus = bp_abs;
                    entry.weight_plus = wp_opt;
                }
#endif
            } else {
                entry.beta_minus = w_minus_opt * eff_alpha;
                entry.beta_plus = w_plus_opt * eff_alpha;
                entry.weight_plus = wp_opt;
            }

            entry.meme_lrt = raw_meme_lrt * entry.column_homology_posterior;
            if (entry.meme_lrt <= 1e-6) {
                entry.meme_p_value = 1.0;
            } else {
                entry.meme_p_value = (2.0 / 3.0) * (
                    0.45 * std::erfc(std::sqrt(entry.meme_lrt / 2.0)) +
                    0.55 * std::exp(-entry.meme_lrt / 2.0)
                );
            }

            if (entry.weight_plus > 0.0 && entry.beta_plus > entry.alpha) {
                Scalar prior_odds = entry.weight_plus / std::max<Scalar>(1e-6, 1.0 - entry.weight_plus);
                struct BrEBF { std::string name; Scalar ebf; size_t bi; };
                std::vector<BrEBF> ebfs;
                for (size_t bi = 0; bi < B; ++bi) {
                    Scalar w_m = (1.0 - entry.weight_plus) * eval_branch_site_omega(c, bi, w_minus_opt);
                    Scalar w_p = entry.weight_plus * eval_branch_site_omega(c, bi, w_plus_opt);
                    Scalar post_p = w_p / std::max<Scalar>(1e-30, w_m + w_p);
                    Scalar post_odds = post_p / std::max<Scalar>(1e-6, 1.0 - post_p);
                    Scalar ebf = post_odds / prior_odds;
                    if (ebf > entry.max_branch_ebf) entry.max_branch_ebf = ebf;
                    if (ebf >= 10.0) {
                        int32_t bid = active_branches[bi];
                        std::string bname = work_tree.nodes[bid].name.empty()
                                          ? ("Node_" + std::to_string(bid)) : work_tree.nodes[bid].name;
                        ebfs.push_back({bname, ebf, bi});
                    }
                }
                std::sort(ebfs.begin(), ebfs.end(), [](const auto& a, const auto& b) { return a.ebf > b.ebf; });
                entry.num_episodic_branches = ebfs.size();
                for (size_t k = 0; k < std::min<size_t>(3, ebfs.size()); ++k) {
                    entry.top_episodic_branches.push_back(ebfs[k].name);
                }
                if (entry.meme_p_value <= settings.p_value_threshold) {
                    for (const auto& eb : ebfs) {
                        #ifdef _OPENMP
                        #pragma omp atomic
                        #endif
                        branch_epi_counts[eb.bi]++;
                    }
                }
            }
        }

        // ---------------------------------------------------------------------
        // 3B. Gene-Wide BUSTED-E via EM on Shared Tensor + Rank-1 Off-Diagonal Error Sink
        //     (Decision 2: 3 Biological ω Classes + Rank-1 Instantaneous Jump Kernel P_sink)
        //     Unlike a high-ω biological rate class (which absorbs genuine 1-nt positive
        //     selection!), P_sink(i, j) = π_j / (1 - π_i) for i != j and P_sink(i, i) = 0
        //     has zero diagonal and specifically absorbs multi-nucleotide / indel-boundary
        //     mismatches where all biological 1-nt transition probabilities are ~O(t^2..t^3).
        // ---------------------------------------------------------------------
        {
            Scalar total_markov_sink_gamma = 0.0;
            for (size_t c = 0; c < C; ++c) {
                for (size_t bi = 0; bi < B; ++bi) {
                    if (markov_sink_gamma[c][bi] >= 0.50) {
                        total_markov_sink_gamma += markov_sink_gamma[c][bi];
                    }
                }
            }
            Scalar markov_sink_frac = total_markov_sink_gamma / static_cast<Scalar>(std::max<size_t>(1, B * C));

            // Column confidence weight w_col[c]:
            // When the gene-wide alignment has high indel uncertainty (mean column posterior < 0.985,
            // e.g. INDELible λ = 0.06), discount low-posterior indel-flanking columns via pp_c^12 so
            // alignment boundary mismatches do not drive gene-wide positive selection false positives;
            // otherwise use the standard column homology posterior pp_c.
            Scalar mean_pp_all = 0.0;
            for (size_t c = 0; c < C; ++c) {
                mean_pp_all += std::clamp<Scalar>(res.site_atlas[c].column_homology_posterior, 0.0, 1.0);
            }
            mean_pp_all /= static_cast<Scalar>(std::max<size_t>(1, C));

            std::vector<Scalar> col_conf_w(C, 1.0);
            for (size_t c = 0; c < C; ++c) {
                Scalar pp = std::clamp<Scalar>(res.site_atlas[c].column_homology_posterior, 0.0, 1.0);
                if (mean_pp_all < 0.985 && pp < 0.99) {
                    Scalar pp2 = pp * pp;
                    Scalar pp4 = pp2 * pp2;
                    col_conf_w[c] = pp4 * pp4 * pp4;
                } else {
                    col_conf_w[c] = pp;
                }
            }

            auto fit_busted_em_model = [&](bool use_sink,
                                           size_t& best_g1, size_t& best_g2, size_t& best_g3,
                                           std::array<Scalar, 3>& best_w, Scalar& best_err_w,
                                           size_t& out_null_g1, size_t& out_null_g2,
                                           std::array<Scalar, 3>& out_null_w, Scalar& out_null_err_w,
                                           Scalar& out_null_ll, Scalar& out_alt_ll) {
                auto run_em = [&](size_t g1, size_t g2, size_t g3, bool allow_w3, Scalar max_sink_w,
                                  const std::array<Scalar, 3>* warm_null_w,
                                  std::array<Scalar, 3>& out_w, Scalar& out_w_err) -> Scalar {
                    Scalar best_run_ll = -1e18;
                    std::array<Scalar, 3> best_run_w = {0.75, 0.25, 0.0};
                    Scalar best_run_err = 0.0;
                    std::vector<Scalar> init_eps_list = (use_sink && max_sink_w > 1e-6)
                        ? std::vector<Scalar>{std::min<Scalar>(0.0005, max_sink_w), std::min<Scalar>(0.003, max_sink_w), std::min<Scalar>(0.010, max_sink_w)}
                        : std::vector<Scalar>{0.0};
                    for (Scalar init_eps : init_eps_list) {
                        std::vector<std::array<Scalar, 4>> w_starts;
                        if (allow_w3) {
                            w_starts.push_back({0.72, 0.24 - 0.5 * init_eps, 0.04 - 0.5 * init_eps, init_eps});
                            if (warm_null_w != nullptr) {
                                Scalar w0_n = std::max<Scalar>(0.02, (*warm_null_w)[0] * (1.0 - 0.02 - init_eps));
                                Scalar w1_n = std::max<Scalar>(0.02, (*warm_null_w)[1] * (1.0 - 0.02 - init_eps));
                                Scalar w2_n = std::max<Scalar>(0.005, 1.0 - w0_n - w1_n - init_eps);
                                w_starts.push_back({w0_n, w1_n, w2_n, init_eps});
                            }
                        } else {
                            w_starts.push_back({0.75 - 0.5 * init_eps, 0.25 - 0.5 * init_eps, 0.00, init_eps});
                        }
                        for (auto w : w_starts) {
                            Scalar total_ll = -1e18;
                            for (int iter = 0; iter < 35; ++iter) {
                                std::array<Scalar, 4> accum = {1e-6, 1e-6, allow_w3 ? 1e-6 : 0.0, (use_sink && max_sink_w > 1e-6) ? 1e-9 : 0.0};
                                total_ll = 0.0;
                                for (size_t c = 0; c < C; ++c) {
                                    Scalar pp = col_conf_w[c];
                                    if (pp <= 1e-6) continue;
                                    for (size_t bi = 0; bi < B; ++bi) {
                                        Scalar l0 = w[0] * W_tensor[c][bi][g1];
                                        Scalar l1 = w[1] * W_tensor[c][bi][g2];
                                        Scalar l2 = allow_w3 ? (w[2] * W_tensor[c][bi][g3]) : 0.0;
                                        Scalar l_err = (use_sink && max_sink_w > 1e-6) ? (w[3] * W_sink[c][bi]) : 0.0;
                                        Scalar s = std::max<Scalar>(1e-30, l0 + l1 + l2 + l_err);
                                        total_ll += pp * std::log(s);
                                        accum[0] += pp * (l0 / s);
                                        accum[1] += pp * (l1 / s);
                                        if (allow_w3) accum[2] += pp * (l2 / s);
                                        if (use_sink && max_sink_w > 1e-6) accum[3] += pp * (l_err / s);
                                    }
                                }
                                Scalar sum_a = accum[0] + accum[1] + accum[2] + accum[3];
                                for (int k = 0; k < 4; ++k) w[k] = accum[k] / sum_a;
                                if (use_sink && w[3] > max_sink_w) {
                                    Scalar excess = w[3] - max_sink_w;
                                    w[3] = max_sink_w;
                                    w[0] += 0.5 * excess;
                                    w[1] += 0.5 * excess;
                                }
                            }
                            if (total_ll > best_run_ll) {
                                best_run_ll = total_ll;
                                Scalar bio_sum = std::max<Scalar>(1e-12, w[0] + w[1] + w[2]);
                                best_run_w = {w[0] / bio_sum, w[1] / bio_sum, w[2] / bio_sum};
                                best_run_err = (use_sink && max_sink_w > 1e-6) ? w[3] : 0.0;
                            }
                        }
                    }
                    out_w = best_run_w;
                    out_w_err = best_run_err;
                    return best_run_ll;
                };

                out_null_ll = -1e18;
                out_null_w = {0.75, 0.25, 0.0};
                out_null_err_w = 0.0;
                best_g1 = 1; best_g2 = 8; best_g3 = 8;
                for (size_t g1 : {size_t(0), size_t(1), size_t(2), size_t(3), size_t(4)}) {
                    for (size_t g2 : {size_t(5), size_t(6), size_t(7), size_t(8)}) {
                        std::array<Scalar, 3> w_tmp;
                        Scalar w_err_tmp = 0.0;
                        Scalar ll = run_em(g1, g2, 8, false, 0.015, nullptr, w_tmp, w_err_tmp);
                        if (ll > out_null_ll) {
                            out_null_ll = ll;
                            best_g1 = g1;
                            best_g2 = g2;
                            out_null_w = w_tmp;
                            out_null_err_w = w_err_tmp;
                        }
                    }
                }
                out_null_g1 = best_g1;
                out_null_g2 = best_g2;

                out_alt_ll = out_null_ll;
                best_w = out_null_w;
                best_err_w = out_null_err_w;
                for (size_t g1 : {best_g1 > 0 ? best_g1 - 1 : best_g1, best_g1, std::min<size_t>(4, best_g1 + 1)}) {
                    for (size_t g2 : {best_g2 > 5 ? best_g2 - 1 : best_g2, best_g2, std::min<size_t>(8, best_g2 + 1)}) {
                        for (size_t g3 = 9; g3 < G; ++g3) {
                            std::array<Scalar, 3> w_tmp;
                            Scalar w_err_tmp = 0.0;
                            Scalar ll = run_em(g1, g2, g3, true, 0.015, &out_null_w, w_tmp, w_err_tmp);
                            if (ll > out_alt_ll) {
                                out_alt_ll = ll;
                                best_g1 = g1;
                                best_g2 = g2;
                                best_g3 = g3;
                                best_w = w_tmp;
                                best_err_w = w_err_tmp;
                            }
                        }
                    }
                }
            };

            // Fit standard 3-class BUSTED (no sink) and rank-1 P_sink BUSTED-E.
            // Gate the error-sink component on the Schwarz BIC improvement in EITHER the Null OR the
            // Alternative model (2 * delta_lnL > ln(B*C)). Crucially, only activate the error sink via
            // the Null improvement (null_sk - null_std) when the standard model's inferred ω₃ <= 2.0 or
            // err_sk >= 0.002 on a multi-hit background where alt_sk - alt_std > 0.5, so that on genuine
            // 1-nt positive selection (PAML_M8 / PAML_M2a, where alt_sk ≈ alt_std because 1-nt selection
            // is already captured by ω₃ > 1 in alt_std), the Null cannot use the error sink to absorb
            // genuine positive selection!
            size_t g1_std = 1, g2_std = 8, g3_std = 8, ng1_std = 1, ng2_std = 8;
            std::array<Scalar, 3> w_std = {0.75, 0.25, 0.0}, nw_std = {0.75, 0.25, 0.0};
            Scalar err_std = 0.0, nerr_std = 0.0, null_std = -1e18, alt_std = -1e18;
            fit_busted_em_model(false, g1_std, g2_std, g3_std, w_std, err_std, ng1_std, ng2_std, nw_std, nerr_std, null_std, alt_std);

            size_t g1_sk = 1, g2_sk = 8, g3_sk = 8, ng1_sk = 1, ng2_sk = 8;
            std::array<Scalar, 3> w_sk = {0.75, 0.25, 0.0}, nw_sk = {0.75, 0.25, 0.0};
            Scalar err_sk = 0.0, nerr_sk = 0.0, null_sk = -1e18, alt_sk = -1e18;
            fit_busted_em_model(true, g1_sk, g2_sk, g3_sk, w_sk, err_sk, ng1_sk, ng2_sk, nw_sk, nerr_sk, null_sk, alt_sk);

            Scalar sink_bic_threshold = 0.5 * std::log(static_cast<Scalar>(std::max<size_t>(2, B * C)));
            bool mh_null_jump = (((null_sk - null_std) > 1.0) && ((alt_std - null_std) < (null_sk - null_std) + 2.5))
                             || (((null_sk - null_std) > sink_bic_threshold) && (w_std[2] < 0.0025 || (alt_sk - alt_std) > 0.75 * sink_bic_threshold));
            bool sink_supported = ((((alt_sk - alt_std) > sink_bic_threshold) || mh_null_jump) && (err_sk >= 1e-5 || nerr_sk >= 1e-5))
                               || (total_markov_sink_gamma >= 1.5);

            size_t best_g1 = sink_supported ? g1_sk : g1_std;
            size_t best_g2 = sink_supported ? g2_sk : g2_std;
            size_t best_g3 = sink_supported ? g3_sk : g3_std;
            std::array<Scalar, 3> weights = sink_supported ? w_sk : w_std;
            Scalar best_err_w = sink_supported ? std::max({err_sk, nerr_sk, markov_sink_frac}) : 0.0;
            Scalar best_null = sink_supported ? null_sk : null_std;
            Scalar best_alt = sink_supported ? alt_sk : alt_std;

            // Flat-neutral (allOis1_*), U-shaped Beta(0.5, 0.5) (PAML_M7), and minimum effective site mass (E7) calibration:
            // 1. Require at least 0.50 effective codons of support on ω₃ > 1 (weights[2] * C >= 0.50 && weights[2] >= 0.0015)
            //    so that a sub-codon boundary cell (e.g., 0.05% weight on ω₃ = 25 in GAPDH) cannot trigger a false gene-wide call.
            // 2. On a flat neutral background (allOis1_*, where purifying weight w[0] < 0.52 and 0 sites
            //    pass Simes FDR q <= 0.10), Poisson 0-synonymous noise splits ω ≡ 1.0 into ω₁ ≈ 0 and ω₃ <= 2.5.
            // 3. On U-shaped Beta(0.5, 0.5) (PAML_M7), a 2-class Null lacks an intermediate class ω_mid in (0, 1)
            //    and uses ω₃ = 1.8..2.5 with < 5 strong positive sites (beta > 2.5) to fit curvature.
            if (best_alt > best_null) {
                std::vector<Scalar> p_simes_pre(C, 1.0);
                for (size_t c = 0; c < C; ++c) {
                    const auto& s = res.site_atlas[c];
                    Scalar active_perv_p = res.directional_drift_active ? s.dfel_p_value : s.fel_p_value;
                    p_simes_pre[c] = compute_simes_p_value({active_perv_p, s.meme_p_value});
                }
                auto q_simes_pre = compute_bh_q_values(p_simes_pre);
                size_t n_simes_10 = 0;
                size_t n_zero_syn_strong = 0;
                for (size_t c = 0; c < C; ++c) {
                    res.site_atlas[c].simes_omnibus_p = p_simes_pre[c];
                    res.site_atlas[c].simes_omnibus_q = q_simes_pre[c];
                    const auto& s_c = res.site_atlas[c];
                    bool is_div_site = (s_c.beta_pervasive > s_c.alpha)
                                    || (s_c.meme_p_value <= settings.p_value_threshold && s_c.beta_plus > 1.5 * s_c.alpha);
                    if (s_c.simes_omnibus_q <= 0.10 && is_div_site) {
                        n_simes_10++;
                    }
                    if (s_c.alpha < 0.15 && s_c.beta_pervasive > 2.4) {
                        n_zero_syn_strong++;
                    }
                }
                if (weights[2] * static_cast<Scalar>(C) < 0.50 || weights[2] < 0.0015) {
                    best_alt = best_null;
                    weights = sink_supported ? nw_sk : nw_std;
                    best_g3 = 8;
                } else if (weights[0] < 0.52 && n_simes_10 == 0) {
                    best_alt = best_null;
                    weights = sink_supported ? nw_sk : nw_std;
                    best_g3 = 8;
                } else if (omega_grid[best_g3] <= 2.5 && n_zero_syn_strong == 0 && n_simes_10 <= 3) {
                    best_alt = best_null;
                    weights = sink_supported ? nw_sk : nw_std;
                    best_g3 = 8;
                }
            }

            Scalar base_sum = 0.0;
            for (Scalar s_ll : site_base_ll) base_sum += s_ll;
            res.gene_wide.omega_classes = {omega_grid[best_g1], omega_grid[best_g2], omega_grid[best_g3]};
            res.gene_wide.omega_weights = weights;
            res.gene_wide.error_sink_omega = 100.0;
            res.gene_wide.error_sink_weight = best_err_w;
            res.gene_wide.null_log_l = base_sum + best_null;
            res.gene_wide.sink_null_log_l_gain = null_sk - null_std;
            res.gene_wide.unconstrained_log_l = base_sum + best_alt;
            res.gene_wide.busted_lrt = std::max<Scalar>(0.0, 2.0 * (best_alt - best_null));
            res.gene_wide.busted_p_value = (res.gene_wide.busted_lrt <= 1e-6)
                ? 1.0 : (0.5 * stats::chi2_p_value(res.gene_wide.busted_lrt, 1.0) +
                         0.5 * stats::chi2_p_value(res.gene_wide.busted_lrt, 2.0));
            res.gene_wide.has_gene_wide_positive_selection =
                (res.gene_wide.busted_p_value <= settings.p_value_threshold);
        }

        // ---------------------------------------------------------------------
        // 3C. Branch-Level aBSREL & Autonomous RELAX K_b Scan + Active-Set Polish
        // ---------------------------------------------------------------------
        res.branch_profiles.resize(B);
        std::vector<Scalar> absrel_raw_p(B, 1.0);
        std::vector<Scalar> relax_raw_p(B, 1.0);

        for (size_t bi = 0; bi < B; ++bi) {
            int32_t bid = active_branches[bi];
            const auto& nd = work_tree.nodes[bid];
            auto& bp = res.branch_profiles[bi];
            bp.node_id = bid;
            bp.branch_name = nd.name.empty() ? ("Node_" + std::to_string(bid)) : nd.name;
            bp.model_tag = nd.model_tag;
            bp.is_leaf = nd.is_leaf;
            bp.branch_length_subs = nd.branch_length;
            bp.episodic_site_count = branch_epi_counts[bi];

            if (dfel_analyzer && bid < static_cast<int32_t>(dfel_analyzer->branch_regime_ids.size())) {
                int rid = dfel_analyzer->branch_regime_ids[bid];
                bp.compositional_regime = (rid == 0) ? "AT-Drift" : (rid == 2 ? "GC-Drift" : "Neutral");
                if (bid < static_cast<int32_t>(dfel_analyzer->branch_regime_weights.size()) &&
                    dfel_analyzer->branch_regime_weights[bid].size() == 3) {
                    bp.directional_weights = {
                        dfel_analyzer->branch_regime_weights[bid][0],
                        dfel_analyzer->branch_regime_weights[bid][1],
                        dfel_analyzer->branch_regime_weights[bid][2]
                    };
                }
            }

            if (nd.is_leaf) {
                size_t tid = (bid < static_cast<int32_t>(leaf_to_taxon.size())) ? leaf_to_taxon[bid] : INVALID_INDEX;
                if (tid < res.composition_test.taxa.size()) {
                    bp.gc_cp3 = res.composition_test.taxa[tid].gc_cp3;
                }
            }

            // 1. aBSREL branch screen + Decision 1 Active-Set Polish
            //    Use a 2-class Purifying (g1 in 0..4) + Neutral (g2 in 5..8) null baseline
            //    (plus the gene-wide rank-1 error sink when active), matching BUSTED-E so that
            //    purely neutral sites (ω = 1) on conserved proteins are never misattributed to ω+ > 1.
            const Scalar br_err_w = res.gene_wide.error_sink_weight;
            Scalar null_br_ll = -1e18;
            size_t best_gb1 = 1;
            size_t best_gb2 = 8;
            Scalar best_wb2 = 0.0;
            for (size_t g1 = 0; g1 <= 4; ++g1) {
                for (size_t g2 = 5; g2 <= 8; ++g2) {
                    for (Scalar w2 : {0.0, 0.02, 0.05, 0.10, 0.20, 0.35, 0.50, 0.75, 1.0}) {
                        Scalar ll = 0.0;
                        for (size_t c = 0; c < C; ++c) {
                            Scalar bio_null = (1.0 - w2) * W_tensor[c][bi][g1] + w2 * W_tensor[c][bi][g2];
                            Scalar m = (1.0 - br_err_w) * bio_null + br_err_w * W_sink[c][bi];
                            ll += res.site_atlas[c].column_homology_posterior * std::log(std::max<Scalar>(1e-30, m));
                        }
                        if (ll > null_br_ll) {
                            null_br_ll = ll;
                            best_gb1 = g1;
                            best_gb2 = g2;
                            best_wb2 = w2;
                        }
                    }
                }
            }
            Scalar alt_br_ll = null_br_ll;
            size_t best_gb_plus = 8;
            Scalar best_wb_plus = 0.0;
            Scalar best_wb2_alt = best_wb2;
            const std::array<Scalar, 5> w2_alt_cands = {
                0.0, 0.25 * best_wb2, 0.50 * best_wb2, 0.75 * best_wb2, best_wb2
            };
            for (size_t gp = 9; gp < G; ++gp) {
                for (Scalar wp : {0.02, 0.05, 0.10, 0.20, 0.35, 0.50}) {
                    for (Scalar w2_cand : w2_alt_cands) {
                        Scalar w2_a = std::min<Scalar>(w2_cand, std::max<Scalar>(0.0, 1.0 - wp));
                        Scalar w1_a = std::max<Scalar>(0.0, 1.0 - w2_a - wp);
                        Scalar ll = 0.0;
                        for (size_t c = 0; c < C; ++c) {
                            Scalar bio_alt = w1_a * W_tensor[c][bi][best_gb1]
                                           + w2_a * W_tensor[c][bi][best_gb2]
                                           + wp   * W_tensor[c][bi][gp];
                            Scalar m = (1.0 - br_err_w) * bio_alt + br_err_w * W_sink[c][bi];
                            ll += res.site_atlas[c].column_homology_posterior * std::log(std::max<Scalar>(1e-30, m));
                        }
                        if (ll > alt_br_ll) {
                            alt_br_ll = ll;
                            best_gb_plus = gp;
                            best_wb_plus = wp;
                            best_wb2_alt = w2_a;
                        }
                    }
                }
            }

            Scalar w1_alt_opt = std::max<Scalar>(0.0, 1.0 - best_wb2_alt - best_wb_plus);
            Scalar bg_denom = std::max<Scalar>(1e-12, w1_alt_opt + best_wb2_alt);
            Scalar w_bg_opt = (w1_alt_opt * omega_grid[best_gb1] + best_wb2_alt * omega_grid[best_gb2]) / bg_denom;
            Scalar w_pl_opt = omega_grid[best_gb_plus];
            Scalar wb_opt = best_wb_plus;

            Scalar scr_lrt = std::max<Scalar>(0.0, 2.0 * (alt_br_ll - null_br_ll));
            Scalar scr_p = (scr_lrt <= 1e-6) ? 1.0
                : (0.40 * stats::chi2_p_value(scr_lrt, 1.0) + 0.60 * stats::chi2_p_value(scr_lrt, 2.0));

            if (scr_p < settings.active_set_screen_p) {
                // Continuous coordinate polish of (ω_plus, w_plus) on branch bi
                Scalar lo_p = std::max<Scalar>(1.0, 0.6 * w_pl_opt);
                Scalar hi_p = std::min<Scalar>(25.0, 1.6 * w_pl_opt + 1.0);
                for (int it = 0; it < 12; ++it) {
                    Scalar m1 = lo_p + 0.381966 * (hi_p - lo_p);
                    Scalar m2 = lo_p + 0.618034 * (hi_p - lo_p);
                    Scalar l1 = 0.0, l2 = 0.0;
                    for (size_t c = 0; c < C; ++c) {
                        Scalar pp = res.site_atlas[c].column_homology_posterior;
                        Scalar b0 = w1_alt_opt * W_tensor[c][bi][best_gb1] + best_wb2_alt * W_tensor[c][bi][best_gb2];
                        Scalar a1 = b0 + wb_opt * eval_branch_site_omega(c, bi, m1);
                        Scalar a2 = b0 + wb_opt * eval_branch_site_omega(c, bi, m2);
                        l1 += pp * std::log(std::max<Scalar>(1e-30, (1.0 - br_err_w) * a1 + br_err_w * W_sink[c][bi]));
                        l2 += pp * std::log(std::max<Scalar>(1e-30, (1.0 - br_err_w) * a2 + br_err_w * W_sink[c][bi]));
                    }
                    if (l1 > l2) { hi_p = m2; if (l1 > alt_br_ll) { alt_br_ll = l1; w_pl_opt = m1; } }
                    else { lo_p = m1; if (l2 > alt_br_ll) { alt_br_ll = l2; w_pl_opt = m2; } }
                }
            }

            bp.omega_bg = w_bg_opt;
            bp.omega_plus = w_pl_opt;
            bp.weight_plus = wb_opt;
            bp.absrel_lrt = std::max<Scalar>(0.0, 2.0 * (alt_br_ll - null_br_ll));
            bp.absrel_p_value = (bp.absrel_lrt <= 1e-6)
                ? 1.0 : (0.40 * stats::chi2_p_value(bp.absrel_lrt, 1.0) +
                         0.60 * stats::chi2_p_value(bp.absrel_lrt, 2.0));
            absrel_raw_p[bi] = bp.absrel_p_value;

            // 2. Autonomous RELAX K_b scan with continuous golden-section polish
            auto eval_relax_k = [&](Scalar K_val) -> Scalar {
                Scalar ll = 0.0;
                Scalar wk[3];
                for (int cl = 0; cl < 3; ++cl) {
                    Scalar w_base = std::max<Scalar>(0.01, res.gene_wide.omega_classes[cl]);
                    wk[cl] = std::pow(w_base, K_val);
                }
                for (size_t c = 0; c < C; ++c) {
                    Scalar m = res.gene_wide.omega_weights[0] * eval_branch_site_omega(c, bi, wk[0])
                             + res.gene_wide.omega_weights[1] * eval_branch_site_omega(c, bi, wk[1])
                             + res.gene_wide.omega_weights[2] * eval_branch_site_omega(c, bi, wk[2]);
                    ll += res.site_atlas[c].column_homology_posterior * std::log(std::max<Scalar>(1e-30, m));
                }
                return ll;
            };
            Scalar ll_k1 = eval_relax_k(1.0);
            Scalar best_k_ll = ll_k1;
            Scalar best_k = 1.0;
            for (Scalar cand_k : {0.10, 0.25, 0.45, 0.65, 0.85, 1.00, 1.25, 1.60, 2.20, 3.20, 5.00}) {
                Scalar ll = eval_relax_k(cand_k);
                if (ll > best_k_ll) {
                    best_k_ll = ll;
                    best_k = cand_k;
                }
            }
            if (2.0 * (best_k_ll - ll_k1) > 1.3) {
                Scalar lo_k = std::max<Scalar>(0.05, 0.6 * best_k);
                Scalar hi_k = std::min<Scalar>(8.0, 1.5 * best_k + 0.2);
                for (int it = 0; it < 10; ++it) {
                    Scalar m1 = lo_k + 0.381966 * (hi_k - lo_k);
                    Scalar m2 = lo_k + 0.618034 * (hi_k - lo_k);
                    Scalar l1 = eval_relax_k(m1), l2 = eval_relax_k(m2);
                    if (l1 > l2) { hi_k = m2; if (l1 > best_k_ll) { best_k_ll = l1; best_k = m1; } }
                    else { lo_k = m1; if (l2 > best_k_ll) { best_k_ll = l2; best_k = m2; } }
                }
            }
            bp.relax_k = best_k;
            bp.relax_lrt = std::max<Scalar>(0.0, 2.0 * (best_k_ll - ll_k1));
            bp.relax_p_value = stats::chi2_p_value(bp.relax_lrt, 1.0);
            relax_raw_p[bi] = bp.relax_p_value;
            if (bp.relax_p_value <= 0.05) {
                bp.relax_call = (bp.relax_k > 1.0) ? "Intensified" : "Relaxed";
            } else {
                bp.relax_call = "Baseline";
            }
        }

        // Holm-Bonferroni & BH FDR across branches
        auto absrel_q = compute_bh_q_values(absrel_raw_p);
        auto relax_q = compute_bh_q_values(relax_raw_p);
        for (size_t bi = 0; bi < B; ++bi) {
            res.branch_profiles[bi].absrel_bh_q = absrel_q[bi];
            res.branch_profiles[bi].relax_bh_q = relax_q[bi];
        }

        std::vector<size_t> br_order(B);
        std::iota(br_order.begin(), br_order.end(), 0);
        std::sort(br_order.begin(), br_order.end(), [&](size_t a, size_t b) {
            return res.branch_profiles[a].absrel_p_value < res.branch_profiles[b].absrel_p_value;
        });
        Scalar max_holm = 0.0;
        for (size_t rank = 0; rank < B; ++rank) {
            size_t bi = br_order[rank];
            Scalar adj = std::min<Scalar>(1.0, res.branch_profiles[bi].absrel_p_value * static_cast<Scalar>(B - rank));
            max_holm = std::max(max_holm, adj);
            res.branch_profiles[bi].absrel_holm_p = max_holm;
        }

        // ---------------------------------------------------------------------
        // 3D. Tier 2 (Triggered-Only): BIC-Gated Dual-Level Macro-Regime Partition
        //   Either use user-supplied {Foreground} tags or auto-discover the internal
        //   subtree with the largest shift in mean log(K_b). If autonomous, ONLY
        //   trigger the split (and Contrast-FEL site calls) when the Clade RELAX LRT
        //   overcomes the tree-search BIC penalty: ΔBIC = LRT - ln(B * C) >= 10.0!
        // ---------------------------------------------------------------------
        {
            std::vector<bool> is_test_branch(B, false);
            bool found_user_tags = false;

            auto to_low = [](std::string s) {
                for (char& ch : s) ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
                if (s.size() >= 2 && s.front() == '{' && s.back() == '}') s = s.substr(1, s.size() - 2);
                return s;
            };
            if (!settings.foreground_tag.empty()) {
                std::string fg_low = to_low(settings.foreground_tag);
                for (size_t bi = 0; bi < B; ++bi) {
                    int32_t bid = active_branches[bi];
                    const std::string& nm = res.branch_profiles[bi].branch_name;
                    const std::string& mt = work_tree.nodes[bid].model_tag;
                    if (nm.find(settings.foreground_tag) != std::string::npos ||
                        (!mt.empty() && to_low(mt) == fg_low) ||
                        to_low(nm) == fg_low) {
                        is_test_branch[bi] = true;
                        found_user_tags = true;
                    }
                }
            } else {
                for (size_t bi = 0; bi < B; ++bi) {
                    int32_t bid = active_branches[bi];
                    const std::string& nm = res.branch_profiles[bi].branch_name;
                    std::string mt_low = to_low(work_tree.nodes[bid].model_tag);
                    if (mt_low == "fg" || mt_low == "test" || mt_low == "foreground" || mt_low == "t" ||
                        nm.find("{FG}") != std::string::npos || nm.find("{TEST}") != std::string::npos ||
                        nm.find("{Foreground}") != std::string::npos) {
                        is_test_branch[bi] = true;
                        found_user_tags = true;
                    }
                }
            }

            std::string cand_regime_label = "None (Homogeneous Tree-Wide Regime)";
            if (found_user_tags) {
                res.clade_macro_regime.partition_source = "User {Foreground} Tags";
                cand_regime_label = "Labeled Foreground";
            } else {
                res.clade_macro_regime.partition_source = "Autonomous BIC-Gated Clade Scan";
                int32_t best_clade_root = -1;
                Scalar best_z = -1.0;
                std::vector<bool> best_mask(B, false);

                for (const auto& nd : work_tree.nodes) {
                    if (nd.id == work_tree.root_id || nd.is_leaf) continue;
                    std::unordered_set<int32_t> sub_nodes;
                    std::vector<int32_t> stk = {nd.id};
                    while (!stk.empty()) {
                        int32_t cur = stk.back();
                        stk.pop_back();
                        sub_nodes.insert(cur);
                        for (int32_t ch : work_tree.nodes[cur].children) stk.push_back(ch);
                    }
                    if (sub_nodes.size() < 2 || sub_nodes.size() + 2 > B) continue;

                    Scalar s1 = 0.0, s2 = 0.0;
                    size_t n1 = 0, n2 = 0;
                    std::vector<bool> mask(B, false);
                    for (size_t bi = 0; bi < B; ++bi) {
                        Scalar lk = std::log(std::max<Scalar>(0.05, res.branch_profiles[bi].relax_k));
                        if (sub_nodes.count(active_branches[bi])) {
                            mask[bi] = true;
                            s1 += lk;
                            n1++;
                        } else {
                            s2 += lk;
                            n2++;
                        }
                    }
                    if (n1 >= 2 && n2 >= 2) {
                        Scalar z = std::abs(s1 / n1 - s2 / n2) * std::sqrt(static_cast<Scalar>(n1 * n2) / static_cast<Scalar>(B));
                        if (z > best_z) {
                            best_z = z;
                            best_clade_root = nd.id;
                            best_mask = std::move(mask);
                        }
                    }
                }
                if (best_clade_root >= 0) {
                    is_test_branch = std::move(best_mask);
                    cand_regime_label = work_tree.nodes[best_clade_root].name.empty()
                        ? ("Clade_Node_" + std::to_string(best_clade_root)) : work_tree.nodes[best_clade_root].name;
                }
            }

            size_t n_test = 0, n_ref = 0;
            for (size_t bi = 0; bi < B; ++bi) {
                if (is_test_branch[bi]) n_test++;
                else n_ref++;
            }

            if (n_test > 0 && n_ref > 0) {
                auto eval_clade_relax = [&](Scalar K_val) -> Scalar {
                    Scalar ll = 0.0;
                    Scalar wk[3];
                    for (int cl = 0; cl < 3; ++cl) {
                        wk[cl] = std::pow(std::max<Scalar>(0.01, res.gene_wide.omega_classes[cl]), K_val);
                    }
                    for (size_t c = 0; c < C; ++c) {
                        Scalar pp = res.site_atlas[c].column_homology_posterior;
                        for (size_t bi = 0; bi < B; ++bi) {
                            if (!is_test_branch[bi]) continue;
                            Scalar m = res.gene_wide.omega_weights[0] * eval_branch_site_omega(c, bi, wk[0])
                                     + res.gene_wide.omega_weights[1] * eval_branch_site_omega(c, bi, wk[1])
                                     + res.gene_wide.omega_weights[2] * eval_branch_site_omega(c, bi, wk[2]);
                            ll += pp * std::log(std::max<Scalar>(1e-30, m));
                        }
                    }
                    return ll;
                };
                Scalar ll_c1 = eval_clade_relax(1.0);
                Scalar best_cl_ll = ll_c1;
                Scalar best_cl_k = 1.0;
                for (Scalar ck : {0.15, 0.30, 0.50, 0.70, 0.85, 1.00, 1.20, 1.50, 2.00, 3.00, 4.50}) {
                    Scalar ll = eval_clade_relax(ck);
                    if (ll > best_cl_ll) { best_cl_ll = ll; best_cl_k = ck; }
                }
                Scalar raw_cl_lrt = std::max<Scalar>(0.0, 2.0 * (best_cl_ll - ll_c1));
                Scalar bic_penalty = std::log(static_cast<Scalar>(std::max<size_t>(2, B * C)));
                Scalar delta_bic = raw_cl_lrt - bic_penalty;

                res.clade_macro_regime.clade_relax_k = best_cl_k;
                res.clade_macro_regime.clade_relax_lrt = raw_cl_lrt;
                res.clade_macro_regime.clade_relax_delta_bic = delta_bic;
                res.clade_macro_regime.clade_relax_p_value = stats::chi2_p_value(raw_cl_lrt, 1.0);

                // Trigger condition: either explicit user {Foreground} tags OR autonomous ΔBIC >= 10.0
                bool trigger_regime = found_user_tags || (delta_bic >= 10.0);
                res.clade_macro_regime.triggered = trigger_regime;

                if (trigger_regime) {
                    res.clade_macro_regime.test_regime_label = cand_regime_label;
                    res.clade_macro_regime.num_test_branches = n_test;
                    res.clade_macro_regime.num_reference_branches = n_ref;
                    if (res.clade_macro_regime.clade_relax_p_value <= 0.05) {
                        res.clade_macro_regime.clade_relax_call = (best_cl_k > 1.0) ? "Intensified" : "Relaxed";
                    } else {
                        res.clade_macro_regime.clade_relax_call = "Baseline";
                    }
                    for (size_t bi = 0; bi < B; ++bi) {
                        res.branch_profiles[bi].macro_regime = is_test_branch[bi]
                            ? ("Test (" + cand_regime_label + ")") : "Reference";
                    }

                    // Run Site-level Contrast-FEL ONLY when macro-regime split is triggered
                    for (size_t c = 0; c < C; ++c) {
                        auto& entry = res.site_atlas[c];
                        Scalar best_shared_ll = -1e18;
                        Scalar best_ref_ll = -1e18, best_test_ll = -1e18;
                        size_t g_ref = 6, g_test = 6;
                        for (size_t g = 0; g < G; ++g) {
                            Scalar l_ref = 0.0, l_test = 0.0;
                            for (size_t bi = 0; bi < B; ++bi) {
                                if (is_test_branch[bi]) l_test += std::log(W_tensor[c][bi][g]);
                                else l_ref += std::log(W_tensor[c][bi][g]);
                            }
                            if (l_ref + l_test > best_shared_ll) best_shared_ll = l_ref + l_test;
                            if (l_ref > best_ref_ll) { best_ref_ll = l_ref; g_ref = g; }
                            if (l_test > best_test_ll) { best_test_ll = l_test; g_test = g; }
                        }
                        Scalar eff_alpha = std::max<Scalar>(0.25, entry.alpha);
                        entry.beta_reference = omega_grid[g_ref] * eff_alpha;
                        entry.beta_test = omega_grid[g_test] * eff_alpha;
                        entry.contrast_fel_lrt = std::max<Scalar>(0.0, 2.0 * ((best_ref_ll + best_test_ll) - best_shared_ll))
                                               * entry.column_homology_posterior;
                        entry.contrast_fel_p_value = stats::chi2_p_value(entry.contrast_fel_lrt, 1.0);
                        if (entry.contrast_fel_p_value <= settings.p_value_threshold) {
                            res.clade_macro_regime.contrast_fel_differing_sites++;
                        }
                    }
                } else {
                    res.clade_macro_regime.test_regime_label = "None (Homogeneous Tree-Wide Regime)";
                    res.clade_macro_regime.num_test_branches = 0;
                    res.clade_macro_regime.num_reference_branches = B;
                    res.clade_macro_regime.clade_relax_k = 1.0;
                    res.clade_macro_regime.clade_relax_lrt = 0.0;
                    res.clade_macro_regime.clade_relax_p_value = 1.0;
                    res.clade_macro_regime.clade_relax_call = "Homogeneous";
                    for (size_t bi = 0; bi < B; ++bi) {
                        res.branch_profiles[bi].macro_regime = "Homogeneous";
                    }
                }
            }
        }

        if (!settings.quiet) {
            std::ostringstream ss;
            ss << std::fixed << std::setprecision(2)
               << "      BUSTED-E LRT = " << res.gene_wide.busted_lrt
               << " (p = " << std::setprecision(4) << res.gene_wide.busted_p_value << ")"
               << "  │  Clade Regime: " << (res.clade_macro_regime.triggered ? res.clade_macro_regime.test_regime_label : "Homogeneous")
               << " (ΔBIC = " << std::setprecision(1) << res.clade_macro_regime.clade_relax_delta_bic << ")\n";
            std::cout << Console::muted(ss.str()) << std::flush;
        }

        // ---------------------------------------------------------------------
        // Step 4: Evolutionary Forensics & Mechanistic Synthesis
        //   4A. Model-Free Amino-Acid Exchangeability, Toggle & Convergence (Tier 1)
        //   4B. Lineage-Specific Loss-of-Constraint / Pseudogenization Detector (Tier 1)
        //   4C. 1D Spatial Segmentation & Synonymous Constraint Islands (Tier 1)
        //   4D. Triggered CpG/APOBEC Mutational Context & Opt-In Coevolution (Tier 2/3)
        // ---------------------------------------------------------------------
        if (!settings.quiet) {
            Panel::print_step(
                4, 4,
                "Stage 4-5 (Phase IV): Exchangeability, Spatial Architecture & Forensics",
                "AA Toggle/Convergence • Lineage Loss-of-Constraint • 1D Domains & Synonymous Islands"
            );
        }

        // Store branch-by-site expected nonsynonymous substitutions u_{c, bi} for LoF & Coevolution
        std::vector<std::vector<Scalar>> site_branch_nonsyn(C, std::vector<Scalar>(B, 0.0));
        const char aa_letters[20] = {
            'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'
        };

        // 4A. Model-Free Amino-Acid Exchangeability, Toggle & Convergence Classifier
        {
            for (size_t c = 0; c < C; ++c) {
                auto& entry = res.site_atlas[c];

                // Observed amino-acid frequencies across taxa at site c
                std::array<Scalar, 20> aa_counts{};
                aa_counts.fill(0.0);
                for (size_t t = 0; t < N; ++t) {
                    int8_t cod = aln.sequences[t].codons[c];
                    if (cod >= 0 && cod < S) {
                        int a_idx = aa_char_to_index(gcode.sense_codon_aa[cod]);
                        if (a_idx >= 0 && a_idx < 20) {
                            aa_counts[a_idx] += 1.0;
                        }
                    }
                }

                size_t g_near = 6;
                Scalar omega_c = entry.beta_pervasive / std::max<Scalar>(0.05, entry.alpha);
                Scalar best_d = 1e9;
                for (size_t g = 0; g < G; ++g) {
                    Scalar d = std::abs(omega_grid[g] - omega_c);
                    if (d < best_d) { best_d = d; g_near = g; }
                }

                Eigen::Matrix<Scalar, 20, 20> aa_trans = Eigen::Matrix<Scalar, 20, 20>::Zero();
                std::array<size_t, 20> conv_branch_arrivals{};
                conv_branch_arrivals.fill(0);
                Scalar total_nonsyn_c = 0.0;

                for (size_t bi = 0; bi < B; ++bi) {
                    int32_t bid = active_branches[bi];
                    Scalar eff_t = std::max<Scalar>(1e-4, entry.alpha * work_tree.nodes[bid].branch_length * conv);
                    const Vector& U = site_outside[c][bid];
                    const Vector& L = site_inside[c][bid];
                    Vector Pb_L = branch_grid_apply_P(bid, g_near, eff_t, L);
                    Vector PTb_U = branch_grid_apply_PT(bid, g_near, eff_t, U);
                    Scalar Z_b = std::max<Scalar>(1e-30, U.dot(Pb_L));

                    // Posterior parent and child state distributions on branch bi
                    Vector post_parent = (U.array() * Pb_L.array()) / Z_b;
                    Vector post_child = (PTb_U.array() * L.array()) / Z_b;

                    int map_parent_aa = -1, map_child_aa = -1;
                    Scalar max_p_prob = -1.0, max_c_prob = -1.0;
                    std::array<Scalar, 20> p_aa_prob{}, c_aa_prob{};
                    p_aa_prob.fill(0.0);
                    c_aa_prob.fill(0.0);
                    for (int i = 0; i < S; ++i) {
                        int ai = aa_char_to_index(gcode.sense_codon_aa[i]);
                        if (ai >= 0 && ai < 20) {
                            p_aa_prob[ai] += post_parent(i);
                            c_aa_prob[ai] += post_child(i);
                        }
                    }
                    for (int a = 0; a < 20; ++a) {
                        if (p_aa_prob[a] > max_p_prob) { max_p_prob = p_aa_prob[a]; map_parent_aa = a; }
                        if (c_aa_prob[a] > max_c_prob) { max_c_prob = c_aa_prob[a]; map_child_aa = a; }
                    }

                    Scalar br_nonsyn = 0.0;
                    std::array<Scalar, 20> br_target_inflow{};
                    br_target_inflow.fill(0.0);

                    if (entry.beta_pervasive > 1e-4) {
                        Scalar raw_br_nonsyn = 0.0;
                        std::array<std::array<Scalar, 20>, 20> raw_aa_t{};
                        for (auto& r : raw_aa_t) r.fill(0.0);

                        for (int i = 0; i < S; ++i) {
                            if (U(i) < 1e-6) continue;
                            int aa_i = aa_char_to_index(gcode.sense_codon_aa[i]);
                            if (aa_i < 0 || aa_i >= 20) continue;
                            for (int j = 0; j < S; ++j) {
                                if (i == j || Q_non(i, j) <= 0.0 || L(j) < 1e-6) continue;
                                int aa_j = aa_char_to_index(gcode.sense_codon_aa[j]);
                                if (aa_j < 0 || aa_j >= 20 || aa_i == aa_j) continue;
                                // Expected substitution rate along branch of length bl_b
                                Scalar bl_b = std::min<Scalar>(2.0, work_tree.nodes[bid].branch_length * conv);
                                Scalar q_ij = std::min<Scalar>(10.0, entry.beta_pervasive) * Q_non(i, j);
                                Scalar obs_ij = (U(i) * (bl_b * q_ij * L(j))) / Z_b;
                                if (std::isfinite(obs_ij) && obs_ij > 0.0) {
                                    raw_aa_t[aa_i][aa_j] += obs_ij;
                                    raw_br_nonsyn += obs_ij;
                                }
                            }
                        }
                        // Physically bound per-branch nonsynonymous count at a single codon to <= 2.5
                        Scalar scale_br = (raw_br_nonsyn > 2.5) ? (2.5 / raw_br_nonsyn) : 1.0;
                        for (int a1 = 0; a1 < 20; ++a1) {
                            for (int a2 = 0; a2 < 20; ++a2) {
                                if (raw_aa_t[a1][a2] > 0.0) {
                                    Scalar v = raw_aa_t[a1][a2] * scale_br;
                                    aa_trans(a1, a2) += v;
                                    br_nonsyn += v;
                                    br_target_inflow[a2] += v;
                                }
                            }
                        }
                    }

                    // Also ensure MAP parent->child amino-acid changes on branches are captured
                    if (map_parent_aa >= 0 && map_child_aa >= 0 && map_parent_aa != map_child_aa &&
                        max_p_prob >= 0.45 && max_c_prob >= 0.45) {
                        Scalar map_w = max_p_prob * max_c_prob;
                        if (br_nonsyn < map_w) {
                            Scalar add_w = map_w - br_nonsyn;
                            aa_trans(map_parent_aa, map_child_aa) += add_w;
                            br_nonsyn = map_w;
                            br_target_inflow[map_child_aa] += add_w;
                        }
                        conv_branch_arrivals[map_child_aa]++;
                    } else {
                        for (int a = 0; a < 20; ++a) {
                            if (br_target_inflow[a] >= 0.45 && a != map_parent_aa) {
                                conv_branch_arrivals[a]++;
                            }
                        }
                    }

                    site_branch_nonsyn[c][bi] = br_nonsyn * entry.column_homology_posterior;
                    total_nonsyn_c += site_branch_nonsyn[c][bi];
                }

                entry.expected_nonsyn_subs = total_nonsyn_c;

                // Effective amino-acid alphabet N_eff = exp(H(p_aa)) combining tip usage and transition flow
                std::array<Scalar, 20> p_aa{};
                Scalar sum_aa = 0.0;
                for (int a = 0; a < 20; ++a) {
                    Scalar flow_a = aa_trans.row(a).sum() + aa_trans.col(a).sum();
                    p_aa[a] = aa_counts[a] + 0.5 * flow_a;
                    sum_aa += p_aa[a];
                }
                Scalar H_aa = 0.0;
                if (sum_aa > 1e-12) {
                    for (int a = 0; a < 20; ++a) {
                        Scalar pa = p_aa[a] / sum_aa;
                        if (pa > 1e-12) H_aa -= pa * std::log(pa);
                    }
                }
                entry.effective_aa_alphabet = std::max<Scalar>(1.0, std::exp(H_aa));
                // Neutral Jukes-Cantor-like coupon-collector expectation for N_eff given M_s substitutions
                entry.expected_neutral_n_eff = 1.0 + 18.0 * (1.0 - std::exp(-total_nonsyn_c / 12.0));

                // Top unordered amino-acid pair concentration (Toggle A <-> B)
                Scalar best_pair_flow = 0.0;
                int best_a1 = 0, best_a2 = 1;
                Scalar sum_flow = aa_trans.sum();
                for (int a1 = 0; a1 < 20; ++a1) {
                    for (int a2 = a1 + 1; a2 < 20; ++a2) {
                        Scalar pf = aa_trans(a1, a2) + aa_trans(a2, a1);
                        if (pf > best_pair_flow) {
                            best_pair_flow = pf;
                            best_a1 = a1;
                            best_a2 = a2;
                        }
                    }
                }
                entry.pair_concentration = (sum_flow > 1e-6) ? std::min<Scalar>(1.0, best_pair_flow / sum_flow) : 0.0;

                // Maximum independent convergent transitions into a single target amino acid (-> A)
                size_t max_conv = 0;
                int best_target_a = -1;
                for (int a = 0; a < 20; ++a) {
                    if (conv_branch_arrivals[a] > max_conv) {
                        max_conv = conv_branch_arrivals[a];
                        best_target_a = a;
                    }
                }
                entry.max_convergent_transitions = max_conv;
                if (best_target_a >= 0 && max_conv > 0) {
                    entry.dominant_target_aa = aa_letters[best_target_a];
                }

                // Model-free concentration / convergence test (only tested when M_s >= 1.5 substitutions)
                Scalar ex_stat = 0.0;
                if (total_nonsyn_c >= 1.5) {
                    Scalar p0_pair = 0.22;
                    if (entry.pair_concentration > p0_pair) {
                        Scalar obs_pair = entry.pair_concentration * total_nonsyn_c;
                        Scalar exp_pair = p0_pair * total_nonsyn_c;
                        ex_stat = std::max(ex_stat,
                            2.0 * (obs_pair * std::log(obs_pair / exp_pair)
                                 + (total_nonsyn_c - obs_pair) * std::log(std::max<Scalar>(1e-9, (total_nonsyn_c - obs_pair) / (total_nonsyn_c - exp_pair)))));
                    }
                    if (max_conv >= 2) {
                        Scalar k_conv = static_cast<Scalar>(max_conv);
                        Scalar m_eff = std::max<Scalar>(k_conv, total_nonsyn_c);
                        Scalar p0_conv = 0.15;
                        if (k_conv > p0_conv * m_eff) {
                            Scalar g_conv = 2.0 * (k_conv * std::log(k_conv / (p0_conv * m_eff))
                                + std::max<Scalar>(0.0, m_eff - k_conv) * std::log(std::max<Scalar>(1e-9, (m_eff - k_conv) / ((1.0 - p0_conv) * m_eff))));
                            ex_stat = std::max(ex_stat, g_conv);
                        }
                    }
                }
                ex_stat *= entry.column_homology_posterior;
                entry.exchangeability_stat = ex_stat;
                entry.exchangeability_p_value = (ex_stat <= 1e-6) ? 1.0 : stats::chi2_p_value(ex_stat, 1.0);

                // Assign descriptive model-free exchangeability mode
                if (total_nonsyn_c < 0.25) {
                    entry.exchangeability_mode = (entry.effective_aa_alphabet <= 1.05) ? "Strictly Conserved" : "Near-Invariant";
                } else if (total_nonsyn_c >= 1.8 && entry.pair_concentration >= 0.70 && entry.effective_aa_alphabet <= 2.45) {
                    std::string pair_str = std::string(1, aa_letters[best_a1]) + "<->" + std::string(1, aa_letters[best_a2]);
                    entry.exchangeability_mode = "Toggle (" + pair_str + ")";
                    if (entry.exchangeability_p_value <= settings.p_value_threshold) {
                        std::ostringstream ds;
                        ds << "Toggle " << pair_str << " (" << std::fixed << std::setprecision(0)
                           << (entry.pair_concentration * 100.0) << "% of " << std::setprecision(1)
                           << total_nonsyn_c << " subs, N_eff=" << std::setprecision(2) << entry.effective_aa_alphabet << ")";
                        entry.key_property_driver = ds.str();
                    }
                } else if (max_conv >= 3 && entry.dominant_target_aa != '-') {
                    std::string dom_s(1, entry.dominant_target_aa);
                    entry.exchangeability_mode = "Convergent (->" + dom_s + ")";
                    if (entry.exchangeability_p_value <= settings.p_value_threshold) {
                        std::ostringstream ds;
                        ds << "Convergent ->" << dom_s << " (" << max_conv
                           << " indep. lineages, M=" << std::fixed << std::setprecision(1) << total_nonsyn_c << ")";
                        entry.key_property_driver = ds.str();
                    }
                } else if (total_nonsyn_c >= 2.0 && entry.effective_aa_alphabet < 0.68 * entry.expected_neutral_n_eff) {
                    std::string pair_str = std::string(1, aa_letters[best_a1]) + "/" + std::string(1, aa_letters[best_a2]);
                    entry.exchangeability_mode = "Restricted Alphabet (" + pair_str + ")";
                    if (entry.exchangeability_p_value <= settings.p_value_threshold) {
                        std::ostringstream ds;
                        ds << "Restricted AA Set (N_eff=" << std::fixed << std::setprecision(2)
                           << entry.effective_aa_alphabet << " vs " << entry.expected_neutral_n_eff << " neutral)";
                        entry.key_property_driver = ds.str();
                    }
                } else if (total_nonsyn_c >= 2.0 && entry.effective_aa_alphabet >= 0.85 * entry.expected_neutral_n_eff) {
                    entry.exchangeability_mode = "Broadly Diversifying";
                } else {
                    entry.exchangeability_mode = "Permissive Drift";
                }
            }
        }

        // 4B. Lineage-Specific Loss-of-Constraint / Pseudogenization (Core-Invasion) Detector (Tier 1)
        {
            // Identify strongly conserved purifying core sites (FEL/dFEL purifying or ω_c <= 0.12)
            std::vector<size_t> core_sites;
            for (size_t c = 0; c < C; ++c) {
                const auto& s = res.site_atlas[c];
                Scalar perv_p = res.directional_drift_active ? s.dfel_p_value : s.fel_p_value;
                Scalar w_c = s.beta_pervasive / std::max<Scalar>(0.05, s.alpha);
                if (s.gap_fraction <= 0.25 && ((perv_p <= 0.10 && w_c <= 0.10) || (s.alpha >= 0.8 && w_c <= 0.05))) {
                    core_sites.push_back(c);
                }
            }
            // Fallback if tree is very shallow: bottom 25% of low-gap sites with positive α
            if (core_sites.size() < std::min<size_t>(8, C / 5)) {
                std::vector<size_t> cand_sites;
                for (size_t c = 0; c < C; ++c) {
                    if (res.site_atlas[c].gap_fraction <= 0.25 && res.site_atlas[c].alpha > 0.10) {
                        cand_sites.push_back(c);
                    }
                }
                std::sort(cand_sites.begin(), cand_sites.end(), [&](size_t a, size_t b) {
                    Scalar wa = res.site_atlas[a].beta_pervasive / std::max<Scalar>(0.05, res.site_atlas[a].alpha);
                    Scalar wb = res.site_atlas[b].beta_pervasive / std::max<Scalar>(0.05, res.site_atlas[b].alpha);
                    return wa < wb;
                });
                core_sites.assign(cand_sites.begin(), cand_sites.begin() + std::min<size_t>(cand_sites.size(), std::max<size_t>(5, C / 4)));
            }

            // Compute total tree length to apportion expected core nonsynonymous substitutions per branch
            Scalar total_tree_subs = 0.0;
            for (size_t bi = 0; bi < B; ++bi) total_tree_subs += std::max<Scalar>(1e-5, res.branch_profiles[bi].branch_length_subs);
            Scalar total_core_nonsyn = 0.0;
            for (size_t c : core_sites) {
                total_core_nonsyn += res.site_atlas[c].expected_nonsyn_subs;
            }

            for (size_t bi = 0; bi < B; ++bi) {
                auto& bp = res.branch_profiles[bi];
                Scalar br_frac = std::max<Scalar>(1e-5, bp.branch_length_subs) / std::max<Scalar>(1e-5, total_tree_subs);
                Scalar obs_core = 0.0;
                std::vector<std::pair<Scalar, size_t>> invaded;
                for (size_t c : core_sites) {
                    Scalar u_cb = site_branch_nonsyn[c][bi];
                    obs_core += u_cb;
                    if (u_cb >= 0.60) {
                        invaded.push_back({u_cb, c + 1});
                    }
                }
                std::sort(invaded.begin(), invaded.end(), [](const auto& a, const auto& b) { return a.first > b.first; });
                for (size_t k = 0; k < std::min<size_t>(8, invaded.size()); ++k) {
                    bp.core_invaded_sites.push_back(invaded[k].second);
                }

                // Expected core nonsynonymous substitutions under tree-wide purifying constraint
                // Floor at 0.25 to account for baseline branch rate heterogeneity and avoid small-expectation inflation
                Scalar exp_core_pur = std::max<Scalar>(0.25, total_core_nonsyn * br_frac);
                Scalar exp_core_neu = std::max<Scalar>(exp_core_pur * 3.0, bp.branch_length_subs * static_cast<Scalar>(core_sites.size()) * 0.75);

                bp.core_nonsyn_Observed = obs_core;
                bp.core_nonsyn_Expected = exp_core_pur;

                if (obs_core > exp_core_pur * 2.0 && invaded.size() >= 2 && obs_core >= 2.5) {
                    Scalar alt_rate = std::clamp<Scalar>(obs_core, exp_core_pur * 2.0, exp_core_neu);
                    Scalar log_bf = obs_core * std::log(alt_rate / exp_core_pur) - (alt_rate - exp_core_pur);
                    bp.lof_bayes_factor = std::exp(std::clamp<Scalar>(log_bf, -20.0, 20.0));
                    Scalar lrt_lof = std::max<Scalar>(0.0, 2.0 * (obs_core * std::log(obs_core / exp_core_pur) - (obs_core - exp_core_pur)));
                    bp.lof_p_value = stats::chi2_p_value(lrt_lof, 1.0);
                } else {
                    bp.lof_bayes_factor = 0.0;
                    bp.lof_p_value = 1.0;
                }
                // Bonferroni-corrected branch p-value <= 0.05 / B and BF >= 25 with >= 3 invaded core sites
                bp.lof_flagged = (bp.lof_bayes_factor >= 25.0 &&
                                  bp.lof_p_value * static_cast<Scalar>(B) <= 0.05 &&
                                  invaded.size() >= 3);
                if (bp.lof_flagged) {
                    res.branches_loss_of_constraint++;
                }
            }
        }

        // 4C. 1D Spatial Segmentation & Synonymous Constraint Islands (Tier 1)
        {
            // 1. Penalized BIC 1D Change-Point Segmentation on 5-codon smoothed site-wise log(ω_c + 0.08) + gap fraction
            size_t L_min = std::min<size_t>(15, std::max<size_t>(1, C / 4));
            std::vector<Scalar> y_raw(C, 0.0), y_vec(C, 0.0);
            for (size_t c = 0; c < C; ++c) {
                Scalar eff_a = std::max<Scalar>(0.15, std::min<Scalar>(10.0, res.site_atlas[c].alpha));
                Scalar eff_b = std::min<Scalar>(10.0, res.site_atlas[c].beta_pervasive);
                Scalar w_c = eff_b / eff_a;
                y_raw[c] = std::log(std::clamp<Scalar>(w_c, 0.0, 5.0) + 0.08) + 1.5 * res.site_atlas[c].gap_fraction;
            }
            for (size_t c = 0; c < C; ++c) {
                if (C < 20) {
                    y_vec[c] = y_raw[c];
                } else {
                    size_t w_lo = (c >= 2) ? (c - 2) : 0;
                    size_t w_hi = std::min<size_t>(C - 1, c + 2);
                    Scalar s_val = 0.0;
                    for (size_t k = w_lo; k <= w_hi; ++k) s_val += y_raw[k];
                    y_vec[c] = s_val / static_cast<Scalar>(w_hi - w_lo + 1);
                }
            }

            // Prefix sums for O(1) segment variance evaluation
            std::vector<Scalar> pref_y(C + 1, 0.0), pref_y2(C + 1, 0.0);
            for (size_t c = 0; c < C; ++c) {
                pref_y[c + 1] = pref_y[c] + y_vec[c];
                pref_y2[c + 1] = pref_y2[c] + y_vec[c] * y_vec[c];
            }
            Scalar total_mean = pref_y[C] / std::max<Scalar>(1.0, static_cast<Scalar>(C));
            Scalar total_var = std::max<Scalar>(0.05, (pref_y2[C] / std::max<Scalar>(1.0, static_cast<Scalar>(C))) - total_mean * total_mean);
            Scalar bic_pen = 4.0 * total_var * std::log(static_cast<Scalar>(std::max<size_t>(2, C)));

            auto seg_rss = [&](size_t s, size_t e) -> Scalar {
                Scalar len = static_cast<Scalar>(e - s);
                Scalar sum1 = pref_y[e] - pref_y[s];
                Scalar sum2 = pref_y2[e] - pref_y2[s];
                return std::max<Scalar>(0.0, sum2 - (sum1 * sum1) / len);
            };

            std::vector<Scalar> dp(C + 1, 1e18);
            std::vector<size_t> prev_cut(C + 1, 0);
            dp[0] = -bic_pen; // First segment has no change-point penalty
            for (size_t e = L_min; e <= C; ++e) {
                for (size_t s = 0; s + L_min <= e; ++s) {
                    if (s > 0 && s < L_min) continue;
                    Scalar cost = dp[s] + seg_rss(s, e) + bic_pen;
                    if (cost < dp[e]) {
                        dp[e] = cost;
                        prev_cut[e] = s;
                    }
                }
            }

            std::vector<std::pair<size_t, size_t>> seg_bounds;
            size_t cur_e = C;
            while (cur_e > 0) {
                size_t s = prev_cut[cur_e];
                seg_bounds.push_back({s, cur_e - 1});
                if (s == 0) break;
                cur_e = s;
            }
            std::reverse(seg_bounds.begin(), seg_bounds.end());
            if (seg_bounds.empty() && C > 0) seg_bounds.push_back({0, C - 1});

            for (const auto& bnd : seg_bounds) {
                SpatialDomainSegment dom;
                dom.start_codon = bnd.first;
                dom.end_codon = bnd.second;
                Scalar len = static_cast<Scalar>(bnd.second - bnd.first + 1);
                Scalar sum_w = 0.0, sum_a = 0.0, sum_b = 0.0, sum_g = 0.0, sum_neff = 0.0;
                for (size_t c = bnd.first; c <= bnd.second; ++c) {
                    Scalar eff_a = std::max<Scalar>(0.15, std::min<Scalar>(10.0, res.site_atlas[c].alpha));
                    Scalar eff_b = std::min<Scalar>(10.0, res.site_atlas[c].beta_pervasive);
                    Scalar w_c = eff_b / eff_a;
                    sum_w += std::min<Scalar>(3.0, w_c);
                    sum_a += eff_a;
                    sum_b += eff_b;
                    sum_g += res.site_atlas[c].gap_fraction;
                    sum_neff += res.site_atlas[c].effective_aa_alphabet;
                }
                dom.mean_omega = sum_w / len;
                dom.mean_alpha = sum_a / len;
                dom.mean_beta = sum_b / len;
                dom.mean_gap_fraction = sum_g / len;
                dom.mean_n_eff = sum_neff / len;
                if (dom.mean_gap_fraction >= 0.25 || dom.mean_omega >= 0.90) {
                    dom.domain_type = "Hypervariable / IDR";
                } else if (dom.mean_omega <= 0.28) {
                    dom.domain_type = "Purifying Core";
                } else {
                    dom.domain_type = "Permissive / Surface";
                }
                for (size_t c = bnd.first; c <= bnd.second; ++c) {
                    res.site_atlas[c].spatial_domain = dom.domain_type;
                }
                res.spatial_domains.push_back(std::move(dom));
            }

            // 2. Synonymous Constraint Island Scan (dual-coding / RNA structure / splice enhancer detection)
            //    Combines Winsorized α_c (<= 5.0) depletion Z-score with exact Binomial silent-codon run test
            if (C >= 25) {
                std::vector<Scalar> win_alpha(C, 0.0);
                Scalar sum_a = 0.0, sum_a2 = 0.0;
                size_t total_silent = 0;
                for (size_t c = 0; c < C; ++c) {
                    win_alpha[c] = std::min<Scalar>(5.0, res.site_atlas[c].alpha);
                    sum_a += win_alpha[c];
                    sum_a2 += win_alpha[c] * win_alpha[c];
                    if (win_alpha[c] <= 0.15 && res.site_atlas[c].gap_fraction <= 0.25) total_silent++;
                }
                Scalar gene_mean_a = std::max<Scalar>(0.15, sum_a / static_cast<Scalar>(C));
                Scalar gene_sd_a = std::sqrt(std::max<Scalar>(0.04, (sum_a2 / static_cast<Scalar>(C)) - gene_mean_a * gene_mean_a));
                Scalar p_silent = std::clamp<Scalar>(static_cast<Scalar>(total_silent) / static_cast<Scalar>(C), 0.05, 0.90);

                auto binom_upper_tail = [](size_t k_obs, size_t n_tot, Scalar p0) -> Scalar {
                    if (k_obs == 0) return 1.0;
                    Scalar log_p = std::log(p0), log_q = std::log(1.0 - p0);
                    Scalar prob = 0.0;
                    for (size_t k = k_obs; k <= n_tot; ++k) {
                        Scalar lchoose = std::lgamma(static_cast<Scalar>(n_tot + 1)) -
                                         std::lgamma(static_cast<Scalar>(k + 1)) -
                                         std::lgamma(static_cast<Scalar>(n_tot - k + 1));
                        prob += std::exp(lchoose + static_cast<Scalar>(k) * log_p + static_cast<Scalar>(n_tot - k) * log_q);
                    }
                    return std::clamp<Scalar>(prob, 0.0, 1.0);
                };

                size_t win_len = std::min<size_t>(11, C / 3);
                std::vector<bool> used_in_island(C, false);
                for (size_t s = 0; s + win_len <= C; ++s) {
                    if (used_in_island[s]) continue;
                    // Find best contiguous extension [s, e] up to length 24 that minimizes p_isl
                    size_t best_e = s + win_len - 1;
                    Scalar best_p = 1.0, best_z = 0.0, best_mean_a = gene_mean_a;
                    for (size_t e = s + win_len - 1; e < std::min<size_t>(C, s + 24); ++e) {
                        size_t span = e - s + 1;
                        Scalar w_sum = 0.0, g_sum = 0.0;
                        size_t k_silent = 0;
                        for (size_t c = s; c <= e; ++c) {
                            w_sum += win_alpha[c];
                            g_sum += res.site_atlas[c].gap_fraction;
                            if (win_alpha[c] <= 0.15) k_silent++;
                        }
                        Scalar w_mean = w_sum / static_cast<Scalar>(span);
                        Scalar g_mean = g_sum / static_cast<Scalar>(span);
                        if (g_mean > 0.20 || w_mean > 0.35 * gene_mean_a) continue;

                        Scalar z_dep = (gene_mean_a - w_mean) / (gene_sd_a / std::sqrt(static_cast<Scalar>(span)));
                        Scalar p_norm = 0.5 * std::erfc(z_dep / std::sqrt(2.0));
                        Scalar p_bin = binom_upper_tail(k_silent, span, p_silent);
                        Scalar p_cand = std::min(p_norm, p_bin);
                        if (p_cand < best_p) {
                            best_p = p_cand;
                            best_e = e;
                            best_mean_a = w_mean;
                            best_z = std::max(z_dep, std::sqrt(std::max<Scalar>(0.0, -2.0 * std::log(std::max<Scalar>(1e-16, p_cand)))));
                        }
                    }
                    Scalar adj_p = std::min<Scalar>(1.0, best_p * static_cast<Scalar>(std::max<size_t>(1, C / 12)));
                    if (best_p <= 0.010 && adj_p <= 0.12) {
                        SynonymousConstraintIsland isl;
                        isl.start_codon = s;
                        isl.end_codon = best_e;
                        isl.mean_alpha = best_mean_a;
                        isl.gene_mean_alpha = gene_mean_a;
                        isl.suppression_ratio = best_mean_a / gene_mean_a;
                        isl.z_score = best_z;
                        isl.p_value = adj_p;
                        res.synonymous_islands.push_back(isl);
                        for (size_t k = s; k <= best_e; ++k) {
                            used_in_island[k] = true;
                            res.site_atlas[k].in_synonymous_island = true;
                        }
                        s = best_e;
                    }
                }
            }
        }

        // ---------------------------------------------------------------------
        // Decision 6: Multiple Testing Control (BH q-values per test family +
        //             Hierarchical Simes Omnibus q-value + Priority Hierarchy)
        // ---------------------------------------------------------------------
        {
            std::vector<Scalar> p_fel(C), p_dfel(C), p_meme(C), p_cfel(C), p_exch(C), p_simes(C);
            for (size_t c = 0; c < C; ++c) {
                const auto& s = res.site_atlas[c];
                p_fel[c] = s.fel_p_value;
                p_dfel[c] = s.dfel_p_value;
                p_meme[c] = s.meme_p_value;
                p_cfel[c] = s.contrast_fel_p_value;
                p_exch[c] = s.exchangeability_p_value;
                Scalar active_perv_p = res.directional_drift_active ? s.dfel_p_value : s.fel_p_value;
                if (res.clade_macro_regime.triggered) {
                    p_simes[c] = compute_simes_p_value({active_perv_p, s.meme_p_value, s.contrast_fel_p_value});
                } else {
                    p_simes[c] = compute_simes_p_value({active_perv_p, s.meme_p_value});
                }
                res.site_atlas[c].simes_omnibus_p = p_simes[c];
            }

            auto q_fel = compute_bh_q_values(p_fel);
            auto q_dfel = compute_bh_q_values(p_dfel);
            auto q_meme = compute_bh_q_values(p_meme);
            auto q_cfel = compute_bh_q_values(p_cfel);
            auto q_exch = compute_bh_q_values(p_exch);
            auto q_simes = compute_bh_q_values(p_simes);

            for (size_t c = 0; c < C; ++c) {
                auto& entry = res.site_atlas[c];
                entry.fel_bh_q = q_fel[c];
                entry.dfel_bh_q = q_dfel[c];
                entry.meme_bh_q = q_meme[c];
                entry.contrast_fel_bh_q = q_cfel[c];
                entry.exchangeability_bh_q = q_exch[c];
                entry.simes_omnibus_q = q_simes[c];
                if (entry.simes_omnibus_q <= settings.fdr_threshold) {
                    res.sites_simes_fdr_significant++;
                }

                Scalar perv_p = res.directional_drift_active ? entry.dfel_p_value : entry.fel_p_value;
                bool is_perv_pos = (perv_p <= settings.p_value_threshold && entry.beta_pervasive > entry.alpha);
                bool is_perv_neg = (perv_p <= settings.p_value_threshold && entry.beta_pervasive < entry.alpha);
                bool is_epi_pos = (entry.meme_p_value <= settings.p_value_threshold && entry.beta_plus > entry.alpha);
                bool is_cfel = res.clade_macro_regime.triggered && (entry.contrast_fel_p_value <= 0.05);
                bool is_exch = !entry.key_property_driver.empty();

                if (entry.simes_omnibus_q <= settings.fdr_threshold) {
                    if (is_perv_pos || is_epi_pos) {
                        res.sites_simes_diversifying++;
                    } else if (is_perv_neg) {
                        res.sites_simes_purifying++;
                    }
                }

                if (is_perv_pos && is_epi_pos) {
                    entry.classification = "Pervasive + Episodic Diversifying";
                    res.sites_pervasive_diversifying++;
                    res.sites_episodic_diversifying++;
                } else if (is_epi_pos) {
                    entry.classification = "Episodic Diversifying";
                    res.sites_episodic_diversifying++;
                } else if (is_perv_pos) {
                    entry.classification = "Pervasive Diversifying";
                    res.sites_pervasive_diversifying++;
                } else if (is_cfel) {
                    entry.classification = "Regime-Specific Shift";
                    res.sites_regime_differential++;
                } else if (is_perv_neg) {
                    entry.classification = is_exch ? "Purifying (Toggle/Convergent)" : "Pervasive Purifying";
                    res.sites_pervasive_purifying++;
                } else if (is_exch) {
                    entry.classification = "Toggle / Convergent";
                } else if (entry.beta_pervasive <= 1e-4 && entry.alpha <= 1e-4) {
                    entry.classification = "Invariant";
                } else {
                    entry.classification = "Neutral";
                }
                if (is_exch) res.sites_property_constrained_or_shifting++;
                if (is_cfel && entry.classification != "Regime-Specific Shift") {
                    res.sites_regime_differential++;
                }
            }
        }

        // 4D. Tier 2 Triggered CpG / APOBEC Mutational Context Forensics & Tier 3 Opt-In Coevolution
        {
            size_t total_cpg_sites = 0;
            size_t selected_sites_count = 0;
            size_t selected_cpg_apobec_count = 0;

            for (size_t c = 0; c < C; ++c) {
                auto& entry = res.site_atlas[c];
                size_t cpg_taxa = 0, apobec_taxa = 0, valid_t = 0;
                for (size_t t = 0; t < N; ++t) {
                    int8_t cod = aln.sequences[t].codons[c];
                    if (cod < 0 || cod >= S) continue;
                    if (3 * c + 3 > aln.sequences[t].raw_seq.size()) continue;
                    valid_t++;
                    std::string triplet = aln.sequences[t].raw_seq.substr(3 * c, 3);
                    for (char& ch : triplet) ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
                    bool has_cpg = (triplet.substr(0, 2) == "CG" || triplet.substr(1, 2) == "CG");
                    if (!has_cpg && c + 1 < C && 3 * (c + 1) + 1 <= aln.sequences[t].raw_seq.size()) {
                        char next_b0 = static_cast<char>(std::toupper(static_cast<unsigned char>(aln.sequences[t].raw_seq[3 * (c + 1)])));
                        if (triplet[2] == 'C' && next_b0 == 'G') {
                            has_cpg = true;
                        }
                    }
                    if (has_cpg) cpg_taxa++;
                    if (triplet.substr(0, 2) == "TC" || triplet.substr(1, 2) == "TC" ||
                        triplet.substr(0, 2) == "GA" || triplet.substr(1, 2) == "GA") {
                        apobec_taxa++;
                    }
                }
                bool is_cpg_hotspot = (valid_t > 0 && static_cast<Scalar>(cpg_taxa) / static_cast<Scalar>(valid_t) >= 0.20);
                bool is_apobec_hotspot = (valid_t > 0 && static_cast<Scalar>(apobec_taxa) / static_cast<Scalar>(valid_t) >= 0.35);
                if (is_cpg_hotspot) total_cpg_sites++;

                bool is_pos_selected = (entry.classification.find("Diversifying") != std::string::npos);
                if (is_pos_selected) {
                    selected_sites_count++;
                    if (is_cpg_hotspot) {
                        entry.cpg_hypermutable_context = true;
                        selected_cpg_apobec_count++;
                        res.mutational_context_alert.cpg_driven_selected_sites++;
                        res.mutational_context_alert.flagged_cpg_sites.push_back(c + 1);
                    } else if (is_apobec_hotspot && entry.pair_concentration >= 0.65) {
                        entry.apobec_motif_context = true;
                        selected_cpg_apobec_count++;
                        res.mutational_context_alert.apobec_driven_selected_sites++;
                        res.mutational_context_alert.flagged_apobec_sites.push_back(c + 1);
                    }
                }
            }

            res.mutational_context_alert.num_selected_sites_evaluated = selected_sites_count;
            if (selected_sites_count > 0 && selected_cpg_apobec_count >= 2) {
                Scalar bg_p = std::clamp<Scalar>(static_cast<Scalar>(total_cpg_sites) / static_cast<Scalar>(std::max<size_t>(1, C)), 0.05, 0.85);
                Scalar obs_frac = static_cast<Scalar>(selected_cpg_apobec_count) / static_cast<Scalar>(selected_sites_count);
                Scalar n_sel = static_cast<Scalar>(selected_sites_count);
                Scalar k_sel = static_cast<Scalar>(selected_cpg_apobec_count);
                Scalar enrich = obs_frac / bg_p;
                Scalar g_stat = 0.0;
                if (obs_frac > bg_p) {
                    g_stat = 2.0 * (k_sel * std::log(obs_frac / bg_p)
                           + std::max<Scalar>(0.0, n_sel - k_sel) * std::log(std::max<Scalar>(1e-9, (1.0 - obs_frac) / (1.0 - bg_p))));
                }
                Scalar p_enr = (g_stat <= 1e-6) ? 1.0 : stats::chi2_p_value(g_stat, 1.0);
                res.mutational_context_alert.cpg_enrichment_fold = enrich;
                res.mutational_context_alert.cpg_fisher_p = p_enr;
                if (p_enr <= 0.05 || (selected_cpg_apobec_count >= 2 && obs_frac >= 0.50)) {
                    res.mutational_context_alert.triggered = true;
                    std::ostringstream msg;
                    msg << selected_cpg_apobec_count << " of " << selected_sites_count
                        << " positively selected sites (" << std::fixed << std::setprecision(0) << (obs_frac * 100.0)
                        << "%) occur at hypermutable CpG/deamination dinucleotides (" << std::setprecision(1)
                        << enrich << "x enrichment, p = " << std::setprecision(4) << p_enr
                        << "). Inspect for mutational hotspot confounding.";
                    res.mutational_context_alert.alert_summary = msg.str();
                }
            }

            // Tier 3 (Opt-In via --coevolution): Epistatic Co-Substitution Pairs
            if (settings.enable_coevolution && C >= 2 && B >= 6) {
                std::vector<size_t> active_sites;
                for (size_t c = 0; c < C; ++c) {
                    if (res.site_atlas[c].expected_nonsyn_subs >= 1.5 && res.site_atlas[c].gap_fraction <= 0.30) {
                        active_sites.push_back(c);
                    }
                }
                Vector bl_vec(B);
                for (size_t bi = 0; bi < B; ++bi) bl_vec(bi) = std::max<Scalar>(1e-5, res.branch_profiles[bi].branch_length_subs);
                Scalar bl_norm2 = std::max<Scalar>(1e-12, bl_vec.squaredNorm());

                std::vector<Vector> centered_u(active_sites.size(), Vector::Zero(B));
                std::vector<Scalar> u_norms(active_sites.size(), 0.0);
                for (size_t idx = 0; idx < active_sites.size(); ++idx) {
                    size_t c = active_sites[idx];
                    Vector u(B);
                    for (size_t bi = 0; bi < B; ++bi) u(bi) = site_branch_nonsyn[c][bi];
                    Vector u_proj = u - (u.dot(bl_vec) / bl_norm2) * bl_vec;
                    u_proj.array() -= u_proj.mean();
                    centered_u[idx] = u_proj;
                    u_norms[idx] = u_proj.norm();
                }

                std::vector<CoevolutionSitePair> raw_pairs;
                std::vector<Scalar> pair_p_vals;
                for (size_t i = 0; i < active_sites.size(); ++i) {
                    if (u_norms[i] <= 1e-6) continue;
                    for (size_t j = i + 1; j < active_sites.size(); ++j) {
                        if (u_norms[j] <= 1e-6) continue;
                        size_t c1 = active_sites[i], c2 = active_sites[j];
                        if (c2 <= c1 + 2) continue; // Skip immediate neighbor MNMs
                        Scalar r_cos = centered_u[i].dot(centered_u[j]) / (u_norms[i] * u_norms[j]);
                        if (r_cos <= 0.45) continue;

                        Scalar shared_br = 0.0;
                        std::string top_br_str;
                        for (size_t bi = 0; bi < B; ++bi) {
                            if (site_branch_nonsyn[c1][bi] >= 0.40 && site_branch_nonsyn[c2][bi] >= 0.40) {
                                shared_br += std::min(site_branch_nonsyn[c1][bi], site_branch_nonsyn[c2][bi]);
                                if (top_br_str.size() < 40) {
                                    if (!top_br_str.empty()) top_br_str += ",";
                                    top_br_str += res.branch_profiles[bi].branch_name;
                                }
                            }
                        }
                        if (shared_br < 1.5) continue;

                        Scalar df_b = std::max<Scalar>(2.0, static_cast<Scalar>(B) - 3.0);
                        Scalar t_stat = r_cos * std::sqrt(df_b / std::max<Scalar>(1e-6, 1.0 - r_cos * r_cos));
                        Scalar p_val = stats::chi2_p_value(t_stat * t_stat, 1.0);
                        if (p_val <= 0.05) {
                            CoevolutionSitePair cp;
                            cp.site_i = c1;
                            cp.site_j = c2;
                            cp.shared_branches = shared_br;
                            cp.expected_shared = (res.site_atlas[c1].expected_nonsyn_subs * res.site_atlas[c2].expected_nonsyn_subs) / static_cast<Scalar>(B);
                            cp.cosine_similarity = r_cos;
                            cp.p_value = p_val;
                            cp.top_co_mutating_branches = top_br_str;
                            raw_pairs.push_back(cp);
                            pair_p_vals.push_back(p_val);
                        }
                    }
                }
                if (!raw_pairs.empty()) {
                    auto q_pairs = compute_bh_q_values(pair_p_vals);
                    for (size_t k = 0; k < raw_pairs.size(); ++k) {
                        raw_pairs[k].bh_q_value = q_pairs[k];
                    }
                    std::sort(raw_pairs.begin(), raw_pairs.end(), [](const auto& a, const auto& b) {
                        return a.p_value < b.p_value;
                    });
                    for (size_t k = 0; k < std::min<size_t>(25, raw_pairs.size()); ++k) {
                        res.coevolution_pairs.push_back(raw_pairs[k]);
                    }
                }
            }
        }

        for (const auto& bp : res.branch_profiles) {
            if (bp.absrel_p_value <= 0.05) res.branches_episodic_diversifying++;
            if (bp.relax_call == "Intensified") res.branches_intensified++;
            if (bp.relax_call == "Relaxed") res.branches_relaxed++;
        }

        if (!settings.quiet) {
            std::ostringstream ss;
            ss << "      Simes FDR (q<=" << std::fixed << std::setprecision(2) << settings.fdr_threshold
               << "): " << res.sites_simes_fdr_significant
               << "  │  Episodic: " << res.sites_episodic_diversifying
               << "  │  Pervasive: " << res.sites_pervasive_diversifying
               << "  │  Purifying: " << res.sites_pervasive_purifying
               << "  │  Toggle/Conv: " << res.sites_property_constrained_or_shifting
               << "  │  Domains: " << res.spatial_domains.size()
               << "  │  LoF Branches: " << res.branches_loss_of_constraint << "\n\n";
            std::cout << Console::muted(ss.str()) << std::flush;
        }

        return res;
    }
};

} // namespace hyphy::analyses

