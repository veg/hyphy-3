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
    Scalar busted_lrt = 0.0;
    Scalar busted_p_value = 1.0;
    bool has_gene_wide_positive_selection = false;
};

struct CharacterizeResult {
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
        std::shared_ptr<const GeneticCode> gcode = nullptr
    ) {
        auto t0 = std::chrono::steady_clock::now();
        if (!gcode) gcode = GeneticCode::universal();

        JointAlignSettings jas;
        jas.p_value_threshold = settings.p_value_threshold;
        jas.enable_coupled_topology_search = settings.enable_coupled_spr;
        jas.run_marginalized_fel = true;
        jas.quiet = settings.quiet;

        JointAlignEngine jae(jas);
        auto ja_res = jae.run_from_unaligned(unaligned_seqs, gcode);

        auto res = run_unified_tensor_characterization(std::move(ja_res));
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
        auto t1 = std::chrono::steady_clock::now();
        res.runtime_seconds = std::chrono::duration<double>(t1 - t0).count();
        return res;
    }

    static std::string to_json(const CharacterizeResult& res) {
        using json = nlohmann::json;
        json j;
        j["analysis"] = {
            {"id", "characterize"},
            {"name", "HyPhy 3 One-Button Evolutionary Characterization"},
            {"version", "3.0.0"},
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
        // 1. Legacy FEL JSON
        {
            json jf;
            jf["analysis"]["info"] = "HyPhy 3 Characterize Legacy-Compatible FEL Export";
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
            jm["analysis"]["info"] = "HyPhy 3 Characterize Legacy-Compatible MEME Export";
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
            jb["analysis"]["info"] = "HyPhy 3 Characterize Legacy-Compatible BUSTED-E Export";
            jb["test results"] = {
                {"LRT", res.gene_wide.busted_lrt},
                {"p-value", res.gene_wide.busted_p_value}
            };
            std::ofstream out(base_prefix + ".BUSTED.json");
            out << jb.dump(2);
        }
    }

    // Decision 7: Unified, self-contained interactive multi-tab SPA (.CHARACTERIZE.html)
    static std::string render_standalone_html(const CharacterizeResult& res, const std::string& dataset_label) {
        std::ostringstream html;
        html << std::fixed << std::setprecision(4);
        html << R"HTML(<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>HyPhy 3 Unified Characterization — )HTML" << dataset_label << R"HTML(</title>
<style>
  :root {
    --bg: #f6f6f2; --card: #ffffff; --ink: #191917; --muted: #666660;
    --border: #d8d8d0; --accent: #eb5e28; --pos: #c92a2a; --neg: #1864ab;
    --epi: #6741d9; --prop: #087f5b; --diff: #d9480f;
  }
  * { box-sizing: border-box; margin: 0; padding: 0; }
  body {
    background: var(--bg); color: var(--ink);
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
    line-height: 1.45; padding: 24px;
  }
  .container { max-width: 1420px; margin: 0 auto; }
  header {
    background: var(--ink); color: #fff; padding: 20px 26px; border-radius: 8px;
    display: flex; justify-content: space-between; align-items: center; margin-bottom: 16px;
  }
  header h1 { font-size: 19px; letter-spacing: 0.5px; }
  header .sub { color: #b8b8b0; font-size: 13px; margin-top: 4px; }
  .badge {
    background: var(--accent); color: #fff; font-weight: 700; font-size: 11px;
    padding: 4px 10px; border-radius: 4px; letter-spacing: 0.8px; text-transform: uppercase;
  }
  .tabs {
    display: flex; gap: 8px; margin-bottom: 18px; border-bottom: 2px solid var(--border); padding-bottom: 8px;
  }
  .tab-btn {
    background: var(--card); border: 1px solid var(--border); border-radius: 6px;
    padding: 8px 16px; font-size: 13px; font-weight: 600; color: var(--muted); cursor: pointer;
  }
  .tab-btn.active {
    background: var(--ink); color: #fff; border-color: var(--ink);
  }
  .tab-panel { display: none; }
  .tab-panel.active { display: block; }
  .grid-4 { display: grid; grid-template-columns: repeat(4, 1fr); gap: 14px; margin-bottom: 18px; }
  .grid-2 { display: grid; grid-template-columns: repeat(2, 1fr); gap: 14px; margin-bottom: 18px; }
  .card {
    background: var(--card); border: 1px solid var(--border); border-radius: 8px; padding: 16px 18px;
  }
  .card h3 { font-size: 11.5px; text-transform: uppercase; letter-spacing: 0.8px; color: var(--muted); margin-bottom: 6px; }
  .card .metric { font-size: 22px; font-weight: 700; color: var(--ink); }
  .card .detail { font-size: 12px; color: var(--muted); margin-top: 5px; }
  .toolbar {
    display: flex; gap: 12px; align-items: center; margin-bottom: 12px; flex-wrap: wrap;
  }
  .toolbar input, .toolbar select {
    padding: 6px 10px; border: 1px solid var(--border); border-radius: 5px; font-size: 12.5px; background: #fff;
  }
  table { width: 100%; border-collapse: collapse; font-size: 12px; }
  th, td { padding: 7px 9px; border-bottom: 1px solid var(--border); text-align: left; }
  th { background: #efefea; font-weight: 600; font-size: 11px; text-transform: uppercase; color: var(--muted); position: sticky; top: 0; }
  td.num, th.num { text-align: right; font-variant-numeric: tabular-nums; }
  .pill {
    display: inline-block; padding: 2px 7px; border-radius: 4px; font-size: 11px; font-weight: 600;
  }
  .pill-pos { background: #ffe3e3; color: var(--pos); }
  .pill-epi { background: #f3f0ff; color: var(--epi); }
  .pill-diff { background: #fff4e6; color: var(--diff); }
  .pill-neg { background: #e7f5ff; color: var(--neg); }
  .pill-prop { background: #e6fcf5; color: var(--prop); }
  .pill-neu { background: #f1f3f5; color: var(--muted); }
  .section-title { font-size: 14.5px; font-weight: 700; margin: 16px 0 8px; display: flex; align-items: center; gap: 8px; }
  .track-bar { display: flex; width: 100%; height: 22px; border-radius: 4px; overflow: hidden; border: 1px solid var(--border); margin: 8px 0; }
  .track-seg { height: 100%; }
  .heat-cell { text-align: center; font-variant-numeric: tabular-nums; font-weight: 600; border-radius: 3px; }
</style>
</head>
<body>
<div class="container">
  <header>
    <div>
      <h1>HYPHY 3 &bull; ONE-BUTTON EVOLUTIONARY CHARACTERIZATION</h1>
      <div class="sub">Dataset: )HTML" << dataset_label << " &bull; "
             << res.joint_align.map_alignment.num_taxa << " taxa &times; "
             << res.joint_align.map_alignment.num_codons << " codons &bull; Completed in "
             << std::setprecision(2) << res.runtime_seconds << R"HTML( s</div>
    </div>
    <span class="badge">Hierarchical 3-Question Forensics</span>
  </header>

  <div class="grid-4">
    <div class="card">
      <h3>Q1. Integrity &amp; Process Audit</h3>
      <div class="metric">PP = )HTML" << std::setprecision(4) << res.alignment_objective.mean_column_posterior << R"HTML(</div>
      <div class="detail">)HTML" << (res.composition_test.is_heterogeneous ? "3-Regime GC3 Drift" : "Stationary")
             << " (CP3 p=" << std::setprecision(3) << res.composition_test.p_value_cp3
             << ") &bull; LoF Lineages: " << res.branches_loss_of_constraint << R"HTML(</div>
    </div>
    <div class="card">
      <h3>Q2. Gene-Wide Selection (BUSTED-E)</h3>
      <div class="metric">)HTML" << (res.gene_wide.has_gene_wide_positive_selection ? "POSITIVE" : "NO SIGNAL")
             << " (p = " << std::setprecision(4) << res.gene_wide.busted_p_value << R"HTML()</div>
      <div class="detail">&omega;<sub>3</sub> = )HTML" << std::setprecision(2) << res.gene_wide.omega_classes[2]
             << " (" << std::setprecision(1) << (res.gene_wide.omega_weights[2] * 100.0)
             << "%) &bull; Domains: " << res.spatial_domains.size()
             << " &bull; Syn Islands: " << res.synonymous_islands.size() << R"HTML(</div>
    </div>
    <div class="card">
      <h3>Q2. Clade &amp; Lineage Architecture</h3>
      <div class="metric">)HTML" << (res.clade_macro_regime.triggered
             ? ("RELAX K = " + std::to_string(res.clade_macro_regime.clade_relax_k).substr(0, 4) + " (" + res.clade_macro_regime.clade_relax_call + ")")
             : "Homogeneous Tree") << R"HTML(</div>
      <div class="detail">aBSREL Pos: )HTML" << res.branches_episodic_diversifying
             << " br &bull; Contrast-FEL: " << res.clade_macro_regime.contrast_fel_differing_sites << R"HTML( sites</div>
    </div>
    <div class="card">
      <h3>Q3. Site &amp; Mechanism Resolution</h3>
      <div class="metric">)HTML" << res.sites_episodic_diversifying << " Epi / "
             << res.sites_pervasive_diversifying << " Perv / "
             << res.sites_pervasive_purifying << R"HTML( Pur</div>
      <div class="detail">Simes FDR q&le;0.10: )HTML" << res.sites_simes_fdr_significant
             << " &bull; Toggles/Conv: " << res.sites_property_constrained_or_shifting << R"HTML(</div>
    </div>
  </div>

  <div class="tabs">
    <button class="tab-btn active" onclick="showTab('tab-overview', this)">1. Executive &amp; 1D Domain Architecture</button>
    <button class="tab-btn" onclick="showTab('tab-sites', this)">2. Site Atlas (FEL / dFEL / MEME / Exchangeability)</button>
    <button class="tab-btn" onclick="showTab('tab-branches', this)">3. Lineage Regimes &amp; LoF Audit (aBSREL / RELAX)</button>
    <button class="tab-btn" onclick="showTab('tab-prime', this)">4. Exchangeability, Toggles &amp; Convergence</button>
  </div>

  <!-- TAB 1: EXECUTIVE NARRATIVE & SYNCHRONIZED MANHATTAN / POA CONFIDENCE TRACK -->
  <div id="tab-overview" class="tab-panel active">
    <div class="card" style="margin-bottom:16px;">
      <h3>Automated 3-Question Evolutionary Narrative</h3>
      <p style="font-size:13.5px; margin-top:6px;">
        <b>1. Data &amp; Process Integrity:</b> Co-estimation across <b>)HTML" << res.joint_align.map_alignment.num_taxa
        << " taxa</b> and <b>" << res.joint_align.map_alignment.num_codons
        << " codons</b> yielded mean column homology posterior <b>"
        << std::setprecision(4) << res.alignment_objective.mean_column_posterior
        << "</b>. Compositional stationarity at 3rd codon positions reported CP3 &chi;&sup2; = <b>"
        << std::setprecision(2) << res.composition_test.chi2_cp3 << "</b> (p = " << std::setprecision(3) << res.composition_test.p_value_cp3
        << "), activating <b>" << (res.directional_drift_active ? "3-Regime Directional Softmax Kernels" : "1-Regime Stationary Kernels")
        << "</b>. Core-invasion screening flagged <b>" << res.branches_loss_of_constraint
        << " lineage(s)</b> for potential loss-of-constraint / pseudogenization.<br>"
        << "<b>2. Gene &amp; Domain Architecture:</b> Error-marginalized <b>BUSTED-E</b> "
        << (res.gene_wide.has_gene_wide_positive_selection ? "<b>detected significant gene-wide episodic diversifying selection</b>" : "did not detect gene-wide positive selection")
        << " (LRT = " << std::setprecision(2) << res.gene_wide.busted_lrt << ", p = " << std::setprecision(4) << res.gene_wide.busted_p_value
        << "). Penalized 1D segmentation partitioned the CDS into <b>" << res.spatial_domains.size()
        << " selective domain(s)</b> and <b>" << res.synonymous_islands.size()
        << " synonymous constraint island(s)</b>.<br>"
        << "<b>3. Site &amp; Mechanism Resolution:</b> Active-set polish identified <b>"
        << res.sites_episodic_diversifying << " episodic diversifying sites</b> (MEME), <b>"
        << res.sites_pervasive_diversifying << " pervasive diversifying sites</b> (FEL/dFEL), and <b>"
        << res.sites_property_constrained_or_shifting << " model-free amino-acid toggle / convergent sites</b>."
        << (res.mutational_context_alert.triggered ? (" <span style=\"color:#c92a2a;font-weight:600;\">Alert: " + res.mutational_context_alert.alert_summary + "</span>") : std::string(""))
        << R"HTML(</p>
    </div>

    <div class="card" style="margin-bottom:16px;">
      <h3>Synchronized CDS Selection, 1D Domain Segmentation &amp; POA Homology Track</h3>
      <div class="detail" style="margin-bottom:6px;">Top bar: Unified Site Classification. Middle bar: 1D Spatial Domain Architecture (Blue = Purifying Core, Gray = Permissive, Red = Hypervariable/IDR). Bottom bar: Pair-HMM Column Homology Posterior PP<sub>c</sub>.</div>
      <div class="track-bar">
)HTML";
        size_t C_tot = std::max<size_t>(1, res.site_atlas.size());
        for (const auto& s : res.site_atlas) {
            std::string col = "#e9ecef";
            if (s.classification.find("Pervasive") != std::string::npos && s.classification.find("Diversifying") != std::string::npos) col = "#c92a2a";
            else if (s.classification.find("Episodic") != std::string::npos) col = "#6741d9";
            else if (s.classification.find("Regime") != std::string::npos) col = "#d9480f";
            else if (s.classification.find("Purifying") != std::string::npos) col = "#1864ab";
            else if (s.classification.find("Toggle") != std::string::npos || s.classification.find("Convergent") != std::string::npos || s.classification.find("Exchange") != std::string::npos) col = "#087f5b";
            html << "<div class=\"track-seg\" style=\"width:" << (100.0 / C_tot) << "%;background:" << col
                 << ";\" title=\"Site " << (s.site_index + 1) << ": " << s.classification
                 << " (Simes p=" << std::setprecision(4) << s.simes_omnibus_p << ", q=" << s.simes_omnibus_q << ")\"></div>";
        }
        html << R"HTML(      </div>
      <div class="track-bar" style="height:14px;">
)HTML";
        for (const auto& s : res.site_atlas) {
            std::string col = "#74c0fc";
            if (s.spatial_domain.find("Hypervariable") != std::string::npos) col = "#ff8787";
            else if (s.spatial_domain.find("Permissive") != std::string::npos) col = "#ced4da";
            html << "<div class=\"track-seg\" style=\"width:" << (100.0 / C_tot) << "%;background:" << col
                 << ";\" title=\"Site " << (s.site_index + 1) << ": " << s.spatial_domain
                 << (s.in_synonymous_island ? " [Synonymous Constraint Island]" : "") << "\"></div>";
        }
        html << R"HTML(      </div>
      <div class="track-bar" style="height:12px;">
)HTML";
        for (const auto& s : res.site_atlas) {
            Scalar pp = std::clamp<Scalar>(s.column_homology_posterior, 0.0, 1.0);
            std::string col = (pp >= 0.99) ? "#2b8a3e" : (pp >= 0.90 ? "#f08c00" : "#e03131");
            html << "<div class=\"track-seg\" style=\"width:" << (100.0 / C_tot) << "%;background:" << col
                 << ";\" title=\"Site " << (s.site_index + 1) << ": Homology PP = " << std::setprecision(3) << pp << "\"></div>";
        }
        html << R"HTML(      </div>
    </div>
  </div>

  <!-- TAB 2: INTERACTIVE SITE EXPLORER -->
  <div id="tab-sites" class="tab-panel">
    <div class="card">
      <div class="toolbar">
        <input type="text" id="siteSearch" placeholder="Filter by site #, classification, branch, toggle/convergence..." oninput="filterSites()" style="width:320px;">
        <select id="siteClassFilter" onchange="filterSites()">
          <option value="ALL">All Site Classes ()HTML" << res.site_atlas.size() << R"HTML( codons)</option>
          <option value="NON_NEUTRAL">All Non-Neutral / Selected Sites</option>
          <option value="Diversifying">Diversifying (Episodic or Pervasive)</option>
          <option value="Regime">Clade-Differential (Contrast-FEL)</option>
          <option value="Purifying">Purifying Selection</option>
          <option value="Toggle">Toggle / Convergent / Exchange</option>
        </select>
      </div>
      <div style="max-height:620px; overflow-y:auto;">
        <table id="siteTable">
          <thead>
            <tr>
              <th>Site</th>
              <th>Unified Call</th>
              <th>1D Domain</th>
              <th class="num">Simes p (BH q)</th>
              <th class="num">&alpha;</th>
              <th class="num">&beta;<sub>perv</sub></th>
              <th class="num">FEL/dFEL p (q)</th>
              <th class="num">&beta;<sup>+</sup> (w<sup>+</sup>)</th>
              <th class="num">MEME p (q)</th>
              <th class="num">N<sub>eff</sub> (Subs)</th>
              <th>Driving Branches (EBF&ge;10)</th>
              <th>Exchangeability / Mechanism</th>
            </tr>
          </thead>
          <tbody>
)HTML";
        for (const auto& s : res.site_atlas) {
            std::string pill_cls = "pill-neu";
            if (s.classification.find("Diversifying") != std::string::npos) {
                pill_cls = (s.classification.find("Episodic") != std::string::npos && s.classification.find("Pervasive") == std::string::npos)
                         ? "pill-epi" : "pill-pos";
            } else if (s.classification.find("Regime") != std::string::npos) {
                pill_cls = "pill-diff";
            } else if (s.classification.find("Purifying") != std::string::npos) {
                pill_cls = "pill-neg";
            } else if (s.classification.find("Toggle") != std::string::npos || s.classification.find("Convergent") != std::string::npos || s.classification.find("Exchange") != std::string::npos) {
                pill_cls = "pill-prop";
            }
            std::string br_str = "-";
            if (!s.top_episodic_branches.empty()) {
                br_str.clear();
                for (size_t k = 0; k < s.top_episodic_branches.size(); ++k) {
                    if (k > 0) br_str += ", ";
                    br_str += s.top_episodic_branches[k];
                }
            }
            Scalar perv_p = res.directional_drift_active ? s.dfel_p_value : s.fel_p_value;
            Scalar perv_q = res.directional_drift_active ? s.dfel_bh_q : s.fel_bh_q;
            html << "            <tr data-class=\"" << s.classification << " " << s.exchangeability_mode << "\">"
                 << "<td class=\"num\">" << (s.site_index + 1) << "</td>"
                 << "<td><span class=\"pill " << pill_cls << "\">" << s.classification << "</span></td>"
                 << "<td>" << s.spatial_domain << (s.in_synonymous_island ? " [SynIsl]" : "") << "</td>"
                 << "<td class=\"num\">" << std::setprecision(4) << s.simes_omnibus_p << " (" << std::setprecision(3) << s.simes_omnibus_q << ")</td>"
                 << "<td class=\"num\">" << std::setprecision(2) << s.alpha << "</td>"
                 << "<td class=\"num\">" << std::setprecision(2) << s.beta_pervasive << "</td>"
                 << "<td class=\"num\">" << std::setprecision(4) << perv_p << " (" << std::setprecision(3) << perv_q << ")</td>"
                 << "<td class=\"num\">" << std::setprecision(2) << s.beta_plus << " (" << std::setprecision(2) << s.weight_plus << ")</td>"
                 << "<td class=\"num\">" << std::setprecision(4) << s.meme_p_value << " (" << std::setprecision(3) << s.meme_bh_q << ")</td>"
                 << "<td class=\"num\">" << std::setprecision(2) << s.effective_aa_alphabet << " (" << std::setprecision(1) << s.expected_nonsyn_subs << ")</td>"
                 << "<td>" << br_str << "</td>"
                 << "<td>" << (s.key_property_driver.empty() ? "-" : s.key_property_driver) << "</td>"
                 << "</tr>\n";
        }
        html << R"HTML(          </tbody>
        </table>
      </div>
    </div>
  </div>

  <!-- TAB 3: LINEAGE & MACRO-REGIME EXPLORER -->
  <div id="tab-branches" class="tab-panel">
    <div class="card" style="margin-bottom:16px;">
      <h3>Macro-Regime Partition ()HTML" << res.clade_macro_regime.partition_source << R"HTML() &amp; Extended Newick</h3>
      <div class="detail" style="margin-bottom:8px;">
        Test Clade (<b>)HTML" << res.clade_macro_regime.test_regime_label << "</b>, " << res.clade_macro_regime.num_test_branches
        << " branches) vs Reference (<b>" << res.clade_macro_regime.reference_regime_label << "</b>, "
        << res.clade_macro_regime.num_reference_branches << " branches): Clade RELAX K = <b>"
        << std::setprecision(2) << res.clade_macro_regime.clade_relax_k << "</b> (LRT = "
        << res.clade_macro_regime.clade_relax_lrt << ", &Delta;BIC = " << res.clade_macro_regime.clade_relax_delta_bic
        << ", p = " << std::setprecision(4) << res.clade_macro_regime.clade_relax_p_value
        << " &rarr; <b>" << res.clade_macro_regime.clade_relax_call << R"HTML(</b>).</div>
      <pre style="background:#f1f3f5; padding:10px; border-radius:6px; font-size:11.5px; overflow-x:auto;">)HTML"
        << res.joint_align.extended_newick << R"HTML(</pre>
    </div>
    <div class="card">
      <table>
        <thead>
          <tr>
            <th>Branch</th>
            <th>Macro-Regime</th>
            <th class="num">Length</th>
            <th>GC3 Drift Regime</th>
            <th class="num">&omega;<sup>+</sup> (aBSREL)</th>
            <th class="num">aBSREL p (Holm / q)</th>
            <th class="num">RELAX K<sub>b</sub></th>
            <th class="num">RELAX p (q)</th>
            <th>Regime Call</th>
            <th class="num">Core Subs (Exp)</th>
            <th class="num">LoF BF (p)</th>
          </tr>
        </thead>
        <tbody>
)HTML";
        std::vector<BranchEvolutionaryProfile> sorted_br = res.branch_profiles;
        std::sort(sorted_br.begin(), sorted_br.end(), [](const auto& a, const auto& b) {
            return std::min({a.absrel_p_value, a.relax_p_value, a.lof_p_value}) <
                   std::min({b.absrel_p_value, b.relax_p_value, b.lof_p_value});
        });
        for (const auto& bp : sorted_br) {
            html << "          <tr>"
                 << "<td>" << bp.branch_name << (bp.lof_flagged ? " <span class=\"pill pill-pos\">LoF Alert</span>" : "") << "</td>"
                 << "<td>" << bp.macro_regime << "</td>"
                 << "<td class=\"num\">" << std::setprecision(3) << bp.branch_length_subs << "</td>"
                 << "<td>" << bp.compositional_regime << " (" << std::setprecision(1) << (bp.gc_cp3 * 100.0) << "%)</td>"
                 << "<td class=\"num\">" << std::setprecision(2) << bp.omega_plus << " (" << std::setprecision(0) << (bp.weight_plus * 100.0) << "%)</td>"
                 << "<td class=\"num\">" << std::setprecision(4) << bp.absrel_p_value
                 << " (" << std::setprecision(3) << bp.absrel_holm_p << " / " << bp.absrel_bh_q << ")</td>"
                 << "<td class=\"num\">" << std::setprecision(2) << bp.relax_k << "</td>"
                 << "<td class=\"num\">" << std::setprecision(4) << bp.relax_p_value << " (" << std::setprecision(3) << bp.relax_bh_q << ")</td>"
                 << "<td>" << bp.relax_call << "</td>"
                 << "<td class=\"num\">" << std::setprecision(1) << bp.core_nonsyn_Observed << " (" << std::setprecision(2) << bp.core_nonsyn_Expected << ")</td>"
                 << "<td class=\"num\">" << std::setprecision(1) << bp.lof_bayes_factor << " (" << std::setprecision(4) << bp.lof_p_value << ")</td>"
                 << "</tr>\n";
        }
        html << R"HTML(        </tbody>
      </table>
    </div>
  </div>

  <!-- TAB 4: MODEL-FREE EXCHANGEABILITY, TOGGLES & CONVERGENCE -->
  <div id="tab-prime" class="tab-panel">
    <div class="card">
      <h3>Model-Free Amino-Acid Exchangeability, Toggle &amp; Convergence Atlas (Replaces Fixed Property Scales)</h3>
      <div style="max-height:600px; overflow-y:auto; margin-top:10px;">
        <table>
          <thead>
            <tr>
              <th>Site</th>
              <th>Mode</th>
              <th class="num">Nonsyn Subs (M<sub>s</sub>)</th>
              <th class="num">Observed N<sub>eff</sub></th>
              <th class="num">Neutral E[N<sub>eff</sub>]</th>
              <th class="num">Pair Concentration</th>
              <th class="num">Max Convergent Branches</th>
              <th class="num">p-value (BH q)</th>
              <th>Mechanism Summary</th>
            </tr>
          </thead>
          <tbody>
)HTML";
        std::vector<UnifiedSiteAtlasEntry> ex_sites = res.site_atlas;
        std::sort(ex_sites.begin(), ex_sites.end(), [](const auto& a, const auto& b) {
            return a.exchangeability_p_value < b.exchangeability_p_value;
        });
        size_t max_pr = std::min<size_t>(ex_sites.size(), 60);
        for (size_t i = 0; i < max_pr; ++i) {
            const auto& s = ex_sites[i];
            if (s.expected_nonsyn_subs < 0.5) continue;
            html << "            <tr>"
                 << "<td class=\"num\">" << (s.site_index + 1) << "</td>"
                 << "<td>" << s.exchangeability_mode << "</td>"
                 << "<td class=\"num\">" << std::setprecision(2) << s.expected_nonsyn_subs << "</td>"
                 << "<td class=\"num\">" << std::setprecision(2) << s.effective_aa_alphabet << "</td>"
                 << "<td class=\"num\">" << std::setprecision(2) << s.expected_neutral_n_eff << "</td>"
                 << "<td class=\"num\">" << std::setprecision(1) << (s.pair_concentration * 100.0) << "%</td>"
                 << "<td class=\"num\">" << s.max_convergent_transitions << " (" << s.dominant_target_aa << ")</td>"
                 << "<td class=\"num\">" << std::setprecision(4) << s.exchangeability_p_value << " (" << std::setprecision(3) << s.exchangeability_bh_q << ")</td>"
                 << "<td>" << (s.key_property_driver.empty() ? "-" : s.key_property_driver) << "</td></tr>\n";
        }
        html << R"HTML(          </tbody>
        </table>
      </div>
    </div>
  </div>
</div>
<script>
function showTab(tabId, btn) {
  document.querySelectorAll('.tab-panel').forEach(el => el.classList.remove('active'));
  document.querySelectorAll('.tab-btn').forEach(el => el.classList.remove('active'));
  document.getElementById(tabId).classList.add('active');
  btn.classList.add('active');
}
function filterSites() {
  const q = document.getElementById('siteSearch').value.toLowerCase();
  const cls = document.getElementById('siteClassFilter').value;
  document.querySelectorAll('#siteTable tbody tr').forEach(tr => {
    const text = tr.innerText.toLowerCase();
    const rowCls = tr.getAttribute('data-class') || '';
    let matchCls = (cls === 'ALL') ||
                   (cls === 'NON_NEUTRAL' && rowCls.indexOf('Neutral') === -1 && rowCls.indexOf('Invariant') === -1) ||
                   (rowCls.indexOf(cls) !== -1);
    tr.style.display = (matchCls && text.indexOf(q) !== -1) ? '' : 'none';
  });
}
</script>
</body>
</html>
)HTML";
        return html.str();
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
                "Phase I: Homology, Indel & Reticulation Audit",
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
        };
        std::vector<RoutedSegmentTree> routed_spr_trees;
        std::vector<int> site_routed_tree_idx(C, -1);

        if (!res.joint_align.ohmic_network.bridges.empty() && C >= 10) {
            std::vector<Matrix> P_base(work_tree.num_nodes());
            for (const auto& nd : work_tree.nodes) {
                if (nd.id != work_tree.root_id && nd.parent_id != INVALID_INDEX) {
                    P_base[nd.id] = mg94_base.transition_matrix(nd.branch_length * conv);
                }
            }
            std::vector<Scalar> base_site_ll(C, 0.0);
            std::vector<Vector> node_L(work_tree.num_nodes(), Vector::Zero(S));
            std::vector<double> node_sc(work_tree.num_nodes(), 0.0);
            for (size_t c = 0; c < C; ++c) {
                SitePattern pat;
                pat.states.resize(N);
                for (size_t t = 0; t < N; ++t) pat.states[t] = aln.sequences[t].codons[c];
                pat.weight = 1;
                base_site_ll[c] = FELAnalyzer::compute_pattern_log_likelihood(
                    work_tree, pat, leaf_to_taxon, mg94_base, P_base, node_L, node_sc
                );
            }

            for (const auto& br : res.joint_align.ohmic_network.bridges) {
                int32_t u_id = br.hub_node_id;
                int32_t v_id = br.bridge_node_id;
                auto spr_opt = TreeEngine::apply_spr_move(work_tree, u_id, v_id);
                if (!spr_opt.has_value()) continue;
                const Tree& br_tree = spr_opt->tree;
                std::vector<Matrix> P_br(br_tree.num_nodes());
                for (const auto& nd : br_tree.nodes) {
                    if (nd.id != br_tree.root_id && nd.parent_id != INVALID_INDEX) {
                        P_br[nd.id] = mg94_base.transition_matrix(nd.branch_length * conv);
                    }
                }
                auto l2t_br = LikelihoodEngine::build_leaf_to_taxon_map(br_tree, aln);
                std::vector<Vector> nL_br(br_tree.num_nodes(), Vector::Zero(S));
                std::vector<double> nS_br(br_tree.num_nodes(), 0.0);

                Scalar best_sum = 0.0, cur_sum = 0.0;
                size_t best_s = 0, best_e = 0, cur_s = 0;
                for (size_t c = 0; c < C; ++c) {
                    SitePattern pat;
                    pat.states.resize(N);
                    for (size_t t = 0; t < N; ++t) pat.states[t] = aln.sequences[t].codons[c];
                    pat.weight = 1;
                    Scalar ll_br = FELAnalyzer::compute_pattern_log_likelihood(
                        br_tree, pat, l2t_br, mg94_base, P_br, nL_br, nS_br
                    );
                    Scalar diff = ll_br - base_site_ll[c];
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
                if (best_sum > 1.5 && best_e >= best_s + 3) {
                    RecombinationBridgeSegment seg;
                    seg.node_u = u_id;
                    seg.node_v = v_id;
                    seg.clade_u_name = (u_id >= 0 && static_cast<size_t>(u_id) < work_tree.nodes.size())
                                     ? work_tree.nodes[u_id].name : ("Node_" + std::to_string(u_id));
                    seg.target_v_name = (v_id >= 0 && static_cast<size_t>(v_id) < work_tree.nodes.size())
                                      ? work_tree.nodes[v_id].name : ("Node_" + std::to_string(v_id));
                    seg.bridge_conductance = br.shunt_stress;
                    seg.peak_start_codon = best_s;
                    seg.peak_end_codon = best_e;
                    seg.segment_delta_log_l = best_sum;
                    seg.p_value = stats::chi2_p_value(2.0 * best_sum, 2.0);
                    if (seg.p_value <= 0.01 && br_tree.num_nodes() == work_tree.num_nodes()) {
                        seg.routed_to_spr_tree = true;
                        int r_idx = static_cast<int>(routed_spr_trees.size());
                        routed_spr_trees.push_back({
                            br_tree, l2t_br,
                            "SPR(" + seg.clade_u_name + "->" + seg.target_v_name + ")",
                            best_s, best_e
                        });
                        for (size_t c = best_s; c <= best_e && c < C; ++c) {
                            site_routed_tree_idx[c] = r_idx;
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
                "Phase II: Compositional Drift & Root Polarity Audit",
                "CP3 χ² Stationarity Test • Unified 3-Regime Softmax Transition Kernel"
            );
        }
        res.composition_test = DirectionalFELAnalyzer::run_composition_test(aln, 0.05);
        res.directional_drift_active = res.composition_test.is_heterogeneous;

        std::unique_ptr<DirectionalFELAnalyzer> dfel_analyzer;
        if (res.directional_drift_active) {
            DirectionalFELSettings dfs;
            dfs.p_value_threshold = settings.p_value_threshold;
            dfs.root_strategy = "auto";
            dfs.num_regimes = 3;
            AutoRootResult root_res = DirectionalFELAnalyzer::auto_select_root(work_tree, aln, gtr_fit.params, dfs);
            res.chosen_root_branch = root_res.chosen_branch_name;
            res.chosen_root_weight = root_res.chosen_weight;
            dfel_analyzer = std::make_unique<DirectionalFELAnalyzer>(
                root_res.best_tree, aln, mg_params, dfs, res.composition_test, root_res
            );
            dfel_analyzer->run();
            // Adopt the directional rooted tree if node count matches
            if (dfel_analyzer->tree.num_nodes() == work_tree.num_nodes()) {
                work_tree = dfel_analyzer->tree;
                leaf_to_taxon = LikelihoodEngine::build_leaf_to_taxon_map(work_tree, aln);
            }
        } else {
            res.chosen_root_branch = "Stationary (Unrooted)";
            res.chosen_root_weight = 1.0;
        }
        if (!settings.quiet) {
            std::ostringstream ss;
            ss << std::fixed << std::setprecision(3)
               << "      CP3 χ² = " << res.composition_test.chi2_cp3
               << " (p = " << res.composition_test.p_value_cp3 << ")"
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
                "Phase III: Shared Inside-Outside Tensor & Active-Set Polish",
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

            MG94Matrix m_syn_k, m_non_k;
            m_syn_k.update(p_syn, pos_freqs, pi_k, gcode);
            m_non_k.update(p_non, pos_freqs, pi_k, gcode);

            Vector sqrt_pi_k = pi_k.array().max(1e-15).sqrt();
            Vector inv_sqrt_pi_k = 1.0 / sqrt_pi_k.array();

            for (size_t g = 0; g < G; ++g) {
                Matrix Q_g = m_syn_k.Q + omega_grid[g] * m_non_k.Q;
                for (int i = 0; i < S; ++i) {
                    Q_g(i, i) = 0.0;
                    Q_g(i, i) = -Q_g.row(i).sum();
                }
                Matrix sym = sqrt_pi_k.asDiagonal() * Q_g * inv_sqrt_pi_k.asDiagonal();
                sym = 0.5 * (sym + sym.transpose());
                Eigen::SelfAdjointEigenSolver<Matrix> es(sym);
                regime_spectral_grids[k_reg][g].evals = es.eigenvalues();
                regime_spectral_grids[k_reg][g].R = inv_sqrt_pi_k.asDiagonal() * es.eigenvectors();
                regime_spectral_grids[k_reg][g].L = es.eigenvectors().transpose() * sqrt_pi_k.asDiagonal();
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

        // Decision 4: Helper to evaluate branch transition matrix for grid index g (or continuous ω via log-linear spectral interpolation)
        auto branch_grid_transition_matrix = [&](int32_t bid, size_t g_idx, Scalar eff_t) -> Matrix {
            Scalar t_pos = std::max<Scalar>(0.0, eff_t);
            if (num_k_regimes == 1 || !dfel_analyzer || bid >= static_cast<int32_t>(dfel_analyzer->branch_regime_weights.size())) {
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

        // Compute Inside L_{v, c} and Outside U_{u, c} for every site c
        // incorporating Decision 2 (K=3 POA path ensemble marginalization for PP_c < 0.99)
        // and Decision 3 (Segment-specific SPR tree routing for recombinant segments p <= 0.01)
        std::vector<std::vector<Vector>> site_inside(C, std::vector<Vector>(num_nodes, Vector::Ones(S)));
        std::vector<std::vector<Vector>> site_outside(C, std::vector<Vector>(num_nodes, Vector::Zero(S)));
        std::vector<Scalar> site_base_ll(C, 0.0);
        std::vector<size_t> site_ensemble_k(C, 1);

        const Alignment& init_aln = res.joint_align.initial_alignment;
        bool has_init_aln = (init_aln.num_taxa == N && init_aln.num_codons == C);

        for (size_t c = 0; c < C; ++c) {
            Scalar alpha_c = (c < res.joint_align.site_results.size())
                           ? std::max<Scalar>(0.05, res.joint_align.site_results[c].marginalized_alpha) : 1.0;
            Scalar beta_c = (c < res.joint_align.site_results.size())
                          ? std::max<Scalar>(0.0, res.joint_align.site_results[c].marginalized_beta) : mg_params.beta;
            Scalar omega_c = beta_c / std::max<Scalar>(0.05, alpha_c);

            size_t g_near = 0;
            Scalar best_d = 1e9;
            for (size_t g = 0; g < G; ++g) {
                Scalar d = std::abs(omega_grid[g] - omega_c);
                if (d < best_d) { best_d = d; g_near = g; }
            }

            // Decision 3: Choose active tree for site c (Backbone vs Routed Recombinant SPR Tree)
            int r_idx = site_routed_tree_idx[c];
            const Tree& active_tree = (r_idx >= 0) ? routed_spr_trees[r_idx].spr_tree : work_tree;
            const auto& active_l2t = (r_idx >= 0) ? routed_spr_trees[r_idx].leaf_to_taxon : leaf_to_taxon;

            std::vector<Matrix> P_site(num_nodes);
            for (const auto& nd : active_tree.nodes) {
                if (nd.id != active_tree.root_id && nd.parent_id != INVALID_INDEX) {
                    P_site[nd.id] = branch_grid_transition_matrix(nd.id, g_near, alpha_c * nd.branch_length * conv);
                }
            }

            // Decision 2: Construct up to K=3 POA path realizations and weights when PP_c < 0.99
            Scalar pp_c = (c < res.joint_align.site_homology_posteriors.size())
                        ? res.joint_align.site_homology_posteriors[c] : 1.0;
            struct ColumnRealization {
                std::vector<int> codons;
                Scalar weight;
            };
            std::vector<ColumnRealization> realizations;
            std::vector<int> map_col(N, -1);
            for (size_t t = 0; t < N; ++t) map_col[t] = aln.sequences[t].codons[c];

            if (pp_c >= 0.99) {
                realizations.push_back({map_col, 1.0});
            } else {
                Scalar w1 = std::clamp<Scalar>(pp_c, 0.50, 0.98);
                Scalar rem = 1.0 - w1;
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
                        alt2_col[t] = -1; // Indel boundary gap-shift realization
                    }
                }
                realizations.push_back({map_col, w1});
                realizations.push_back({alt1_col, 0.65 * rem});
                realizations.push_back({alt2_col, 0.35 * rem});
            }
            site_ensemble_k[c] = realizations.size();

            // Accumulate posterior-weighted Inside and Outside cavity vectors across the K realizations
            std::vector<Vector> acc_inside(num_nodes, Vector::Zero(S));
            std::vector<Vector> acc_outside(num_nodes, Vector::Zero(S));
            Scalar acc_ll = 0.0;

            for (const auto& realiz : realizations) {
                std::vector<Vector> loc_in(num_nodes, Vector::Ones(S));
                std::vector<Vector> loc_out(num_nodes, Vector::Zero(S));
                std::vector<Scalar> log_sc(num_nodes, 0.0);

                for (int32_t nid : active_tree.post_order) {
                    const auto& nd = active_tree.nodes[nid];
                    if (nd.is_leaf) {
                        size_t tid = (nid < static_cast<int32_t>(active_l2t.size())) ? active_l2t[nid] : INVALID_INDEX;
                        if (tid < N && realiz.codons[tid] >= 0) {
                            loc_in[nid].setZero();
                            loc_in[nid](realiz.codons[tid]) = 1.0;
                        } else {
                            loc_in[nid].setOnes();
                        }
                    } else {
                        Vector prod = Vector::Ones(S);
                        Scalar sc = 0.0;
                        for (int32_t ch : nd.children) {
                            prod.array() *= (P_site[ch] * loc_in[ch]).array();
                            sc += log_sc[ch];
                        }
                        Scalar mx = prod.maxCoeff();
                        if (mx > 0.0) {
                            prod /= mx;
                            sc += std::log(mx);
                        }
                        loc_in[nid] = prod;
                        log_sc[nid] = sc;
                    }
                }
                Scalar root_dot = std::max<Scalar>(1e-300, pi_codon.dot(loc_in[active_tree.root_id]));
                acc_ll += realiz.weight * (log_sc[active_tree.root_id] + std::log(root_dot));

                std::vector<Vector> node_post(num_nodes, Vector::Zero(S));
                node_post[active_tree.root_id] = pi_codon;
                for (auto it = active_tree.post_order.rbegin(); it != active_tree.post_order.rend(); ++it) {
                    int32_t uid = *it;
                    const auto& u_nd = active_tree.nodes[uid];
                    if (u_nd.children.empty()) continue;

                    std::vector<Vector> child_msgs(u_nd.children.size());
                    for (size_t k = 0; k < u_nd.children.size(); ++k) {
                        int32_t vid = u_nd.children[k];
                        child_msgs[k] = P_site[vid] * loc_in[vid];
                    }
                    for (size_t k = 0; k < u_nd.children.size(); ++k) {
                        int32_t vid = u_nd.children[k];
                        Vector cav_u = node_post[uid];
                        for (size_t j = 0; j < u_nd.children.size(); ++j) {
                            if (j != k) cav_u.array() *= child_msgs[j].array();
                        }
                        Scalar s_sum = cav_u.sum();
                        if (s_sum > 0.0) cav_u /= s_sum;
                        else cav_u = pi_codon;
                        loc_out[vid] = cav_u;

                        Vector post_v = P_site[vid].transpose() * cav_u;
                        Scalar pv_sum = post_v.sum();
                        if (pv_sum > 0.0) post_v /= pv_sum;
                        else post_v = pi_codon;
                        node_post[vid] = post_v;
                    }
                }

                for (size_t nid = 0; nid < num_nodes; ++nid) {
                    Scalar in_sum = loc_in[nid].sum();
                    Vector norm_in = (in_sum > 0.0) ? (loc_in[nid] / in_sum).eval() : loc_in[nid];
                    acc_inside[nid].noalias() += realiz.weight * norm_in;
                    acc_outside[nid].noalias() += realiz.weight * loc_out[nid];
                }
            }
            site_inside[c] = std::move(acc_inside);
            site_outside[c] = std::move(acc_outside);
            site_base_ll[c] = acc_ll;
        }

        // Build the 3D Branch-Site-Omega Relative Likelihood Sufficiency Tensor W[c][bi][g]
        // paired with Principles 2 + 3: 2-State Spatial Markov Error-Sink HMM along each branch (bi, c = 0..C-1)
        // and Sink-Quarantined Branch-Length & Tensor Sanitization.
        std::vector<std::vector<std::array<Scalar, 16>>> W_tensor(
            C, std::vector<std::array<Scalar, 16>>(B)
        );
        std::vector<std::vector<Scalar>> W_sink(C, std::vector<Scalar>(B, 0.0));
        std::vector<std::vector<Scalar>> markov_sink_gamma(C, std::vector<Scalar>(B, 0.0));

        Vector inv_one_minus_pi = Vector::Zero(S);
        for (int i = 0; i < S; ++i) {
            inv_one_minus_pi(i) = 1.0 / std::max<Scalar>(1e-12, 1.0 - pi_codon(i));
        }

        std::vector<Scalar> pos_bls;
        for (size_t bi = 0; bi < B; ++bi) {
            int32_t bid = active_branches[bi];
            if (work_tree.nodes[bid].branch_length > 1e-6) {
                pos_bls.push_back(work_tree.nodes[bid].branch_length * conv);
            }
        }
        std::sort(pos_bls.begin(), pos_bls.end());
        Scalar med_bl_conv = pos_bls.empty() ? 0.05 : pos_bls[pos_bls.size() / 2];

        for (size_t bi = 0; bi < B; ++bi) {
            int32_t bid = active_branches[bi];
            Scalar raw_bl = work_tree.nodes[bid].branch_length * conv;
            Scalar ref_bl = raw_bl;

            // Principle 3 (Pass 0): If a terminal branch is >2.5x longer than its sibling, evaluate R_sink at 1.5x sibling length
            // so a 10-codon contamination stretch that stretched t_b 3x in Stage 1 doesn't mask its own Bayes Factor!
            const auto& nd_b = work_tree.nodes[bid];
            if (nd_b.is_leaf && nd_b.parent_id >= 0) {
                const auto& par = work_tree.nodes[nd_b.parent_id];
                Scalar min_sib = 1e9;
                for (int32_t ch : par.children) {
                    if (ch != bid && work_tree.nodes[ch].branch_length > 1e-6) {
                        min_sib = std::min(min_sib, work_tree.nodes[ch].branch_length * conv);
                    }
                }
                if (min_sib < 1e8 && raw_bl > 2.5 * std::max(min_sib, 0.5 * med_bl_conv)) {
                    ref_bl = 1.5 * std::max(min_sib, 0.5 * med_bl_conv);
                }
            }

            for (int pass = 0; pass < 2; ++pass) {
                Scalar eval_bl = (pass == 0) ? ref_bl : (work_tree.nodes[bid].branch_length * conv);
                Matrix P_ref_sink = branch_grid_transition_matrix(bid, 6, eval_bl);
                for (size_t c = 0; c < C; ++c) {
                    const Vector& U = site_outside[c][bid];
                    const Vector& L = site_inside[c][bid];
                    Scalar ref_val = std::max<Scalar>(1e-30, U.dot(P_ref_sink * L));
                    Scalar pi_dot_L = pi_codon.dot(L);
                    Scalar sink_val = 0.0;
                    for (int i = 0; i < S; ++i) {
                        Scalar p_sink_L_i = (pi_dot_L - pi_codon(i) * L(i)) * inv_one_minus_pi(i);
                        sink_val += U(i) * std::max<Scalar>(0.0, p_sink_L_i);
                    }
                    W_sink[c][bi] = std::max<Scalar>(0.0, sink_val / ref_val);
                }

                // Principle 2: 2-State Forward-Backward Markov HMM along sites c = 0..C-1 on branch bi
                if (C >= 2) {
                    const Scalar eps_entry = 0.0005;
                    const Scalar rho_stay = 0.80;
                    const Scalar log_A00 = std::log(1.0 - eps_entry);
                    const Scalar log_A01 = std::log(eps_entry);
                    const Scalar log_A10 = std::log(1.0 - rho_stay);
                    const Scalar log_A11 = std::log(rho_stay);
                    auto lse = [](Scalar a, Scalar b) -> Scalar {
                        Scalar mx = std::max(a, b);
                        return mx + std::log(std::exp(a - mx) + std::exp(b - mx));
                    };
                    std::vector<Scalar> f0(C, 0.0), f1(C, 0.0), b0(C, 0.0), b1(C, 0.0), log_em1(C, -30.0);
                    for (size_t c = 0; c < C; ++c) {
                        log_em1[c] = std::min<Scalar>(8.0, std::log(std::max<Scalar>(1e-15, W_sink[c][bi])));
                    }
                    f0[0] = log_A00;
                    f1[0] = log_A01 + log_em1[0];
                    for (size_t c = 1; c < C; ++c) {
                        f0[c] = lse(f0[c - 1] + log_A00, f1[c - 1] + log_A10);
                        f1[c] = log_em1[c] + lse(f0[c - 1] + log_A01, f1[c - 1] + log_A11);
                    }
                    for (int c = static_cast<int>(C) - 2; c >= 0; --c) {
                        b0[c] = lse(log_A00 + b0[c + 1], log_A01 + log_em1[c + 1] + b1[c + 1]);
                        b1[c] = lse(log_A10 + b0[c + 1], log_A11 + log_em1[c + 1] + b1[c + 1]);
                    }
                    Scalar log_Z = lse(f0[C - 1], f1[C - 1]);
                    for (size_t c = 0; c < C; ++c) {
                        Scalar g1 = std::exp(std::clamp<Scalar>(f1[c] + b1[c] - log_Z, -50.0, 0.0));
                        markov_sink_gamma[c][bi] = std::clamp<Scalar>(g1, 0.0, 1.0);
                    }

                    // Retain only contiguous Markov sink contamination runs of >= 3 codons (run_sum >= 2.6)
                    std::vector<Scalar> filtered_gamma(C, 0.0);
                    size_t c_idx = 0;
                    while (c_idx < C) {
                        if (markov_sink_gamma[c_idx][bi] >= 0.50) {
                            size_t run_start = c_idx;
                            Scalar run_sum = 0.0;
                            while (c_idx < C && (markov_sink_gamma[c_idx][bi] >= 0.50 ||
                                   (c_idx + 1 < C && markov_sink_gamma[c_idx + 1][bi] >= 0.50))) {
                                run_sum += markov_sink_gamma[c_idx][bi];
                                ++c_idx;
                            }
                            if (c_idx - run_start >= 3 && run_sum >= 2.6) {
                                for (size_t k = run_start; k < c_idx; ++k) {
                                    filtered_gamma[k] = markov_sink_gamma[k][bi];
                                }
                            }
                        } else {
                            ++c_idx;
                        }
                    }
                    Scalar sum_gamma = 0.0;
                    for (size_t c = 0; c < C; ++c) {
                        markov_sink_gamma[c][bi] = filtered_gamma[c];
                        sum_gamma += filtered_gamma[c];
                    }

                    // Principle 3 (Pass 0 -> Pass 1): If branch bi has a contiguous Markov error-sink tract (sum_gamma >= 1.8 codons),
                    // sanitize its inflated branch length before Pass 1 and before building W_tensor!
                    if (pass == 0 && sum_gamma >= 1.8) {
                        Scalar frac_contam = sum_gamma / static_cast<Scalar>(C);
                        Scalar est_spurious_bl = 2.2 * frac_contam;
                        Scalar clean_bl = raw_bl;
                        if (raw_bl > est_spurious_bl + 0.005) {
                            clean_bl = std::min(clean_bl, std::max<Scalar>(0.005, raw_bl - est_spurious_bl));
                        }
                        if (ref_bl < raw_bl) {
                            clean_bl = std::min(clean_bl, ref_bl);
                        }
                        work_tree.nodes[bid].branch_length = clean_bl / conv;
                    }
                }
            }

            Scalar bl = work_tree.nodes[bid].branch_length * conv;
            Matrix U_mat(S, C), L_mat(S, C);
            for (size_t c = 0; c < C; ++c) {
                U_mat.col(static_cast<Eigen::Index>(c)) = site_outside[c][bid];
                L_mat.col(static_cast<Eigen::Index>(c)) = site_inside[c][bid];
            }

            for (size_t g = 0; g < G; ++g) {
                Matrix P_bg = branch_grid_transition_matrix(bid, g, bl);
                Matrix PL = P_bg * L_mat; // Single BLAS GEMM across all C sites
                for (size_t c = 0; c < C; ++c) {
                    Scalar l_val = U_mat.col(static_cast<Eigen::Index>(c)).dot(PL.col(static_cast<Eigen::Index>(c)));
                    W_tensor[c][bi][g] = std::max<Scalar>(1e-30, l_val);
                }
            }
            for (size_t c = 0; c < C; ++c) {
                Scalar ref = std::max<Scalar>(1e-30, W_tensor[c][bi][6]); // ω = 0.50 reference
                for (size_t g = 0; g < G; ++g) {
                    W_tensor[c][bi][g] /= ref;
                }
                Scalar g_mkt = markov_sink_gamma[c][bi];
                if (g_mkt >= 0.50) {
                    for (size_t g = 0; g < G; ++g) {
                        W_tensor[c][bi][g] = (1.0 - g_mkt) * W_tensor[c][bi][g] + g_mkt * 1.0;
                    }
                    W_sink[c][bi] = 1.0;
                }
            }
        }

        // Post-pass: if any branch bi had a contiguous Markov error-sink tract (sum_gamma >= 1.8 codons),
        // also quarantine its sibling branches at those exact codons so the contaminated leaf does not
        // leak spurious substitutions into its sister branch via the parent node!
        for (size_t bi = 0; bi < B; ++bi) {
            Scalar branch_gamma_sum = 0.0;
            for (size_t c = 0; c < C; ++c) branch_gamma_sum += markov_sink_gamma[c][bi];
            if (branch_gamma_sum < 1.8) continue;
            int32_t bid = active_branches[bi];
            const auto& nd_b = work_tree.nodes[bid];
            if (nd_b.parent_id < 0) continue;
            const auto& par = work_tree.nodes[nd_b.parent_id];
            for (size_t c = 0; c < C; ++c) {
                Scalar g_mkt = markov_sink_gamma[c][bi];
                if (g_mkt >= 0.50) {
                    for (int32_t sib_id : par.children) {
                        if (sib_id != bid) {
                            for (size_t s_bi = 0; s_bi < B; ++s_bi) {
                                if (active_branches[s_bi] == sib_id) {
                                    for (size_t g = 0; g < G; ++g) {
                                        W_tensor[c][s_bi][g] = (1.0 - g_mkt) * W_tensor[c][s_bi][g] + g_mkt * 1.0;
                                    }
                                    W_sink[c][s_bi] *= (1.0 - g_mkt);
                                }
                            }
                        }
                    }
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

        auto spectral_decompose_omega = [&](Scalar omega_val) -> SpectralDecomp {
            Matrix Q_w = m_syn_base.Q + omega_val * m_non_base.Q;
            for (int i = 0; i < S; ++i) {
                Q_w(i, i) = 0.0;
                Q_w(i, i) = -Q_w.row(i).sum();
            }
            Matrix sym = sqrt_pi_base.asDiagonal() * Q_w * inv_sqrt_pi_base.asDiagonal();
            sym = 0.5 * (sym + sym.transpose());
            Eigen::SelfAdjointEigenSolver<Matrix> es(sym);
            SpectralDecomp sd;
            sd.evals = es.eigenvalues();
            sd.R = inv_sqrt_pi_base.asDiagonal() * es.eigenvectors();
            sd.L = es.eigenvectors().transpose() * sqrt_pi_base.asDiagonal();
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
            } else {
                entry.dfel_lrt = entry.fel_lrt;
                entry.dfel_p_value = entry.fel_p_value;
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
            // Diagonalizes Q(1, w1) and Q(1, w+) ONCE per Nelder-Mead evaluation and computes all B branch
            // matrices via R * exp(alpha * t_b * evals) * L!
            if (screen_p < settings.active_set_screen_p || entry.beta_pervasive > entry.alpha) {
                entry.meme_polished = true;
                int r_idx = site_routed_tree_idx[c];
                const Tree& active_tree = (r_idx >= 0) ? routed_spr_trees[r_idx].spr_tree : work_tree;
                const auto& active_l2t = (r_idx >= 0) ? routed_spr_trees[r_idx].leaf_to_taxon : leaf_to_taxon;

                SitePattern pat;
                pat.states.resize(N);
                bool is_invar = true;
                int8_t first_st = -1;
                for (size_t t = 0; t < N; ++t) {
                    int8_t st = aln.sequences[t].codons[c];
                    pat.states[t] = st;
                    if (st >= 0) {
                        if (first_st == -1) first_st = st;
                        else if (st != first_st) is_invar = false;
                    }
                }
                pat.weight = 1;

                if (!is_invar) {
                    std::vector<Vector> t_node_L(active_tree.num_nodes(), Vector::Zero(S));
                    std::vector<double> t_node_sc(active_tree.num_nodes(), 0.0);
                    std::vector<Matrix> P_br(active_tree.num_nodes());

                    auto eval_exact_meme_ll = [&](Scalar a_val, Scalar w1_val, Scalar p1_val, Scalar b_plus_val) -> Scalar {
                        Scalar a_pos = std::max<Scalar>(1e-5, a_val);
                        Scalar w_plus_val = b_plus_val / a_pos;
                        SpectralDecomp sd1 = spectral_decompose_omega(w1_val);
                        SpectralDecomp sdp = spectral_decompose_omega(w_plus_val);
                        Scalar p_plus_wt = 1.0 - p1_val;
                        for (const auto& nd : active_tree.nodes) {
                            if (nd.id != active_tree.root_id && nd.parent_id != INVALID_INDEX) {
                                Scalar eff_t = a_pos * nd.branch_length * conv;
                                Vector e1 = (sd1.evals.array() * eff_t).exp();
                                Vector ep = (sdp.evals.array() * eff_t).exp();
                                Matrix m1 = sd1.R * e1.asDiagonal() * sd1.L;
                                Matrix mp = sdp.R * ep.asDiagonal() * sdp.L;
                                P_br[nd.id] = (p1_val * m1 + p_plus_wt * mp).cwiseMax(1e-15);
                            }
                        }
                        return FELAnalyzer::compute_pattern_log_likelihood(
                            active_tree, pat, active_l2t, mg94_base, P_br, t_node_L, t_node_sc
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
            } else {
                entry.beta_minus = w_minus_opt * eff_alpha;
                entry.beta_plus = w_plus_opt * eff_alpha;
                entry.weight_plus = wp_opt;
            }

            Scalar max_site_sink_gamma = 0.0;
            for (size_t bi = 0; bi < B; ++bi) {
                max_site_sink_gamma = std::max(max_site_sink_gamma, markov_sink_gamma[c][bi]);
            }
            Scalar clean_site_factor = (max_site_sink_gamma >= 0.50) ? (1.0 - max_site_sink_gamma) : 1.0;
            entry.meme_lrt = raw_meme_lrt * entry.column_homology_posterior * clean_site_factor;
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

            auto fit_busted_em_model = [&](bool use_sink,
                                           size_t& best_g1, size_t& best_g2, size_t& best_g3,
                                           std::array<Scalar, 3>& best_w, Scalar& best_err_w,
                                           Scalar& out_null_ll, Scalar& out_alt_ll) {
                auto run_em = [&](size_t g1, size_t g2, size_t g3, bool allow_w3,
                                  std::array<Scalar, 3>& out_w, Scalar& out_w_err) -> Scalar {
                    Scalar best_run_ll = -1e18;
                    std::array<Scalar, 3> best_run_w = {0.75, 0.25, 0.0};
                    Scalar best_run_err = 0.0;
                    std::vector<Scalar> init_eps_list = use_sink
                        ? std::vector<Scalar>{0.0005, 0.0025, 0.008}
                        : std::vector<Scalar>{0.0};
                    for (Scalar init_eps : init_eps_list) {
                        std::array<Scalar, 4> w = allow_w3
                            ? std::array<Scalar, 4>{0.72, 0.24 - 0.5 * init_eps, 0.04 - 0.5 * init_eps, init_eps}
                            : std::array<Scalar, 4>{0.75, 0.25 - init_eps, 0.00, init_eps};
                        Scalar total_ll = -1e18;
                        for (int iter = 0; iter < 25; ++iter) {
                            std::array<Scalar, 4> accum = {1e-6, 1e-6, allow_w3 ? 1e-6 : 0.0, use_sink ? 1e-9 : 0.0};
                            total_ll = 0.0;
                            for (size_t c = 0; c < C; ++c) {
                                Scalar pp = res.site_atlas[c].column_homology_posterior;
                                for (size_t bi = 0; bi < B; ++bi) {
                                    Scalar l0 = w[0] * W_tensor[c][bi][g1];
                                    Scalar l1 = w[1] * W_tensor[c][bi][g2];
                                    Scalar l2 = allow_w3 ? (w[2] * W_tensor[c][bi][g3]) : 0.0;
                                    Scalar l_err = use_sink ? (w[3] * W_sink[c][bi]) : 0.0;
                                    Scalar s = std::max<Scalar>(1e-30, l0 + l1 + l2 + l_err);
                                    total_ll += pp * std::log(s);
                                    accum[0] += pp * (l0 / s);
                                    accum[1] += pp * (l1 / s);
                                    if (allow_w3) accum[2] += pp * (l2 / s);
                                    if (use_sink) accum[3] += pp * (l_err / s);
                                }
                            }
                            Scalar sum_a = accum[0] + accum[1] + accum[2] + accum[3];
                            for (int k = 0; k < 4; ++k) w[k] = accum[k] / sum_a;
                            if (use_sink && w[3] > 0.015) {
                                Scalar excess = w[3] - 0.015;
                                w[3] = 0.015;
                                w[0] += 0.5 * excess;
                                w[1] += 0.5 * excess;
                            }
                        }
                        if (total_ll > best_run_ll) {
                            best_run_ll = total_ll;
                            Scalar bio_sum = std::max<Scalar>(1e-12, w[0] + w[1] + w[2]);
                            best_run_w = {w[0] / bio_sum, w[1] / bio_sum, w[2] / bio_sum};
                            best_run_err = use_sink ? w[3] : 0.0;
                        }
                    }
                    out_w = best_run_w;
                    out_w_err = best_run_err;
                    return best_run_ll;
                };

                out_null_ll = -1e18;
                std::array<Scalar, 3> null_w = {0.75, 0.25, 0.0};
                Scalar null_err_w = 0.0;
                best_g1 = 1; best_g2 = 6; best_g3 = 8;
                for (size_t g1 : {size_t(0), size_t(1), size_t(2), size_t(3), size_t(4)}) {
                    for (size_t g2 : {size_t(5), size_t(6), size_t(7), size_t(8)}) {
                        std::array<Scalar, 3> w_tmp;
                        Scalar w_err_tmp = 0.0;
                        Scalar ll = run_em(g1, g2, 8, false, w_tmp, w_err_tmp);
                        if (ll > out_null_ll) {
                            out_null_ll = ll;
                            best_g1 = g1;
                            best_g2 = g2;
                            null_w = w_tmp;
                            null_err_w = w_err_tmp;
                        }
                    }
                }

                out_alt_ll = out_null_ll;
                best_w = null_w;
                best_err_w = null_err_w;
                for (size_t g1 : {best_g1 > 0 ? best_g1 - 1 : best_g1, best_g1, std::min<size_t>(4, best_g1 + 1)}) {
                    for (size_t g2 : {best_g2 > 5 ? best_g2 - 1 : best_g2, best_g2, std::min<size_t>(8, best_g2 + 1)}) {
                        for (size_t g3 = 9; g3 < G; ++g3) {
                            std::array<Scalar, 3> w_tmp;
                            Scalar w_err_tmp = 0.0;
                            Scalar ll = run_em(g1, g2, g3, true, w_tmp, w_err_tmp);
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

            // Fit standard 3-class BUSTED (no sink) and rank-1 P_sink BUSTED-E;
            // only retain the error sink if it improves unconstrained log-likelihood by >= 0.5 * ln(B*C) (BIC check)
            // or if a contiguous Markov error-sink tract was quarantined.
            size_t g1_std = 1, g2_std = 6, g3_std = 8;
            std::array<Scalar, 3> w_std = {0.75, 0.25, 0.0};
            Scalar err_std = 0.0, null_std = -1e18, alt_std = -1e18;
            fit_busted_em_model(false, g1_std, g2_std, g3_std, w_std, err_std, null_std, alt_std);

            size_t g1_sk = 1, g2_sk = 6, g3_sk = 8;
            std::array<Scalar, 3> w_sk = {0.75, 0.25, 0.0};
            Scalar err_sk = 0.0, null_sk = -1e18, alt_sk = -1e18;
            fit_busted_em_model(true, g1_sk, g2_sk, g3_sk, w_sk, err_sk, null_sk, alt_sk);

            Scalar sink_bic_threshold = 0.5 * std::log(static_cast<Scalar>(std::max<size_t>(2, B * C)));
            bool sink_supported = ((2.0 * (alt_sk - alt_std) > 2.0 * sink_bic_threshold) && (err_sk >= 1e-4))
                                  || (total_markov_sink_gamma >= 1.5);

            size_t best_g1 = sink_supported ? g1_sk : g1_std;
            size_t best_g2 = sink_supported ? g2_sk : g2_std;
            size_t best_g3 = sink_supported ? g3_sk : g3_std;
            std::array<Scalar, 3> weights = sink_supported ? w_sk : w_std;
            Scalar best_err_w = sink_supported ? std::max(err_sk, markov_sink_frac) : 0.0;
            Scalar best_null = sink_supported ? null_sk : null_std;
            Scalar best_alt = sink_supported ? alt_sk : alt_std;

            Scalar base_sum = 0.0;
            for (Scalar s_ll : site_base_ll) base_sum += s_ll;
            res.gene_wide.omega_classes = {omega_grid[best_g1], omega_grid[best_g2], omega_grid[best_g3]};
            res.gene_wide.omega_weights = weights;
            res.gene_wide.error_sink_omega = 100.0;
            res.gene_wide.error_sink_weight = best_err_w;
            res.gene_wide.null_log_l = base_sum + best_null;
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
            Scalar null_br_ll = -1e18;
            size_t best_gb0 = 6;
            for (size_t g0 = 0; g0 <= 8; ++g0) {
                Scalar ll = 0.0;
                for (size_t c = 0; c < C; ++c) {
                    ll += res.site_atlas[c].column_homology_posterior * std::log(W_tensor[c][bi][g0]);
                }
                if (ll > null_br_ll) {
                    null_br_ll = ll;
                    best_gb0 = g0;
                }
            }
            Scalar alt_br_ll = null_br_ll;
            size_t best_gb_plus = 8;
            Scalar best_wb_plus = 0.0;
            for (size_t gp = 9; gp < G; ++gp) {
                for (Scalar wp : {0.02, 0.05, 0.10, 0.20, 0.35, 0.50}) {
                    Scalar ll = 0.0;
                    for (size_t c = 0; c < C; ++c) {
                        Scalar m = (1.0 - wp) * W_tensor[c][bi][best_gb0] + wp * W_tensor[c][bi][gp];
                        ll += res.site_atlas[c].column_homology_posterior * std::log(std::max<Scalar>(1e-30, m));
                    }
                    if (ll > alt_br_ll) {
                        alt_br_ll = ll;
                        best_gb_plus = gp;
                        best_wb_plus = wp;
                    }
                }
            }

            Scalar w_bg_opt = omega_grid[best_gb0];
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
                        Scalar b0 = eval_branch_site_omega(c, bi, w_bg_opt);
                        l1 += pp * std::log(std::max<Scalar>(1e-30, (1.0 - wb_opt) * b0 + wb_opt * eval_branch_site_omega(c, bi, m1)));
                        l2 += pp * std::log(std::max<Scalar>(1e-30, (1.0 - wb_opt) * b0 + wb_opt * eval_branch_site_omega(c, bi, m2)));
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

            if (!settings.foreground_tag.empty()) {
                for (size_t bi = 0; bi < B; ++bi) {
                    if (res.branch_profiles[bi].branch_name.find(settings.foreground_tag) != std::string::npos) {
                        is_test_branch[bi] = true;
                        found_user_tags = true;
                    }
                }
            } else {
                for (size_t bi = 0; bi < B; ++bi) {
                    const std::string& nm = res.branch_profiles[bi].branch_name;
                    if (nm.find("{FG}") != std::string::npos || nm.find("{TEST}") != std::string::npos ||
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
                "Phase IV: Model-Free Exchangeability, Spatial Architecture & Forensics",
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
                    Matrix P_b = branch_grid_transition_matrix(bid, g_near, eff_t);
                    const Vector& U = site_outside[c][bid];
                    const Vector& L = site_inside[c][bid];
                    Scalar Z_b = std::max<Scalar>(1e-30, U.dot(P_b * L));

                    // Posterior parent and child state distributions on branch bi
                    Vector post_parent = (U.array() * (P_b * L).array()) / Z_b;
                    Vector post_child = ((P_b.transpose() * U).array() * L.array()) / Z_b;

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

