// Out-of-line JSON emitter for MSS (Multi-class Synonymous Substitution) analyses.
#include "hyphy/analyses/mss.hpp"
#include "hyphy/core/json_output.hpp"
#include "hyphy/core/provenance.hpp"

namespace hyphy::analyses {

using namespace hyphy::core;

nlohmann::json MSSAnalyzer::to_json(const MSSResult& res) {
    nlohmann::json j;

    std::string norm_str = "flux";
    if (res.settings.normalization == MSSNormalization::NeutralRef) norm_str = "neutral_ref";
    else if (res.settings.normalization == MSSNormalization::Unweighted) norm_str = "unweighted";

    j["analysis"] = {
        {"name", "MSS"},
        {"description", "MG94xREV Codon Model with Multi-class Synonymous Substitutions"},
        {"authors", "Sergei L. Kosakovsky Pond"},
        {"contact", "spond@temple.edu"},
        {"version", "3.0.0"},
        {"settings", {
            {"model", res.partition.name},
            {"normalization", norm_str},
            {"neutral_reference", (res.partition.neutral_class_id >= 0 && static_cast<size_t>(res.partition.neutral_class_id) < res.partition.class_names.size())
                                      ? res.partition.class_names[res.partition.neutral_class_id]
                                      : ""},
            {"fix_omega", res.settings.fix_omega},
            {"f81", res.settings.use_f81},
            {"num_classes", res.partition.num_classes()},
            {"num_synonymous_pairs", res.partition.num_pairs()}
        }}
    };

    j["baseline"] = {
        {"Log Likelihood", res.baseline_log_l},
        {"AIC-c", res.baseline_aicc},
        {"BIC", res.baseline_bic}
    };

    j["mss_fit"] = {
        {"Log Likelihood", res.mss_log_l},
        {"AIC-c", res.mss_aicc},
        {"BIC", res.mss_bic},
        {"parameters", res.total_parameters},
        {"sample_size", res.total_sample_size},
        {"omnibus_test", {
            {"LRT", res.omnibus_lrt},
            {"df", res.omnibus_df},
            {"p-value", res.omnibus_p_value}
        }}
    };

    nlohmann::json rates_json = nlohmann::json::object();
    nlohmann::json tests_json = nlohmann::json::object();
    for (const auto& r : res.rates) {
        rates_json[r.class_name] = {
            {"MLE", r.rate},
            {"SE", r.std_error},
            {"CI_low", r.ci_low},
            {"CI_high", r.ci_high},
            {"num_pairs", r.num_pairs},
            {"equilibrium_flux_share", r.equilibrium_flux_share},
            {"is_reference", r.is_reference}
        };
        if (!r.is_reference) {
            tests_json[r.class_name] = {
                {"LRT", r.lrt_vs_one},
                {"uncorrected_p_value", r.p_value_uncorrected},
                {"corrected_p_value", r.p_value_holm}
            };
        }
    }
    j["synonymous_rates"] = rates_json;
    j["test results"] = tests_json;

    nlohmann::json genes_json = nlohmann::json::array();
    for (const auto& g : res.genes) {
        genes_json.push_back({
            {"filepath", g.filepath},
            {"sequences", g.num_taxa},
            {"sites", g.num_codons},
            {"baseline_log_l", g.baseline_log_l},
            {"mss_log_l", g.mss_log_l},
            {"baseline_omega", g.baseline_omega},
            {"mss_omega", g.mss_omega},
            {"tree_scaler", g.tree_scaler},
            {"baseline_tree_length", g.baseline_tree_length},
            {"synonymous_tree_length", g.synonymous_tree_length},
            {"nonsynonymous_tree_length", g.nonsynonymous_tree_length}
        });
    }
    j["genes"] = genes_json;

    if (res.selection.enabled) {
        nlohmann::json sel_json;
        sel_json["ic_type"] = res.selection.ic_type;
        sel_json["num_classes"] = res.selection.num_classes;
        sel_json["null_mg94_ic"] = res.selection.null_mg94_ic;
        sel_json["synrev_codon_log_l"] = res.selection.synrev_codon_log_l;
        sel_json["synrev_codon_ic"] = res.selection.synrev_codon_ic;
        sel_json["best_partition_log_l"] = res.selection.best_partition_log_l;
        sel_json["best_partition_ic"] = res.selection.best_partition_ic;
        sel_json["models_evaluated"] = res.selection.models_evaluated;
        sel_json["class_rates"] = res.selection.class_rates;
        sel_json["class_labels"] = res.selection.class_labels;

        nlohmann::json pairs_arr = nlohmann::json::array();
        for (const auto& pe : res.selection.pair_entries) {
            pairs_arr.push_back({
                {"aa", std::string(1, pe.aa)},
                {"codon1", pe.codon1},
                {"codon2", pe.codon2},
                {"mle_rate", pe.mle_rate},
                {"curvature", pe.curvature},
                {"assigned_class", pe.assigned_class},
                {"group", pe.assigned_label},
                {"support", pe.max_support},
                {"model_averaged_rate", pe.model_averaged_rate},
                {"class_support", pe.class_support}
            });
        }
        sel_json["pairs"] = pairs_arr;
        j["partition_selection"] = sel_json;
    }

    j["timers"] = {
        {"Overall", res.runtime_seconds}
    };

    return j;
}

} // namespace hyphy::analyses
