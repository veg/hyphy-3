// Out-of-line JSON & TSV emitters for ForecastAnalyzer.
#include "hyphy/analyses/forecast.hpp"
#include "hyphy/core/json_output.hpp"
#include "hyphy/core/provenance.hpp"

namespace hyphy::analyses {

using namespace hyphy::core;

namespace {

inline nlohmann::json matrix_to_json_2d(const Matrix& M) {
    nlohmann::json rows = nlohmann::json::array();
    rows.get_ref<nlohmann::json::array_t&>().reserve(static_cast<size_t>(M.rows()));
    for (Eigen::Index r = 0; r < M.rows(); ++r) {
        nlohmann::json row = nlohmann::json::array();
        row.get_ref<nlohmann::json::array_t&>().reserve(static_cast<size_t>(M.cols()));
        for (Eigen::Index c = 0; c < M.cols(); ++c) {
            row.push_back(finite_or_null(M(r, c)));
        }
        rows.push_back(std::move(row));
    }
    return rows;
}

} // namespace

nlohmann::json ForecastResult::to_modern_json(
    const Tree& tree,
    const Alignment& aln,
    const Provenance& prov,
    const std::string& input_filepath
) const {
    AnalysisMetadata meta;
    meta.id = "forecast";
    meta.name = "Mutational Fitness & Evolutionary Trajectory Forecasting (FORECAST)";
    meta.version = HYPHY3_VERSION;
    meta.category = "evolutionary_forecasting";
    meta.description =
        "Reversible site-specific codon generator combining lineage-partitioned FEL-SRV, "
        "hierarchical G+S-PRIME Empirical Bayes biophysical exchangeabilities, "
        "multi-class synonymous substitution (SynREVCodon/MSS) codon-pair rates, and "
        "exact Inside-Outside marginal internal-node amino-acid preferences.";
    meta.citations = {
        citation("Kosakovsky Pond SL, Frost SDW (2005) Not So Different After All: A Comparison of Methods for Detecting Amino Acid Sites Under Selection. Mol Biol Evol 22(5):1208-1222.", "10.1093/molbev/msi105", "15703242"),
        citation("Conant GC, Stadler PF (2009) Solvent exposure imparts similar selective pressures across a range of yeast proteins. Mol Biol Evol 26(5):1155-1161.", "10.1093/molbev/msp031", "19233963"),
        citation("Lucaci AG et al. (2023) Extra base hits: widespread empirical support for instantaneous multiple-nucleotide changes and synonymous rate variation.", "10.1371/journal.pone.0248337", "33720978")
    };
    meta.settings = {
        {"branches", effective_branches},
        {"property_set", settings.custom_property_file.empty() ? settings.property_set : settings.custom_property_file},
        {"shrinkage_lambda", finite_or_null(settings.shrinkage_lambda)},
        {"shrinkage_prior", settings.shrinkage_prior},
        {"max_bound", finite_or_null(settings.max_bound)},
        {"syn_model", effective_syn_model},
        {"syn_shrinkage", finite_or_null(settings.syn_shrinkage)},
        {"ancestral_exponent", finite_or_null(settings.ancestral_exponent)},
        {"ancestral_pseudocount", finite_or_null(settings.ancestral_pseudocount)},
        {"omega_prior_weight", finite_or_null(settings.omega_prior_weight)},
        {"rate_prior_weight", finite_or_null(settings.rate_prior_weight)},
        {"horizon", finite_or_null(settings.horizon)},
        {"reference_taxon", resolved_reference_taxon}
    };

    nlohmann::json j = new_analysis_json(meta, prov, input_filepath);
    j["dataset"] = dataset_block(aln);
    j["phylogeny"] = phylogeny_block(
        tree, "",
        [](const TreeNode& node) {
            return std::map<std::string, Scalar>{{"codon_prime", node.branch_length}};
        }
    );

    // Model fits block
    nlohmann::json gprime_fit = {
        {"description", "Global property-informed codon baseline (G-PRIME)"},
        {"log_likelihood", finite_or_null(prime_global_log_l)},
        {"aicc", finite_or_null(prime_global_aicc)},
        {"global_omega", finite_or_null(global_omega)},
        {"global_lambdas", sanitized_doubles(global_lambdas)},
        {"property_names", property_names}
    };

    nlohmann::json mss_fit = {
        {"description", "Multi-class synonymous codon-pair exchangeability fit (" + effective_syn_model + ")"},
        {"model", effective_syn_model},
        {"baseline_log_likelihood", finite_or_null(mss_baseline_log_l)},
        {"log_likelihood", finite_or_null(mss_log_l)},
        {"classes_count", static_cast<int>(mss_num_classes)},
        {"class_names", mss_class_names},
        {"class_rates", sanitized_doubles(mss_class_rates)}
    };

    j["model_fits"] = {
        {"nucleotide_gtr", {
            {"description", "General Time Reversible nucleotide baseline"},
            {"log_likelihood", finite_or_null(gtr_log_l)}
        }},
        {"global_prime", gprime_fit},
        {"synonymous_mss", mss_fit}
    };

    j["statistical_tests"] = {
        {"synonymous_selection", {
            {"model", effective_syn_model},
            {"lrt", finite_or_null(mss_omnibus_lrt)},
            {"p_value", finite_or_null(mss_omnibus_p_value)},
            {"df", static_cast<int>((mss_num_classes > 0) ? (mss_num_classes - 1) : 0)}
        }}
    };

    // Columnar site summary table
    ColumnarResults table;
    std::vector<int> col_site;
    std::vector<std::string> col_focal_cod, col_focal_aa, col_root_cod, col_root_aa;
    std::vector<double> col_alpha, col_beta, col_beta_bg, col_omega_eff,
                        col_olrt, col_op, col_oq, col_exp_sub, col_ent, col_neff;
    std::map<std::string, std::vector<double>> col_lambdas;

    for (const auto& pname : property_names) {
        col_lambdas[pname] = {};
        col_lambdas[pname].reserve(site_results.size());
    }

    col_site.reserve(site_results.size());
    col_focal_cod.reserve(site_results.size());
    col_focal_aa.reserve(site_results.size());
    col_root_cod.reserve(site_results.size());
    col_root_aa.reserve(site_results.size());
    col_alpha.reserve(site_results.size());
    col_beta.reserve(site_results.size());
    col_beta_bg.reserve(site_results.size());
    col_omega_eff.reserve(site_results.size());
    col_olrt.reserve(site_results.size());
    col_op.reserve(site_results.size());
    col_oq.reserve(site_results.size());
    col_exp_sub.reserve(site_results.size());
    col_ent.reserve(site_results.size());
    col_neff.reserve(site_results.size());

    for (const auto& sr : site_results) {
        col_site.push_back(static_cast<int>(sr.site_index) + 1);
        col_focal_cod.push_back(sr.focal_codon);
        col_focal_aa.push_back(std::string(1, sr.focal_aa));
        col_root_cod.push_back(sr.root_codon);
        col_root_aa.push_back(std::string(1, sr.root_aa));
        col_alpha.push_back(sr.alpha);
        col_beta.push_back(sr.beta);
        col_beta_bg.push_back(sr.beta_background);
        col_omega_eff.push_back(sr.effective_omega);
        col_olrt.push_back(sr.omnibus_lrt);
        col_op.push_back(sr.omnibus_p_value);
        col_oq.push_back(sr.omnibus_q_value);
        col_exp_sub.push_back(sr.expected_substitutions_horizon);
        col_ent.push_back(sr.stationary_aa_entropy);
        col_neff.push_back(sr.effective_aa_count);

        for (size_t k = 0; k < property_names.size(); ++k) {
            Scalar v = (k < sr.lambdas.size()) ? sr.lambdas[k] : 0.0;
            col_lambdas[property_names[k]].push_back(v);
        }
    }

    table.add_column("site", {"integer", "1-based codon site index"}, std::move(col_site));
    table.add_column("focal_codon", {"string", "Focal / wildtype sense codon"}, std::move(col_focal_cod));
    table.add_column("focal_aa", {"string", "Focal / wildtype amino acid"}, std::move(col_focal_aa));
    table.add_column("root_codon", {"string", "ML marginal root sense codon"}, std::move(col_root_cod));
    table.add_column("root_aa", {"string", "ML marginal root amino acid"}, std::move(col_root_aa));
    table.add_column("alpha", {"number", "Site synonymous substitution rate (SRV)"}, std::move(col_alpha));
    table.add_column("beta", {"number", "Site nonsynonymous substitution rate on tested branches"}, std::move(col_beta));
    table.add_column("beta_background", {"number", "Site nonsynonymous substitution rate on background branches"}, std::move(col_beta_bg));
    table.add_column("effective_omega", {"number", "Empirical Bayes smoothed site dN/dS ratio"}, std::move(col_omega_eff));
    table.add_column("omnibus_lrt", {"number", "S-PRIME omnibus LRT statistic"}, std::move(col_olrt));
    table.add_column("omnibus_p_value", {"number", "S-PRIME raw omnibus p-value"}, std::move(col_op));
    table.add_column("omnibus_q_value", {"number", "S-PRIME Benjamini-Hochberg FDR q-value"}, std::move(col_oq));
    table.add_column("expected_substitutions_horizon", {"number", "Probability of leaving focal codon over forecast horizon"}, std::move(col_exp_sub));
    table.add_column("stationary_aa_entropy", {"number", "Shannon entropy (nats) of stationary amino-acid distribution"}, std::move(col_ent));
    table.add_column("effective_aa_count", {"number", "Effective number of tolerated amino acids exp(H)"}, std::move(col_neff));

    for (const auto& pname : property_names) {
        table.add_column("lambda_" + pname, {"number", "Hierarchical Empirical Bayes property importance for " + pname}, std::move(col_lambdas[pname]));
    }

    j["site_results"] = table.to_json();
    j["timers"] = {
        {"total_seconds", finite_or_null(runtime_seconds)}
    };

    // Method-specific dense forecasting tensors
    j["method_results"] = {
        {"aa_alphabet", aa_alphabet},
        {"sense_codons", sense_codons},
        {"sense_codon_aa", sense_codon_aa},
        {"property_names", property_names},
        {"horizon", finite_or_null(settings.horizon)},
        {"horizon_tau", finite_or_null(horizon_tau)},
        {"baseline_scale_factor", finite_or_null(baseline_scale_factor)},
        {"aa_dms_matrix", matrix_to_json_2d(aa_dms_matrix)},
        {"aa_stationary_matrix", matrix_to_json_2d(aa_stationary_matrix)},
        {"ancestral_aa_matrix", matrix_to_json_2d(ancestral_aa_matrix)},
        {"codon_rate_matrix", matrix_to_json_2d(codon_rate_matrix)},
        {"codon_log_prob_matrix", matrix_to_json_2d(codon_log_prob_matrix)}
    };

    return j;
}

void ForecastResult::write_dms_tsv(const std::string& filepath) const {
    std::ofstream out(filepath);
    if (!out.is_open()) {
        throw std::runtime_error("Cannot open DMS TSV output file for writing: " + filepath);
    }

    out << "site\tfocal_codon\tfocal_aa\troot_codon\troot_aa\talpha\tbeta\teffective_omega\teffective_aa_count";
    for (const auto& aa : aa_alphabet) out << "\tdms_" << aa;
    for (const auto& aa : aa_alphabet) out << "\tpi_" << aa;
    for (const auto& aa : aa_alphabet) out << "\tanc_" << aa;
    out << "\n";

    out << std::fixed << std::setprecision(6);
    for (size_t s = 0; s < site_results.size(); ++s) {
        const auto& sr = site_results[s];
        out << (sr.site_index + 1) << "\t"
            << sr.focal_codon << "\t"
            << sr.focal_aa << "\t"
            << sr.root_codon << "\t"
            << sr.root_aa << "\t"
            << sr.alpha << "\t"
            << sr.beta << "\t"
            << sr.effective_omega << "\t"
            << sr.effective_aa_count;
        for (int a = 0; a < 20; ++a) out << "\t" << aa_dms_matrix(s, a);
        for (int a = 0; a < 20; ++a) out << "\t" << aa_stationary_matrix(s, a);
        for (int a = 0; a < 20; ++a) out << "\t" << ancestral_aa_matrix(s, a);
        out << "\n";
    }

    close_or_throw(out, filepath);
}

void ForecastResult::write_codon_tsv(const std::string& filepath) const {
    std::ofstream out(filepath);
    if (!out.is_open()) {
        throw std::runtime_error("Cannot open codon forecast TSV output file for writing: " + filepath);
    }

    out << "site\tfocal_codon\tfocal_aa\texpected_subs_horizon";
    for (const auto& cod : sense_codons) out << "\trate_" << cod;
    for (const auto& cod : sense_codons) out << "\tlogprob_" << cod;
    out << "\n";

    out << std::fixed << std::setprecision(6);
    const int S = static_cast<int>(sense_codons.size());
    for (size_t s = 0; s < site_results.size(); ++s) {
        const auto& sr = site_results[s];
        out << (sr.site_index + 1) << "\t"
            << sr.focal_codon << "\t"
            << sr.focal_aa << "\t"
            << sr.expected_substitutions_horizon;
        for (int c = 0; c < S; ++c) out << "\t" << codon_rate_matrix(s, c);
        for (int c = 0; c < S; ++c) out << "\t" << codon_log_prob_matrix(s, c);
        out << "\n";
    }

    close_or_throw(out, filepath);
}

} // namespace hyphy::analyses
