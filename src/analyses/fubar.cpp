// Out-of-line JSON emitter definitions for fubar.
// Moved from include/hyphy/analyses/fubar.hpp so that emitter
// edits recompile this single TU instead of every consumer
// (tests, apps, Python bindings).
#include "hyphy/analyses/fubar.hpp"
#include "hyphy/core/json_output.hpp"
#include "hyphy/core/provenance.hpp"

namespace hyphy::analyses {

using namespace hyphy::core;
using namespace hyphy::opt;

std::string FUBARAnalyzer::to_modern_json(
        const FUBARResult& res,
        const std::string& cli_cmd,
        const std::string& input_filepath) const {
        AnalysisMetadata meta{
            "fubar",
            "Fast, Unconstrained Bayesian AppRoximation",
            "3.0.0",
            "site-level selection",
            "Site-by-site detection of positive and purifying selection via Bayesian rate grid approximation",
            {citation("Murrell B, et al. FUBAR: A Fast, Unconstrained Bayesian AppRoximation for Inferring Selection. Mol Biol Evol. 30(5):1196-1205 (2013).",
                      "10.1093/molbev/mst030", "23420840")},
            {
                {"grid_points", settings.grid_points},
                {"total_grid_points", settings.grid_points * settings.grid_points},
                {"method", fubar_method_to_string(settings.method)},
                {"concentration", settings.concentration},
                {"posterior_threshold", settings.posterior_threshold},
                {"non_zero", settings.non_zero},
                {"tree_scaling", settings.tree_scaling}
            }
        };

        nlohmann::json j;
        j["$schema"] = "https://raw.githubusercontent.com/veg/hyphy-3/main/schemas/v3/fubar.v3.schema.json";
        j["schema_version"] = meta.version;
        j["analysis"] = meta.to_json();
        j["software"] = software_block();
        // FUBAR's signature carries cli_cmd (not a Provenance), so build the
        // auto block inline and attach the real runtime.
        {
            Provenance prov;
            prov.invocation.cli_command = cli_cmd;
            prov.invocation.working_directory = Provenance::get_cwd();
            if (!input_filepath.empty()) {
                size_t sz = 0;
                std::string hash = crypto::SHA256::hash_file(input_filepath, &sz);
                if (!hash.empty()) {
                    prov.inputs["alignment"] = {input_filepath, "FASTA/NEXUS", hash, sz};
                }
            }
            prov.execution.start_time = Provenance::current_iso8601();
            prov.execution.end_time = prov.execution.start_time;
            prov.execution.wall_time_seconds = res.time_total_sec;
            prov.execution.cpu_threads = 1;
            prov.execution.hostname = Provenance::get_hostname();
            prov.execution.os = Provenance::detect_os();
            prov.execution.compiler = Provenance::detect_compiler();
            j["provenance"] = prov.to_json();
        }

        // Dataset — shape pinned by fubar.v3.schema.json (file/taxa/partitions
        // with coverage), distinct from the shared dataset_block contract.
        j["dataset"] = {
            {"file", input_filepath},
            {"taxa", aln.taxon_names},
            {"taxa_count", aln.num_taxa},
            {"codon_sites", aln.num_codons},
            {"nucleotide_sites", aln.num_nucleotides > 0 ? aln.num_nucleotides : aln.num_codons * 3},
            {"partitions", nlohmann::json::array({
                {
                    {"name", "default"},
                    {"coverage", nlohmann::json::array()}
                }
            })}
        };
        for (size_t s = 0; s < aln.num_codons; ++s) {
            j["dataset"]["partitions"][0]["coverage"].push_back(s);
        }

        // Phylogeny — shape pinned by fubar.v3.schema.json
        j["phylogeny"] = {
            {"tree", tree.to_newick(true)},
            {"total_tree_length", res.baseline_tree_length},
            {"nodes_count", tree.num_nodes()}
        };

        // Model fits
        if (res.has_gtr_fit) {
            j["model_fits"]["nucleotide_gtr"] = {
                {"log_likelihood", finite_or_null(res.gtr_log_l)},
                {"aicc", finite_or_null(res.gtr_aicc)},
                {"rates", {
                    {"theta_AC", finite_or_null(res.gtr_rates.theta_AC)},
                    {"theta_AT", finite_or_null(res.gtr_rates.theta_AT)},
                    {"theta_CG", finite_or_null(res.gtr_rates.theta_CG)},
                    {"theta_CT", finite_or_null(res.gtr_rates.theta_CT)},
                    {"theta_GT", finite_or_null(res.gtr_rates.theta_GT)}
                }},
                {"branch_lengths", res.gtr_branch_lengths}
            };
        }

        // Rate grid
        j["rate_grid"] = {
            {"dimensions", {settings.grid_points, settings.grid_points}},
            {"total_points", res.grid.size()},
            {"prior_weight_non_positive", res.weight_non_positive},
            {"best_scaler_alpha", res.best_scaler_alpha},
            {"best_scaler_beta", res.best_scaler_beta},
            {"data", {
                {"alpha", nlohmann::json::array()},
                {"beta", nlohmann::json::array()},
                {"posterior_weight", nlohmann::json::array()}
            }}
        };
        for (size_t g = 0; g < res.grid.size(); ++g) {
            j["rate_grid"]["data"]["alpha"].push_back(res.grid[g].first);
            j["rate_grid"]["data"]["beta"].push_back(res.grid[g].second);
            j["rate_grid"]["data"]["posterior_weight"].push_back(res.posterior_grid_weights[g]);
        }

        // Summary
        j["summary"] = {
            {"positive_sites_count", res.positive_sites.size()},
            {"purifying_sites_count", res.purifying_sites.size()},
            {"posterior_threshold", settings.posterior_threshold},
            {"expected_false_positives", res.expected_false_positives},
            {"false_positives_ci_95", {res.false_positives_ci_95.first, res.false_positives_ci_95.second}},
            {"positive_sites", res.positive_sites},
            {"purifying_sites", res.purifying_sites}
        };

        // Site results (columnar; ColumnarResults enforces equal-length)
        ColumnarResults table;
        std::vector<int> col_site;
        std::vector<double> col_alpha, col_beta, col_bma, col_pneg, col_ppos, col_bf, col_psrf, col_neff;
        std::vector<std::string> col_class;
        const size_t n = res.site_results.size();
        col_site.reserve(n); col_alpha.reserve(n); col_beta.reserve(n);
        col_bma.reserve(n); col_pneg.reserve(n); col_ppos.reserve(n);
        col_bf.reserve(n); col_psrf.reserve(n); col_neff.reserve(n);
        col_class.reserve(n);
        for (const auto& sr : res.site_results) {
            col_site.push_back(static_cast<int>(sr.site));
            col_alpha.push_back(sr.alpha);
            col_beta.push_back(sr.beta);
            col_bma.push_back(sr.beta_minus_alpha);
            col_pneg.push_back(sr.prob_neg);
            col_ppos.push_back(sr.prob_pos);
            col_bf.push_back(sr.bayes_factor);
            col_psrf.push_back(sr.psrf);
            col_neff.push_back(sr.neff);
            col_class.push_back(sr.classification);
        }
        table.add_column("site", {"integer", "1-based codon site index"}, col_site);
        table.add_column("alpha", {"number", "Mean posterior synonymous rate", "substitutions/site"}, col_alpha);
        table.add_column("beta", {"number", "Mean posterior non-synonymous rate", "substitutions/site"}, col_beta);
        table.add_column("beta_minus_alpha", {"number", "Mean posterior beta - alpha", "substitutions/site"}, col_bma);
        table.add_column("prob_neg", {"number", "Posterior probability of negative selection P(alpha > beta | Data)"}, col_pneg);
        table.add_column("prob_pos", {"number", "Posterior probability of positive selection P(beta > alpha | Data)"}, col_ppos);
        table.add_column("bayes_factor", {"number", "Empirical Bayes Factor for positive selection"}, col_bf);
        table.add_column("psrf", {"number", "Potential Scale Reduction Factor (Gelman-Rubin R-hat)"}, col_psrf);
        table.add_column("neff", {"number", "Effective sample size (N_eff)"}, col_neff);
        table.add_column("classification", {"string", "Selection classification: diversifying, purifying, or neutral"}, col_class);
        j["site_results"] = table.to_json();

        return j.dump(2);
    }

std::string FUBARAnalyzer::to_legacy_json(const FUBARResult& res, const std::string& input_filepath) const {
        nlohmann::json j;

        j["analysis"] = {
            {"info", "Perform a Fast Unbiased AppRoximate Bayesian (FUBAR) analysis of a coding sequence alignment."},
            {"version", "2.2"},
            {"reference", "FUBAR: a fast, unconstrained bayesian approximation for inferring selection (2013), Mol Biol Evol. 30(5):1196-205"},
            {"authors", "Sergei L Kosakovsky Pond"},
            {"contact", "spond@temple.edu"}
        };

        j["input"] = {
            {"file", input_filepath},
            {"sequences", aln.num_taxa},
            {"sites", aln.num_codons},
            {"partition count", 1}
        };

        // GTR Fit
        j["fits"] = nlohmann::json::object();
        if (res.has_gtr_fit) {
            j["fits"]["Nucleotide GTR"] = {
                {"Log Likelihood", res.gtr_log_l},
                {"AIC-c", res.gtr_aicc},
                {"Rate Distributions", {
                    {"Substitution rate from nucleotide A to nucleotide C", res.gtr_rates.theta_AC},
                    {"Substitution rate from nucleotide A to nucleotide G", 1.0},
                    {"Substitution rate from nucleotide A to nucleotide T", res.gtr_rates.theta_AT},
                    {"Substitution rate from nucleotide C to nucleotide G", res.gtr_rates.theta_CG},
                    {"Substitution rate from nucleotide C to nucleotide T", res.gtr_rates.theta_CT},
                    {"Substitution rate from nucleotide G to nucleotide T", res.gtr_rates.theta_GT}
                }},
                {"display order", 0},
                {"estimated parameters", 24}
            };
        }

        // Settings
        j["settings"] = {
            {"grid size", settings.grid_points},
            {"chains", settings.mcmc_chains},
            {"chain-length", settings.mcmc_chain_length},
            {"burn-in", settings.mcmc_burn_in},
            {"samples", settings.mcmc_samples},
            {"concentration", settings.concentration},
            {"posterior", settings.posterior_threshold},
            {"non-zero", settings.non_zero ? 1 : 0},
            {"method", fubar_method_to_string(settings.method)}
        };

        // Grid: [G x 3] matrix
        j["grid"] = nlohmann::json::array();
        for (size_t g = 0; g < res.grid.size(); ++g) {
            j["grid"].push_back({res.grid[g].first, res.grid[g].second, res.posterior_grid_weights[g]});
        }

        // Posterior: {"0": {"0": [weights...], ...}}
        j["posterior"] = nlohmann::json::object();
        j["posterior"]["0"] = nlohmann::json::object();
        for (size_t s = 0; s < res.site_results.size(); ++s) {
            j["posterior"]["0"][std::to_string(s)] = nlohmann::json::array({res.site_results[s].posterior_weights});
        }

        // MLE headers and content
        j["MLE"] = {
            {"headers", nlohmann::json::array({
                nlohmann::json::array({"alpha", "Mean posterior synonymous substitution rate at a site"}),
                nlohmann::json::array({"beta", "Mean posterior non-synonymous substitution rate at a site"}),
                nlohmann::json::array({"beta-alpha", "Mean posterior beta-alpha"}),
                nlohmann::json::array({"Prob[alpha>beta]", "Posterior probability of negative selection at a site"}),
                nlohmann::json::array({"Prob[alpha<beta]", "Posterior probability of positive selection at a site"}),
                nlohmann::json::array({"BayesFactor[alpha<beta]", "Empirical Bayes Factor for positive selection at a site"}),
                nlohmann::json::array({"PSRF", "Potential scale reduction factor - an MCMC mixing measure"}),
                nlohmann::json::array({"Neff", "Estimated effective sample size for Prob[alpha<beta]"})
            })},
            {"content", {
                {"0", nlohmann::json::array()}
            }}
        };

        for (const auto& sr : res.site_results) {
            j["MLE"]["content"]["0"].push_back({
                sr.alpha,
                sr.beta,
                sr.beta_minus_alpha,
                sr.prob_neg,
                sr.prob_pos,
                sr.bayes_factor,
                sr.psrf,
                sr.neff
            });
        }

        // Timers
        j["timers"] = {
            {"Nucleotide Fit", {{"timer", res.time_gtr_sec}}},
            {"Grid Calculations", {{"timer", res.time_grid_sec}}},
            {"Posterior estimation", {{"timer", res.time_posterior_sec}}},
            {"Overall", {{"timer", res.time_total_sec}}}
        };

        j["runtime"] = "3.0.0";
        return j.dump(2);
    }

} // namespace hyphy::analyses
