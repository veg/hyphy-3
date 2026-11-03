// Out-of-line JSON emitter definitions for slac.
// Moved from include/hyphy/analyses/slac.hpp so that emitter
// edits recompile this single TU instead of every consumer
// (tests, apps, Python bindings).
#include "hyphy/analyses/slac.hpp"
#include "hyphy/core/json_output.hpp"
#include "hyphy/core/provenance.hpp"

namespace hyphy::analyses {

using namespace hyphy::core;
using namespace hyphy::opt;

nlohmann::json SLACResult::to_legacy_json(
        const std::string& input_filepath,
        const std::string& tree_string) const {
        nlohmann::json j;

        // Analysis description
        j["analysis"]["authors"] = "Sergei L Kosakovsky Pond and Simon DW Frost";
        j["analysis"]["citation"] = "Not So Different After All: A Comparison of Methods for Detecting Amino Acid Sites Under Selection (2005). _Mol Biol Evol_ 22 (5): 1208-1222";
        j["analysis"]["contact"] = "spond@temple.edu";
        j["analysis"]["info"] = "SLAC (Single Likelihood Ancestor Counting) uses a maximum likelihood ancestral state reconstruction and minimum path substitution counting to estimate site-level dS and dN, and applies a simple binomial-based test to test if dS differs from dN.";
        j["analysis"]["requirements"] = "in-frame codon alignment and a phylogenetic tree";
        j["analysis"]["settings"]["branches"] = settings.test_branches;
        j["analysis"]["settings"]["pvalue"] = settings.p_value_threshold;
        j["analysis"]["settings"]["samples"] = settings.samples;
        j["analysis"]["version"] = "2.00";

        // Input
        j["input"]["file name"] = input_filepath;
        j["input"]["genetic code"] = aln.code ? aln.code->name : "Universal";
        j["input"]["number of sequences"] = aln.num_taxa;
        j["input"]["number of sites"] = aln.num_codons;
        j["input"]["partition count"] = 1;
        if (!tree_string.empty()) {
            j["input"]["trees"]["0"] = tree_string;
        }

        // Tested branches
        nlohmann::json tested_obj;
        for (const auto& b : tested_branch_names) {
            tested_obj[b] = "test";
        }
        for (const auto& b : background_branch_names) {
            tested_obj[b] = "background";
        }
        j["tested"]["0"] = tested_obj;

        // Data partitions
        nlohmann::json cov = nlohmann::json::array();
        for (size_t i = 0; i < total_sites; ++i) {
            cov.push_back(i);
        }
        j["data partitions"]["0"]["coverage"] = nlohmann::json::array({cov});
        j["data partitions"]["0"]["name"] = "slac.filter.default";

        // Fits
        // GTR
        j["fits"]["Nucleotide GTR"]["Log Likelihood"] = gtr_log_l;
        j["fits"]["Nucleotide GTR"]["AIC-c"] = gtr_aicc;
        j["fits"]["Nucleotide GTR"]["display order"] = 0;
        nlohmann::json gtr_ef = nlohmann::json::array();
        for (int i = 0; i < 4; ++i) {
            gtr_ef.push_back({aln.nuc_frequencies(i)});
        }
        j["fits"]["Nucleotide GTR"]["Equilibrium frequencies"] = gtr_ef;
        j["fits"]["Nucleotide GTR"]["Rate Distributions"] = {
            {"Substitution rate from nucleotide A to nucleotide C", gtr_rates.theta_AC},
            {"Substitution rate from nucleotide A to nucleotide G", 1.0},
            {"Substitution rate from nucleotide A to nucleotide T", gtr_rates.theta_AT},
            {"Substitution rate from nucleotide C to nucleotide G", gtr_rates.theta_CG},
            {"Substitution rate from nucleotide C to nucleotide T", gtr_rates.theta_CT},
            {"Substitution rate from nucleotide G to nucleotide T", gtr_rates.theta_GT}
        };

        // MG94
        j["fits"]["Global MG94xREV"]["Log Likelihood"] = baseline_log_l;
        j["fits"]["Global MG94xREV"]["AIC-c"] = baseline_aicc;
        j["fits"]["Global MG94xREV"]["display order"] = 1;
        nlohmann::json mg_ef = nlohmann::json::array();
        for (int i = 0; i < aln.codon_frequencies_f3x4.size(); ++i) {
            mg_ef.push_back({aln.codon_frequencies_f3x4(i)});
        }
        j["fits"]["Global MG94xREV"]["Equilibrium frequencies"] = mg_ef;
        j["fits"]["Global MG94xREV"]["Rate Distributions"]["non-synonymous/synonymous rate ratio for *background*"] = global_omega;

        // Branch attributes
        nlohmann::json branch_attr;
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                branch_attr[node.name]["Global MG94xREV"] = node.branch_length;
                branch_attr[node.name]["Nucleotide GTR"] = gtr_branch_lengths.count(node.name) ? gtr_branch_lengths.at(node.name) : node.branch_length;
                branch_attr[node.name]["original name"] = node.name;

                if (branch_site_attrs.count(node.name)) {
                    const auto& bsa = branch_site_attrs.at(node.name);
                    branch_attr[node.name]["codon"] = nlohmann::json::array({bsa.codons});
                    branch_attr[node.name]["amino-acid"] = nlohmann::json::array({bsa.amino_acids});
                    branch_attr[node.name]["synonymous substitution count"] = nlohmann::json::array({bsa.syn_subs});
                    branch_attr[node.name]["nonsynonymous substitution count"] = nlohmann::json::array({bsa.nonsyn_subs});
                }
            }
        }
        j["branch attributes"]["0"] = branch_attr;
        j["branch attributes"]["attributes"]["Global MG94xREV"]["attribute type"] = "branch length";
        j["branch attributes"]["attributes"]["Global MG94xREV"]["display order"] = 1;
        j["branch attributes"]["attributes"]["Nucleotide GTR"]["attribute type"] = "branch length";
        j["branch attributes"]["attributes"]["Nucleotide GTR"]["display order"] = 0;
        j["branch attributes"]["attributes"]["original name"]["attribute type"] = "node label";
        j["branch attributes"]["attributes"]["original name"]["display order"] = -1;
        j["branch attributes"]["attributes"]["codon"]["attribute type"] = "node label";
        j["branch attributes"]["attributes"]["codon"]["display order"] = 0;
        j["branch attributes"]["attributes"]["amino-acid"]["attribute type"] = "node label";
        j["branch attributes"]["attributes"]["amino-acid"]["display order"] = 1;
        j["branch attributes"]["attributes"]["synonymous substitution count"]["attribute type"] = "branch label";
        j["branch attributes"]["attributes"]["synonymous substitution count"]["display order"] = 0;
        j["branch attributes"]["attributes"]["nonsynonymous substitution count"]["attribute type"] = "branch label";
        j["branch attributes"]["attributes"]["nonsynonymous substitution count"]["display order"] = 1;

        // MLE headers & content
        j["MLE"]["headers"] = nlohmann::json::array({
            nlohmann::json::array({"ES", "Expected synonymous sites"}),
            nlohmann::json::array({"EN", "Expected non-synonymous sites"}),
            nlohmann::json::array({"S", "Inferred synonymous substitutions"}),
            nlohmann::json::array({"N", "Inferred non-synonymous substitutions"}),
            nlohmann::json::array({"P[S]", "Expected proportion of synonymous sites"}),
            nlohmann::json::array({"dS", "Inferred synonymous susbsitution rate"}),
            nlohmann::json::array({"dN", "Inferred non-synonymous susbsitution rate"}),
            nlohmann::json::array({"dN-dS", "Scaled by the length of the tested branches"}),
            nlohmann::json::array({"P [dN/dS > 1]", "Binomial probability that S is no greater than the observed value, with P<sub>s</sub> probability of success"}),
            nlohmann::json::array({"P [dN/dS < 1]", "Binomial probability that S is no less than the observed value, with P<sub>s</sub> probability of success"}),
            nlohmann::json::array({"Total branch length", "The total length of branches contributing to inference at this site, and used to scale dN-dS"})
        });

        nlohmann::json res_by_site = nlohmann::json::array();
        nlohmann::json avg_by_site = nlohmann::json::array();
        for (size_t s = 0; s < total_sites; ++s) {
            res_by_site.push_back(resolved_sites[s].to_legacy_row());
            avg_by_site.push_back(averaged_sites[s].to_legacy_row());
        }
        j["MLE"]["content"]["0"]["by-site"]["RESOLVED"] = res_by_site;
        j["MLE"]["content"]["0"]["by-site"]["AVERAGED"] = avg_by_site;

        nlohmann::json b_names = nlohmann::json::array();
        nlohmann::json res_by_branch = nlohmann::json::array();
        nlohmann::json avg_by_branch = nlohmann::json::array();
        for (size_t b = 0; b < resolved_branches.size(); ++b) {
            b_names.push_back(nlohmann::json::array({resolved_branches[b].branch_name}));
            res_by_branch.push_back(resolved_branches[b].to_legacy_row());
            avg_by_branch.push_back(averaged_branches[b].to_legacy_row());
        }
        j["MLE"]["content"]["0"]["by-branch"]["NAMES"] = b_names;
        j["MLE"]["content"]["0"]["by-branch"]["RESOLVED"] = res_by_branch;
        j["MLE"]["content"]["0"]["by-branch"]["AVERAGED"] = avg_by_branch;

        // Ancestral sample quantiles if available
        if (settings.samples > 0 && !resolved_sites.empty() && !resolved_sites[0].quantiles_median.empty()) {
            nlohmann::json med_res = nlohmann::json::array();
            nlohmann::json ci25_res = nlohmann::json::array();
            nlohmann::json ci975_res = nlohmann::json::array();
            nlohmann::json med_avg = nlohmann::json::array();
            nlohmann::json ci25_avg = nlohmann::json::array();
            nlohmann::json ci975_avg = nlohmann::json::array();

            for (size_t s = 0; s < total_sites; ++s) {
                med_res.push_back(resolved_sites[s].quantiles_median);
                ci25_res.push_back(resolved_sites[s].quantiles_2_5);
                ci975_res.push_back(resolved_sites[s].quantiles_97_5);

                med_avg.push_back(averaged_sites[s].quantiles_median);
                ci25_avg.push_back(averaged_sites[s].quantiles_2_5);
                ci975_avg.push_back(averaged_sites[s].quantiles_97_5);
            }

            j["sample-median"]["0"]["RESOLVED"] = med_res;
            j["sample-median"]["0"]["AVERAGED"] = med_avg;
            j["sample-2.5"]["0"]["RESOLVED"] = ci25_res;
            j["sample-2.5"]["0"]["AVERAGED"] = ci25_avg;
            j["sample-97.5"]["0"]["RESOLVED"] = ci975_res;
            j["sample-97.5"]["0"]["AVERAGED"] = ci975_avg;
        }

        j["timers"]["Total time"] = runtime_seconds;
        j["runtime"] = "3.0.0";

        return j;
    }

nlohmann::json SLACResult::to_modern_json(
        const std::string& input_filepath,
        const std::string& tree_string,
        const Provenance& prov) const {
        AnalysisMetadata meta{
            "slac",
            "Single-Likelihood Ancestor Counting",
            "3.0.0",
            "site_selection",
            "SLAC uses maximum likelihood ancestral state reconstruction and minimum path substitution counting to estimate site-level dS and dN, with exact binomial hypothesis testing.",
            {citation("Kosakovsky Pond SL, Frost SDW (2005). Not So Different After All: A Comparison of Methods for Detecting Amino Acid Sites Under Selection. Mol Biol Evol 22(5): 1208-1222.",
                      "10.1093/molbev/msi105", "15703277")},
            {
                {"pvalue_threshold", settings.p_value_threshold},
                {"samples", settings.samples},
                {"test_branches", settings.test_branches}
            }
        };

        nlohmann::json j = new_analysis_json(meta, prov, input_filepath);
        if (prov.invocation.cli_command.empty() && prov.inputs.empty() &&
            j["provenance"].contains("execution")) {
            j["provenance"]["execution"]["wall_time_seconds"] = runtime_seconds;
        }
        j["dataset"] = dataset_block(aln);
        j["phylogeny"] = phylogeny_block(
            tree, tree_string,
            [](const TreeNode& node) {
                return std::map<std::string, Scalar>{{"mg94_baseline", node.branch_length}};
            },
            [this](const TreeNode& node) -> nlohmann::json {
                if (node.id == tree.root_id) return nlohmann::json::object();
                const bool is_test = std::find(
                    tested_branch_names.begin(), tested_branch_names.end(),
                    node.name) != tested_branch_names.end();
                return {{"length", finite_or_null(node.branch_length)},
                        {"tested", is_test}};
            });
        // SLAC also reports the legacy scalar summary for downstream parity.
        j["phylogeny"]["tree_length"] = finite_or_null(total_tree_length);

        // Model fits
        if (gtr_log_l != 0.0 || gtr_aicc != 0.0) {
            j["model_fits"]["Nucleotide GTR"] = nucleotide_gtr_fit(
                aln, gtr_rates, gtr_log_l, gtr_aicc,
                5 + static_cast<int>(tree.num_nodes()) - 1);
        }
        nlohmann::json mg94 = fit_entry(
            "Muse-Gaut 94 codon model with GTR core",
            baseline_log_l, baseline_aicc,
            6 + static_cast<int>(tree.num_nodes()) - 1);
        mg94["global_parameters"] = {{"omega", finite_or_null(global_omega)}};
        j["model_fits"]["Global MG94xREV"] = mg94;

        // Statistical tests
        j["statistical_tests"]["site_level_summary"] = {
            {"method", "Single-Likelihood Ancestor Counting (SLAC)"},
            {"test_statistic", "Exact Binomial Distribution"},
            {"threshold", settings.p_value_threshold},
            {"total_sites", static_cast<int>(total_sites)},
            {"diversifying_count", static_cast<int>(diversifying_sites)},
            {"purifying_count", static_cast<int>(purifying_sites)},
            {"neutral_count", static_cast<int>(neutral_sites)},
            {"tested_branches_count", static_cast<int>(num_tested_branches)},
            {"background_branches_count", static_cast<int>(num_background_branches)},
            {"total_tested_tree_length", finite_or_null(total_tested_tree_length)}
        };

        // Site results (columnar; ColumnarResults enforces equal-length).
        // dS/dN/dN-dS are nullable: sites with zero nonsynonymous expectation
        // have no rate estimate, so the column carries null per upstream.
        ColumnarResults table;
        std::vector<int> col_site;
        std::vector<double> col_ES, col_EN, col_S, col_N, col_P_S, col_p_pos, col_p_neg, col_tbl;
        nlohmann::json col_dS = nlohmann::json::array(),
                       col_dN = nlohmann::json::array(),
                       col_diff = nlohmann::json::array();
        std::vector<std::string> col_class;

        for (size_t s = 0; s < total_sites; ++s) {
            const auto& r = resolved_sites[s];
            col_site.push_back(static_cast<int>(s + 1));
            col_ES.push_back(r.ES);
            col_EN.push_back(r.EN);
            col_S.push_back(r.S);
            col_N.push_back(r.N);
            col_P_S.push_back(r.P_S);
            col_dS.push_back(r.has_dS ? nlohmann::json(r.dS) : nlohmann::json(nullptr));
            col_dN.push_back(r.has_dN ? nlohmann::json(r.dN) : nlohmann::json(nullptr));
            col_diff.push_back(r.has_dN_minus_dS ? nlohmann::json(r.dN_minus_dS) : nlohmann::json(nullptr));
            col_p_pos.push_back(r.p_positive);
            col_p_neg.push_back(r.p_negative);
            col_tbl.push_back(r.total_branch_length);

            if (r.p_positive <= settings.p_value_threshold) {
                col_class.push_back("Diversifying");
            } else if (r.p_negative <= settings.p_value_threshold) {
                col_class.push_back("Purifying");
            } else {
                col_class.push_back("Neutral");
            }
        }

        table.add_column("site", {"integer", "1-based codon position", "1-based codon position"}, col_site);
        table.add_column("ES", {"number", "Expected synonymous sites"}, col_ES);
        table.add_column("EN", {"number", "Expected non-synonymous sites"}, col_EN);
        table.add_column("S", {"number", "Inferred synonymous substitutions"}, col_S);
        table.add_column("N", {"number", "Inferred non-synonymous substitutions"}, col_N);
        table.add_column("P_S", {"number", "Expected proportion of synonymous sites"}, col_P_S);
        table.add_column("dS", {"number", "Inferred synonymous substitution rate"}, col_dS);
        table.add_column("dN", {"number", "Inferred non-synonymous substitution rate"}, col_dN);
        table.add_column("dN_minus_dS", {"number", "Scaled substitution rate difference"}, col_diff);
        table.add_column("p_positive", {"number", "Binomial p-value for positive (diversifying) selection"}, col_p_pos);
        table.add_column("p_negative", {"number", "Binomial p-value for negative (purifying) selection"}, col_p_neg);
        table.add_column("total_branch_length", {"number", "Effective branch length contributing to inference"}, col_tbl);
        table.add_column("classification", {"string", "Selection status: Diversifying, Purifying, or Neutral"}, col_class);

        if (settings.samples > 0 && !resolved_sites.empty() && !resolved_sites[0].quantiles_median.empty()) {
            std::vector<double> S_med, S_lo, S_hi, N_med, N_lo, N_hi;
            nlohmann::json dS_med = nlohmann::json::array(), dS_lo = nlohmann::json::array(),
                           dS_hi = nlohmann::json::array(), dN_med = nlohmann::json::array(),
                           dN_lo = nlohmann::json::array(), dN_hi = nlohmann::json::array();

            for (size_t s = 0; s < total_sites; ++s) {
                const auto& r = resolved_sites[s];
                S_med.push_back(r.quantiles_median[2]);
                S_lo.push_back(r.quantiles_2_5[2]);
                S_hi.push_back(r.quantiles_97_5[2]);
                N_med.push_back(r.quantiles_median[3]);
                N_lo.push_back(r.quantiles_2_5[3]);
                N_hi.push_back(r.quantiles_97_5[3]);
                dS_med.push_back(r.has_dS ? nlohmann::json(r.quantiles_median[5]) : nlohmann::json(nullptr));
                dS_lo.push_back(r.has_dS ? nlohmann::json(r.quantiles_2_5[5]) : nlohmann::json(nullptr));
                dS_hi.push_back(r.has_dS ? nlohmann::json(r.quantiles_97_5[5]) : nlohmann::json(nullptr));
                dN_med.push_back(r.has_dN ? nlohmann::json(r.quantiles_median[6]) : nlohmann::json(nullptr));
                dN_lo.push_back(r.has_dN ? nlohmann::json(r.quantiles_2_5[6]) : nlohmann::json(nullptr));
                dN_hi.push_back(r.has_dN ? nlohmann::json(r.quantiles_97_5[6]) : nlohmann::json(nullptr));
            }

            table.add_column("S_median", {"number", "Median synonymous substitutions across samples"}, S_med);
            table.add_column("S_ci_lower", {"number", "2.5% CI for synonymous substitutions"}, S_lo);
            table.add_column("S_ci_upper", {"number", "97.5% CI for synonymous substitutions"}, S_hi);
            table.add_column("N_median", {"number", "Median non-synonymous substitutions across samples"}, N_med);
            table.add_column("N_ci_lower", {"number", "2.5% CI for non-synonymous substitutions"}, N_lo);
            table.add_column("N_ci_upper", {"number", "97.5% CI for non-synonymous substitutions"}, N_hi);
            table.add_column("dS_median", {"number", "Median dS across samples"}, dS_med);
            table.add_column("dS_ci_lower", {"number", "2.5% CI for dS"}, dS_lo);
            table.add_column("dS_ci_upper", {"number", "97.5% CI for dS"}, dS_hi);
            table.add_column("dN_median", {"number", "Median dN across samples"}, dN_med);
            table.add_column("dN_ci_lower", {"number", "2.5% CI for dN"}, dN_lo);
            table.add_column("dN_ci_upper", {"number", "97.5% CI for dN"}, dN_hi);
        }
        j["site_results"] = table.to_json();

        return j;
    }

} // namespace hyphy::analyses
