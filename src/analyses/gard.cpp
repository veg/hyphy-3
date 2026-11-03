// Out-of-line JSON emitter definitions for gard.
// Moved from include/hyphy/analyses/gard.hpp so that emitter
// edits recompile this single TU instead of every consumer
// (tests, apps, Python bindings).
#include "hyphy/analyses/gard.hpp"
#include "hyphy/core/json_output.hpp"
#include "hyphy/core/provenance.hpp"

namespace hyphy::analyses {

using namespace hyphy::core;
using namespace hyphy::opt;

nlohmann::json GARDResult::to_json() const {
        nlohmann::json j;
        j["$schema"] = "https://raw.githubusercontent.com/veg/hyphy-3/main/schemas/v3/gard.v3.schema.json";
        j["schema_version"] = "3.0.0";
        j["analysis"] = {
            {"id", "gard"},
            {"name", "GARD"},
            {"version", "3.0.0"},
            {"citation", "Mol Biol Evol 23(10), 1891–1901"},
            {"authors", "Sergei L Kosakovsky Pond"},
            {"contact", "spond@temple.edu"}
        };
        j["software"] = software_block();

        // Provenance — to_json() carries no cli_cmd, so build the block
        // inline with a fixed invocation label.
        {
            Provenance prov;
            prov.invocation.cli_command = "hyphy3 gard";
            prov.invocation.working_directory = Provenance::get_cwd();
            if (!input_file.empty()) {
                size_t sz = 0;
                std::string hash = crypto::SHA256::hash_file(input_file, &sz);
                if (!hash.empty()) {
                    prov.inputs["alignment"] = {input_file, "FASTA/NEXUS", hash, sz};
                }
            }
            prov.execution.start_time = Provenance::current_iso8601();
            prov.execution.end_time = prov.execution.start_time;
            prov.execution.wall_time_seconds = time_elapsed;
            prov.execution.cpu_threads = 1;
            prov.execution.hostname = Provenance::get_hostname();
            prov.execution.os = Provenance::detect_os();
            prov.execution.compiler = Provenance::detect_compiler();
            j["provenance"] = prov.to_json();
        }

        // GARD's schema pins this dataset shape (sites_count rather than the
        // shared codon_sites/unique_patterns contract), so it stays custom.
        j["dataset"] = {
            {"file", input_file},
            {"taxa_count", static_cast<int>(num_taxa)},
            {"sites_count", static_cast<int>(num_sites)},
            {"variable_sites_count", static_cast<int>(variable_sites.size())},
            {"potential_breakpoints", static_cast<int>(potential_breakpoints)},
            {"minimum_partition_size", static_cast<int>(min_partition_size)}
        };
        j["baseline"] = {
            {"log_likelihood", finite_or_null(baseline_log_l)},
            {"aicc", finite_or_null(baseline_aicc)},
            {"tree", baseline_tree.to_newick(true)},
            {"rates", {
                {"AC", finite_or_null(baseline_rates.theta_AC)},
                {"AT", finite_or_null(baseline_rates.theta_AT)},
                {"CG", finite_or_null(baseline_rates.theta_CG)},
                {"CT", finite_or_null(baseline_rates.theta_CT)},
                {"GT", finite_or_null(baseline_rates.theta_GT)}
            }},
            {"frequencies", {
                {"A", finite_or_null(baseline_frequencies(0))},
                {"C", finite_or_null(baseline_frequencies(1))},
                {"G", finite_or_null(baseline_frequencies(2))},
                {"T", finite_or_null(baseline_frequencies(3))}
            }}
        };

        j["rate_variation"] = {
            {"mode", (srv_mode == GARDSRVMode::DiscreteGamma) ? "discrete_gamma" : "none"},
            {"classes", (srv_mode == GARDSRVMode::DiscreteGamma) ? static_cast<int>(gamma_rates.size()) : 1},
            {"shape_alpha", finite_or_null(gamma_alpha)},
            {"rates", gamma_rates},
            {"weights", gamma_weights}
        };

        j["tree_engine"] = {
            {"method", (tree_method == GARDTreeMethod::Ohmic) ? "ohmic_spectral" : "neighbor_joining"},
            {"complexity", (tree_method == GARDTreeMethod::Ohmic) ? "O(N log N)" : "O(N^3)"}
        };

        nlohmann::json bridges_arr = nlohmann::json::array();
        for (const auto& b : bridge_summaries) {
            bridges_arr.push_back({
                {"bridge_index", b.bridge_index},
                {"hub_node", b.hub_name},
                {"bridge_node", b.bridge_name},
                {"initial_weights", {b.initial_weights[0], b.initial_weights[1], b.initial_weights[2]}},
                {"fitted_weights", {b.fitted_weights[0], b.fitted_weights[1], b.fitted_weights[2]}},
                {"ohmic_support", finite_or_null(b.ohmic_support)},
                {"shunt_stress", finite_or_null(b.shunt_stress)}
            });
        }

        j["ohmic_network"] = {
            {"bridges_detected", ohmic_bridges_detected},
            {"num_bridges", static_cast<int>(num_bridges)},
            {"ohmic_log_likelihood", finite_or_null(ohmic_log_l)},
            {"ohmic_aicc", finite_or_null(ohmic_aicc)},
            {"delta_aicc_vs_baseline", finite_or_null(ohmic_delta_aicc)},
            {"bridge_candidate_breakpoints", bridge_candidate_breakpoints},
            {"bridges", bridges_arr}
        };

        nlohmann::json improv_arr = nlohmann::json::array();
        for (size_t i = 0; i < improvements.size(); ++i) {
            const auto& [bps, delta] = improvements[i];
            const auto& step_bps = (i < improvement_breakpoints.size()) ? improvement_breakpoints[i] : best_model.breakpoints;
            improv_arr.push_back({
                {"breakpoints_count", bps},
                {"breakpoints", step_bps},
                {"delta_aicc", delta}
            });
        }

        j["model_selection"] = {
            {"best_model_breakpoints", best_model.breakpoints},
            {"best_model_aicc", finite_or_null(best_model.aicc)},
            {"best_model_log_likelihood", finite_or_null(best_model.log_likelihood)},
            {"improvements", improv_arr}
        };

        nlohmann::json part_arr = nlohmann::json::array();
        for (size_t p = 0; p < best_model.partitions.size(); ++p) {
            const auto& part = best_model.partitions[p];
            part_arr.push_back({
                {"partition_index", p},
                {"start_site", part.start_site},
                {"end_site", part.end_site},
                {"log_likelihood", finite_or_null(part.log_likelihood)},
                {"tree", part.tree.to_newick(true)}
            });
        }
        j["partitions"] = part_arr;

        nlohmann::json bp_tests_arr = nlohmann::json::array();
        for (const auto& t : topology_test_result.breakpoint_tests) {
            bp_tests_arr.push_back({
                {"breakpoint", t.breakpoint_site},
                {"left_partition", {t.left_interval.first, t.left_interval.second}},
                {"right_partition", {t.right_interval.first, t.right_interval.second}},
                {"delta_log_l_left", finite_or_null(t.delta_lnl_left)},
                {"delta_log_l_right", finite_or_null(t.delta_lnl_right)},
                {"kh_p_value_left", finite_or_null(t.kh_p_value_left)},
                {"kh_p_value_right", finite_or_null(t.kh_p_value_right)},
                {"sh_p_value_left", finite_or_null(t.sh_p_value_left)},
                {"sh_p_value_right", finite_or_null(t.sh_p_value_right)},
                {"statistically_significant", t.statistically_significant}
            });
        }

        nlohmann::json part_matrix_arr = nlohmann::json::array();
        for (const auto& pm : topology_test_result.partition_matrix) {
            part_matrix_arr.push_back({
                {"partition_index", pm.partition_index},
                {"interval", {pm.interval.first, pm.interval.second}},
                {"best_tree_index", pm.best_tree_index},
                {"tree_log_likelihoods", pm.tree_log_likelihoods},
                {"kh_p_values", pm.kh_p_values},
                {"sh_p_values", pm.sh_p_values}
            });
        }

        j["topology_test"] = {
            {"single_tree_aicc", finite_or_null(single_tree_aicc)},
            {"topological_incongruence_supported", topological_incongruence},
            {"delta_aicc_topology", finite_or_null(single_tree_aicc - best_model.aicc)},
            {"test_type", topology_test_result.test_type},
            {"bootstrap_replicates", topology_test_result.bootstrap_replicates},
            {"any_incongruence_significant", topology_test_result.any_incongruence_significant},
            {"breakpoints", bp_tests_arr},
            {"partition_tree_matrix", part_matrix_arr}
        };

        j["site_breakpoint_support"] = site_breakpoint_support;
        j["total_models_evaluated"] = total_models_evaluated;
        j["time_elapsed"] = time_elapsed;

        return j;
    }

nlohmann::json GARDResult::to_legacy_json() const {
        nlohmann::json j;
        j["analysis"] = {
            {"authors", "Sergei L Kosakovsky Pond"},
            {"citation", "**Automated Phylogenetic Detection of Recombination Using a Genetic Algorithm**, _Mol Biol Evol 23(10), 1891–1901"},
            {"contact", "spond@temple.edu"},
            {"info", "GARD : Genetic Algorithms for Recombination Detection. Implements a heuristic approach to screening alignments of sequences for recombination, by using the CHC genetic algorithm to search for phylogenetic incongruence among different partitions of the data."},
            {"requirements", "A sequence alignment."},
            {"version", "0.2"}
        };

        j["baselineScore"] = baseline_aicc;
        j["bestModelAICc"] = best_model.aicc;

        nlohmann::json bp_data = nlohmann::json::object();
        nlohmann::json trees_obj = nlohmann::json::object();

        if (best_model.breakpoints.empty()) {
            bp_data["0"] = {
                {"bps", nlohmann::json::array({{std::to_string(1), std::to_string(num_sites)}})},
                {"tree", baseline_tree.to_newick(true)}
            };
            trees_obj["newickString"] = baseline_tree.to_newick(true);
        } else {
            for (size_t p = 0; p < best_model.partitions.size(); ++p) {
                const auto& part = best_model.partitions[p];
                bp_data[std::to_string(p)] = {
                    {"bps", nlohmann::json::array({{std::to_string(part.start_site), std::to_string(part.end_site)}})},
                    {"tree", part.tree.to_newick(true)}
                };
                trees_obj[std::to_string(p)] = {
                    {"newickString", part.tree.to_newick(true)}
                };
            }
        }
        j["breakpointData"] = bp_data;
        j["trees"] = trees_obj;

        nlohmann::json improv_obj = nlohmann::json::object();
        if (improvements.empty() || best_model.breakpoints.empty()) {
            improv_obj["0"] = {
                {"breakpoints", nullptr},
                {"deltaAICc", 0.0}
            };
        } else {
            for (size_t i = 0; i < improvements.size(); ++i) {
                const auto& step_bps = (i < improvement_breakpoints.size()) ? improvement_breakpoints[i] : best_model.breakpoints;
                improv_obj[std::to_string(i)] = {
                    {"breakpoints", step_bps},
                    {"deltaAICc", improvements[i].second}
                };
            }
        }
        j["improvements"] = improv_obj;

        j["input"] = {
            {"file name", input_file},
            {"number of sequences", num_taxa},
            {"number of sites", num_sites},
            {"partition count", best_model.partitions.size()}
        };
        j["potentialBreakpoints"] = potential_breakpoints;
        j["singleTreeAICc"] = single_tree_aicc;
        j["siteBreakPointSupport"] = site_breakpoint_support;
        j["timeElapsed"] = static_cast<size_t>(time_elapsed + 0.5);
        j["totalModelCount"] = total_models_evaluated;

        nlohmann::json legacy_bp_tests = nlohmann::json::array();
        for (const auto& t : topology_test_result.breakpoint_tests) {
            legacy_bp_tests.push_back({
                {"breakpoint", t.breakpoint_site},
                {"left_partition", {t.left_interval.first, t.left_interval.second}},
                {"right_partition", {t.right_interval.first, t.right_interval.second}},
                {"delta_log_l_left", t.delta_lnl_left},
                {"delta_log_l_right", t.delta_lnl_right},
                {"kh_p_value_left", t.kh_p_value_left},
                {"kh_p_value_right", t.kh_p_value_right},
                {"sh_p_value_left", t.sh_p_value_left},
                {"sh_p_value_right", t.sh_p_value_right},
                {"statistically_significant", t.statistically_significant}
            });
        }
        j["topologyTests"] = {
            {"testType", topology_test_result.test_type},
            {"replicates", topology_test_result.bootstrap_replicates},
            {"anySignificant", topology_test_result.any_incongruence_significant},
            {"breakpoints", legacy_bp_tests}
        };

        return j;
    }

} // namespace hyphy::analyses
