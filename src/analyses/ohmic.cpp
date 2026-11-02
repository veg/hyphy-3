// Out-of-line JSON emitter definitions for ohmic.
// Moved from include/hyphy/analyses/ohmic.hpp so that emitter
// edits recompile this single TU instead of every consumer
// (tests, apps, Python bindings).
#include "hyphy/analyses/ohmic.hpp"
#include "hyphy/core/json_output.hpp"
#include "hyphy/core/provenance.hpp"

namespace hyphy::analyses {

using namespace hyphy::core;
using namespace hyphy::opt;

nlohmann::json OhmicResult::to_modern_json(const std::string& filepath) const {
        nlohmann::json j;
        j["$schema"] = "https://raw.githubusercontent.com/veg/hyphy-3/main/schemas/v3/ohmic.v3.schema.json";
        j["schema_version"] = "3.0.0";
        j["analysis"] = {
            {"id", "ohmic"},
            {"name", "OHMIC"},
            {"description", "Ohmic Phylogenetic Networks & Linear-Time Wheatstone Pruning"},
            {"version", "3.0.0"},
            {"authors", "Sergei L Kosakovsky Pond"},
            {"contact", "spond@temple.edu"}
        };
        j["software"] = software_block();

        // Provenance — to_modern_json carries no cli_cmd, so build the
        // block inline with a fixed invocation label.
        {
            const std::string in_path = filepath.empty() ? settings.input_file : filepath;
            Provenance prov;
            prov.invocation.cli_command = "hyphy3 ohmic";
            prov.invocation.working_directory = Provenance::get_cwd();
            if (!in_path.empty()) {
                size_t sz = 0;
                std::string hash = crypto::SHA256::hash_file(in_path, &sz);
                if (!hash.empty()) {
                    prov.inputs["alignment"] = {in_path, "FASTA/NEXUS", hash, sz};
                }
            }
            prov.execution.end_time = end_time_iso.empty() ? Provenance::current_iso8601() : end_time_iso;
            prov.execution.start_time = start_time_iso.empty() ? prov.execution.end_time : start_time_iso;
            prov.execution.wall_time_seconds = runtime_seconds;
            prov.execution.cpu_threads = cpu_threads > 0 ? cpu_threads : 1;
            prov.execution.hostname = Provenance::get_hostname();
            prov.execution.os = Provenance::detect_os();
            prov.execution.compiler = Provenance::detect_compiler();
            j["provenance"] = prov.to_json();
        }

        j["dataset"] = {
            {"file", filepath.empty() ? settings.input_file : filepath},
            {"datatype", settings.use_codon_model ? "codon" : "nucleotide"},
            {"taxa_count", static_cast<int>(num_taxa)},
            {"sites_count", static_cast<int>(num_sites)},
            {"patterns_count", static_cast<int>(num_patterns)},
            {"super_nodes_count", static_cast<int>(num_super_nodes)},
            {"taxa", taxon_names}
        };
        std::string mode_str = (network.active_mode == OhmicScaleMode::ExactDense)
            ? "ExactDense" : "ScalableWoodbury";
        j["electrical_network"] = {
            {"scale_mode", mode_str},
            {"used_hyphaeon_embeddings", network.used_hyphaeon_embeddings},
            {"embedding_dim", network.used_hyphaeon_embeddings ? static_cast<int>(network.hyphaeon_embeddings.cols()) : 0},
            {"foster_sum", foster_sum},
            {"num_bridges", static_cast<int>(num_bridges)},
            {"num_collapsed_polytomies", static_cast<int>(num_collapsed_polytomies)},
            {"backbone_tree", network.backbone_tree.to_newick(true)},
            {"collapsed_tree", network.collapsed_tree.to_newick(true)}
        };

        // Export pairwise matrices in ExactDense mode for interactive visualization (when N <= 128)
        if (network.P_edge.rows() > 0 && network.P_edge.rows() <= 128) {
            auto mat_to_json = [](const Matrix& m) {
                nlohmann::json rows = nlohmann::json::array();
                for (Eigen::Index r = 0; r < m.rows(); ++r) {
                    nlohmann::json row = nlohmann::json::array();
                    for (Eigen::Index c = 0; c < m.cols(); ++c) {
                        row.push_back(m(r, c));
                    }
                    rows.push_back(std::move(row));
                }
                return rows;
            };
            j["electrical_network"]["pairwise_foster_probabilities"] = mat_to_json(network.P_edge);
            j["electrical_network"]["pairwise_effective_resistances"] = mat_to_json(network.Omega);
            j["electrical_network"]["pairwise_tn93_distances"] = mat_to_json(network.D_tn93);
        }

        j["fits"] = {
            {"single_tree", {
                {"log_likelihood", finite_or_null(single_tree_log_l)},
                {"aicc", finite_or_null(single_tree_aicc)}
            }},
            {"ohmic_network", {
                {"log_likelihood", finite_or_null(ohmic_log_l)},
                {"aicc", finite_or_null(ohmic_aicc)},
                {"delta_log_likelihood", finite_or_null(delta_log_l)},
                {"delta_aicc", finite_or_null(delta_aicc)}
            }}
        };
        if (settings.use_codon_model) {
            j["fits"]["omega"] = mg94_params.beta;
        }
        nlohmann::json bridges_arr = nlohmann::json::array();
        for (const auto& b : bridge_summaries) {
            bridges_arr.push_back({
                {"bridge_index", b.bridge_index},
                {"hub_node", b.hub_name},
                {"bridge_node", b.bridge_name},
                {"children", {b.child_names[0], b.child_names[1], b.child_names[2]}},
                {"quartet_clades", {b.quartet_clades[0], b.quartet_clades[1], b.quartet_clades[2], b.quartet_clades[3]}},
                {"ohmic_support", b.ohmic_support},
                {"shunt_stress", b.shunt_stress},
                {"initial_weights", {b.initial_weights[0], b.initial_weights[1], b.initial_weights[2]}},
                {"fitted_weights", {b.fitted_weights[0], b.fitted_weights[1], b.fitted_weights[2]}},
                {"fitted_branch_lengths", {b.fitted_branch_lengths[0], b.fitted_branch_lengths[1], b.fitted_branch_lengths[2]}}
            });
        }
        j["wheatstone_bridges"] = bridges_arr;
        if (!site_delta_log_l.empty()) {
            j["site_results"] = {
                {"delta_log_l", site_delta_log_l},
                {"single_tree_log_l", site_single_tree_log_l},
                {"ohmic_log_l", site_ohmic_log_l}
            };
        }
        j["timers"] = {
            {"total_seconds", runtime_seconds}
        };
        return j;
    }

nlohmann::json OhmicResult::to_legacy_json(const std::string& filepath) const {
        nlohmann::json j = to_modern_json(filepath);
        j["input"] = {
            {"file name", filepath.empty() ? settings.input_file : filepath},
            {"number of sequences", static_cast<int>(num_taxa)},
            {"number of sites", static_cast<int>(num_sites)}
        };
        return j;
    }

} // namespace hyphy::analyses
