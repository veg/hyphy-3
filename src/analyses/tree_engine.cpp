// Out-of-line JSON emitter definitions for TreeEngine.
// Moved from include/hyphy/analyses/tree_engine.hpp so that emitter
// edits recompile this single TU instead of every consumer.
#include "hyphy/analyses/tree_engine.hpp"
#include "hyphy/core/json_output.hpp"
#include "hyphy/core/provenance.hpp"

namespace hyphy::analyses {

using namespace hyphy::core;

nlohmann::json TreeEngine::to_json(const TreeEngineResult& res, const Provenance& prov) const {
    nlohmann::json j;
    j["$schema"] = "https://raw.githubusercontent.com/veg/hyphy-3/main/schemas/v3/tree.v3.schema.json";
    j["schema_version"] = "3.0.0";
    j["analysis"] = {
        {"id", "tree"},
        {"name", "HyPhy-3 ML/Ohmic Tree Engine"},
        {"description", "Fast Spectral Scaffolding, Ohmic Divide-and-Conquer, and Native Polytomy Resolution"},
        {"version", "3.0.0"},
        {"authors", "Sergei L Kosakovsky Pond"},
        {"contact", "spond@temple.edu"},
        {"model", res.active_model},
        {"taxa", res.num_taxa},
        {"patterns", res.num_patterns},
        {"runtime_seconds", res.runtime_seconds}
    };
    j["software"] = software_block();

    if (!prov.invocation.cli_command.empty()) {
        j["provenance"] = prov.to_json();
    } else {
        Provenance auto_prov;
        auto_prov.invocation.cli_command = "hyphy3 tree";
        auto_prov.invocation.working_directory = Provenance::get_cwd();
        auto_prov.execution.start_time = Provenance::current_iso8601();
        auto_prov.execution.end_time = auto_prov.execution.start_time;
        auto_prov.execution.wall_time_seconds = res.runtime_seconds;
        auto_prov.execution.cpu_threads = 1;
        auto_prov.execution.hostname = Provenance::get_hostname();
        auto_prov.execution.os = Provenance::detect_os();
        auto_prov.execution.compiler = Provenance::detect_compiler();
        j["provenance"] = auto_prov.to_json();
    }

    j["likelihoods"] = {
        {"initial_gtr", finite_or_null(res.initial_gtr_log_l)},
        {"optimized_gtr", finite_or_null(res.optimized_gtr_log_l)},
        {"final", finite_or_null(res.final_log_l)}
    };
    j["topology"] = {
        {"initial_nodes", res.num_initial_nodes},
        {"final_nodes", res.num_final_nodes},
        {"polytomies_retained", res.num_polytomies_retained},
        {"branches_collapsed", res.num_branches_collapsed},
        {"polytomies_resolved", res.num_polytomies_resolved},
        {"spr_moves_accepted", res.num_spr_moves_accepted},
        {"newick", res.optimized_tree.to_newick(true)}
    };
    j["divide_and_conquer"] = {
        {"enabled", res.dnc_used},
        {"num_islands", res.num_islands},
        {"island_sizes", res.island_sizes},
        {"island_spr_moves", res.dnc_spr_moves}
    };
    j["parameters"] = {
        {"gtr_rates", {
            {"AC", finite_or_null(res.gtr_rates[0])},
            {"AT", finite_or_null(res.gtr_rates[1])},
            {"CG", finite_or_null(res.gtr_rates[2])},
            {"CT", finite_or_null(res.gtr_rates[3])},
            {"GT", finite_or_null(res.gtr_rates[4])}
        }}
    };
    if (res.active_model == "MG94") {
        j["parameters"]["mg94_omega"] = finite_or_null(res.mg94_omega);
    }
    return j;
}

} // namespace hyphy::analyses
