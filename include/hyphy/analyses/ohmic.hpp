#pragma once

#include "hyphy/core/types.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/core/rate_matrix.hpp"
#include "hyphy/core/likelihood.hpp"
#include "hyphy/core/ohmic_network.hpp"
#include "hyphy/core/wheatstone_likelihood.hpp"
#include "hyphy/core/provenance.hpp"
#include "hyphy/core/json_output.hpp"
#include "hyphy/core/stats.hpp"
#include "hyphy/opt/optimizer.hpp"
#include "LBFGSB.h"
#include <nlohmann/json.hpp>

#include <chrono>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>

namespace hyphy::analyses {

using namespace hyphy::core;

struct OhmicAnalysisSettings {
    OhmicSettings network_settings;
    bool use_codon_model = false;          // true: MG94xREV codon model; false: GTR nucleotide model
    bool optimize_bridge_weights = true;   // Refine shunt weights w_{k,r} during ML optimization
    bool optimize_branch_lengths = true;   // Refine backbone & bridge branch lengths
    bool quiet = false;
    std::string input_file;
    std::string output_file;
};

struct OhmicBridgeSummary {
    size_t bridge_index = 0;
    std::string hub_name;
    std::string bridge_name;
    std::array<std::string, 3> child_names;
    std::array<std::vector<std::string>, 4> quartet_clades; // [A, B, C, D] leaf taxon names
    std::array<Scalar, 3> initial_weights{1.0, 0.0, 0.0};
    std::array<Scalar, 3> fitted_weights{1.0, 0.0, 0.0};
    std::array<Scalar, 3> fitted_branch_lengths{0.01, 0.01, 0.01};
    Scalar ohmic_support = 1.0;
    Scalar shunt_stress = 0.0;
};

struct OhmicResult {
    OhmicAnalysisSettings settings;
    OhmicNetwork network;
    std::vector<std::string> taxon_names;
    size_t num_taxa = 0;
    size_t num_sites = 0;
    size_t num_patterns = 0;
    size_t num_super_nodes = 0;
    size_t num_bridges = 0;
    size_t num_collapsed_polytomies = 0;
    Scalar foster_sum = 0.0;

    Scalar single_tree_log_l = 0.0;
    Scalar single_tree_aicc = 0.0;
    Scalar ohmic_log_l = 0.0;
    Scalar ohmic_aicc = 0.0;
    Scalar delta_log_l = 0.0;
    Scalar delta_aicc = 0.0;

    GTRParameters gtr_params;
    MG94Parameters mg94_params;
    std::vector<OhmicBridgeSummary> bridge_summaries;
    std::vector<Scalar> site_single_tree_log_l;
    std::vector<Scalar> site_ohmic_log_l;
    std::vector<Scalar> site_delta_log_l;
    double runtime_seconds = 0.0;
    // Execution provenance captured around run(); empty/0 for hand-built results.
    std::string start_time_iso;
    std::string end_time_iso;
    int cpu_threads = 0;

nlohmann::json to_modern_json(const std::string& filepath = "") const;

nlohmann::json to_legacy_json(const std::string& filepath = "") const;
};

class OhmicAnalyzer {
public:
    Alignment aln;
    OhmicAnalysisSettings settings;

    OhmicAnalyzer(Alignment a, OhmicAnalysisSettings s = {})
        : aln(std::move(a)), settings(std::move(s)) {}

    OhmicResult run(const Matrix* external_embeddings = nullptr) {
        auto t_start = std::chrono::steady_clock::now();

        OhmicResult result;
        result.start_time_iso = Provenance::current_iso8601();
#ifdef _OPENMP
        result.cpu_threads = omp_get_max_threads();
#else
        result.cpu_threads = 1;
#endif
        result.settings = settings;
        result.taxon_names = aln.taxon_names;
        result.num_taxa = aln.num_taxa;
        result.num_sites = settings.use_codon_model ? aln.num_codons : aln.num_nucleotides;
        if (result.num_sites == 0) {
            result.num_sites = aln.num_nucleotides > 0 ? aln.num_nucleotides : aln.num_codons;
        }
        result.num_patterns = aln.patterns.size();

        // 1. Infer Ohmic electrical network and Wheatstone bridges
        result.network = OhmicNetwork::infer(aln, settings.network_settings, external_embeddings);
        result.num_super_nodes = result.network.num_super_nodes;
        result.num_bridges = result.network.bridges.size();
        result.num_collapsed_polytomies = result.network.num_collapsed_polytomies;
        result.foster_sum = result.network.foster_sum;

        init_bridge_summaries(result);

        if (settings.use_codon_model && aln.num_codons > 0) {
            fit_codon_network(result);
        } else {
            fit_nucleotide_network(result);
        }

        for (size_t k = 0; k < result.num_bridges; ++k) {
            result.bridge_summaries[k].fitted_weights = result.network.bridges[k].weights;
            result.bridge_summaries[k].fitted_branch_lengths = result.network.bridges[k].branch_lengths;
        }
        result.network.collapsed_tree.bridges = result.network.bridges;
        result.network.collapsed_tree.node_to_bridge_idx = result.network.node_to_bridge_idx;

        auto t_end = std::chrono::steady_clock::now();
        result.runtime_seconds = std::chrono::duration<double>(t_end - t_start).count();
        result.end_time_iso = Provenance::current_iso8601();
        return result;
    }

private:
    void init_bridge_summaries(OhmicResult& result) const {
        auto collect_clade_leaves = [&](int32_t root_nid) -> std::vector<std::string> {
            std::vector<std::string> leaves;
            if (root_nid < 0 || static_cast<size_t>(root_nid) >= result.network.collapsed_tree.nodes.size()) return leaves;
            std::vector<int32_t> stack = {root_nid};
            while (!stack.empty()) {
                int32_t curr = stack.back();
                stack.pop_back();
                const auto& node = result.network.collapsed_tree.nodes[curr];
                if (node.is_leaf) {
                    leaves.push_back(node.name);
                } else {
                    for (int32_t c : node.children) stack.push_back(c);
                }
            }
            return leaves;
        };

        result.bridge_summaries.clear();
        result.bridge_summaries.resize(result.num_bridges);
        for (size_t k = 0; k < result.num_bridges; ++k) {
            const auto& wb = result.network.bridges[k];
            auto& bs = result.bridge_summaries[k];
            bs.bridge_index = k;
            bs.hub_name = (wb.hub_node_id >= 0) ? result.network.collapsed_tree.nodes[wb.hub_node_id].name : "";
            bs.bridge_name = (wb.bridge_node_id >= 0) ? result.network.backbone_tree.nodes[wb.bridge_node_id].name : "";
            std::unordered_set<std::string> in_abc;
            for (int r = 0; r < 3; ++r) {
                int32_t cid = wb.children[r];
                bs.child_names[r] = (cid >= 0) ? result.network.collapsed_tree.nodes[cid].name : "";
                bs.quartet_clades[r] = collect_clade_leaves(cid);
                for (const auto& nm : bs.quartet_clades[r]) in_abc.insert(nm);
            }
            for (const auto& nm : aln.taxon_names) {
                if (in_abc.find(nm) == in_abc.end()) {
                    bs.quartet_clades[3].push_back(nm);
                }
            }
            bs.initial_weights = wb.weights;
            bs.fitted_weights = wb.weights;
            bs.fitted_branch_lengths = wb.branch_lengths;
            bs.ohmic_support = wb.ohmic_support;
            bs.shunt_stress = wb.shunt_stress;
        }
    }

    static Scalar compute_aicc(Scalar log_l, size_t k_params, size_t n_obs) {
        size_t n_clamped = std::max<size_t>(n_obs, k_params + 2);
        return hyphy::core::stats::compute_aicc(log_l, k_params, n_clamped);
    }

    void fit_nucleotide_network(OhmicResult& result) {
        Alignment nuc_aln = (aln.num_nucleotides > 0 && aln.num_codons == 0)
            ? aln
            : aln.to_nucleotide_alignment();

        // 1. Fit baseline single tree (GTR) + Ohmic-regularized Spectral-Cavity NNI
        opt::GTRFitter gtr_fitter(result.network.backbone_tree, nuc_aln);
        auto base_fit = gtr_fitter.fit();
        result.gtr_params = base_fit.params;
        result.network.backbone_tree = base_fit.tree;

        GTRMatrix gtr;
        gtr.update(result.gtr_params, nuc_aln.nuc_frequencies);
        auto bb_leaf_map = LikelihoodEngine::build_leaf_to_taxon_map(result.network.backbone_tree, nuc_aln);

        int total_nni_swaps = 0;
        const Matrix* d_syn_ptr = (result.network.D_backbone.rows() > 0) ? &result.network.D_backbone : nullptr;
        for (int nni_pass = 0; nni_pass < 4; ++nni_pass) {
            int sw = opt::GTRFitter::sweep_gtr_nni_cavity(
                result.network.backbone_tree, nuc_aln, gtr, bb_leaf_map, d_syn_ptr
            );
            total_nni_swaps += sw;
            if (sw == 0) break;
            opt::GTRFitter::sweep_gtr_branches_cavity(
                result.network.backbone_tree, nuc_aln, gtr, bb_leaf_map, 1e-3, 10
            );
            if (sw <= 1) break;
        }

        if (total_nni_swaps > 0) {
            result.single_tree_log_l = LikelihoodEngine::compute_gtr_log_likelihood(
                result.network.backbone_tree, nuc_aln, gtr, &bb_leaf_map
            );
            OhmicNetwork::rebuild_bridges(result.network, aln, settings.network_settings);
            result.num_bridges = result.network.bridges.size();
            result.num_collapsed_polytomies = result.network.num_collapsed_polytomies;
            init_bridge_summaries(result);
        } else {
            result.single_tree_log_l = base_fit.log_likelihood;
        }

        // Transfer optimized backbone branch lengths onto collapsed_tree and bridges
        for (size_t nid = 0; nid < result.network.collapsed_tree.nodes.size(); ++nid) {
            if (nid < result.network.backbone_tree.nodes.size() &&
                (result.network.collapsed_tree.nodes[nid].is_leaf ||
                 !result.network.collapsed_tree.nodes[nid].children.empty())) {
                // Keep 0.0 on root-child bridge hub if hub_node_id == bridge_node_id
                int32_t b_idx = result.network.node_to_bridge_idx[nid];
                if (b_idx >= 0 && result.network.bridges[b_idx].hub_node_id == result.network.bridges[b_idx].bridge_node_id) {
                    result.network.collapsed_tree.nodes[nid].branch_length = 0.0;
                } else {
                    result.network.collapsed_tree.nodes[nid].branch_length = result.network.backbone_tree.nodes[nid].branch_length;
                }
            }
        }
        for (auto& wb : result.network.bridges) {
            if (wb.bridge_node_id >= 0 && static_cast<size_t>(wb.bridge_node_id) < result.network.backbone_tree.nodes.size()) {
                wb.branch_lengths[0] = result.network.backbone_tree.nodes[wb.bridge_node_id].branch_length;
                wb.branch_lengths[1] = std::max<Scalar>(1e-5, wb.branch_lengths[0] * wb.shunt_stress);
            }
        }

        auto leaf_map = LikelihoodEngine::build_leaf_to_taxon_map(result.network.collapsed_tree, nuc_aln);

        // 2. Refine Wheatstone bridge branch lengths and shunt weights under GTR via exact O(N * L) spectral-cavity DFS
        if (!result.network.bridges.empty() && (settings.optimize_branch_lengths || settings.optimize_bridge_weights)) {
            Tree& tree = result.network.collapsed_tree;
            size_t num_nodes = tree.num_nodes();
            size_t num_bridges = result.network.bridges.size();
            size_t num_patterns = nuc_aln.patterns.size();

            std::vector<std::vector<Vector4>> D(num_nodes, std::vector<Vector4>(num_patterns, Vector4::Ones()));
            std::vector<std::vector<Vector4>> P_D(num_nodes, std::vector<Vector4>(num_patterns, Vector4::Ones()));
            std::vector<std::vector<Vector4>> U(num_nodes, std::vector<Vector4>(num_patterns, Vector4::Ones()));
            std::vector<Matrix4> P_branch(num_nodes, Matrix4::Identity());
            std::vector<std::array<Matrix4, 3>> P_bridges(num_bridges);

            std::vector<Vector4> c_spec(num_patterns);
            std::vector<Scalar> base_vec(num_patterns, 0.0);
            std::vector<Scalar> p1_vec(num_patterns, 0.0);
            std::vector<Scalar> p2_vec(num_patterns, 0.0);
            std::vector<Vector4> u_below(num_patterns);

            for (int pass = 0; pass < 2; ++pass) {
                for (const auto& node : tree.nodes) {
                    if (node.id != tree.root_id && (node.is_leaf || !node.children.empty())) {
                        P_branch[node.id] = gtr.transition_matrix(node.branch_length);
                    }
                }
                for (size_t k = 0; k < num_bridges; ++k) {
                    for (int r = 0; r < 3; ++r) {
                        P_bridges[k][r] = (result.network.bridges[k].weights[r] > 0.0 || r < 2)
                            ? gtr.transition_matrix(result.network.bridges[k].branch_lengths[r])
                            : Matrix4::Identity();
                    }
                }

                // Bottom-up postorder pass
                for (int32_t nid : tree.post_order) {
                    const auto& node = tree.nodes[nid];
                    if (node.is_leaf) {
                        size_t t_idx = leaf_map[nid];
                        for (size_t p = 0; p < num_patterns; ++p) {
                            int8_t state = (t_idx != static_cast<size_t>(-1)) ? nuc_aln.patterns[p].states[t_idx] : -1;
                            if (state >= 0 && state < 4) {
                                D[nid][p].setZero();
                                D[nid][p](state) = 1.0;
                            } else {
                                D[nid][p].setOnes();
                            }
                            if (nid != tree.root_id) {
                                P_D[nid][p] = P_branch[nid] * D[nid][p];
                            }
                        }
                    } else {
                        int32_t b_idx = (static_cast<size_t>(nid) < result.network.node_to_bridge_idx.size())
                            ? result.network.node_to_bridge_idx[nid] : -1;
                        for (size_t p = 0; p < num_patterns; ++p) {
                            Vector4 d_vec = Vector4::Zero();
                            if (b_idx >= 0) {
                                const auto& wb = result.network.bridges[static_cast<size_t>(b_idx)];
                                const Vector4& M_A = P_D[wb.children[0]][p];
                                const Vector4& M_B = P_D[wb.children[1]][p];
                                const Vector4& M_C = P_D[wb.children[2]][p];
                                if (wb.weights[0] > 0.0) {
                                    d_vec.noalias() += wb.weights[0] * (P_bridges[b_idx][0] * M_A.cwiseProduct(M_B)).cwiseProduct(M_C);
                                }
                                if (wb.weights[1] > 0.0) {
                                    d_vec.noalias() += wb.weights[1] * (P_bridges[b_idx][1] * M_A.cwiseProduct(M_C)).cwiseProduct(M_B);
                                }
                            } else {
                                d_vec.setOnes();
                                for (int32_t cid : node.children) {
                                    d_vec.array() *= P_D[cid][p].array();
                                    Scalar mx_c = d_vec.maxCoeff();
                                    if (mx_c > 0.0 && mx_c < 1e-100) d_vec /= mx_c;
                                }
                            }
                            Scalar mx = d_vec.maxCoeff();
                            if (mx > 0.0 && mx < 1e-100) d_vec /= mx;
                            D[nid][p] = d_vec;
                            if (nid != tree.root_id) {
                                P_D[nid][p] = P_branch[nid] * d_vec;
                            }
                        }
                    }
                }

                for (size_t p = 0; p < num_patterns; ++p) {
                    U[tree.root_id][p] = gtr.pi;
                }

                std::function<void(int32_t)> dfs_bridge = [&](int32_t uid) {
                    const auto& node = tree.nodes[uid];
                    if (node.is_leaf || node.children.empty()) return;

                    for (size_t p = 0; p < num_patterns; ++p) {
                        u_below[p] = (uid == tree.root_id)
                            ? gtr.pi
                            : (P_branch[uid].transpose() * U[uid][p]);
                        Scalar mx = u_below[p].maxCoeff();
                        if (mx > 0.0 && mx < 1e-100) u_below[p] /= mx;
                    }

                    int32_t b_idx = (static_cast<size_t>(uid) < result.network.node_to_bridge_idx.size())
                        ? result.network.node_to_bridge_idx[uid] : -1;

                    if (b_idx >= 0) {
                        size_t k = static_cast<size_t>(b_idx);
                        auto& wb = result.network.bridges[k];
                        int32_t id_A = wb.children[0];
                        int32_t id_B = wb.children[1];
                        int32_t id_C = wb.children[2];

                        // 1. Optimize secondary bridge length t_{k,2}
                        if (settings.optimize_branch_lengths && wb.weights[1] > 0.0) {
                            for (size_t p = 0; p < num_patterns; ++p) {
                                const Vector4& M_A = P_D[id_A][p];
                                const Vector4& M_B = P_D[id_B][p];
                                const Vector4& M_C = P_D[id_C][p];
                                Vector4 Y1 = M_A.cwiseProduct(M_B);
                                base_vec[p] = wb.weights[0] * u_below[p].dot((P_bridges[k][0] * Y1).cwiseProduct(M_C));
                                Vector4 Y2 = M_A.cwiseProduct(M_C);
                                Vector4 W2 = wb.weights[1] * u_below[p].cwiseProduct(M_B);
                                c_spec[p] = (gtr.U_pi.transpose() * W2).cwiseProduct(gtr.V_pi * Y2);
                            }
                            Scalar init_t2 = wb.branch_lengths[1];
                            auto t2_obj = [&](Scalar t2) -> Scalar {
                                Vector4 exp_lam = (gtr.eigenvalues * t2).array().exp();
                                Scalar lnl = 0.0;
                                for (size_t p = 0; p < num_patterns; ++p) {
                                    Scalar lk = base_vec[p] + c_spec[p].dot(exp_lam);
                                    lnl += nuc_aln.patterns[p].weight * ((lk > 1e-300) ? std::log(lk) : -690.7755);
                                }
                                return -lnl;
                            };
                            auto [best_t2, _] = opt::Brent1D::minimize(t2_obj, 1e-6, std::max<Scalar>(init_t2, 1e-4), 2.0, 1e-4, 20);
                            wb.branch_lengths[1] = best_t2;
                            P_bridges[k][1] = gtr.transition_matrix(best_t2);
                        }

                        // 2. Optimize primary bridge length t_{k,1}
                        if (settings.optimize_branch_lengths) {
                            for (size_t p = 0; p < num_patterns; ++p) {
                                const Vector4& M_A = P_D[id_A][p];
                                const Vector4& M_B = P_D[id_B][p];
                                const Vector4& M_C = P_D[id_C][p];
                                Vector4 Y2 = M_A.cwiseProduct(M_C);
                                base_vec[p] = wb.weights[1] * u_below[p].dot((P_bridges[k][1] * Y2).cwiseProduct(M_B));
                                Vector4 Y1 = M_A.cwiseProduct(M_B);
                                Vector4 W1 = wb.weights[0] * u_below[p].cwiseProduct(M_C);
                                c_spec[p] = (gtr.U_pi.transpose() * W1).cwiseProduct(gtr.V_pi * Y1);
                            }
                            Scalar init_t1 = wb.branch_lengths[0];
                            auto t1_obj = [&](Scalar t1) -> Scalar {
                                Vector4 exp_lam = (gtr.eigenvalues * t1).array().exp();
                                Scalar lnl = 0.0;
                                for (size_t p = 0; p < num_patterns; ++p) {
                                    Scalar lk = base_vec[p] + c_spec[p].dot(exp_lam);
                                    lnl += nuc_aln.patterns[p].weight * ((lk > 1e-300) ? std::log(lk) : -690.7755);
                                }
                                return -lnl;
                            };
                            auto [best_t1, _] = opt::Brent1D::minimize(t1_obj, 1e-6, std::max<Scalar>(init_t1, 1e-4), 2.0, 1e-4, 20);
                            wb.branch_lengths[0] = best_t1;
                            P_bridges[k][0] = gtr.transition_matrix(best_t1);
                        }

                        // 3. Optimize shunt weight w_{k,2} in [0, 0.99] via pure scalar mixture (zero matrix multiplies)
                        if (settings.optimize_bridge_weights) {
                            for (size_t p = 0; p < num_patterns; ++p) {
                                const Vector4& M_A = P_D[id_A][p];
                                const Vector4& M_B = P_D[id_B][p];
                                const Vector4& M_C = P_D[id_C][p];
                                p1_vec[p] = u_below[p].dot((P_bridges[k][0] * M_A.cwiseProduct(M_B)).cwiseProduct(M_C));
                                p2_vec[p] = u_below[p].dot((P_bridges[k][1] * M_A.cwiseProduct(M_C)).cwiseProduct(M_B));
                            }
                            Scalar init_w2 = wb.weights[1];
                            auto w2_obj = [&](Scalar w2) -> Scalar {
                                Scalar w1 = 1.0 - w2;
                                Scalar lnl = 0.0;
                                for (size_t p = 0; p < num_patterns; ++p) {
                                    Scalar lk = w1 * p1_vec[p] + w2 * p2_vec[p];
                                    lnl += nuc_aln.patterns[p].weight * ((lk > 1e-300) ? std::log(lk) : -690.7755);
                                }
                                return -lnl;
                            };
                            auto [best_w2, _] = opt::Brent1D::minimize(w2_obj, 0.0, std::clamp<Scalar>(init_w2, 0.01, 0.95), 0.99, 1e-4, 20);
                            wb.weights = {1.0 - best_w2, best_w2, 0.0};
                        }

                        // 4. Refine the 3 incident child branches {id_A, id_B, id_C} via spectral dot product
                        if (settings.optimize_branch_lengths && wb.weights[1] > 1e-4) {
                            for (int c_idx = 0; c_idx < 3; ++c_idx) {
                                int32_t cid = wb.children[c_idx];
                                if (cid < 0) continue;
                                for (size_t p = 0; p < num_patterns; ++p) {
                                    const Vector4& M_A = P_D[id_A][p];
                                    const Vector4& M_B = P_D[id_B][p];
                                    const Vector4& M_C = P_D[id_C][p];
                                    Vector4 v_c;
                                    if (c_idx == 0) {
                                        Vector4 Z1 = P_bridges[k][0].transpose() * u_below[p].cwiseProduct(M_C);
                                        Vector4 Z2 = P_bridges[k][1].transpose() * u_below[p].cwiseProduct(M_B);
                                        v_c = wb.weights[0] * Z1.cwiseProduct(M_B) + wb.weights[1] * Z2.cwiseProduct(M_C);
                                    } else if (c_idx == 1) {
                                        Vector4 Z1 = P_bridges[k][0].transpose() * u_below[p].cwiseProduct(M_C);
                                        Vector4 PY2 = P_bridges[k][1] * M_A.cwiseProduct(M_C);
                                        v_c = wb.weights[0] * Z1.cwiseProduct(M_A) + wb.weights[1] * u_below[p].cwiseProduct(PY2);
                                    } else {
                                        Vector4 PY1 = P_bridges[k][0] * M_A.cwiseProduct(M_B);
                                        Vector4 Z2 = P_bridges[k][1].transpose() * u_below[p].cwiseProduct(M_B);
                                        v_c = wb.weights[0] * u_below[p].cwiseProduct(PY1) + wb.weights[1] * Z2.cwiseProduct(M_A);
                                    }
                                    c_spec[p] = (gtr.U_pi.transpose() * v_c).cwiseProduct(gtr.V_pi * D[cid][p]);
                                }
                                Scalar init_bl = tree.nodes[cid].branch_length;
                                auto c_obj = [&](Scalar bl) -> Scalar {
                                    Vector4 exp_lam = (gtr.eigenvalues * bl).array().exp();
                                    Scalar lnl = 0.0;
                                    for (size_t p = 0; p < num_patterns; ++p) {
                                        Scalar lk = c_spec[p].dot(exp_lam);
                                        lnl += nuc_aln.patterns[p].weight * ((lk > 1e-300) ? std::log(lk) : -690.7755);
                                    }
                                    return -lnl;
                                };
                                auto [best_c, _] = opt::Brent1D::minimize(c_obj, 1e-6, std::max<Scalar>(init_bl, 1e-4), 2.0, 1e-4, 15);
                                tree.nodes[cid].branch_length = best_c;
                                P_branch[cid] = gtr.transition_matrix(best_c);
                                for (size_t p = 0; p < num_patterns; ++p) {
                                    P_D[cid][p] = P_branch[cid] * D[cid][p];
                                }
                            }
                        }

                        // Propagate upward cavity vectors U[id_A], U[id_B], U[id_C] and recurse into child subtrees
                        for (size_t p = 0; p < num_patterns; ++p) {
                            const Vector4& M_A = P_D[id_A][p];
                            const Vector4& M_B = P_D[id_B][p];
                            const Vector4& M_C = P_D[id_C][p];
                            Vector4 Z1 = P_bridges[k][0].transpose() * u_below[p].cwiseProduct(M_C);
                            Vector4 Z2 = P_bridges[k][1].transpose() * u_below[p].cwiseProduct(M_B);
                            Vector4 PY1 = P_bridges[k][0] * M_A.cwiseProduct(M_B);
                            Vector4 PY2 = P_bridges[k][1] * M_A.cwiseProduct(M_C);

                            Vector4 v_a = wb.weights[0] * Z1.cwiseProduct(M_B) + wb.weights[1] * Z2.cwiseProduct(M_C);
                            Vector4 v_b = wb.weights[0] * Z1.cwiseProduct(M_A) + wb.weights[1] * u_below[p].cwiseProduct(PY2);
                            Vector4 v_c = wb.weights[0] * u_below[p].cwiseProduct(PY1) + wb.weights[1] * Z2.cwiseProduct(M_A);

                            Scalar mxa = v_a.maxCoeff(); if (mxa > 0.0 && mxa < 1e-100) v_a /= mxa;
                            Scalar mxb = v_b.maxCoeff(); if (mxb > 0.0 && mxb < 1e-100) v_b /= mxb;
                            Scalar mxc = v_c.maxCoeff(); if (mxc > 0.0 && mxc < 1e-100) v_c /= mxc;

                            U[id_A][p] = v_a;
                            U[id_B][p] = v_b;
                            U[id_C][p] = v_c;
                        }
                        dfs_bridge(id_A);
                        dfs_bridge(id_B);
                        dfs_bridge(id_C);

                        // Refresh D[uid] and P_D[uid] after child subtrees return
                        for (size_t p = 0; p < num_patterns; ++p) {
                            const Vector4& M_A = P_D[id_A][p];
                            const Vector4& M_B = P_D[id_B][p];
                            const Vector4& M_C = P_D[id_C][p];
                            Vector4 d_vec = wb.weights[0] * (P_bridges[k][0] * M_A.cwiseProduct(M_B)).cwiseProduct(M_C)
                                          + wb.weights[1] * (P_bridges[k][1] * M_A.cwiseProduct(M_C)).cwiseProduct(M_B);
                            Scalar mx = d_vec.maxCoeff();
                            if (mx > 0.0 && mx < 1e-100) d_vec /= mx;
                            D[uid][p] = d_vec;
                            if (uid != tree.root_id) {
                                P_D[uid][p] = P_branch[uid] * d_vec;
                            }
                        }
                    } else {
                        size_t n_ch = node.children.size();
                        for (size_t i = 0; i < n_ch; ++i) {
                            int32_t cid = node.children[i];
                            for (size_t p = 0; p < num_patterns; ++p) {
                                Vector4 u_c = (uid == tree.root_id)
                                    ? gtr.pi
                                    : (P_branch[uid].transpose() * U[uid][p]);
                                for (size_t j = 0; j < n_ch; ++j) {
                                    if (j != i) {
                                        u_c.array() *= P_D[node.children[j]][p].array();
                                        Scalar mx_j = u_c.maxCoeff();
                                        if (mx_j > 0.0 && mx_j < 1e-100) u_c /= mx_j;
                                    }
                                }
                                Scalar mx = u_c.maxCoeff();
                                if (mx > 0.0 && mx < 1e-100) u_c /= mx;
                                U[cid][p] = u_c;
                            }
                            dfs_bridge(cid);
                        }
                        for (size_t p = 0; p < num_patterns; ++p) {
                            Vector4 d_vec = Vector4::Ones();
                            for (int32_t cid : node.children) {
                                d_vec.array() *= P_D[cid][p].array();
                                Scalar mx_c = d_vec.maxCoeff();
                                if (mx_c > 0.0 && mx_c < 1e-100) d_vec /= mx_c;
                            }
                            D[uid][p] = d_vec;
                            if (uid != tree.root_id) {
                                P_D[uid][p] = P_branch[uid] * d_vec;
                            }
                        }
                    }
                };

                dfs_bridge(tree.root_id);
            }
        }

        result.ohmic_log_l = WheatstoneLikelihood::compute_gtr_log_likelihood(
            result.network, nuc_aln, gtr, 1.0, &leaf_map
        );
        // Ensure Ohmic network is at least as good as the nested single-tree special case (w2 = 0)
        if (result.ohmic_log_l < result.single_tree_log_l) {
            for (auto& wb : result.network.bridges) {
                wb.weights = {1.0, 0.0, 0.0};
            }
            result.ohmic_log_l = WheatstoneLikelihood::compute_gtr_log_likelihood(
                result.network, nuc_aln, gtr, 1.0, &leaf_map
            );
        }

        size_t n_obs = nuc_aln.num_nucleotides * nuc_aln.num_taxa;
        size_t k_single = 5 + (2 * nuc_aln.num_taxa - 3);
        size_t k_active_bridges = 0;
        for (const auto& wb : result.network.bridges) {
            if (wb.weights[1] > 1e-4) k_active_bridges++;
        }
        size_t k_ohmic = k_single + k_active_bridges - result.num_collapsed_polytomies;

        result.single_tree_aicc = compute_aicc(result.single_tree_log_l, k_single, n_obs);
        result.ohmic_aicc = compute_aicc(result.ohmic_log_l, k_ohmic, n_obs);
        result.delta_log_l = result.ohmic_log_l - result.single_tree_log_l;
        result.delta_aicc = result.single_tree_aicc - result.ohmic_aicc;

        // Compute per-site log-likelihood profile (Ohmic vs Single-Tree)
        std::vector<Scalar> pat_ohm, pat_single;
        WheatstoneLikelihood::compute_gtr_log_likelihood(result.network, nuc_aln, gtr, 1.0, &leaf_map, &pat_ohm);
        OhmicNetwork single_net = result.network;
        for (auto& wb : single_net.bridges) wb.weights = {1.0, 0.0, 0.0};
        WheatstoneLikelihood::compute_gtr_log_likelihood(single_net, nuc_aln, gtr, 1.0, &leaf_map, &pat_single);

        size_t S_sites = nuc_aln.site_to_pattern.size();
        result.site_single_tree_log_l.resize(S_sites, 0.0);
        result.site_ohmic_log_l.resize(S_sites, 0.0);
        result.site_delta_log_l.resize(S_sites, 0.0);
        for (size_t s = 0; s < S_sites; ++s) {
            size_t p = nuc_aln.site_to_pattern[s];
            result.site_single_tree_log_l[s] = pat_single[p];
            result.site_ohmic_log_l[s] = pat_ohm[p];
            result.site_delta_log_l[s] = pat_ohm[p] - pat_single[p];
        }
    }

    void fit_codon_network(OhmicResult& result) {
        // 1. Fit baseline single tree under MG94xREV
        opt::MG94Fitter mg_fitter(result.network.backbone_tree, aln);
        auto base_fit = mg_fitter.fit_full_model(1.0, MG94Parameters{}, nullptr, 40, 1e-4);
        result.single_tree_log_l = base_fit.log_likelihood;
        result.mg94_params = base_fit.params;
        result.network.backbone_tree = base_fit.tree;

        for (size_t nid = 0; nid < result.network.collapsed_tree.nodes.size(); ++nid) {
            if (nid < base_fit.tree.nodes.size() &&
                (result.network.collapsed_tree.nodes[nid].is_leaf ||
                 !result.network.collapsed_tree.nodes[nid].children.empty())) {
                int32_t b_idx = result.network.node_to_bridge_idx[nid];
                if (b_idx >= 0 && result.network.bridges[b_idx].hub_node_id == result.network.bridges[b_idx].bridge_node_id) {
                    result.network.collapsed_tree.nodes[nid].branch_length = 0.0;
                } else {
                    result.network.collapsed_tree.nodes[nid].branch_length = base_fit.tree.nodes[nid].branch_length;
                }
            }
        }
        for (auto& wb : result.network.bridges) {
            if (wb.bridge_node_id >= 0 && static_cast<size_t>(wb.bridge_node_id) < base_fit.tree.nodes.size()) {
                wb.branch_lengths[0] = base_fit.tree.nodes[wb.bridge_node_id].branch_length;
                wb.branch_lengths[1] = std::max<Scalar>(1e-5, wb.branch_lengths[0] * wb.shunt_stress);
            }
        }

        // 2. Joint Fréchet quasi-Newton optimization over Wheatstone network
        auto leaf_map = LikelihoodEngine::build_leaf_to_taxon_map(result.network.collapsed_tree, aln);
        size_t K = result.network.bridges.size();
        if (K > 0 && (settings.optimize_branch_lengths || settings.optimize_bridge_weights)) {
            // Optimize bridge lengths t_{k,1}, t_{k,2} and logit weights w_{k,2} + omega via L-BFGS-B
            size_t dim = 3 * K + 1;
            Vector x(dim), lb(dim), ub(dim);
            for (size_t k = 0; k < K; ++k) {
                x(3 * k + 0) = std::log(std::clamp<Scalar>(result.network.bridges[k].branch_lengths[0], 1e-5, 5.0));
                lb(3 * k + 0) = -12.0; ub(3 * k + 0) = 2.0;
                x(3 * k + 1) = std::log(std::clamp<Scalar>(result.network.bridges[k].branch_lengths[1], 1e-5, 5.0));
                lb(3 * k + 1) = -12.0; ub(3 * k + 1) = 2.0;
                Scalar w2 = std::clamp<Scalar>(result.network.bridges[k].weights[1], 1e-4, 0.99);
                x(3 * k + 2) = std::log(w2 / (1.0 - w2));
                lb(3 * k + 2) = -10.0; ub(3 * k + 2) = 4.0;
            }
            x(3 * K) = std::log(std::clamp<Scalar>(result.mg94_params.beta, 1e-3, 10.0));
            lb(3 * K) = -8.0; ub(3 * K) = 4.0;

            auto obj = [&](const Vector& y, Vector& grad) -> Scalar {
                for (size_t k = 0; k < K; ++k) {
                    result.network.bridges[k].branch_lengths[0] = std::exp(y(3 * k + 0));
                    result.network.bridges[k].branch_lengths[1] = std::exp(y(3 * k + 1));
                    Scalar w2 = 1.0 / (1.0 + std::exp(-y(3 * k + 2)));
                    result.network.bridges[k].weights = {1.0 - w2, w2, 0.0};
                }
                result.mg94_params.beta = std::exp(y(3 * K));

                auto gres = WheatstoneLikelihood::compute_joint_frechet_gradients(
                    result.network, aln, result.mg94_params, &leaf_map
                );

                grad.resize(dim);
                for (size_t k = 0; k < K; ++k) {
                    Scalar t1 = result.network.bridges[k].branch_lengths[0];
                    Scalar t2 = result.network.bridges[k].branch_lengths[1];
                    Scalar w2 = result.network.bridges[k].weights[1];
                    grad(3 * k + 0) = -t1 * gres.grad_bridge_lengths[k][0];
                    grad(3 * k + 1) = -t2 * gres.grad_bridge_lengths[k][1];
                    Scalar dL_dw2 = gres.grad_bridge_weights[k][1] - gres.grad_bridge_weights[k][0];
                    grad(3 * k + 2) = -w2 * (1.0 - w2) * dL_dw2;
                }
                grad(3 * K) = -result.mg94_params.beta * gres.grad_omega;
                return -gres.log_likelihood;
            };

            LBFGSpp::LBFGSBParam<Scalar> opt_param;
            opt_param.m = 6;
            opt_param.epsilon = 1e-4;
            opt_param.max_iterations = 30;
            LBFGSpp::LBFGSBSolver<Scalar> solver(opt_param);
            Scalar fx = 0.0;
            try {
                solver.minimize(obj, x, fx, lb, ub);
            } catch (...) {}

            for (size_t k = 0; k < K; ++k) {
                result.network.bridges[k].branch_lengths[0] = std::exp(x(3 * k + 0));
                result.network.bridges[k].branch_lengths[1] = std::exp(x(3 * k + 1));
                Scalar w2 = 1.0 / (1.0 + std::exp(-x(3 * k + 2)));
                result.network.bridges[k].weights = {1.0 - w2, w2, 0.0};
            }
            result.mg94_params.beta = std::exp(x(3 * K));
        }

        auto final_g = WheatstoneLikelihood::compute_joint_frechet_gradients(
            result.network, aln, result.mg94_params, &leaf_map
        );
        result.ohmic_log_l = final_g.log_likelihood;

        size_t n_obs = aln.num_codons * aln.num_taxa;
        size_t k_single = 6 + (2 * aln.num_taxa - 3);
        size_t k_active_bridges = 0;
        for (const auto& wb : result.network.bridges) {
            if (wb.weights[1] > 1e-4) k_active_bridges++;
        }
        size_t k_ohmic = k_single + k_active_bridges - result.num_collapsed_polytomies;

        result.single_tree_aicc = compute_aicc(result.single_tree_log_l, k_single, n_obs);
        result.ohmic_aicc = compute_aicc(result.ohmic_log_l, k_ohmic, n_obs);
        result.delta_log_l = result.ohmic_log_l - result.single_tree_log_l;
        result.delta_aicc = result.single_tree_aicc - result.ohmic_aicc;

        // Compute per-site log-likelihood profile (Ohmic vs Single-Tree)
        const GeneticCode& gcode = *(aln.code ? aln.code : GeneticCode::universal());
        MG94Matrix mg94;
        mg94.update(result.mg94_params, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, gcode);
        std::vector<Scalar> pat_ohm, pat_single;
        WheatstoneLikelihood::compute_mg94_log_likelihood(result.network, aln, mg94, 1.0, &leaf_map, &pat_ohm);
        OhmicNetwork single_net = result.network;
        for (auto& wb : single_net.bridges) wb.weights = {1.0, 0.0, 0.0};
        WheatstoneLikelihood::compute_mg94_log_likelihood(single_net, aln, mg94, 1.0, &leaf_map, &pat_single);

        size_t S_sites = aln.site_to_pattern.size();
        result.site_single_tree_log_l.resize(S_sites, 0.0);
        result.site_ohmic_log_l.resize(S_sites, 0.0);
        result.site_delta_log_l.resize(S_sites, 0.0);
        for (size_t s = 0; s < S_sites; ++s) {
            size_t p = aln.site_to_pattern[s];
            result.site_single_tree_log_l[s] = pat_single[p];
            result.site_ohmic_log_l[s] = pat_ohm[p];
            result.site_delta_log_l[s] = pat_ohm[p] - pat_single[p];
        }
    }
};

} // namespace hyphy::analyses
