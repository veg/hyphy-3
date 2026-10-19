#pragma once

#include "hyphy/core/types.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/core/rate_matrix.hpp"
#include "hyphy/core/likelihood.hpp"
#include "hyphy/core/param_map.hpp"
#include "hyphy/opt/nelder_mead.hpp"
#include "LBFGSB.h"
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <functional>

namespace hyphy::opt {

using namespace hyphy::core;

struct GTRFitResult {
    Scalar log_likelihood = 0.0;
    Scalar aicc = 0.0;
    GTRParameters params;
    Tree tree;
    int iterations = 0;
};

// Fits nucleotide GTR model (substitution rates + branch lengths) on an alignment
class GTRFitter {
public:
    Tree tree;
    Alignment aln;
    Alignment nuc_aln;

    GTRFitter(Tree t, const Alignment& a)
        : tree(std::move(t)), aln(a), nuc_aln(a.to_nucleotide_alignment()) {}

    GTRFitResult fit() {
        // Collect active non-root branch node indices (skipping collapsed polytomy tombstones)
        std::vector<int32_t> branch_node_ids;
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id && node.parent_id != INVALID_INDEX) {
                branch_node_ids.push_back(node.id);
            }
        }
        size_t num_branches = branch_node_ids.size();

        GTRParameters params;
        params.theta_AC = 0.5;
        params.theta_AT = 0.3;
        params.theta_CG = 0.5;
        params.theta_CT = 1.0;
        params.theta_GT = 0.3;

        // Initialize branch lengths to 0.05 if not set or zero
        Tree curr_tree = tree;
        for (size_t b = 0; b < num_branches; ++b) {
            int32_t nid = branch_node_ids[b];
            if (curr_tree.nodes[nid].branch_length < 1e-4) {
                curr_tree.nodes[nid].branch_length = 0.05;
            }
        }

        // Precompute leaf-to-taxon map ONCE for the entire GTR fit
        std::vector<size_t> leaf_to_taxon = LikelihoodEngine::build_leaf_to_taxon_map(curr_tree, nuc_aln);

        auto eval_tree_lnl = [&](const Tree& t, const GTRParameters& p) -> Scalar {
            GTRMatrix gtr;
            gtr.update(p, nuc_aln.nuc_frequencies);
            return LikelihoodEngine::compute_gtr_log_likelihood(t, nuc_aln, gtr, &leaf_to_taxon);
        };

        // Stage 1: Optimize global tree scale factor s
        auto scale_obj = [&](Scalar s) -> Scalar {
            Tree t_scale = curr_tree;
            for (size_t b = 0; b < num_branches; ++b) {
                int32_t nid = branch_node_ids[b];
                t_scale.nodes[nid].branch_length *= s;
            }
            t_scale.scale_bridge_lengths(s);
            return -eval_tree_lnl(t_scale, params);
        };

        auto [opt_s, opt_s_nll] = Brent1D::minimize(scale_obj, 0.01, 1.0, 50.0, 1e-2, 12);
        for (size_t b = 0; b < num_branches; ++b) {
            int32_t nid = branch_node_ids[b];
            curr_tree.nodes[nid].branch_length *= opt_s;
        }
        curr_tree.scale_bridge_lengths(opt_s);

        // Stage 2 & 3: Alternating coordinate-wise optimization of GTR rates and individual branch lengths
        int total_evals = 0;
        for (int cycle = 0; cycle < 2; ++cycle) {
            // Optimize 5 GTR rate parameters coordinate-wise. Thetas must be
            // re-fit every cycle: the branch-length sweep below changes the tree,
            // so single-pass theta optimization leaves rates stale at the joint
            // optimum (observed as ~13% theta drift and ~14 log units short of
            // the converged likelihood on benchmark datasets).
            auto opt_rate = [&](Scalar& r_ref) {
                auto rate_obj = [&](Scalar val) -> Scalar {
                    Scalar old_val = r_ref;
                    r_ref = val;
                    Scalar lnl = eval_tree_lnl(curr_tree, params);
                    r_ref = old_val;
                    return -lnl;
                };
                auto [best_r, _] = Brent1D::minimize(rate_obj, 1e-4, r_ref, 20.0, 1e-3, 20);
                r_ref = best_r;
                total_evals += 20;
            };

            opt_rate(params.theta_AC);
            opt_rate(params.theta_AT);
            opt_rate(params.theta_CG);
            opt_rate(params.theta_CT);
            opt_rate(params.theta_GT);

            // Optimize individual branch lengths via exact O(N * L) spectral-cavity DFS sweep
            if (curr_tree.bridges.empty()) {
                GTRMatrix gtr_cur;
                gtr_cur.update(params, nuc_aln.nuc_frequencies);
                sweep_gtr_branches_cavity(curr_tree, nuc_aln, gtr_cur, leaf_to_taxon, 1e-3, 12);
                total_evals += static_cast<int>(num_branches);
            } else {
                for (size_t b = 0; b < num_branches; ++b) {
                    int32_t nid = branch_node_ids[b];
                    Scalar init_bl = curr_tree.nodes[nid].branch_length;
                    auto bl_obj = [&](Scalar bl) -> Scalar {
                        Scalar old_bl = curr_tree.nodes[nid].branch_length;
                        curr_tree.nodes[nid].branch_length = bl;
                        Scalar lnl = eval_tree_lnl(curr_tree, params);
                        curr_tree.nodes[nid].branch_length = old_bl;
                        return -lnl;
                    };
                    auto [best_bl, _] = Brent1D::minimize(bl_obj, 1e-6, std::max(init_bl, 1e-4), 5.0, 1e-3, 20);
                    curr_tree.nodes[nid].branch_length = best_bl;
                    total_evals += 20;
                }
            }
        }
        for (auto& br : curr_tree.bridges) {
            if (br.bridge_node_id >= 0 && static_cast<size_t>(br.bridge_node_id) < curr_tree.nodes.size() &&
                br.bridge_node_id != br.hub_node_id) {
                br.branch_lengths[0] = curr_tree.nodes[br.bridge_node_id].branch_length;
            }
        }

        Scalar final_lnl = eval_tree_lnl(curr_tree, params);

        GTRFitResult res;
        res.log_likelihood = final_lnl;
        res.params = params;
        res.tree = curr_tree;
        res.iterations = 3;

        size_t K = 5 + num_branches;
        size_t N = nuc_aln.num_nucleotides;
        res.aicc = -2.0 * res.log_likelihood + 2.0 * K * N / (N - K - 1.0);
        return res;
    }

    static void sweep_gtr_branches_cavity(
        Tree& curr_tree,
        const Alignment& nuc_aln,
        const GTRMatrix& gtr,
        const std::vector<size_t>& leaf_to_taxon,
        Scalar tol = 1e-3,
        int max_iter = 20
    ) {
        size_t num_nodes = curr_tree.num_nodes();
        size_t num_patterns = nuc_aln.patterns.size();
        if (num_nodes == 0 || num_patterns == 0) return;

        std::vector<std::vector<Vector4>> D(num_nodes, std::vector<Vector4>(num_patterns, Vector4::Ones()));
        std::vector<std::vector<Vector4>> P_D(num_nodes, std::vector<Vector4>(num_patterns, Vector4::Ones()));
        std::vector<std::vector<Vector4>> U(num_nodes, std::vector<Vector4>(num_patterns, Vector4::Ones()));
        std::vector<Matrix4> P_branch(num_nodes, Matrix4::Identity());

        for (const auto& node : curr_tree.nodes) {
            if (node.id != curr_tree.root_id && (node.is_leaf || !node.children.empty())) {
                P_branch[node.id] = gtr.transition_matrix(node.branch_length);
            }
        }

        for (int32_t nid : curr_tree.post_order) {
            const auto& node = curr_tree.nodes[nid];
            if (node.is_leaf) {
                size_t t_idx = leaf_to_taxon[nid];
                for (size_t p = 0; p < num_patterns; ++p) {
                    int8_t state = (t_idx != static_cast<size_t>(-1)) ? nuc_aln.patterns[p].states[t_idx] : -1;
                    if (state >= 0 && state < 4) {
                        D[nid][p].setZero();
                        D[nid][p](state) = 1.0;
                    } else {
                        D[nid][p].setOnes();
                    }
                    if (nid != curr_tree.root_id) {
                        P_D[nid][p] = P_branch[nid] * D[nid][p];
                    }
                }
            } else {
                for (size_t p = 0; p < num_patterns; ++p) {
                    Vector4 d_vec = Vector4::Ones();
                    for (int32_t cid : node.children) {
                        d_vec.array() *= P_D[cid][p].array();
                    }
                    Scalar mx = d_vec.maxCoeff();
                    if (mx > 0.0 && mx < 1e-100) d_vec /= mx;
                    D[nid][p] = d_vec;
                    if (nid != curr_tree.root_id) {
                        P_D[nid][p] = P_branch[nid] * d_vec;
                    }
                }
            }
        }

        for (size_t p = 0; p < num_patterns; ++p) {
            U[curr_tree.root_id][p] = gtr.pi;
        }

        std::vector<Vector4> c_spec(num_patterns);

        std::function<void(int32_t)> dfs = [&](int32_t uid) {
            const auto& node = curr_tree.nodes[uid];
            if (uid != curr_tree.root_id && (node.is_leaf || !node.children.empty())) {
                for (size_t p = 0; p < num_patterns; ++p) {
                    Vector4 u_hat = gtr.U_pi.transpose() * U[uid][p];
                    Vector4 d_hat = gtr.V_pi * D[uid][p];
                    c_spec[p] = u_hat.cwiseProduct(d_hat);
                }
                Scalar init_bl = curr_tree.nodes[uid].branch_length;
                auto bl_obj = [&](Scalar bl) -> Scalar {
                    Vector4 exp_lam = (gtr.eigenvalues * bl).array().exp();
                    Scalar lnl = 0.0;
                    for (size_t p = 0; p < num_patterns; ++p) {
                        Scalar lk = c_spec[p].dot(exp_lam);
                        lnl += nuc_aln.patterns[p].weight * ((lk > 1e-300) ? std::log(lk) : -690.7755);
                    }
                    return -lnl;
                };
                auto [best_bl, _] = Brent1D::minimize(bl_obj, 1e-6, std::max<Scalar>(init_bl, 1e-4), 5.0, tol, max_iter);
                curr_tree.nodes[uid].branch_length = best_bl;
                P_branch[uid] = gtr.transition_matrix(best_bl);
                for (size_t p = 0; p < num_patterns; ++p) {
                    P_D[uid][p] = P_branch[uid] * D[uid][p];
                }
            }

            if (!node.is_leaf && !node.children.empty()) {
                size_t n_ch = node.children.size();
                for (size_t i = 0; i < n_ch; ++i) {
                    int32_t cid = node.children[i];
                    for (size_t p = 0; p < num_patterns; ++p) {
                        Vector4 u_c = (uid == curr_tree.root_id)
                            ? gtr.pi
                            : (P_branch[uid].transpose() * U[uid][p]);
                        for (size_t j = 0; j < n_ch; ++j) {
                            if (j != i) {
                                u_c.array() *= P_D[node.children[j]][p].array();
                            }
                        }
                        Scalar mx = u_c.maxCoeff();
                        if (mx > 0.0 && mx < 1e-100) u_c /= mx;
                        U[cid][p] = u_c;
                    }
                    dfs(cid);
                }
                for (size_t p = 0; p < num_patterns; ++p) {
                    Vector4 d_vec = Vector4::Ones();
                    for (int32_t cid : node.children) {
                        d_vec.array() *= P_D[cid][p].array();
                    }
                    Scalar mx = d_vec.maxCoeff();
                    if (mx > 0.0 && mx < 1e-100) d_vec /= mx;
                    D[uid][p] = d_vec;
                    if (uid != curr_tree.root_id) {
                        P_D[uid][p] = P_branch[uid] * d_vec;
                    }
                }
            }
        };

        dfs(curr_tree.root_id);
    }

    static int sweep_gtr_nni_cavity(
        Tree& curr_tree,
        const Alignment& nuc_aln,
        const GTRMatrix& gtr,
        const std::vector<size_t>& leaf_to_taxon,
        const Matrix* D_backbone = nullptr
    ) {
        size_t num_nodes = curr_tree.num_nodes();
        size_t num_patterns = nuc_aln.patterns.size();
        size_t N_taxa = nuc_aln.num_taxa;
        if (num_nodes == 0 || num_patterns == 0 || N_taxa < 4) return 0;

        curr_tree.compute_post_order();

        (void)D_backbone;
        std::vector<std::vector<Vector4>> D(num_nodes, std::vector<Vector4>(num_patterns, Vector4::Ones()));
        std::vector<std::vector<Vector4>> P_D(num_nodes, std::vector<Vector4>(num_patterns, Vector4::Ones()));
        std::vector<std::vector<Vector4>> U(num_nodes, std::vector<Vector4>(num_patterns, Vector4::Ones()));
        std::vector<Matrix4> P_branch(num_nodes, Matrix4::Identity());

        for (const auto& node : curr_tree.nodes) {
            if (node.id != curr_tree.root_id && (node.is_leaf || !node.children.empty())) {
                P_branch[node.id] = gtr.transition_matrix(node.branch_length);
            }
        }

        for (int32_t nid : curr_tree.post_order) {
            const auto& node = curr_tree.nodes[nid];
            if (node.is_leaf) {
                size_t t_idx = leaf_to_taxon[nid];
                for (size_t p = 0; p < num_patterns; ++p) {
                    int8_t state = (t_idx != static_cast<size_t>(-1)) ? nuc_aln.patterns[p].states[t_idx] : -1;
                    if (state >= 0 && state < 4) {
                        D[nid][p].setZero();
                        D[nid][p](state) = 1.0;
                    } else {
                        D[nid][p].setOnes();
                    }
                    if (nid != curr_tree.root_id) {
                        P_D[nid][p] = P_branch[nid] * D[nid][p];
                    }
                }
            } else {
                for (size_t p = 0; p < num_patterns; ++p) {
                    Vector4 d_vec = Vector4::Ones();
                    for (int32_t cid : node.children) {
                        d_vec.array() *= P_D[cid][p].array();
                    }
                    Scalar mx = d_vec.maxCoeff();
                    if (mx > 0.0 && mx < 1e-100) d_vec /= mx;
                    D[nid][p] = d_vec;
                    if (nid != curr_tree.root_id) {
                        P_D[nid][p] = P_branch[nid] * d_vec;
                    }
                }
            }
        }

        for (size_t p = 0; p < num_patterns; ++p) {
            U[curr_tree.root_id][p] = gtr.pi;
        }

        std::vector<Vector4> c_spec1(num_patterns);
        std::vector<Vector4> c_spec2(num_patterns);
        std::vector<Vector4> c_spec3(num_patterns);
        std::vector<Vector4> M_D(num_patterns);
        std::vector<Scalar> dp_vec(num_patterns, 0.0);
        int accepted_swaps = 0;
        size_t L_sites = nuc_aln.site_to_pattern.size();
        bool is_coding = (L_sites >= 30 && (L_sites % 3) == 0);
        size_t n_codons = L_sites / 3;
        size_t k_trim_codons = std::max<size_t>(2, static_cast<size_t>(std::round(0.025 * n_codons)));
        Scalar trim_budget_sites = std::max<Scalar>(2.0, 0.025 * static_cast<Scalar>(std::max<size_t>(1, nuc_aln.num_nucleotides)));

        // Stratify patterns into 4 empirical Gamma rate quartiles by site Shannon entropy
        std::vector<std::pair<Scalar, size_t>> pat_ent(num_patterns);
        for (size_t p = 0; p < num_patterns; ++p) {
            std::array<int, 4> cnt{0, 0, 0, 0};
            int tot = 0;
            for (int8_t st : nuc_aln.patterns[p].states) {
                if (st >= 0 && st < 4) {
                    cnt[st]++;
                    tot++;
                }
            }
            Scalar ent = 0.0;
            if (tot > 0) {
                for (int b = 0; b < 4; ++b) {
                    if (cnt[b] > 0) {
                        Scalar f = static_cast<Scalar>(cnt[b]) / static_cast<Scalar>(tot);
                        ent -= f * std::log(f);
                    }
                }
            }
            pat_ent[p] = {ent, p};
        }
        std::sort(pat_ent.begin(), pat_ent.end());
        std::vector<uint8_t> pat_bin(num_patterns, 1);
        for (size_t rank = 0; rank < num_patterns; ++rank) {
            uint8_t bin = static_cast<uint8_t>(std::min<size_t>(3, (rank * 4) / num_patterns));
            pat_bin[pat_ent[rank].second] = bin;
        }
        const std::array<Scalar, 4> gamma_rates{0.22, 0.60, 1.08, 2.10};

        auto eval_nni_at = [&](int32_t pid) {
            size_t n_ch = curr_tree.nodes[pid].children.size();
            for (size_t c_slot = 0; c_slot < n_ch; ++c_slot) {
                int32_t uid = curr_tree.nodes[pid].children[c_slot];
                if (curr_tree.nodes[uid].is_leaf || curr_tree.nodes[uid].children.size() != 2) continue;

                // Test against each sibling slot under pid
                for (size_t sib_slot = 0; sib_slot < n_ch; ++sib_slot) {
                    if (sib_slot == c_slot) continue;
                    int32_t id_C = curr_tree.nodes[pid].children[sib_slot];
                    int32_t id_A = curr_tree.nodes[uid].children[0];
                    int32_t id_B = curr_tree.nodes[uid].children[1];

                    for (size_t p = 0; p < num_patterns; ++p) {
                        Vector4 md = (pid == curr_tree.root_id)
                            ? gtr.pi
                            : (P_branch[pid].transpose() * U[pid][p]);
                        for (size_t k = 0; k < n_ch; ++k) {
                            if (k != c_slot && k != sib_slot) {
                                md.array() *= P_D[curr_tree.nodes[pid].children[k]][p].array();
                            }
                        }
                        Scalar mx = md.maxCoeff();
                        if (mx > 0.0 && mx < 1e-100) md /= mx;
                        M_D[p] = md;

                        const Vector4& M_A = P_D[id_A][p];
                        const Vector4& M_B = P_D[id_B][p];
                        const Vector4& M_C = P_D[id_C][p];

                        Vector4 Y1 = M_A.cwiseProduct(M_B);
                        Vector4 W1 = md.cwiseProduct(M_C);
                        c_spec1[p] = (gtr.U_pi.transpose() * W1).cwiseProduct(gtr.V_pi * Y1);

                        Vector4 Y2 = M_A.cwiseProduct(M_C);
                        Vector4 W2 = md.cwiseProduct(M_B);
                        c_spec2[p] = (gtr.U_pi.transpose() * W2).cwiseProduct(gtr.V_pi * Y2);

                        Vector4 Y3 = M_B.cwiseProduct(M_C);
                        Vector4 W3 = md.cwiseProduct(M_A);
                        c_spec3[p] = (gtr.U_pi.transpose() * W3).cwiseProduct(gtr.V_pi * Y3);
                    }

                    auto make_obj = [&](const std::vector<Vector4>& spec) {
                        return [&](Scalar bl) -> Scalar {
                            std::array<Vector4, 4> exp_lam;
                            for (int b = 0; b < 4; ++b) {
                                exp_lam[b] = (gtr.eigenvalues * (gamma_rates[b] * bl)).array().exp();
                            }
                            Scalar lnl = 0.0;
                            for (size_t p = 0; p < num_patterns; ++p) {
                                Scalar lk = spec[p].dot(exp_lam[pat_bin[p]]);
                                lnl += nuc_aln.patterns[p].weight * ((lk > 1e-300) ? std::log(lk) : -690.7755);
                            }
                            return -lnl;
                        };
                    };

                    Scalar init_bl = std::clamp<Scalar>(curr_tree.nodes[uid].branch_length, 1e-4, 1.0);
                    auto [t1, nll1] = Brent1D::minimize(make_obj(c_spec1), 1e-6, init_bl, 2.0, 2e-3, 10);
                    auto [t2, nll2] = Brent1D::minimize(make_obj(c_spec2), 1e-6, init_bl, 2.0, 2e-3, 10);
                    auto [t3, nll3] = Brent1D::minimize(make_obj(c_spec3), 1e-6, init_bl, 2.0, 2e-3, 10);

                    Scalar lnl1 = -nll1;
                    Scalar lnl2 = -nll2;
                    Scalar lnl3 = -nll3;

                    int best_r = (lnl2 >= lnl3) ? 2 : 3;
                    Scalar best_lnl = (best_r == 2) ? lnl2 : lnl3;
                    Scalar best_t = (best_r == 2) ? t2 : t3;
                    const auto& best_spec = (best_r == 2) ? c_spec2 : c_spec3;

                    if (best_lnl > lnl1 + 0.25) {
                        std::array<Vector4, 4> exp1, expr;
                        for (int b = 0; b < 4; ++b) {
                            exp1[b] = (gtr.eigenvalues * (gamma_rates[b] * t1)).array().exp();
                            expr[b] = (gtr.eigenvalues * (gamma_rates[b] * best_t)).array().exp();
                        }
                        Scalar d_pos = 0.0, d_neg = 0.0;

                        for (size_t p = 0; p < num_patterns; ++p) {
                            Scalar lk1 = c_spec1[p].dot(exp1[pat_bin[p]]);
                            Scalar lkr = best_spec[p].dot(expr[pat_bin[p]]);
                            Scalar l1 = (lk1 > 1e-300) ? std::log(lk1) : -690.7755;
                            Scalar lr = (lkr > 1e-300) ? std::log(lkr) : -690.7755;
                            Scalar dp = lr - l1;
                            dp_vec[p] = dp;
                            Scalar w = nuc_aln.patterns[p].weight;
                            if (dp > 0.0) {
                                d_pos += w * dp;
                            } else if (dp < 0.0) {
                                d_neg += w * (-dp);
                            }
                        }

                        bool accept_move = false;
                        if (is_coding) {
                            Scalar tot_p12 = 0.0, tot_p3 = 0.0;
                            std::vector<Scalar> pos_codons;
                            pos_codons.reserve(n_codons);
                            int n_pos_cod = 0, n_neg_cod = 0;

                            for (size_t c = 0; c < n_codons; ++c) {
                                Scalar d12 = dp_vec[nuc_aln.site_to_pattern[3 * c]] +
                                             dp_vec[nuc_aln.site_to_pattern[3 * c + 1]];
                                Scalar d3 = dp_vec[nuc_aln.site_to_pattern[3 * c + 2]];
                                tot_p12 += d12;
                                tot_p3 += d3;
                                Scalar dc = d12 + d3;
                                if (dc > 0.05) {
                                    pos_codons.push_back(dc);
                                    n_pos_cod++;
                                } else if (dc < -0.05) {
                                    n_neg_cod++;
                                }
                            }
                            std::sort(pos_codons.begin(), pos_codons.end(), std::greater<Scalar>());
                            Scalar top_outliers = 0.0;
                            for (size_t k = 0; k < std::min(k_trim_codons, pos_codons.size()); ++k) {
                                top_outliers += pos_codons[k];
                            }
                            Scalar d_trimmed_codon = (tot_p12 + tot_p3) - top_outliers;

                            if (tot_p3 >= -0.25 && (d_neg <= 1.5 || d_trimmed_codon > 0.1)) {
                                accept_move = true;
                            } else if (d_trimmed_codon > 1.5 && n_pos_cod >= n_neg_cod + 2) {
                                accept_move = true;
                            }
                        } else {
                            std::vector<std::pair<Scalar, Scalar>> pos_pats;
                            for (size_t p = 0; p < num_patterns; ++p) {
                                if (dp_vec[p] > 0.0) pos_pats.push_back({dp_vec[p], nuc_aln.patterns[p].weight});
                            }
                            std::sort(pos_pats.begin(), pos_pats.end(), [](const auto& a, const auto& b) {
                                return a.first > b.first;
                            });
                            Scalar rem_trim = trim_budget_sites;
                            Scalar trimmed_pos = d_pos;
                            for (const auto& [dp, w] : pos_pats) {
                                if (rem_trim <= 0.0) break;
                                Scalar take = std::min(w, rem_trim);
                                trimmed_pos -= take * dp;
                                rem_trim -= take;
                            }
                            accept_move = (d_neg <= 1.5 || (trimmed_pos - d_neg) > 0.2);
                        }

                        if (accept_move) {
                            if (best_r == 2) {
                                curr_tree.nodes[uid].children[1] = id_C;
                                curr_tree.nodes[pid].children[sib_slot] = id_B;
                                curr_tree.nodes[id_C].parent_id = uid;
                                curr_tree.nodes[id_B].parent_id = pid;
                            } else {
                                curr_tree.nodes[uid].children[0] = id_C;
                                curr_tree.nodes[pid].children[sib_slot] = id_A;
                                curr_tree.nodes[id_C].parent_id = uid;
                                curr_tree.nodes[id_A].parent_id = pid;
                            }
                            curr_tree.nodes[uid].branch_length = best_t;
                            P_branch[uid] = gtr.transition_matrix(best_t);

                            int32_t new_A = curr_tree.nodes[uid].children[0];
                            int32_t new_B = curr_tree.nodes[uid].children[1];
                            for (size_t p = 0; p < num_patterns; ++p) {
                                Vector4 d_vec = P_D[new_A][p].cwiseProduct(P_D[new_B][p]);
                                Scalar mx = d_vec.maxCoeff();
                                if (mx > 0.0 && mx < 1e-100) d_vec /= mx;
                                D[uid][p] = d_vec;
                                P_D[uid][p] = P_branch[uid] * d_vec;
                            }
                            accepted_swaps++;
                        }
                    }
                }
            }
        };

        std::function<void(int32_t)> dfs_nni = [&](int32_t pid) {
            if (curr_tree.nodes[pid].is_leaf || curr_tree.nodes[pid].children.size() < 2) return;

            // 1. Pre-order NNI check at pid
            eval_nni_at(pid);

            // 2. Propagate U[cid] and recurse into child subtrees
            size_t n_ch = curr_tree.nodes[pid].children.size();
            for (size_t i = 0; i < n_ch; ++i) {
                int32_t cid = curr_tree.nodes[pid].children[i];
                for (size_t p = 0; p < num_patterns; ++p) {
                    Vector4 u_c = (pid == curr_tree.root_id)
                        ? gtr.pi
                        : (P_branch[pid].transpose() * U[pid][p]);
                    for (size_t j = 0; j < n_ch; ++j) {
                        if (j != i) {
                            u_c.array() *= P_D[curr_tree.nodes[pid].children[j]][p].array();
                        }
                    }
                    Scalar mx = u_c.maxCoeff();
                    if (mx > 0.0 && mx < 1e-100) u_c /= mx;
                    U[cid][p] = u_c;
                }
                dfs_nni(cid);
            }

            // 3. Post-order NNI check at pid using freshly refined child subtrees
            eval_nni_at(pid);

            for (size_t p = 0; p < num_patterns; ++p) {
                Vector4 d_vec = Vector4::Ones();
                for (int32_t cid : curr_tree.nodes[pid].children) {
                    d_vec.array() *= P_D[cid][p].array();
                }
                Scalar mx = d_vec.maxCoeff();
                if (mx > 0.0 && mx < 1e-100) d_vec /= mx;
                D[pid][p] = d_vec;
                if (pid != curr_tree.root_id) {
                    P_D[pid][p] = P_branch[pid] * d_vec;
                }
            }
        };

        dfs_nni(curr_tree.root_id);
        if (accepted_swaps > 0) {
            curr_tree.compute_post_order();
        }
        return accepted_swaps;
    }
};

// Optimizes branch lengths on a fixed GTR model
class GTRBranchLengthFitter {
public:
    static std::pair<Tree, Scalar> fit(
        Tree tree,
        const Alignment& aln,
        const GTRParameters& params,
        const Vector4& freqs,
        int passes = 1
    ) {
        std::vector<int32_t> branch_node_ids;
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id && node.parent_id != INVALID_INDEX) {
                branch_node_ids.push_back(node.id);
            }
        }
        size_t num_branches = branch_node_ids.size();

        Tree curr_tree = std::move(tree);
        for (size_t b = 0; b < num_branches; ++b) {
            int32_t nid = branch_node_ids[b];
            if (curr_tree.nodes[nid].branch_length < 1e-6) {
                curr_tree.nodes[nid].branch_length = 0.01;
            }
        }

        std::vector<size_t> leaf_to_taxon = LikelihoodEngine::build_leaf_to_taxon_map(curr_tree, aln);

        GTRMatrix gtr;
        gtr.update(params, freqs);

        auto eval_tree_lnl = [&](const Tree& t) -> Scalar {
            return LikelihoodEngine::compute_gtr_log_likelihood(t, aln, gtr, &leaf_to_taxon);
        };

        // Stage 1: Scale factor optimization
        auto scale_obj = [&](Scalar s) -> Scalar {
            Tree t_scale = curr_tree;
            for (size_t b = 0; b < num_branches; ++b) {
                int32_t nid = branch_node_ids[b];
                t_scale.nodes[nid].branch_length *= s;
            }
            t_scale.scale_bridge_lengths(s);
            return -eval_tree_lnl(t_scale);
        };

        auto [opt_s, _] = Brent1D::minimize(scale_obj, 0.01, 1.0, 50.0, 1e-3, 20);
        for (size_t b = 0; b < num_branches; ++b) {
            int32_t nid = branch_node_ids[b];
            curr_tree.nodes[nid].branch_length *= opt_s;
        }
        curr_tree.scale_bridge_lengths(opt_s);

        // Stage 2: Coordinate sweeps
        for (int p = 0; p < passes; ++p) {
            if (curr_tree.bridges.empty()) {
                GTRFitter::sweep_gtr_branches_cavity(curr_tree, aln, gtr, leaf_to_taxon, 1e-3, 15);
            } else {
                for (size_t b = 0; b < num_branches; ++b) {
                    int32_t nid = branch_node_ids[b];
                    Scalar init_bl = curr_tree.nodes[nid].branch_length;
                    auto bl_obj = [&](Scalar bl) -> Scalar {
                        Scalar old_bl = curr_tree.nodes[nid].branch_length;
                        curr_tree.nodes[nid].branch_length = bl;
                        Scalar lnl = eval_tree_lnl(curr_tree);
                        curr_tree.nodes[nid].branch_length = old_bl;
                        return -lnl;
                    };
                    auto [best_bl, _] = Brent1D::minimize(bl_obj, 1e-6, std::max(init_bl, static_cast<Scalar>(1e-4)), 5.0, 1e-3, 15);
                    curr_tree.nodes[nid].branch_length = best_bl;
                }
            }
        }
        for (auto& br : curr_tree.bridges) {
            if (br.bridge_node_id >= 0 && static_cast<size_t>(br.bridge_node_id) < curr_tree.nodes.size() &&
                br.bridge_node_id != br.hub_node_id) {
                br.branch_lengths[0] = curr_tree.nodes[br.bridge_node_id].branch_length;
            }
        }

        Scalar final_lnl = eval_tree_lnl(curr_tree);
        return {curr_tree, final_lnl};
    }
};

struct FitResult {
    Scalar log_likelihood = 0.0;
    Vector x_opt;
    std::vector<std::string> param_names;
    int iterations = 0;
    bool converged = false;
    Tree tree;
    MG94Parameters params;
};

// Fits global parameters (e.g. omega, transition rates) on a fixed tree
class MG94Fitter {
public:
    Tree tree;
    Alignment aln;

    MG94Fitter(Tree t, Alignment a) : tree(std::move(t)), aln(std::move(a)) {}

    // Functor for L-BFGS-B: minimizes negative log-likelihood
    struct Objective {
        const Tree& tree;
        const Alignment& aln;
        MG94Parameters base_params;
        bool fit_rates = false;

        Scalar operator()(const Vector& x, Vector& grad) {
            Scalar omega = x(0);

            MG94Parameters p = base_params;
            p.alpha = 1.0;
            p.beta = omega;
            if (fit_rates && x.size() >= 6) {
                p.theta_AC = x(1);
                p.theta_AT = x(2);
                p.theta_CG = x(3);
                p.theta_CT = x(4);
                p.theta_GT = x(5);
            }

            auto eval_lnl = [&](const MG94Parameters& params) -> Scalar {
                MG94Matrix mg;
                mg.update(params, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *(aln.code ? aln.code : GeneticCode::universal()));
                return LikelihoodEngine::compute_mg94_log_likelihood(tree, aln, mg);
            };

            Scalar f_val = -eval_lnl(p);

            grad.resize(x.size());
            for (int i = 0; i < x.size(); ++i) {
                Scalar h = std::max(1e-6, 1e-5 * std::abs(x(i)));
                Vector x_plus = x;
                Vector x_minus = x;
                x_plus(i) += h;
                x_minus(i) -= h;

                MG94Parameters p_plus = p;
                MG94Parameters p_minus = p;

                if (i == 0) {
                    p_plus.beta = x_plus(0);
                    p_minus.beta = x_minus(0);
                } else if (fit_rates) {
                    if (i == 1) { p_plus.theta_AC = x_plus(1); p_minus.theta_AC = x_minus(1); }
                    if (i == 2) { p_plus.theta_AT = x_plus(2); p_minus.theta_AT = x_minus(2); }
                    if (i == 3) { p_plus.theta_CG = x_plus(3); p_minus.theta_CG = x_minus(3); }
                    if (i == 4) { p_plus.theta_CT = x_plus(4); p_minus.theta_CT = x_minus(4); }
                    if (i == 5) { p_plus.theta_GT = x_plus(5); p_minus.theta_GT = x_minus(5); }
                }

                Scalar f_plus = -eval_lnl(p_plus);
                Scalar f_minus = -eval_lnl(p_minus);
                grad(i) = (f_plus - f_minus) / (2.0 * h);
            }

            return f_val;
        }
    };

    FitResult fit_omega(Scalar init_omega = 1.0, MG94Parameters base_params = MG94Parameters{}) {
        LBFGSpp::LBFGSBParam<Scalar> param;
        param.m = 6;
        param.epsilon = 1e-5;
        param.max_iterations = 50;

        LBFGSpp::LBFGSBSolver<Scalar> solver(param);

        Objective obj{tree, aln, base_params, false};

        Vector x(1);
        x(0) = init_omega;

        Vector lb(1), ub(1);
        lb(0) = 1e-4;
        ub(0) = 100.0;

        Scalar fx = 0.0;
        int niter = 0;
        try {
            niter = solver.minimize(obj, x, fx, lb, ub);
        } catch (const std::exception& e) {
            std::cerr << "L-BFGS-B warning: " << e.what() << "\n";
        }

        FitResult res;
        res.log_likelihood = -fx;
        res.x_opt = x;
        res.param_names = {"omega"};
        res.iterations = niter;
        res.converged = true;
        res.tree = tree;
        return res;
    }

    // Jointly fits omega and global tree branch length scaling factor s
    FitResult fit_omega_and_scale(Scalar init_omega = 1.0, MG94Parameters base_params = MG94Parameters{}) {
        auto eval_at = [&](Scalar omega, Scalar s) -> Scalar {
            MG94Parameters p = base_params;
            p.alpha = 1.0;
            p.beta = omega;

            MG94Matrix mg;
            mg.update(p, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *(aln.code ? aln.code : GeneticCode::universal()));
            Scalar conv = (mg.scale_factor > 1e-12) ? (3.0 / mg.scale_factor) : 1.0;

            Tree t_eval = tree;
            for (auto& node : t_eval.nodes) {
                if (node.id != t_eval.root_id) {
                    node.branch_length *= (s * conv);
                }
            }
            t_eval.scale_bridge_lengths(s * conv);

            return -LikelihoodEngine::compute_mg94_log_likelihood(t_eval, aln, mg);
        };

        auto opt_pt = NelderMead2D::minimize(eval_at, init_omega, 1.0, 1e-4, 50.0, 1e-4, 80);

        FitResult res;
        res.log_likelihood = -opt_pt.f;
        res.x_opt.resize(2);
        res.x_opt(0) = opt_pt.a; // omega
        res.x_opt(1) = opt_pt.b; // scale s
        res.param_names = {"omega", "tree_scale"};
        res.iterations = 50;
        res.converged = true;

        res.tree = tree;
        for (auto& node : res.tree.nodes) {
            if (node.id != res.tree.root_id) {
                node.branch_length *= opt_pt.b;
            }
        }
        res.tree.scale_bridge_lengths(opt_pt.b);
        res.params = base_params;
        res.params.alpha = 1.0;
        res.params.beta = opt_pt.a;
        return res;
    }

    // Jointly fits all branch lengths, global omega, and GTR exchangeability rates
    // under unconstrained MG94 in a unified single-stage quasi-Newton solve via
    // spectral Fréchet matrix exponential derivatives.
    FitResult fit_full_model(
        Scalar init_omega = 1.0,
        MG94Parameters base_params = MG94Parameters{},
        std::function<void(const std::string&, double)> progress_cb = nullptr,
        int max_iters = 60,
        Scalar tolerance = 1e-5
    ) {
        if (progress_cb) progress_cb("Phase 2: Initializing branch lengths and omega", 0.30);

        // 1. Initial quick fit of omega & global branch scale (warm-up)
        auto init_res = fit_omega_and_scale(init_omega, base_params);
        Tree cur_tree = init_res.tree;
        MG94Parameters p = base_params;
        p.alpha = 1.0;
        p.beta = init_res.x_opt(0);

        std::vector<size_t> leaf_to_taxon = LikelihoodEngine::build_leaf_to_taxon_map(cur_tree, aln);

        // Collect active non-root branches (skipping collapsed polytomy tombstones)
        std::vector<int32_t> branch_node_ids;
        for (const auto& node : cur_tree.nodes) {
            if (node.id != cur_tree.root_id && node.parent_id != INVALID_INDEX) {
                branch_node_ids.push_back(node.id);
            }
        }
        size_t B = branch_node_ids.size();
        size_t dim = B + 6;
        size_t num_nodes = cur_tree.num_nodes();

        Vector x(dim), lb(dim), ub(dim);
        for (size_t i = 0; i < B; ++i) {
            x(i) = std::log(std::clamp(cur_tree.nodes[branch_node_ids[i]].branch_length, 1e-7, 10.0));
            lb(i) = -18.420680744;
            ub(i) = 3.0;
        }

        x(B) = std::log(std::clamp(p.beta, 1e-3, 10.0));
        lb(B) = -8.0; ub(B) = 4.0;

        x(B + 1) = std::log(std::clamp(p.theta_AC, 1e-3, 20.0)); lb(B + 1) = -8.0; ub(B + 1) = 4.0;
        x(B + 2) = std::log(std::clamp(p.theta_AT, 1e-3, 20.0)); lb(B + 2) = -8.0; ub(B + 2) = 4.0;
        x(B + 3) = std::log(std::clamp(p.theta_CG, 1e-3, 20.0)); lb(B + 3) = -8.0; ub(B + 3) = 4.0;
        x(B + 4) = std::log(std::clamp(p.theta_CT, 1e-3, 20.0)); lb(B + 4) = -8.0; ub(B + 4) = 4.0;
        x(B + 5) = std::log(std::clamp(p.theta_GT, 1e-3, 20.0)); lb(B + 5) = -8.0; ub(B + 5) = 4.0;

        struct JointMG94Objective {
            const Tree& tree;
            const Alignment& aln;
            const std::vector<int32_t>& nodes;
            size_t num_nodes;
            const std::vector<size_t>& leaf_map;

            Scalar operator()(const Vector& y, Vector& grad) {
                size_t B = nodes.size();
                std::vector<Scalar> bls(num_nodes, 0.0);
                for (size_t i = 0; i < B; ++i) {
                    bls[nodes[i]] = std::exp(y(i));
                }

                MG94Parameters params;
                params.alpha = 1.0;
                params.beta = std::exp(y(B));
                params.theta_AC = std::exp(y(B + 1));
                params.theta_AT = std::exp(y(B + 2));
                params.theta_CG = std::exp(y(B + 3));
                params.theta_CT = std::exp(y(B + 4));
                params.theta_GT = std::exp(y(B + 5));

                auto res = LikelihoodEngine::compute_joint_frechet_gradients(
                    tree, aln, params, bls, &leaf_map
                );

                grad.resize(y.size());
                for (size_t i = 0; i < B; ++i) {
                    int32_t nid = nodes[i];
                    grad(i) = -bls[nid] * res.grad_branches[nid];
                }
                grad(B)     = -params.beta * res.grad_omega;
                grad(B + 1) = -params.theta_AC * res.grad_theta_AC;
                grad(B + 2) = -params.theta_AT * res.grad_theta_AT;
                grad(B + 3) = -params.theta_CG * res.grad_theta_CG;
                grad(B + 4) = -params.theta_CT * res.grad_theta_CT;
                grad(B + 5) = -params.theta_GT * res.grad_theta_GT;

                return -res.log_likelihood;
            }
        };

        if (progress_cb) {
            progress_cb("Phase 2: Unified joint Fréchet quasi-Newton optimization", 0.45);
        }

        JointMG94Objective obj{cur_tree, aln, branch_node_ids, num_nodes, leaf_to_taxon};
        LBFGSpp::LBFGSBParam<Scalar> opt_param;
        opt_param.m = 8;
        opt_param.epsilon = tolerance;
        opt_param.max_iterations = max_iters;
        LBFGSpp::LBFGSBSolver<Scalar> solver(opt_param);
        Scalar fx = 0.0;
        int iters = 0;
        try {
            iters = solver.minimize(obj, x, fx, lb, ub);
        } catch (...) {}

        for (size_t i = 0; i < B; ++i) {
            cur_tree.nodes[branch_node_ids[i]].branch_length = std::exp(x(i));
        }
        p.beta = std::exp(x(B));
        p.theta_AC = std::exp(x(B + 1));
        p.theta_AT = std::exp(x(B + 2));
        p.theta_CG = std::exp(x(B + 3));
        p.theta_CT = std::exp(x(B + 4));
        p.theta_GT = std::exp(x(B + 5));

        FitResult res;
        res.log_likelihood = -fx;
        res.x_opt.resize(6);
        res.x_opt(0) = p.beta;
        res.x_opt(1) = p.theta_AC;
        res.x_opt(2) = p.theta_AT;
        res.x_opt(3) = p.theta_CG;
        res.x_opt(4) = p.theta_CT;
        res.x_opt(5) = p.theta_GT;
        res.param_names = {"omega", "theta_AC", "theta_AT", "theta_CG", "theta_CT", "theta_GT"};
        res.iterations = iters;
        res.converged = true;
        res.tree = cur_tree;
        res.params = p;
        return res;
    }
};

} // namespace hyphy::opt
