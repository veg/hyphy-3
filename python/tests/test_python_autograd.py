import os
import sys
import math
import unittest

# Ensure hyphy3 package is discoverable
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

import hyphy3

try:
    import torch
    from hyphy3.autograd import tree_log_likelihood
    TORCH_AVAILABLE = True
except ImportError:
    TORCH_AVAILABLE = False


class TestHyPhy3Python(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        candidate_paths = [
            os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "benchmarks", "data", "adh.nex")),
            os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", "tests", "data", "adh.nex")),
        ]
        cls.adh_path = next((p for p in candidate_paths if os.path.exists(p)), candidate_paths[0])
        cls.assertTrue(os.path.exists(cls.adh_path), f"Test data not found at {cls.adh_path}")

    def test_version_and_imports(self):
        self.assertEqual(hyphy3.__version__, "3.0.0")
        code = hyphy3.GeneticCode.universal()
        self.assertEqual(code.name, "Universal")
        self.assertEqual(code.num_sense_codons, 61)

    def test_alignment_and_tree(self):
        aln = hyphy3.Alignment.load(self.adh_path)
        self.assertEqual(aln.num_taxa, 23)
        self.assertEqual(aln.num_codons, 254)
        self.assertTrue(aln.num_patterns > 0)
        self.assertTrue(len(aln.embedded_tree_newick) > 0)

        tree = hyphy3.Tree.from_newick(aln.embedded_tree_newick)
        self.assertEqual(tree.num_leaves(), 23)
        self.assertEqual(tree.num_nodes(), 44)
        self.assertEqual(len(tree.get_branch_lengths()), 44)

    def test_fel_analyzer(self):
        aln = hyphy3.Alignment.load(self.adh_path)
        tree = hyphy3.Tree.from_newick(aln.embedded_tree_newick)

        # Full unconstrained MG94 fit (default, matching HyPhy 2.5 standard fit)
        fel = hyphy3.FELAnalyzer.create_and_fit(tree, aln, 0.1)
        fel.run()
        self.assertAlmostEqual(fel.global_log_l, -4686.73, delta=5.0)
        self.assertEqual(len(fel.site_results), 254)

        # Quick scaled MG94 fit (full_model=False)
        tree2 = hyphy3.Tree.from_newick(aln.embedded_tree_newick)
        fel_quick = hyphy3.FELAnalyzer.create_and_fit(tree2, aln, 0.1, full_model=False)
        self.assertAlmostEqual(fel_quick.global_log_l, -4769.79, delta=5.0)

    def test_busted_analyzer(self):
        aln = hyphy3.Alignment.load(self.adh_path)
        tree = hyphy3.Tree.from_newick(aln.embedded_tree_newick)

        busted = hyphy3.BUSTEDAnalyzer.create_and_fit(tree, aln)
        settings = hyphy3.BUSTEDSettings()
        settings.auto_select_k = True
        settings.refine_branch_lengths = True
        settings.max_k = 2

        res = busted.run(settings)
        self.assertIn(res.optimal_k, [1, 2])
        self.assertTrue(res.unconstrained.log_likelihood >= res.constrained.log_likelihood - 1e-4)
        self.assertTrue(res.runtime_seconds > 0.0)

    def test_absrel_analyzer(self):
        bglobin_candidates = [
            os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "benchmarks", "data", "bglobin.nex")),
            os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", "tests", "data", "bglobin.nex")),
        ]
        bglobin_path = next((p for p in bglobin_candidates if os.path.exists(p)), None)
        if bglobin_path is None:
            self.skipTest("bglobin.nex not found")

        aln = hyphy3.Alignment.load(bglobin_path)
        tree = hyphy3.Tree.from_newick(aln.embedded_tree_newick)

        settings = hyphy3.ABSRELSettings()
        settings.max_rate_classes = 2
        settings.p_threshold = 0.05

        absrel = hyphy3.ABSRELAnalyzer.create(tree, aln, settings)
        res = absrel.run()

        self.assertAlmostEqual(res.gtr_fit.log_likelihood, -3926.03, delta=5.0)
        self.assertAlmostEqual(res.baseline_fit.log_likelihood, -3765.93, delta=15.0)
        self.assertTrue(res.full_adaptive_fit.log_likelihood >= res.baseline_fit.log_likelihood)
        self.assertEqual(len(res.tested_branches), 31)
        self.assertTrue(len(res.positive_branches) > 0)
        self.assertTrue(res.runtime_seconds > 0.0)

        legacy_str = res.to_legacy_json(tree, aln)
        self.assertIn("Baseline MG94xREV", legacy_str)
        self.assertIn("Full adaptive model", legacy_str)

        modern_str = res.to_modern_json(tree, aln)
        self.assertIn("baseline_mg94", modern_str)
        self.assertIn("full_adaptive", modern_str)
        self.assertIn("branch_results", modern_str)

    @unittest.skipUnless(TORCH_AVAILABLE, "PyTorch required for autograd tests")
    def test_autograd_gradients_vs_finite_differences(self):
        aln = hyphy3.Alignment.load(self.adh_path)
        tree = hyphy3.Tree.from_newick(aln.embedded_tree_newick)
        params = hyphy3.MG94Parameters()
        params.alpha = 1.0
        params.beta = 0.4

        # Non-zero realistic branch lengths for gradient evaluation
        init_lengths = [0.05 if i != tree.root_id else 0.0 for i in range(tree.num_nodes())]
        bl = torch.tensor(init_lengths, dtype=torch.float64, requires_grad=True)

        # Forward pass
        ll = tree_log_likelihood(bl, tree, aln, params)
        self.assertTrue(torch.isfinite(ll))

        # Backward pass
        ll.backward()
        ana_grad = bl.grad.clone()

        # Compare with finite differences on all non-root branches
        eps = 1e-6
        for nid in range(tree.num_nodes()):
            if nid == tree.root_id:
                continue

            bl_plus = list(init_lengths)
            bl_plus[nid] += eps
            ll_plus, _ = hyphy3.compute_branch_length_gradients(tree, aln, params, bl_plus)

            bl_minus = list(init_lengths)
            bl_minus[nid] -= eps
            ll_minus, _ = hyphy3.compute_branch_length_gradients(tree, aln, params, bl_minus)

            num_grad = (ll_plus - ll_minus) / (2.0 * eps)
            ana_val = ana_grad[nid].item()

            denom = max(abs(ana_val), abs(num_grad), 1e-4)
            rel_err = abs(ana_val - num_grad) / denom
            self.assertLess(rel_err, 1e-5, f"Gradient mismatch at node {nid}: ana={ana_val}, num={num_grad}, err={rel_err}")

    @unittest.skipUnless(TORCH_AVAILABLE, "PyTorch required for autograd tests")
    def test_pytorch_adam_branch_optimization(self):
        aln = hyphy3.Alignment.load(self.adh_path)
        tree = hyphy3.Tree.from_newick(aln.embedded_tree_newick)
        params = hyphy3.MG94Parameters()
        params.alpha = 1.0
        params.beta = 0.5

        # Parameterize in log space
        init_t = [0.05 if i != tree.root_id else 0.0 for i in range(tree.num_nodes())]
        log_t = torch.tensor([math.log(max(x, 1e-6)) for x in init_t], dtype=torch.float64, requires_grad=True)

        optimizer = torch.optim.Adam([log_t], lr=0.05)
        initial_ll = tree_log_likelihood(torch.exp(log_t), tree, aln, params).item()

        for step in range(5):
            optimizer.zero_grad()
            loss = -tree_log_likelihood(torch.exp(log_t), tree, aln, params)
            loss.backward()
            with torch.no_grad():
                log_t.grad[tree.root_id] = 0.0
            optimizer.step()

        final_ll = tree_log_likelihood(torch.exp(log_t), tree, aln, params).item()
        # Likelihood must strictly improve
        self.assertGreater(final_ll, initial_ll + 50.0)


if __name__ == "__main__":
    unittest.main()
