"""
PyTorch Autograd Integration for HyPhy 3
Enables end-to-end backpropagation through phylogenetic tree likelihoods using
exact analytical Inside-Outside adjoint states.
"""

import os
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

try:
    import torch
    TORCH_AVAILABLE = True
except ImportError:
    TORCH_AVAILABLE = False

if TORCH_AVAILABLE:
    import _hyphy3

    class PhylogeneticTreeLogLikelihood(torch.autograd.Function):
        """
        Differentiable phylogenetic tree log-likelihood function.
        Forward: Evaluates exact log-likelihood ln L(t) under MG94.
        Backward: Evaluates exact analytical gradients d ln L / d t_b via Inside-Outside adjoints.
        """
        @staticmethod
        def forward(ctx, branch_lengths, tree, aln, params):
            """
            branch_lengths: torch.Tensor of shape (num_nodes,) representing branch lengths.
            tree: _hyphy3.Tree
            aln: _hyphy3.Alignment
            params: _hyphy3.MG94Parameters
            """
            ctx.tree = tree
            ctx.aln = aln
            ctx.params = params

            bl_list = branch_lengths.detach().cpu().tolist()
            ll, grads = _hyphy3.compute_branch_length_gradients(tree, aln, params, bl_list)
            ctx.grads = torch.tensor(grads, dtype=branch_lengths.dtype, device=branch_lengths.device)

            return torch.tensor(ll, dtype=branch_lengths.dtype, device=branch_lengths.device)

        @staticmethod
        def backward(ctx, grad_output):
            # Chain rule: grad_output * (d ln L / d t)
            grad_branch_lengths = grad_output * ctx.grads
            return grad_branch_lengths, None, None, None

    def tree_log_likelihood(branch_lengths, tree, aln, params=None):
        """
        Convenience wrapper returning a differentiable scalar torch.Tensor.
        """
        if params is None:
            params = _hyphy3.MG94Parameters()
        return PhylogeneticTreeLogLikelihood.apply(branch_lengths, tree, aln, params)
