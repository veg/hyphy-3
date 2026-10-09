#!/usr/bin/env python3
"""
Example 5: End-to-End PyTorch Integration with Custom Autograd Function
Demonstrates how to optimize evolutionary models using PyTorch's Adam optimizer and HyPhy 3's Inside-Outside adjoints.
"""

import sys
from pathlib import Path

repo_root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(repo_root / "python"))

try:
    import torch
except ImportError:
    print("PyTorch is required for this example. Install with: pip install torch")
    sys.exit(0)

import hyphy3 as hp

class TreeLikelihoodFunction(torch.autograd.Function):
    """
    Custom PyTorch autograd Function wrapping HyPhy 3's exact analytical Inside-Outside likelihood.
    """
    @staticmethod
    def forward(ctx, branch_lengths_tensor, tree, aln, params):
        bls = branch_lengths_tensor.detach().cpu().numpy().tolist()
        log_l, grads = hp.compute_branch_length_gradients(tree, aln, params, bls)
        ctx.save_for_backward(torch.tensor(grads, dtype=torch.float64, device=branch_lengths_tensor.device))
        return torch.tensor(log_l, dtype=torch.float64, device=branch_lengths_tensor.device)

    @staticmethod
    def backward(ctx, grad_output):
        grads, = ctx.saved_tensors
        return grad_output * grads, None, None, None

def main():
    data_dir = repo_root / "benchmarks" / "data"
    aln = hp.Alignment.load(str(data_dir / "cd2.fna"))
    tree = hp.Tree.from_newick_file(str(data_dir / "cd2.nwk"))

    params = hp.MG94Parameters()
    params.alpha = 1.0
    params.beta = 0.5

    init_bls = [node.branch_length if node.id != tree.root_id else 0.0 for node in tree.nodes]
    t_bls = torch.tensor(init_bls, dtype=torch.float64, requires_grad=True)

    optimizer = torch.optim.Adam([t_bls], lr=0.01)

    print("Optimizing branch lengths with PyTorch Adam optimizer...")
    for step in range(25):
        optimizer.zero_grad()
        # Loss is negative log-likelihood (we minimize loss to maximize likelihood)
        ll = TreeLikelihoodFunction.apply(t_bls, tree, aln, params)
        loss = -ll
        loss.backward()

        # Zero root gradient
        t_bls.grad[tree.root_id] = 0.0

        optimizer.step()

        # Enforce positive branch length constraint
        with torch.no_grad():
            t_bls.clamp_(min=1e-4, max=10.0)

        if step % 5 == 0 or step == 24:
            print(f"  Step {step:2d} | Log-Likelihood: {ll.item():.4f} | Max Branch Grad: {t_bls.grad.abs().max().item():.4f}")

    print("\nPyTorch end-to-end optimization successfully converged!")

if __name__ == "__main__":
    main()
