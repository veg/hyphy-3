"""
HyPhy 3: Next-Generation High-Performance Molecular Evolution & Phylogenetics
"""
from _hyphy3 import (
    GeneticCode,
    Tree,
    TreeNode,
    Alignment,
    MG94Parameters,
    FELSiteResult,
    FELAnalyzer,
    MEMESiteResult,
    MEMEAnalyzer,
    BUSTEDSettings,
    BUSTEDRateDistribution,
    BUSTEDFit,
    BUSTEDResult,
    BUSTEDAnalyzer,
    ABSRELSettings,
    ABSRELRateDistribution,
    ABSRELBranchResult,
    ABSRELFitSummary,
    ABSRELResult,
    ABSRELAnalyzer,
    compute_branch_length_gradients
)

from . import autograd

__version__ = "3.0.0"
__all__ = [
    "GeneticCode",
    "Tree",
    "TreeNode",
    "Alignment",
    "MG94Parameters",
    "FELSiteResult",
    "FELAnalyzer",
    "MEMESiteResult",
    "MEMEAnalyzer",
    "BUSTEDSettings",
    "BUSTEDRateDistribution",
    "BUSTEDFit",
    "BUSTEDResult",
    "BUSTEDAnalyzer",
    "ABSRELSettings",
    "ABSRELRateDistribution",
    "ABSRELBranchResult",
    "ABSRELFitSummary",
    "ABSRELResult",
    "ABSRELAnalyzer",
    "compute_branch_length_gradients",
    "autograd"
]
