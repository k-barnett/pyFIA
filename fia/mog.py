"""
Backward-compatible path; MOG engine implementation lives in ``fia_mog.engine``.
"""

from fia_mog.engine import ConditionContext, MOGEngine, TreeMetrics, compute_tree_metrics

__all__ = ["ConditionContext", "MOGEngine", "TreeMetrics", "compute_tree_metrics"]
