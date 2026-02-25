"""
Optiland bridge – builds optiland Optic objects from thin-lens Candidates
and evaluates their optical performance.

Stage A: placeholder thicknesses  (fast, 100 % trace success)
Stage B: thickened prescription   (future)
"""

from .builder import build_optic
from .evaluator import evaluate, batch_evaluate

__all__ = [
    "build_optic",
    "evaluate",
    "batch_evaluate",
]
