"""
optimizer.py – Stage B: thick-lens local optimisation via optiland's LeastSquares.

Refines a Stage-A starting design (minimum-thickness thickened prescription)
by simultaneously optimising surface radii, centre thicknesses, and air gap
while preserving the target EFL and respecting manufacturing minimum thicknesses.

The module is intentionally thin: it maps the AutoAchromat domain model onto
optiland's OptimizationProblem API, then hands off to scipy's Trust-Region
Reflective (TRF) solver.  No optical logic lives here.

Abstraction note
----------------
All optimisation work is encapsulated in ``optimize_optic()``.  The caller in
``builder.py`` uses this function via the ``OpticOptimizer`` protocol so the
backend can be swapped (e.g. to a torch-based solver) without touching the
rest of the pipeline.
"""

from __future__ import annotations

import math
import logging
from typing import Optional

import numpy as np
from optiland import optic as optic_module
from optiland.optimization.problem import OptimizationProblem
from optiland.optimization.optimizer.scipy import LeastSquares
from optiland.optimization.operand.operand import operand_registry

from ..models import ThickPrescription, Inputs
from ..thickening import lookup_t_edge_min, lookup_t_center_min

logger = logging.getLogger(__name__)


def _get_radius(op: optic_module.Optic, surf_idx: int) -> float:
    """Read a surface radius from the optic (handles numpy scalars/arrays)."""
    return float(np.asarray(op.surface_group.radii[surf_idx]).ravel()[0])


def _sag_local(R: float, a: float) -> float:
    """Signed sagitta at semi-aperture *a*; returns 0 for flat/degenerate."""
    if not math.isfinite(R) or abs(R) < 1e-12 or abs(R) <= a + 1e-6:
        return 0.0
    return R - math.copysign(math.sqrt(R * R - a * a), R)


# ---------------------------------------------------------------------------
# Custom edge-thickness operand
#
# optiland's built-in "edge_thickness" operand reads `surface.semi_aperture`,
# which is only populated after a full ray trace.  Because we call image_solve()
# (paraxial only) before optimisation, semi_aperture is None → TypeError.
# This custom operand avoids that by accepting the semi-aperture as an argument.
# ---------------------------------------------------------------------------

_CUSTOM_ET_OPERAND = "edge_thickness_fixed_a"


def _edge_thickness_fixed_a(optic, surface_number: int, semi_aperture: float):
    """Edge thickness between surface *surface_number* and the next surface.

    Uses the caller-supplied *semi_aperture* instead of the (potentially None)
    surface.semi_aperture attribute.
    """
    surf1 = optic.surface_group.surfaces[surface_number]
    surf2 = optic.surface_group.surfaces[surface_number + 1]
    sag1 = float(np.asarray(surf1.geometry.sag(y=semi_aperture)).ravel()[0])
    sag2 = float(np.asarray(surf2.geometry.sag(y=semi_aperture)).ravel()[0])
    thickness = float(np.asarray(
        optic.surface_group.surfaces[surface_number].thickness
    ).ravel()[0])
    return thickness - sag1 + sag2


# Register once (idempotent if already registered with the same function)
if _CUSTOM_ET_OPERAND not in operand_registry:
    operand_registry.register(_CUSTOM_ET_OPERAND, _edge_thickness_fixed_a)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def optimize_optic(
    op: optic_module.Optic,
    rx: ThickPrescription,
    inputs: Inputs,
) -> Optional[optic_module.Optic]:
    """Optimise a Stage-A optic via optiland's LeastSquares (TRF) optimiser.

    The optic is modified **in-place** and also returned for convenience.
    On failure, ``None`` is returned and the caller should fall back to the
    unoptimised Stage-A optic.

    Parameters
    ----------
    op :
        A fully configured Stage-A optiland ``Optic`` (aperture, fields,
        wavelengths set; ``image_solve`` already called).
    rx :
        The ``ThickPrescription`` used to build ``op``.  Used to read element
        thicknesses (for manufacturing bounds) and the air gap value.
    inputs :
        System specification supplying ``fprime``, ``D``, and ``air_gap``.

    Returns
    -------
    optiland.optic.Optic | None
        The same ``op`` object, modified in-place with optimised parameters.
        ``None`` on any failure.
    """
    try:
        return _run_optimization(op, rx, inputs)
    except Exception as exc:
        logger.debug(
            "Stage B optimization failed (%s): %s",
            rx.system_type,
            exc,
            exc_info=True,
        )
        return None


# ---------------------------------------------------------------------------
# Internal implementation
# ---------------------------------------------------------------------------


def _run_optimization(
    op: optic_module.Optic,
    rx: ThickPrescription,
    inputs: Inputs,
) -> Optional[optic_module.Optic]:
    """Core optimisation logic — raises on any error."""
    problem = OptimizationProblem()
    is_spaced = rx.system_type == "spaced"

    # ------------------------------------------------------------------
    # Surface layout (0-indexed, surface 0 = object at infinity):
    #
    #   Cemented: [0:obj] [1:R1 stop] [2:R2 cement] [3:R3] [4:image]
    #   Spaced:   [0:obj] [1:R1 stop] [2:R2] [3:R3]  [4:R4] [5:image]
    #
    # Thickness of surface k is the gap *after* surface k before surface k+1.
    # ------------------------------------------------------------------
    if is_spaced:
        radius_surfaces = [1, 2, 3, 4]
        # (element_index_in_rx, surface_index_for_thickness)
        glass_surfaces = [(0, 1), (1, 3)]
        air_gap_surface = 2
    else:  # cemented
        radius_surfaces = [1, 2, 3]
        glass_surfaces = [(0, 1), (1, 2)]
        air_gap_surface = None

    margin = inputs.D / 2.0 + 1.0  # minimum |R|: semi-aperture + 1 mm guard

    # ------------------------------------------------------------------
    # Radius variables — sign-preserving bounds prevent topological change
    # (e.g. biconvex → meniscus).  Flat surfaces (|R| ≈ ∞) are skipped.
    # Upper bound is at least 10× f' AND at least 5× the starting value.
    # ------------------------------------------------------------------
    n_radius_vars = 0
    for surf_idx in radius_surfaces:
        R = op.surface_group.radii[surf_idx]
        if not math.isfinite(R) or abs(R) < 1e-10:
            continue  # flat surface — do not optimise
        R_max = max(10.0 * inputs.fprime, abs(R) * 5.0)
        if R > 0:
            problem.add_variable(
                op, "radius",
                surface_number=surf_idx,
                min_val=+margin,
                max_val=+R_max,
            )
        else:
            problem.add_variable(
                op, "radius",
                surface_number=surf_idx,
                min_val=-R_max,
                max_val=-margin,
            )
        n_radius_vars += 1

    if n_radius_vars == 0:
        logger.debug("Stage B: no optimisable radius variables — returning as-is")
        return op

    # ------------------------------------------------------------------
    # Centre-thickness variables
    #
    # For a positive lens the manufacturing limit is on t_edge:
    #   t_edge = t_center - sag(R_front) + sag(R_back) >= te_min
    #   → t_center_min = te_min + sag(R_front) - sag(R_back)
    #
    # The sag values are evaluated at the Stage-A (starting) radii so that
    # the lower bound is geometrically correct at the initial design point.
    # An edge_thickness inequality operand (added below) provides a dynamic
    # soft constraint that keeps t_edge valid as radii change during optimisation.
    #
    # For a negative lens the limit is directly on t_center.
    # ------------------------------------------------------------------
    D = inputs.D
    a = D / 2.0
    te_min = lookup_t_edge_min(D)
    tc_min_neg = lookup_t_center_min(D)

    for elem_idx, surf_idx in glass_surfaces:
        elem = rx.elements[elem_idx]
        is_positive = elem.t_center >= elem.t_edge

        if is_positive:
            # Correct minimum t_center from sag geometry at Stage-A radii
            R_f = _get_radius(op, surf_idx)
            R_b = _get_radius(op, surf_idx + 1)   # back surface radius
            sag_f = _sag_local(R_f, a)
            sag_b = _sag_local(R_b, a)
            t_min = max(te_min + sag_f - sag_b, 0.5)
        else:
            t_min = max(tc_min_neg, 0.5)

        t_max = max(t_min * 5.0, 20.0)
        problem.add_variable(
            op, "thickness",
            surface_number=surf_idx,
            min_val=t_min,
            max_val=t_max,
        )

    # ------------------------------------------------------------------
    # Air-gap variable (spaced doublet only)
    # ------------------------------------------------------------------
    if is_spaced and air_gap_surface is not None:
        gap_current = rx.air_gap if rx.air_gap is not None else inputs.air_gap
        gap_min = max(gap_current * 0.1, 1.0)
        gap_max = gap_current * 20.0
        problem.add_variable(
            op, "thickness",
            surface_number=air_gap_surface,
            min_val=gap_min,
            max_val=gap_max,
        )

    # ------------------------------------------------------------------
    # Operands
    # ------------------------------------------------------------------
    # Image-plane surface index (last surface of the system)
    image_surf = op.surface_group.num_surfaces - 1

    # ------------------------------------------------------------------
    # Edge-thickness soft constraints — custom operand that avoids the
    # semi_aperture=None bug in optiland's built-in edge_thickness.
    # Applied to ALL elements (not just positive ones) because radii can
    # change sign of (sag_front - sag_back) during optimisation.
    # ------------------------------------------------------------------
    for elem_idx, surf_idx in glass_surfaces:
        problem.add_operand(
            _CUSTOM_ET_OPERAND,
            min_val=te_min,
            weight=5,
            input_data={
                "optic": op,
                "surface_number": surf_idx,
                "semi_aperture": a,
            },
        )

    # ------------------------------------------------------------------
    # Snapshot Stage-A variable values AND on-axis RMS so we can revert
    # if the optimiser worsens the primary performance metric.
    # The on-axis RMS operand is the first "rms_spot_size" we add below.
    # ------------------------------------------------------------------
    x0_saved = [var.value for var in problem.variables]

    # Index of the on-axis RMS operand (added next — track before we add it)
    _on_axis_rms_idx = len(problem.operands)

    # EFL constraint — highest priority, preserves the target focal length.
    problem.add_operand(
        "f2",
        target=inputs.fprime,
        weight=10,
        input_data={"optic": op},
    )

    # RMS spot radius — on-axis, polychromatic.  Primary aberration objective.
    _on_axis_rms_idx = len(problem.operands)
    problem.add_operand(
        "rms_spot_size",
        target=0.0,
        weight=4,
        input_data={
            "optic": op,
            "surface_number": image_surf,
            "Hx": 0.0,
            "Hy": 0.0,
            "num_rays": 12,
            "wavelength": "all",
        },
    )

    # RMS spot radius — off-axis (1° field), polychromatic.
    # Lower weight than on-axis so the optimiser does not sacrifice on-axis
    # performance for off-axis gains.
    problem.add_operand(
        "rms_spot_size",
        target=0.0,
        weight=1,
        input_data={
            "optic": op,
            "surface_number": image_surf,
            "Hx": 0.0,
            "Hy": 1.0,
            "num_rays": 12,
            "wavelength": "all",
        },
    )

    # Record Stage-A on-axis RMS BEFORE running the optimiser
    try:
        rms_initial = float(problem.operands[_on_axis_rms_idx].value)
    except Exception:
        rms_initial = float("inf")

    # ------------------------------------------------------------------
    # Optimise (TRF: Trust-Region Reflective — supports variable bounds)
    # ------------------------------------------------------------------
    optimizer = LeastSquares(problem)
    result = optimizer.optimize(maxiter=500, tol=1e-6, method_choice="trf")

    # Check on-axis RMS after optimisation
    try:
        rms_final = float(problem.operands[_on_axis_rms_idx].value)
    except Exception:
        rms_final = float("inf")

    logger.debug(
        "Stage B (%s): status=%s nfev=%d rms %.4f→%.4f mm message=%s",
        rx.system_type,
        result.status,
        result.nfev,
        rms_initial,
        rms_final,
        result.message,
    )

    # Revert to Stage-A if on-axis RMS got worse (>2 % tolerance)
    if rms_final > rms_initial * 1.02:
        logger.debug("Stage B: on-axis RMS worsened (%.4f > %.4f) — restoring Stage-A",
                     rms_final, rms_initial)
        for var, x in zip(problem.variables, x0_saved):
            var.update(x)
        problem.update_optics()

    # Recompute image-plane position (back focal distance) after optimisation
    op.image_solve()

    return op
