"""
builder.py – Construct an optiland ``Optic`` from a thin-lens ``Candidate``.

Stage A  (current):  placeholder thicknesses, fast trace.
Stage B  (future) :  physically thickened prescription (same interface).
"""

from __future__ import annotations

import math
import logging
from typing import Optional

from optiland import optic
from optiland.materials import AbbeMaterial

from ..models import Candidate, Inputs

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Thickness helpers
# ---------------------------------------------------------------------------

_MIN_CENTER = 2.0  # mm – absolute minimum centre thickness
_MIN_EDGE = 0.5  # mm – absolute minimum edge thickness


def _sag(R: float, h: float) -> float:
    """Signed sagitta of a sphere of radius *R* at height *h*.

    Returns the z-offset of the surface at height *h* relative to the
    vertex (z = 0).  Positive for R > 0, negative for R < 0.
    If |R| < h the surface cannot accommodate the aperture; returns the
    full sag (R itself) as a fallback.
    """
    if not math.isfinite(R) or abs(R) < 1e-12:
        return 0.0
    if abs(R) < h:
        # Surface too curved – can't even reach this height; return max sag
        return R  # sign preserving; caller will get a large thickness
    return R - math.copysign(math.sqrt(R * R - h * h), R)


def _element_center_thickness(R_front: float, R_back: float, h: float) -> float:
    """Minimum safe centre thickness for a single lens element.

    Ensures both positive centre *and* positive edge thickness, with
    margins ``_MIN_CENTER`` and ``_MIN_EDGE``.

    ``h`` is the semi-aperture (D / 2).
    """
    sag_f = _sag(R_front, h)  # positive for convex-toward-object
    sag_b = _sag(R_back, h)  # negative for concave-toward-image

    # Edge thickness = t_center + sag_back - sag_front
    # We need  t_center >= sag_front - sag_back + _MIN_EDGE
    t_from_edge = (sag_f - sag_b) + _MIN_EDGE

    return max(t_from_edge, _MIN_CENTER)


def _placeholder_air_gap(D: float) -> float:
    """Air gap between two singlets (spaced doublet only).

    Rule: max(D / 50, 0.5) mm.
    """
    return max(D / 50.0, 0.5)


def _initial_back_focus(fprime: float) -> float:
    """Initial guess for the back focal distance before ``image_solve``.

    ``image_solve()`` will correct this, but a reasonable starting value
    helps avoid occasional tracing failures during the solve.
    """
    return abs(fprime) * 0.9


# Optiland's ray-surface intersection has a numerical precision issue
# that causes NaN when radius values carry too many decimal places.
# Rounding to 6 decimal places (nanometre precision for mm-scale radii)
# avoids this while being optically insignificant.
_RADIUS_DP = 6


def _safe_radius(r: float) -> float:
    """Round a radius value to ``_RADIUS_DP`` decimal places."""
    if not math.isfinite(r):
        return r
    return round(r, _RADIUS_DP)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def build_optic(
    candidate: Candidate,
    inputs: Inputs,
    *,
    stage: str = "A",
) -> Optional[optic.Optic]:
    """Build an optiland ``Optic`` from a thin-lens *candidate*.

    Parameters
    ----------
    candidate : Candidate
        Output of ``run_cemented`` or ``run_spaced``.
    inputs : Inputs
        Same config that drove the thin-lens synthesis.
    stage : str
        ``"A"`` for placeholder thicknesses (default).
        ``"B"`` reserved for future thickened prescription.

    Returns
    -------
    optiland.optic.Optic | None
        The constructed system, or *None* if the model could not be built
        (e.g. trace failure, degenerate radii).
    """
    if stage == "B":
        raise NotImplementedError("Stage B (thickened) not yet implemented")

    if candidate.system_type == "cemented":
        return _build_cemented_a(candidate, inputs)
    elif candidate.system_type == "spaced":
        return _build_spaced_a(candidate, inputs)
    else:
        logger.warning("Unknown system_type=%s", candidate.system_type)
        return None


# ---------------------------------------------------------------------------
# Cemented doublet  (3 radii: R1, R2, R3)
# ---------------------------------------------------------------------------
#   Surface 0 : Object     (flat, thickness → ∞)
#   Surface 1 : R1         thickness = t1, material = glass1, aperture stop
#   Surface 2 : R2         thickness = t2, material = glass2
#   Surface 3 : R3         thickness = back-focus  (air)
#   Surface 4 : Image      (flat)
# ---------------------------------------------------------------------------


def _build_cemented_a(cand: Candidate, inp: Inputs) -> Optional[optic.Optic]:
    R1 = _safe_radius(cand.radii[0])
    R2 = _safe_radius(cand.radii[1])
    R3 = _safe_radius(cand.radii[2])
    h = inp.D / 2.0
    t1 = _element_center_thickness(R1, R2, h)
    t2 = _element_center_thickness(R2, R3, h)
    bf = _initial_back_focus(inp.fprime)

    mat1 = AbbeMaterial(cand.n1, cand.nu1)
    mat2 = AbbeMaterial(cand.n2, cand.nu2)

    try:
        op = optic.Optic()

        # Surface 0 – Object (at infinity)
        op.add_surface(index=0, radius=float("inf"), thickness=1e10)

        # Surface 1 – first element front face (aperture stop)
        op.add_surface(
            index=1,
            radius=R1,
            thickness=t1,
            material=mat1,
            is_stop=True,
        )

        # Surface 2 – cemented interface
        op.add_surface(
            index=2,
            radius=R2,
            thickness=t2,
            material=mat2,
        )

        # Surface 3 – last surface of the doublet
        op.add_surface(
            index=3,
            radius=R3,
            thickness=bf,
        )

        # Surface 4 – Image plane
        op.add_surface(index=4, radius=float("inf"), thickness=0.0)

        _configure_system(op, inp)
        return op

    except Exception:
        logger.debug(
            "Failed to build cemented optic for %s:%s + %s:%s",
            cand.catalog1,
            cand.glass1,
            cand.catalog2,
            cand.glass2,
            exc_info=True,
        )
        return None


# ---------------------------------------------------------------------------
# Air-spaced doublet  (4 radii: R1, R2, R3, R4)
# ---------------------------------------------------------------------------
#   Surface 0 : Object     (flat, thickness → ∞)
#   Surface 1 : R1         thickness = t1, material = glass1, aperture stop
#   Surface 2 : R2         thickness = air_gap  (air)
#   Surface 3 : R3         thickness = t2, material = glass2
#   Surface 4 : R4         thickness = back-focus  (air)
#   Surface 5 : Image      (flat)
# ---------------------------------------------------------------------------


def _build_spaced_a(cand: Candidate, inp: Inputs) -> Optional[optic.Optic]:
    R1 = _safe_radius(cand.radii[0])
    R2 = _safe_radius(cand.radii[1])
    R3 = _safe_radius(cand.radii[2])
    R4 = _safe_radius(cand.radii[3])
    h = inp.D / 2.0
    t1 = _element_center_thickness(R1, R2, h)
    t2 = _element_center_thickness(R3, R4, h)
    gap = _placeholder_air_gap(inp.D)
    bf = _initial_back_focus(inp.fprime)

    mat1 = AbbeMaterial(cand.n1, cand.nu1)
    mat2 = AbbeMaterial(cand.n2, cand.nu2)

    try:
        op = optic.Optic()

        # Surface 0 – Object (at infinity)
        op.add_surface(index=0, radius=float("inf"), thickness=1e10)

        # Surface 1 – lens 1 front face (aperture stop)
        op.add_surface(
            index=1,
            radius=R1,
            thickness=t1,
            material=mat1,
            is_stop=True,
        )

        # Surface 2 – lens 1 rear face → air gap
        op.add_surface(
            index=2,
            radius=R2,
            thickness=gap,
        )

        # Surface 3 – lens 2 front face
        op.add_surface(
            index=3,
            radius=R3,
            thickness=t2,
            material=mat2,
        )

        # Surface 4 – lens 2 rear face → back-focus
        op.add_surface(
            index=4,
            radius=R4,
            thickness=bf,
        )

        # Surface 5 – Image plane
        op.add_surface(index=5, radius=float("inf"), thickness=0.0)

        _configure_system(op, inp)
        return op

    except Exception:
        logger.debug(
            "Failed to build spaced optic for %s:%s + %s:%s",
            cand.catalog1,
            cand.glass1,
            cand.catalog2,
            cand.glass2,
            exc_info=True,
        )
        return None


# ---------------------------------------------------------------------------
# Shared configuration (aperture, field, wavelengths, image solve)
# ---------------------------------------------------------------------------


def _configure_system(op: optic.Optic, inp: Inputs) -> None:
    """Set aperture, field, wavelengths, and run image_solve."""

    # Aperture  –  entrance-pupil diameter = D (from inputs)
    op.set_aperture(aperture_type="EPD", value=inp.D)

    # Field  –  on-axis only for Stage A (fast)
    op.set_field_type(field_type="angle")
    op.add_field(y=0.0)

    # Wavelengths  –  three spectral lines from the synthesis
    #   lam1 (short), lam0 (centre / primary), lam2 (long)
    op.add_wavelength(value=inp.lam1)  # short wavelength
    op.add_wavelength(value=inp.lam0, is_primary=True)  # primary / design
    op.add_wavelength(value=inp.lam2)  # long wavelength

    # Move image plane to paraxial focus
    op.image_solve()
