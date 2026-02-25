"""
builder.py – Construct an optiland ``Optic`` from a ``ThickPrescription``.

This module is the thin adapter between the pure-math ``thickening`` module
and the optiland ray-tracing library.  It contains **no** optical design
logic – only the translation of a ``ThickPrescription`` into optiland API
calls.
"""

from __future__ import annotations

import math
import logging
from typing import Optional

from optiland import optic
from optiland.materials import AbbeMaterial

from ..models import Candidate, Inputs, ThickPrescription
from ..thickening import thicken

logger = logging.getLogger(__name__)

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

    rx = thicken(candidate, inputs)
    if rx is None:
        return None
    return _build_from_prescription(rx)


def build_optic_from_prescription(
    rx: ThickPrescription,
) -> Optional[optic.Optic]:
    """Build an optiland ``Optic`` directly from a ``ThickPrescription``.

    Useful when the caller has already computed the prescription via
    ``thickening.thicken()`` and wants to avoid recomputing it.
    """
    return _build_from_prescription(rx)


# ---------------------------------------------------------------------------
# Internal: ThickPrescription → optiland Optic
# ---------------------------------------------------------------------------


def _build_from_prescription(
    rx: ThickPrescription,
) -> Optional[optic.Optic]:
    """Translate a ``ThickPrescription`` into an optiland ``Optic``.

    Cemented doublet  (5 surfaces):
        Object | R1(stop) | R2(cemented) | R3 | Image

    Air-spaced doublet  (6 surfaces):
        Object | R1(stop) | R2 | R3 | R4 | Image
    """
    try:
        if rx.system_type == "cemented":
            return _build_cemented(rx)
        elif rx.system_type == "spaced":
            return _build_spaced(rx)
        else:
            logger.warning("Unknown system_type=%s in prescription", rx.system_type)
            return None
    except Exception:
        logger.debug(
            "Failed to build optic from prescription (type=%s)",
            rx.system_type,
            exc_info=True,
        )
        return None


def _build_cemented(rx: ThickPrescription) -> optic.Optic:
    """Build cemented doublet: 5 surfaces."""
    e1, e2 = rx.elements[0], rx.elements[1]

    R1 = _safe_radius(e1.R_front)
    R2 = _safe_radius(e1.R_back)  # == e2.R_front (cemented)
    R3 = _safe_radius(e2.R_back)

    mat1 = AbbeMaterial(e1.nd, e1.vd)
    mat2 = AbbeMaterial(e2.nd, e2.vd)

    op = optic.Optic()

    # Surface 0 – Object (at infinity)
    op.add_surface(index=0, radius=float("inf"), thickness=1e10)

    # Surface 1 – first element front face (aperture stop)
    op.add_surface(
        index=1,
        radius=R1,
        thickness=e1.t_center,
        material=mat1,
        is_stop=True,
    )

    # Surface 2 – cemented interface
    op.add_surface(
        index=2,
        radius=R2,
        thickness=e2.t_center,
        material=mat2,
    )

    # Surface 3 – last surface of the doublet
    op.add_surface(
        index=3,
        radius=R3,
        thickness=rx.back_focus_guess,
    )

    # Surface 4 – Image plane
    op.add_surface(index=4, radius=float("inf"), thickness=0.0)

    _configure_system(op, rx)
    return op


def _build_spaced(rx: ThickPrescription) -> optic.Optic:
    """Build air-spaced doublet: 6 surfaces."""
    e1, e2 = rx.elements[0], rx.elements[1]

    R1 = _safe_radius(e1.R_front)
    R2 = _safe_radius(e1.R_back)
    R3 = _safe_radius(e2.R_front)
    R4 = _safe_radius(e2.R_back)

    mat1 = AbbeMaterial(e1.nd, e1.vd)
    mat2 = AbbeMaterial(e2.nd, e2.vd)

    assert rx.air_gap is not None, "Spaced doublet must have an air_gap"

    op = optic.Optic()

    # Surface 0 – Object (at infinity)
    op.add_surface(index=0, radius=float("inf"), thickness=1e10)

    # Surface 1 – lens 1 front face (aperture stop)
    op.add_surface(
        index=1,
        radius=R1,
        thickness=e1.t_center,
        material=mat1,
        is_stop=True,
    )

    # Surface 2 – lens 1 rear face → air gap
    op.add_surface(
        index=2,
        radius=R2,
        thickness=rx.air_gap,
    )

    # Surface 3 – lens 2 front face
    op.add_surface(
        index=3,
        radius=R3,
        thickness=e2.t_center,
        material=mat2,
    )

    # Surface 4 – lens 2 rear face → back-focus
    op.add_surface(
        index=4,
        radius=R4,
        thickness=rx.back_focus_guess,
    )

    # Surface 5 – Image plane
    op.add_surface(index=5, radius=float("inf"), thickness=0.0)

    _configure_system(op, rx)
    return op


# ---------------------------------------------------------------------------
# Shared configuration (aperture, field, wavelengths, image solve)
# ---------------------------------------------------------------------------


def _configure_system(op: optic.Optic, rx: ThickPrescription) -> None:
    """Set aperture, field, wavelengths, and run image_solve."""

    # Aperture  –  entrance-pupil diameter = D
    op.set_aperture(aperture_type="EPD", value=rx.D)

    # Field  –  on-axis only for Stage A (fast)
    op.set_field_type(field_type="angle")
    op.add_field(y=0.0)

    # Wavelengths  –  three spectral lines
    #   wavelengths = (lam1_short, lam0_primary, lam2_long)
    lam1, lam0, lam2 = rx.wavelengths
    op.add_wavelength(value=lam1)  # short wavelength
    op.add_wavelength(value=lam0, is_primary=True)  # primary / design
    op.add_wavelength(value=lam2)  # long wavelength

    # Move image plane to paraxial focus
    op.image_solve()
