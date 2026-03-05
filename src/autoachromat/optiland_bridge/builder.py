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
from dataclasses import replace
from typing import Optional

import numpy as np
from optiland import optic
from optiland.materials import AbbeMaterial
from optiland.materials.base import BaseMaterial

from ..models import Candidate, Inputs, ThickPrescription, ElementRx
from ..optics import refractive_index
from ..glass_reader import Glass
from ..thickening import thicken
from .optimizer import optimize_optic

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Custom material: uses AGF dispersion formula for arbitrary wavelengths
# ---------------------------------------------------------------------------


class AGFMaterial(BaseMaterial):
    """Material that uses AGF catalog dispersion coefficients.

    Unlike ``AbbeMaterial`` (limited to 0.38–0.75 µm), this class can
    compute *n(λ)* at any wavelength supported by the dispersion formula,
    making it suitable for IR / UV designs.
    """

    def __init__(self, formula_id: int, cd: list[float], name: str = ""):
        super().__init__()
        self._formula_id = formula_id
        self._cd = list(cd)
        self._glass = Glass(name=name, catalog="", formula_id=formula_id, cd=list(cd))

    def _calculate_n(self, wavelength, **kwargs):
        """Compute n(λ) using the AGF dispersion formula."""
        wavelength = np.atleast_1d(np.asarray(wavelength, dtype=float))
        result = np.array(
            [refractive_index(self._glass, float(w)) for w in wavelength.ravel()]
        ).reshape(wavelength.shape)
        if result.size == 1:
            return float(result.ravel()[0])
        return result

    def _calculate_k(self, wavelength, **kwargs):
        """Extinction coefficient – always 0 for catalog glasses."""
        return np.zeros_like(np.asarray(wavelength, dtype=float)) * 0.0

    def to_dict(self):
        d = super().to_dict()
        d.update(
            {
                "formula_id": self._formula_id,
                "cd": self._cd,
            }
        )
        return d

    @classmethod
    def from_dict(cls, data):
        return cls(data["formula_id"], data["cd"])


def _make_material(elem: ElementRx) -> BaseMaterial:
    """Create the best available optiland material for an element.

    If AGF dispersion data is present, use ``AGFMaterial`` (works at any
    wavelength).  Otherwise, fall back to ``AbbeMaterial`` (visible-only).
    """
    if elem.formula_id is not None and elem.cd:
        return AGFMaterial(elem.formula_id, elem.cd)
    return AbbeMaterial(elem.nd, elem.vd)


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
    rx = thicken(candidate, inputs)
    if rx is None:
        return None

    if stage == "B":
        op_a = _build_from_prescription(rx)
        if op_a is None:
            return None
        return optimize_optic(op_a, rx, inputs)

    return _build_from_prescription(rx)


def build_optic_from_prescription(
    rx: ThickPrescription,
) -> Optional[optic.Optic]:
    """Build an optiland ``Optic`` directly from a ``ThickPrescription``.

    Useful when the caller has already computed the prescription via
    ``thickening.thicken()`` and wants to avoid recomputing it.
    """
    return _build_from_prescription(rx)


def rx_from_optic(
    op: optic.Optic,
    rx_template: ThickPrescription,
) -> ThickPrescription:
    """Read surface parameters back from an optimised optic.

    Creates a new ``ThickPrescription`` by reading the current radii and
    thicknesses from *op*'s ``surface_group``.  Glass properties (nd, vd,
    formula_id, cd) are preserved from *rx_template*.  Edge thickness is
    recomputed from sag geometry.

    Used after Stage-B optimisation to ensure ``PipelineResult.rx`` reflects
    the optimised geometry rather than the original Stage-A prescription.
    """

    def _get_t(surf_idx: int) -> float:
        t = op.surface_group.surfaces[surf_idx].thickness
        return float(np.asarray(t).ravel()[0])

    def _get_R(surf_idx: int) -> float:
        r = op.surface_group.radii[surf_idx]
        return float(np.asarray(r).ravel()[0])

    def _sag(R: float, a: float) -> float:
        """Sag at semi-aperture *a*; returns 0 for flat/degenerate surface."""
        if not math.isfinite(R) or abs(R) < 1e-12:
            return 0.0
        if abs(R) <= a + 1e-6:
            return 0.0
        return R - math.copysign(math.sqrt(R * R - a * a), R)

    a = rx_template.D / 2.0

    if rx_template.system_type == "cemented":
        R1, R2, R3 = _get_R(1), _get_R(2), _get_R(3)
        t1, t2 = _get_t(1), _get_t(2)
        bfd = _get_t(3)
        te1 = max(t1 - _sag(R1, a) + _sag(R2, a), 0.1)
        te2 = max(t2 - _sag(R2, a) + _sag(R3, a), 0.1)
        e1 = replace(rx_template.elements[0],
                     R_front=R1, R_back=R2, t_center=t1, t_edge=te1)
        e2 = replace(rx_template.elements[1],
                     R_front=R2, R_back=R3, t_center=t2, t_edge=te2)
        return replace(rx_template, elements=[e1, e2], back_focus_guess=bfd)

    else:  # spaced
        R1, R2, R3, R4 = _get_R(1), _get_R(2), _get_R(3), _get_R(4)
        t1, gap, t2 = _get_t(1), _get_t(2), _get_t(3)
        bfd = _get_t(4)
        te1 = max(t1 - _sag(R1, a) + _sag(R2, a), 0.1)
        te2 = max(t2 - _sag(R3, a) + _sag(R4, a), 0.1)
        e1 = replace(rx_template.elements[0],
                     R_front=R1, R_back=R2, t_center=t1, t_edge=te1)
        e2 = replace(rx_template.elements[1],
                     R_front=R3, R_back=R4, t_center=t2, t_edge=te2)
        return replace(rx_template, elements=[e1, e2],
                       air_gap=gap, back_focus_guess=bfd)


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

    mat1 = _make_material(e1)
    mat2 = _make_material(e2)

    op = optic.Optic()

    # Surface 0 – Object (at infinity)
    op.add_surface(index=0, radius=float("inf"), thickness=float("inf"))

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

    mat1 = _make_material(e1)
    mat2 = _make_material(e2)

    assert rx.air_gap is not None, "Spaced doublet must have an air_gap"

    op = optic.Optic()

    # Surface 0 – Object (at infinity)
    op.add_surface(index=0, radius=float("inf"), thickness=float("inf"))

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

    # Field  –  on-axis + off-axis for proper Seidel aberration evaluation.
    # Seidel CC, AC, PC, DC, TchC all depend on chief-ray height, which is
    # zero for on-axis only, making those aberrations identically 0.
    op.set_field_type(field_type="angle")
    op.add_field(y=0.0)
    op.add_field(y=1.0)  # 1° off-axis for aberration evaluation

    # Wavelengths  –  three spectral lines
    #   wavelengths = (lam1_short, lam0_primary, lam2_long)
    lam1, lam0, lam2 = rx.wavelengths
    op.add_wavelength(value=lam1)  # short wavelength
    op.add_wavelength(value=lam0, is_primary=True)  # primary / design
    op.add_wavelength(value=lam2)  # long wavelength

    # Move image plane to paraxial focus
    op.image_solve()
