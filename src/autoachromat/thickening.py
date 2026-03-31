"""
thickening.py – Thin-lens → thick-lens prescription (pure math, no optical library).

Converts a thin-lens ``Candidate`` into a ``ThickPrescription`` containing
corrected radii, centre/edge thicknesses, and outside diameters, following
the textbook procedure (Table 10-2 / 10-3 + sag geometry + power-preserving
radius correction).

Zero external dependencies beyond the Python standard library.
"""

from __future__ import annotations

import math
import logging
from collections.abc import Sequence

from .models import Candidate, Inputs, ElementRx, ThickPrescription

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_MIN_CENTER_ABSOLUTE = 0.5  # mm – always enforced (any lens)
_SAG_MARGIN = 1e-6  # mm – reject if |R| < a + margin
_CORRECTION_ITER = 4  # thickness → radii → re-thickness iterations

# ---------------------------------------------------------------------------
# Table 10-2: Δ(D) for outside diameter (压圈法固定 / retaining ring)
# ---------------------------------------------------------------------------

_DELTA_RETAINING_RING: list[tuple[float, float | None]] = [
    # (D_upper_bound, Delta)  — first match wins
    (6.0, None),  # "—" in table → invalid for retaining ring
    (10.0, 1.0),
    (18.0, 1.5),
    (30.0, 2.0),
    (50.0, 2.5),
    (80.0, 3.0),
    (120.0, 3.5),
    (math.inf, 4.5),
]

# ---------------------------------------------------------------------------
# Table 10-3: min edge thickness (positive) / min centre thickness (negative)
# ---------------------------------------------------------------------------

_T_EDGE_MIN_TABLE: list[tuple[float, float]] = [
    (6.0, 0.4),
    (10.0, 0.6),
    (18.0, 0.8),
    (30.0, 1.2),
    (50.0, 1.8),
    (80.0, 2.4),
    (120.0, 3.0),
    (math.inf, 4.0),
]

_T_CENTER_MIN_TABLE: list[tuple[float, float]] = [
    (6.0, 0.6),
    (10.0, 0.8),
    (18.0, 1.0),
    (30.0, 1.5),
    (50.0, 2.2),
    (80.0, 3.5),
    (120.0, 5.0),
    (math.inf, 8.0),
]

# ---------------------------------------------------------------------------
# Piecewise lookup
# ---------------------------------------------------------------------------


def _piecewise_lookup(
    D: float,
    table: Sequence[tuple[float, float | None]],
    *,
    D_min: float = 3.0,
) -> float | None:
    """Return the first table value whose upper bound >= *D*.

    Returns ``None`` when *D* < *D_min* or the matched entry is ``None``.
    """
    if D < D_min:
        return None
    for D_upper, value in table:
        if D <= D_upper:
            return value
    return None


# ---------------------------------------------------------------------------
# Public lookup helpers
# ---------------------------------------------------------------------------


def lookup_delta(D: float, mount: str = "retaining_ring") -> float | None:
    """Outside-diameter increment Δ(D) [mm]  (Table 10-2).

    Returns ``None`` if *D* is out of range for *mount*.
    Default mount = 压圈法固定 (retaining ring).
    """
    if mount != "retaining_ring":
        raise ValueError(f"Unsupported mount method: {mount!r}")
    return _piecewise_lookup(D, _DELTA_RETAINING_RING, D_min=0.0)


def lookup_t_edge_min(D: float) -> float:
    """Min edge thickness for a **positive** lens (Table 10-3)."""
    val = _piecewise_lookup(D, _T_EDGE_MIN_TABLE)
    return val if val is not None else _T_EDGE_MIN_TABLE[0][1]  # 0.4


def lookup_t_center_min(D: float) -> float:
    """Min centre thickness for a **negative** lens (Table 10-3)."""
    val = _piecewise_lookup(D, _T_CENTER_MIN_TABLE)
    return val if val is not None else _T_CENTER_MIN_TABLE[0][1]  # 0.6


def outside_diameter(D: float, mount: str = "retaining_ring") -> float | None:
    """Outside diameter φ = D + Δ  (Table 10-2).

    Returns ``None`` when *D* is out of range for *mount*.
    """
    delta = lookup_delta(D, mount)
    if delta is None:
        return None
    return D + delta


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------


def _sag(R: float, a: float) -> float:
    """Signed sagitta at semi-aperture *a*.

    Formula:  ``sag = R − sign(R) · √(R² − a²)``

    Raises ``ValueError`` when the surface cannot accommodate the clear
    aperture (|R| < a + margin) or when *R* is degenerate.
    """
    if not math.isfinite(R) or abs(R) < 1e-12:
        return 0.0  # flat surface → sag = 0
    if abs(R) < a + _SAG_MARGIN:
        raise ValueError(
            f"|R|={abs(R):.6f} mm < semi-aperture+margin="
            f"{a + _SAG_MARGIN:.6f} mm — surface cannot exist"
        )
    return R - math.copysign(math.sqrt(R * R - a * a), R)


def _lens_power(R_front: float, R_back: float, n: float) -> float:
    """Thin-lens paraxial power Φ = (n−1)(1/R₁ − 1/R₂)."""

    def _curv(R: float) -> float:
        if not math.isfinite(R) or abs(R) < 1e-12:
            return 0.0
        return 1.0 / R

    return (n - 1.0) * (_curv(R_front) - _curv(R_back))


def element_thickness(
    R_front: float,
    R_back: float,
    D: float,
    n: float,
) -> tuple[float, float] | None:
    """Compute *(t_center, t_edge)* for a single lens element.

    Uses the textbook procedure (Table 10-3 + sag geometry):

    * **Positive lens** (Φ ≥ 0) → enforce min *edge* thickness.
    * **Negative lens** (Φ < 0) → enforce min *centre* thickness.
    * Always enforce absolute min centre thickness (0.5 mm).

    Returns ``None`` when geometry is invalid (|R| too small, non-finite, etc.).
    """
    a = D / 2.0

    # -- sag validity --
    try:
        sag_f = _sag(R_front, a)
        sag_b = _sag(R_back, a)
    except ValueError as exc:
        logger.debug("element_thickness rejected: %s", exc)
        return None

    # -- lens sign --
    phi = _lens_power(R_front, R_back, n)

    # -- lookup limits --
    te_min = lookup_t_edge_min(D)
    tc_min_neg = lookup_t_center_min(D)

    # -- centre thickness from edge-thickness constraint --
    #    t_edge = t_center − sag_f + sag_b  ≥  te_min
    #  → t_center ≥ te_min + sag_f − sag_b
    t_from_edge = te_min + sag_f - sag_b

    if phi >= 0:
        # Positive / plano lens: edge constraint dominates
        t_center = max(t_from_edge, _MIN_CENTER_ABSOLUTE)
    else:
        # Negative lens: centre constraint dominates
        t_center = max(t_from_edge, tc_min_neg, _MIN_CENTER_ABSOLUTE)

    t_edge = t_center - sag_f + sag_b

    # -- manufacturability guard --
    if not math.isfinite(t_center) or t_center <= 0:
        logger.debug("element_thickness: invalid t_center=%.4f", t_center)
        return None
    if not math.isfinite(t_edge) or t_edge <= 0:
        logger.debug("element_thickness: invalid t_edge=%.4f", t_edge)
        return None

    return (t_center, t_edge)


# ---------------------------------------------------------------------------
# System EFL via ABCD (paraxial ray-transfer) matrix
# ---------------------------------------------------------------------------

_EFL_TOL = 1e-9  # relative convergence tolerance for EFL correction


def _mat2_mul(
    a: tuple[tuple[float, float], tuple[float, float]],
    b: tuple[tuple[float, float], tuple[float, float]],
) -> tuple[tuple[float, float], tuple[float, float]]:
    """Multiply two 2\u00d72 matrices ``a @ b``."""
    return (
        (
            a[0][0] * b[0][0] + a[0][1] * b[1][0],
            a[0][0] * b[0][1] + a[0][1] * b[1][1],
        ),
        (
            a[1][0] * b[0][0] + a[1][1] * b[1][0],
            a[1][0] * b[0][1] + a[1][1] * b[1][1],
        ),
    )


def _refraction_matrix(
    R: float, n_before: float, n_after: float
) -> tuple[tuple[float, float], tuple[float, float]]:
    """Refraction matrix at a spherical surface of radius *R*."""
    if not math.isfinite(R) or abs(R) < 1e-12:
        phi = 0.0
    else:
        phi = (n_after - n_before) / R
    return ((1.0, 0.0), (-phi, 1.0))


def _transfer_matrix(
    t: float, n: float
) -> tuple[tuple[float, float], tuple[float, float]]:
    """Transfer (propagation) matrix for thickness *t* in medium index *n*."""
    return ((1.0, t / n), (0.0, 1.0))


def _system_efl_cemented(
    R1: float,
    R2: float,
    R3: float,
    t1: float,
    t2: float,
    n1: float,
    n2: float,
) -> float:
    """Paraxial system EFL for a cemented doublet (ABCD matrix method).

    Surface model::

        air(1) \u2192 R1 \u2192 glass(n1), thickness t1
               \u2192 R2 \u2192 glass(n2), thickness t2
               \u2192 R3 \u2192 air(1)

    Returns ``f' = -1/C`` where the system matrix is ``((A,B),(C,D))``.
    """
    M = _refraction_matrix(R1, 1.0, n1)
    M = _mat2_mul(_transfer_matrix(t1, n1), M)
    M = _mat2_mul(_refraction_matrix(R2, n1, n2), M)
    M = _mat2_mul(_transfer_matrix(t2, n2), M)
    M = _mat2_mul(_refraction_matrix(R3, n2, 1.0), M)
    C = M[1][0]
    if abs(C) < 1e-15:
        return float("inf")
    return -1.0 / C


def _system_efl_spaced(
    R1: float,
    R2: float,
    R3: float,
    R4: float,
    t1: float,
    t2: float,
    air_gap: float,
    n1: float,
    n2: float,
) -> float:
    """Paraxial system EFL for an air-spaced doublet (ABCD matrix method).

    Surface model::

        air(1) \u2192 R1 \u2192 glass(n1), thickness t1
               \u2192 R2 \u2192 air(1),    thickness air_gap
               \u2192 R3 \u2192 glass(n2), thickness t2
               \u2192 R4 \u2192 air(1)
    """
    M = _refraction_matrix(R1, 1.0, n1)
    M = _mat2_mul(_transfer_matrix(t1, n1), M)
    M = _mat2_mul(_refraction_matrix(R2, n1, 1.0), M)
    M = _mat2_mul(_transfer_matrix(air_gap, 1.0), M)
    M = _mat2_mul(_refraction_matrix(R3, 1.0, n2), M)
    M = _mat2_mul(_transfer_matrix(t2, n2), M)
    M = _mat2_mul(_refraction_matrix(R4, n2, 1.0), M)
    C = M[1][0]
    if abs(C) < 1e-15:
        return float("inf")
    return -1.0 / C


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _initial_back_focus(fprime: float) -> float:
    """Initial guess for the back focal distance before ``image_solve``.

    ``image_solve()`` will correct this, but a reasonable starting value
    helps avoid occasional tracing failures during the solve.
    """
    return abs(fprime) * 0.9


# ---------------------------------------------------------------------------
# Prescription assembly  –  Candidate + Inputs → ThickPrescription
# ---------------------------------------------------------------------------


def thicken_cemented(
    cand: Candidate,
    inp: Inputs,
) -> ThickPrescription | None:
    """Compute a thick-lens prescription for a cemented doublet.

    Assigns physical thicknesses and iteratively scales the radii so that
    the thick-lens EFL matches the target ``inp.fprime``.  Scaling all three
    radii by the same factor ``k = fprime / actual_efl`` preserves the
    bending (shape factor Q) and the achromatic power split while correcting
    for the EFL shortening caused by thick-lens effects.  ``_CORRECTION_ITER``
    rounds are applied; the residual is reported in ``efl_deviation``.

    Returns ``None`` if any element has invalid geometry.
    """
    # Mutable local copies of radii (corrected in-place across iterations)
    R1, R2, R3 = cand.radii[0], cand.radii[1], cand.radii[2]

    t1 = te1 = t2 = te2 = 0.0
    actual_efl: float = math.nan

    for _iter in range(_CORRECTION_ITER + 1):
        # --- assign physical thicknesses ---
        result1 = element_thickness(R1, R2, inp.D, cand.n1)
        if result1 is None:
            return None
        t1, te1 = result1

        result2 = element_thickness(R2, R3, inp.D, cand.n2)
        if result2 is None:
            return None
        t2, te2 = result2

        # --- actual thick-lens EFL ---
        actual_efl = _system_efl_cemented(R1, R2, R3, t1, t2, cand.n1, cand.n2)

        # --- EFL correction: scale all radii uniformly ---
        if _iter < _CORRECTION_ITER:
            if not math.isfinite(actual_efl) or abs(actual_efl) < 1e-9:
                break
            k = inp.fprime / actual_efl
            if not math.isfinite(k) or k <= 0.0 or k > 5.0:
                break  # pathological case – stop correcting
            if abs(k - 1.0) < _EFL_TOL:
                break  # converged
            R1 *= k
            R2 *= k
            R3 *= k

    efl_dev = (
        (actual_efl - inp.fprime) / inp.fprime
        if math.isfinite(actual_efl) and abs(inp.fprime) > 1e-12
        else None
    )

    elem1 = ElementRx(
        R_front=R1,
        R_back=R2,
        t_center=t1,
        t_edge=te1,
        nd=cand.n1,
        vd=cand.nu1,
        formula_id=cand.formula_id1,
        cd=list(cand.cd1),
    )
    elem2 = ElementRx(
        R_front=R2,
        R_back=R3,
        t_center=t2,
        t_edge=te2,
        nd=cand.n2,
        vd=cand.nu2,
        formula_id=cand.formula_id2,
        cd=list(cand.cd2),
    )

    return ThickPrescription(
        system_type="cemented",
        elements=[elem1, elem2],
        air_gap=None,
        back_focus_guess=_initial_back_focus(inp.fprime),
        D=inp.D,
        wavelengths=(inp.lam1, inp.lam0, inp.lam2),
        half_field_angle=inp.half_field_angle,
        actual_efl=actual_efl,
        efl_deviation=efl_dev,
    )


def thicken_spaced(
    cand: Candidate,
    inp: Inputs,
) -> ThickPrescription | None:
    """Compute a thick-lens prescription for an air-spaced doublet.

    Assigns physical thicknesses and iteratively scales the radii so that
    the thick-lens EFL matches the target ``inp.fprime``.  Scaling all four
    radii by the same factor ``k = fprime / actual_efl`` preserves the
    bending (shape factors Q1, Q2) and the achromatic power split while
    correcting the EFL shortening caused by thick-lens effects (typically
    5–15 % for f/4 systems).  ``_CORRECTION_ITER`` rounds are applied; the
    residual is reported in ``efl_deviation``.

    The air gap is intentionally NOT scaled: it is a user-specified design
    constraint and the achromatic condition is only weakly sensitive to the
    gap for air_gap « EFL.

    Returns ``None`` if any element has invalid geometry.
    """
    R1, R2 = cand.radii[0], cand.radii[1]
    R3, R4 = cand.radii[2], cand.radii[3]
    air_gap = inp.air_gap
    a = inp.D / 2.0

    t1 = te1 = t2 = te2 = 0.0
    actual_efl: float = math.nan

    for _iter in range(_CORRECTION_ITER + 1):
        # --- assign physical thicknesses ---
        result1 = element_thickness(R1, R2, inp.D, cand.n1)
        if result1 is None:
            return None
        t1, te1 = result1

        result2 = element_thickness(R3, R4, inp.D, cand.n2)
        if result2 is None:
            return None
        t2, te2 = result2

        # --- Reject if inner surfaces physically overlap at the aperture edge ---
        # No-overlap condition: air_gap + sag(R3, a) − sag(R2, a) ≥ 0
        try:
            sag_back1 = _sag(R2, a)
            sag_front2 = _sag(R3, a)
        except ValueError:
            return None
        edge_clearance = air_gap + sag_front2 - sag_back1
        if edge_clearance < 0:
            logger.debug(
                "thicken_spaced rejected: inner surfaces overlap at aperture edge "
                "(air_gap=%.4f, sag_R2=%.4f, sag_R3=%.4f, clearance=%.4f mm)",
                air_gap,
                sag_back1,
                sag_front2,
                edge_clearance,
            )
            return None

        # --- actual thick-lens EFL ---
        actual_efl = _system_efl_spaced(
            R1, R2, R3, R4, t1, t2, air_gap, cand.n1, cand.n2
        )

        # --- EFL correction: scale all radii uniformly ---
        if _iter < _CORRECTION_ITER:
            if not math.isfinite(actual_efl) or abs(actual_efl) < 1e-9:
                break
            k = inp.fprime / actual_efl
            if not math.isfinite(k) or k <= 0.0 or k > 5.0:
                break  # pathological – stop
            if abs(k - 1.0) < _EFL_TOL:
                break  # converged
            R1 *= k
            R2 *= k
            R3 *= k
            R4 *= k

    efl_dev = (
        (actual_efl - inp.fprime) / inp.fprime
        if math.isfinite(actual_efl) and abs(inp.fprime) > 1e-12
        else None
    )

    elem1 = ElementRx(
        R_front=R1,
        R_back=R2,
        t_center=t1,
        t_edge=te1,
        nd=cand.n1,
        vd=cand.nu1,
        formula_id=cand.formula_id1,
        cd=list(cand.cd1),
    )
    elem2 = ElementRx(
        R_front=R3,
        R_back=R4,
        t_center=t2,
        t_edge=te2,
        nd=cand.n2,
        vd=cand.nu2,
        formula_id=cand.formula_id2,
        cd=list(cand.cd2),
    )

    return ThickPrescription(
        system_type="spaced",
        elements=[elem1, elem2],
        air_gap=air_gap,
        back_focus_guess=_initial_back_focus(inp.fprime),
        D=inp.D,
        wavelengths=(inp.lam1, inp.lam0, inp.lam2),
        half_field_angle=inp.half_field_angle,
        actual_efl=actual_efl,
        efl_deviation=efl_dev,
    )


# ---------------------------------------------------------------------------
# Public dispatch
# ---------------------------------------------------------------------------


def thicken(
    candidate: Candidate,
    inputs: Inputs,
) -> ThickPrescription | None:
    """Convert a thin-lens *candidate* to a thick-lens prescription.

    Dispatches to ``thicken_cemented`` or ``thicken_spaced`` based on
    ``candidate.system_type``.

    Returns ``None`` if the candidate cannot be thickened (invalid geometry).
    """
    if candidate.system_type == "cemented":
        return thicken_cemented(candidate, inputs)
    elif candidate.system_type == "spaced":
        return thicken_spaced(candidate, inputs)
    else:
        logger.warning("Unknown system_type=%s", candidate.system_type)
        return None
