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
_CORRECTION_ITER = 2  # thickness → radii → re-thickness (2 rounds)
_SCALE_EPS = 1e-12  # treat B_prod ≈ 0 as no correction needed

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
# Radius correction  –  keep Φ_thick = Φ_thin (Strategy A)
# ---------------------------------------------------------------------------


def correct_radii_for_thickness(
    R_front: float,
    R_back: float,
    t: float,
    n: float,
) -> tuple[float, float] | None:
    """Scale both curvatures so thick-lens power equals thin-lens power.

    Keeps the *shape factor* (bending) invariant.  Both curvatures are
    multiplied by the **same** scale *s* so that:

        Φ_thick(s·c₁, s·c₂, t, n)  =  Φ_thin(c₁, c₂, n)

    Returns ``(R_front', R_back')`` or ``None`` if correction fails
    (e.g. discriminant < 0, s non-positive, or degenerate curvatures).
    """

    def _curv(R: float) -> float:
        if not math.isfinite(R) or abs(R) < 1e-12:
            return 0.0
        return 1.0 / R

    c1 = _curv(R_front)
    c2 = _curv(R_back)

    # If either surface is flat, Φ_thick = (n-1)·s·(c1-c2) → s = 1.
    # Similarly if c1·c2 ≈ 0 the thickness term vanishes.
    A = c1 - c2  # proportional to Φ_thin
    B_prod = (n - 1.0) * t * c1 * c2 / n  # thickness coupling term

    if abs(B_prod) < _SCALE_EPS:
        # No correction needed (plano or negligible coupling)
        return (R_front, R_back)

    if abs(A) < 1e-15:
        # Zero power lens – nothing to preserve
        return (R_front, R_back)

    # Solve:  B_prod · s² + A · s − A = 0
    discriminant = A * A + 4.0 * A * B_prod
    if discriminant < 0:
        logger.debug(
            "correct_radii: negative discriminant=%.6e for "
            "R_f=%.4f R_b=%.4f t=%.4f n=%.4f",
            discriminant,
            R_front,
            R_back,
            t,
            n,
        )
        return None

    sqrt_disc = math.sqrt(discriminant)

    # Two roots: pick the positive one closest to 1.0
    s1 = (-A + sqrt_disc) / (2.0 * B_prod)
    s2 = (-A - sqrt_disc) / (2.0 * B_prod)

    candidates = [x for x in (s1, s2) if x > 0 and math.isfinite(x)]
    if not candidates:
        logger.debug(
            "correct_radii: no positive root (s1=%.6e, s2=%.6e) for "
            "R_f=%.4f R_b=%.4f t=%.4f n=%.4f",
            s1,
            s2,
            R_front,
            R_back,
            t,
            n,
        )
        return None

    s = min(candidates, key=lambda x: abs(x - 1.0))

    R_front_new = R_front / s if abs(c1) > 1e-12 else R_front
    R_back_new = R_back / s if abs(c2) > 1e-12 else R_back

    return (R_front_new, R_back_new)


def _reconcile_cemented_radius(
    R_from_elem1: float,
    R_from_elem2: float,
) -> float:
    """Weighted average of two corrections of the shared cemented surface."""
    return (R_from_elem1 + R_from_elem2) / 2.0


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


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


# ---------------------------------------------------------------------------
# Prescription assembly  –  Candidate + Inputs → ThickPrescription
# ---------------------------------------------------------------------------


def _round6(r: float) -> float:
    """Round radius to 6 decimal places (optiland precision workaround)."""
    if not math.isfinite(r):
        return r
    return round(r, 6)


def thicken_cemented(
    cand: Candidate,
    inp: Inputs,
) -> ThickPrescription | None:
    """Compute a thick-lens prescription for a cemented doublet.

    Iteratively determines thicknesses and corrects radii to preserve
    each element's thin-lens power.

    Returns ``None`` if any element has invalid geometry.
    """
    R1, R2, R3 = cand.radii[0], cand.radii[1], cand.radii[2]

    # --- iterative thickness ↔ radius correction ---
    for _iter in range(_CORRECTION_ITER):
        result1 = element_thickness(R1, R2, inp.D, cand.n1)
        if result1 is None:
            return None
        t1, _ = result1

        result2 = element_thickness(R2, R3, inp.D, cand.n2)
        if result2 is None:
            return None
        t2, _ = result2

        # correct radii to keep Φ_thick = Φ_thin
        corr1 = correct_radii_for_thickness(R1, R2, t1, cand.n1)
        if corr1 is None:
            return None
        R1_new, R2_from_e1 = corr1

        corr2 = correct_radii_for_thickness(R2, R3, t2, cand.n2)
        if corr2 is None:
            return None
        R2_from_e2, R3_new = corr2

        # reconcile the shared cemented surface
        R2_new = _reconcile_cemented_radius(R2_from_e1, R2_from_e2)

        R1 = _round6(R1_new)
        R2 = _round6(R2_new)
        R3 = _round6(R3_new)

    # final thickness with corrected radii
    result1 = element_thickness(R1, R2, inp.D, cand.n1)
    if result1 is None:
        return None
    t1, te1 = result1

    result2 = element_thickness(R2, R3, inp.D, cand.n2)
    if result2 is None:
        return None
    t2, te2 = result2

    elem1 = ElementRx(
        R_front=R1,
        R_back=R2,
        t_center=t1,
        t_edge=te1,
        nd=cand.n1,
        vd=cand.nu1,
    )
    elem2 = ElementRx(
        R_front=R2,
        R_back=R3,
        t_center=t2,
        t_edge=te2,
        nd=cand.n2,
        vd=cand.nu2,
    )

    return ThickPrescription(
        system_type="cemented",
        elements=[elem1, elem2],
        air_gap=None,
        back_focus_guess=_initial_back_focus(inp.fprime),
        D=inp.D,
        wavelengths=(inp.lam1, inp.lam0, inp.lam2),
    )


def thicken_spaced(
    cand: Candidate,
    inp: Inputs,
) -> ThickPrescription | None:
    """Compute a thick-lens prescription for an air-spaced doublet.

    Iteratively determines thicknesses and corrects radii to preserve
    each element's thin-lens power.

    Returns ``None`` if any element has invalid geometry.
    """
    R1, R2 = cand.radii[0], cand.radii[1]
    R3, R4 = cand.radii[2], cand.radii[3]

    # --- iterative thickness ↔ radius correction ---
    for _iter in range(_CORRECTION_ITER):
        result1 = element_thickness(R1, R2, inp.D, cand.n1)
        if result1 is None:
            return None
        t1, _ = result1

        result2 = element_thickness(R3, R4, inp.D, cand.n2)
        if result2 is None:
            return None
        t2, _ = result2

        # correct radii (independent elements – no shared surface)
        corr1 = correct_radii_for_thickness(R1, R2, t1, cand.n1)
        if corr1 is None:
            return None
        R1, R2 = _round6(corr1[0]), _round6(corr1[1])

        corr2 = correct_radii_for_thickness(R3, R4, t2, cand.n2)
        if corr2 is None:
            return None
        R3, R4 = _round6(corr2[0]), _round6(corr2[1])

    # final thickness with corrected radii
    result1 = element_thickness(R1, R2, inp.D, cand.n1)
    if result1 is None:
        return None
    t1, te1 = result1

    result2 = element_thickness(R3, R4, inp.D, cand.n2)
    if result2 is None:
        return None
    t2, te2 = result2

    elem1 = ElementRx(
        R_front=R1,
        R_back=R2,
        t_center=t1,
        t_edge=te1,
        nd=cand.n1,
        vd=cand.nu1,
    )
    elem2 = ElementRx(
        R_front=R3,
        R_back=R4,
        t_center=t2,
        t_edge=te2,
        nd=cand.n2,
        vd=cand.nu2,
    )

    return ThickPrescription(
        system_type="spaced",
        elements=[elem1, elem2],
        air_gap=_placeholder_air_gap(inp.D),
        back_focus_guess=_initial_back_focus(inp.fprime),
        D=inp.D,
        wavelengths=(inp.lam1, inp.lam0, inp.lam2),
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
