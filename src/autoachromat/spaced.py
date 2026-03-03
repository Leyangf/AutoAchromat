from __future__ import annotations

import logging
import math
from typing import List, Tuple

import numpy as np

from .glass_reader import Glass
from .models import Inputs, Candidate
from .optics import (
    achromat_power,
    check_min_radius,
    prepare_glass_data,
)
from .thermal import compute_thermal_metrics

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Sweep parameters
# ---------------------------------------------------------------------------

# Number of Q1 grid points scanned along the coma-zero line.
# 61 points → step ≈ 0.1 over [-3, 3].
_Q1_SWEEP_N = 61

# Treat two Q1 values as duplicates when they are closer than this threshold.
# Used to avoid adding a sweep solution that is nearly identical to an
# algebraic solution.
_Q1_DEDUP_TOL = 0.05


# ---------------------------------------------------------------------------
# Inner-surface sagitta (local copy; avoids import from thickening)
# ---------------------------------------------------------------------------


def _sag(R: float, a: float) -> float:
    """Signed sagitta of a spherical surface at semi-aperture *a*."""
    if not math.isfinite(R) or abs(R) < 1e-12:
        return 0.0
    inner = R * R - a * a
    if inner < 0.0:
        raise ValueError(f"|R|={abs(R):.4f} < semi-aperture a={a:.4f}")
    return R - math.copysign(math.sqrt(inner), R)


# ---------------------------------------------------------------------------
# Seidel coefficient helpers
# ---------------------------------------------------------------------------


def _coeffs(n1: float, n2: float, phi1: float, phi2: float) -> dict:
    A1 = (phi1**3) * (1.0 + 2.0 / n1)
    B1 = (phi1**3) * (3.0 / (n1 - 1.0))

    A2 = (phi2**3) * (1.0 + 2.0 / n2)
    B2 = (phi2**3) * (3.0 / (n2 - 1.0)) - 4.0 * phi1 * (phi2**2) * (1.0 + 1.0 / n2)

    C = (
        (phi1**3) * (n1 / ((n1 - 1.0) ** 2))
        + (phi2**3) * (n2 / ((n2 - 1.0) ** 2))
        - (phi1 * (phi2**2)) * ((4.0 / (n2 - 1.0)) + 1.0)
        + ((phi1**2) * phi2) * (3.0 + (2.0 / n2))
    )

    K1 = (phi1**2) * (1.0 + 1.0 / n1)
    K2 = (phi2**2) * (1.0 + 1.0 / n2)

    L = (
        (phi1**2) * (1.0 / (n1 - 1.0))
        + (phi2**2) * (1.0 / (n2 - 1.0))
        - phi1 * phi2 * (2.0 + 1.0 / n2)
    )

    return dict(A1=A1, B1=B1, A2=A2, B2=B2, C=C, K1=K1, K2=K2, L=L)


# ---------------------------------------------------------------------------
# Q-pair solvers
# ---------------------------------------------------------------------------


def _solve_Q_pairs(inputs: Inputs, c: dict) -> List[Tuple[float, float]]:
    """Return up to 2 (Q1, Q2) pairs that minimise spherical aberration on
    the coma-zero line  K1·Q1 + K2·Q2 + (L − W0) = 0.

    These are the algebraic solutions of the Seidel thin-lens equations for a
    spaced doublet.  For fast systems (f/4 – f/6) they often fall in the
    geometrically forbidden zone (inner-surface overlap); in that case the
    coma-line sweep (``_sweep_Q1_coma_line``) finds the nearest valid design.
    """
    A1, B1, A2, B2, C, K1, K2, L = (
        c["A1"],
        c["B1"],
        c["A2"],
        c["B2"],
        c["C"],
        c["K1"],
        c["K2"],
        c["L"],
    )

    if abs(K2) < inputs.eps:
        return []

    # Q2 = a·Q1 + b  (coma-zero line)
    a = -K1 / K2
    b = (inputs.W0 - L) / K2

    # Substitute → quadratic in Q1
    q2 = A1 + A2 * (a * a)
    q1 = B1 + A2 * (2.0 * a * b) + B2 * a
    q0 = A2 * (b * b) + B2 * b + (C - inputs.P0)

    roots = np.roots(np.array([q2, q1, q0], dtype=float))

    pairs: list[Tuple[float, float]] = []
    for r in roots:
        if abs(r.imag) > inputs.root_imag_tol:
            continue
        Q1 = float(r.real)
        Q2 = float(a * Q1 + b)
        pairs.append((Q1, Q2))

    # Return both solutions; the outer loop applies geometric filters and
    # PE-based ranking.  Sorting by minimum norm here would bias all designs
    # toward minimum bending, making all results look structurally identical.
    return pairs


def _seidel_SA_residual(c: dict, inputs: Inputs, Q1: float) -> float:
    """Seidel spherical aberration residual on the coma-zero line.

    Returns  S_I(Q1) − P0  where S_I is the Seidel third-order SA coefficient
    of the spaced doublet evaluated at (Q1, Q2=a·Q1+b) on the coma-zero line.

    The algebraic optimum has residual = 0 by construction (it is the root of
    the quadratic).  Sweep designs have |residual| > 0, proportional to how
    far they are from the aberration-minimum bending.

    This is the correct ranking metric for spaced doublets – it replaces the
    ad-hoc per-surface P-sum used previously.
    """
    A1, B1, A2, B2, C, K1, K2, L = (
        c["A1"],
        c["B1"],
        c["A2"],
        c["B2"],
        c["C"],
        c["K1"],
        c["K2"],
        c["L"],
    )
    if abs(K2) < inputs.eps:
        return math.inf
    a = -K1 / K2
    b = (inputs.W0 - L) / K2
    q2 = A1 + A2 * (a * a)
    q1 = B1 + A2 * (2.0 * a * b) + B2 * a
    q0 = A2 * (b * b) + B2 * b + (C - inputs.P0)
    return q2 * Q1**2 + q1 * Q1 + q0


def _sweep_Q1_coma_line(
    inputs: Inputs,
    c: dict,
    n1: float,
    n2: float,
    phi1: float,
    phi2: float,
) -> List[Tuple[float, float]]:
    """Sweep Q1 ∈ [−3, 3] along the coma-zero line and return the single
    (Q1, Q2) pair that is geometrically valid and has minimum |S_I residual|.

    Rationale
    ---------
    The algebraic solver gives the 2 Q1 values that minimise spherical
    aberration.  For fast systems (f/4–f/6 with a small air gap) those Q1
    values fall in the geometrically forbidden zone: the inner surfaces of the
    two elements overlap at the aperture edge.  This sweep finds the
    feasible design closest to the aberration minimum for each glass pair.

    Because the geometrically feasible Q1 range is glass-pair–dependent
    (it depends on n1, n2, φ1, φ2 and the air gap), different glass pairs
    yield different optimal Q1 values → structurally diverse outputs.
    """
    K1 = c["K1"]
    K2 = c["K2"]
    L = c["L"]

    if abs(K2) < inputs.eps:
        return []

    a_coma = -K1 / K2
    b_coma = (inputs.W0 - L) / K2
    semi_ap = inputs.D / 2.0

    best_SI: float = math.inf  # minimum |S_I residual| seen so far
    best_pair: Tuple[float, float] | None = None
    low_pair: Tuple[float, float] | None = None  # most negative valid Q1
    high_pair: Tuple[float, float] | None = None  # most positive valid Q1

    for Q1 in np.linspace(-3.0, 3.0, _Q1_SWEEP_N):
        Q2 = a_coma * Q1 + b_coma
        try:
            R1, R2, R3, R4 = _radii(inputs, n1, n2, phi1, phi2, Q1, Q2)
        except (ZeroDivisionError, ValueError):
            continue

        if not check_min_radius((R1, R2, R3, R4), inputs.D):
            continue

        try:
            sag_R2 = _sag(R2, semi_ap)
            sag_R3 = _sag(R3, semi_ap)
        except ValueError:
            continue

        if inputs.air_gap + sag_R3 - sag_R2 < 0.0:
            continue

        SI_resid = abs(_seidel_SA_residual(c, inputs, Q1))

        if SI_resid < best_SI:
            best_SI = SI_resid
            best_pair = (Q1, Q2)

        # linspace is sorted ascending, so first valid = most negative (low),
        # each update overwrites high with the most positive seen so far.
        if low_pair is None:
            low_pair = (Q1, Q2)
        high_pair = (Q1, Q2)

    # Return up to 3 structurally distinct designs: best + two boundary extremes.
    # Collapse near-duplicates within this set using _Q1_DEDUP_TOL.
    candidates: list[Tuple[float, float]] = []
    for pair in (best_pair, low_pair, high_pair):
        if pair is None:
            continue
        if not any(abs(pair[0] - ex[0]) < _Q1_DEDUP_TOL for ex in candidates):
            candidates.append(pair)

    return candidates


# ---------------------------------------------------------------------------
# Radii from (Q1, Q2)
# ---------------------------------------------------------------------------


def _radii(
    inputs: Inputs, n1: float, n2: float, phi1: float, phi2: float, Q1: float, Q2: float
) -> Tuple[float, float, float, float]:
    f = inputs.fprime
    if abs(phi1) < inputs.eps or abs(phi2) < inputs.eps:
        raise ZeroDivisionError("phi1 or phi2 near zero in spaced _radii")
    f1 = f / phi1
    f2 = f / phi2

    denom_R1 = n1 / (n1 - 1.0) + Q1
    denom_R2 = 1.0 + Q1
    denom_R3 = n2 / (n2 - 1.0) + Q2
    denom_R4 = 1.0 + Q2

    if (
        abs(denom_R1) < 1e-12
        or abs(denom_R2) < 1e-12
        or abs(denom_R3) < 1e-12
        or abs(denom_R4) < 1e-12
    ):
        raise ZeroDivisionError("Degenerate Q value: denominator near zero")

    R1 = f1 / denom_R1
    R2 = f1 / denom_R2
    R3 = f2 / denom_R3
    R4 = f2 / denom_R4
    return R1, R2, R3, R4


# ---------------------------------------------------------------------------
# PE (primary aberration) evaluator
# ---------------------------------------------------------------------------


def _Ps_and_PE(
    n1: float, n2: float, phi1: float, phi2: float, Q1: float, Q2: float
) -> Tuple[List[float], float]:
    u1 = 0.0
    u2 = Q1 * (1.0 - 1.0 / n1)
    u3 = phi1
    u4 = phi1 + Q2 * (1.0 - 1.0 / n2)
    u5 = 1.0

    P1 = ((u2 - u1) / (1.0 / n1 - 1.0)) ** 2 * (u2 / n1 - u1)
    P2 = ((u3 - u2) / (1.0 - 1.0 / n1)) ** 2 * (u3 - u2 / n1)
    P3 = ((u4 - u3) / (1.0 / n2 - 1.0)) ** 2 * (u4 / n2 - u3)
    P4 = ((u5 - u4) / (1.0 - 1.0 / n2)) ** 2 * (u5 - u4 / n2)

    Ps = [P1, P2, P3, P4]
    PE = sum(abs(x) for x in Ps) / 4.0
    return Ps, PE


# ---------------------------------------------------------------------------
# Main synthesis entry point
# ---------------------------------------------------------------------------


def run_spaced(inputs: Inputs, glasses: list[Glass]) -> list[Candidate]:
    """Enumerate glass pairs and return achromatic spaced-doublet candidates.

    Strategy
    --------
    For each ordered glass pair (g1, g2):

    1. **Algebraic solutions** – solve the thin-lens Seidel equations for the
       two Q1 values that minimise spherical aberration subject to zero coma
       (``_solve_Q_pairs``).  These are optimal but may fail the geometric
       constraints (min-radius, inner-surface non-overlap) for fast systems.

    2. **Coma-line sweep** – scan Q1 ∈ [−3, 3] along the coma-zero line and
       find the geometrically feasible design with minimum PE
       (``_sweep_Q1_coma_line``).  This finds a valid design even when both
       algebraic solutions are in the forbidden zone, and produces structurally
       diverse results because the feasible Q1 range is glass-pair–specific.

    The two sources are deduplicated (sweep result dropped when it is within
    ``_Q1_DEDUP_TOL`` of an algebraic solution), then evaluated together.
    All valid candidates are ranked by PE.
    """
    gdata = prepare_glass_data(glasses, inputs.lam0, inputs.lam1, inputs.lam2)
    semi_ap = inputs.D / 2.0

    out: list[Candidate] = []

    for g1, n1, nu1 in gdata:
        for g2, n2, nu2 in gdata:
            if (g1.name == g2.name) and (g1.catalog == g2.catalog):
                continue
            if abs(nu1 - nu2) < inputs.min_delta_nu:
                continue

            try:
                phi1, phi2 = achromat_power(nu1, nu2, inputs.C0)
                c = _coeffs(n1, n2, phi1, phi2)
            except (ValueError, ZeroDivisionError):
                continue

            # --- Collect (Q1, Q2) candidates from both sources ---
            try:
                pairs_alg = _solve_Q_pairs(inputs, c)
            except Exception:
                pairs_alg = []

            try:
                pairs_sweep = _sweep_Q1_coma_line(inputs, c, n1, n2, phi1, phi2)
            except Exception:
                pairs_sweep = []

            # Deduplicate: drop sweep result if within _Q1_DEDUP_TOL of any
            # algebraic Q1 (the algebraic solution is already the min-SA point).
            alg_Q1s = [p[0] for p in pairs_alg]
            deduped_sweep = [
                (Q1, Q2)
                for (Q1, Q2) in pairs_sweep
                if all(abs(Q1 - aq1) > _Q1_DEDUP_TOL for aq1 in alg_Q1s)
            ]

            all_pairs = pairs_alg + deduped_sweep

            for Q1, Q2 in all_pairs:
                try:
                    R1, R2, R3, R4 = _radii(inputs, n1, n2, phi1, phi2, Q1, Q2)
                    if not check_min_radius((R1, R2, R3, R4), inputs.D):
                        continue

                    # Non-overlap: air_gap + sag(R3, a) − sag(R2, a) ≥ 0
                    if inputs.air_gap + _sag(R3, semi_ap) - _sag(R2, semi_ap) < 0.0:
                        continue

                    Ps, _ = _Ps_and_PE(n1, n2, phi1, phi2, Q1, Q2)
                    PE = abs(_seidel_SA_residual(c, inputs, Q1))
                    if PE > inputs.max_PE:
                        continue

                    thermal = compute_thermal_metrics(
                        g1,
                        g2,
                        n1,
                        n2,
                        phi1,
                        phi2,
                        wavelength_um=inputs.lam0,
                    )

                    out.append(
                        Candidate(
                            system_type="spaced",
                            glass1=g1.name,
                            catalog1=g1.catalog,
                            glass2=g2.name,
                            catalog2=g2.catalog,
                            n1=n1,
                            n2=n2,
                            nu1=nu1,
                            nu2=nu2,
                            phi1=phi1,
                            phi2=phi2,
                            Q1=Q1,
                            Q2=Q2,
                            P_surfaces=Ps,
                            radii=[R1, R2, R3, R4],
                            PE=PE,
                            cost1=g1.relative_cost,
                            cost2=g2.relative_cost,
                            formula_id1=g1.formula_id,
                            cd1=list(g1.cd),
                            formula_id2=g2.formula_id,
                            cd2=list(g2.cd),
                            notes={"C_const_used": c["C"]},
                            thermal=thermal,
                        )
                    )
                except (ValueError, ZeroDivisionError) as exc:
                    logger.debug(
                        "Skipping candidate %s+%s Q1=%.4g Q2=%.4g – %s",
                        g1.name,
                        g2.name,
                        Q1,
                        Q2,
                        exc,
                    )
                    continue

    # Primary sort: |S_I residual| (algebraic solutions = 0, sweep solutions > 0).
    # Secondary sort: sum(|P_surfaces|) — measures per-surface aberration load;
    #   lower load means less aggressive correction, better tolerance sensitivity.
    #   This acts as a tiebreaker among all PE=0 algebraic designs.
    def _sort_key(cand: Candidate) -> tuple:
        primary = cand.PE if cand.PE is not None else float("inf")
        secondary = (
            sum(abs(x) for x in cand.P_surfaces) if cand.P_surfaces else float("inf")
        )
        return (primary, secondary)

    out.sort(key=_sort_key)
    return out
