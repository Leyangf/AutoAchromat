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
    geometrically forbidden zone (inner-surface overlap); in that case both
    roots are rejected and no numerical sweep fallback is performed.
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

    return pairs


# ---------------------------------------------------------------------------
# Radii from (Q1, Q2)
# ---------------------------------------------------------------------------


def radii_from_Q_pair(
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


def _Ps_and_PE(n1, n2, phi1, phi2, Q1, Q2):
    """Seidel surface contributions P1‒P4 and PE = Σ|Pᵢ|/4.

    Spaced doublet with 4 surfaces:
      Surface 1:  air (1)  → glass1 (n1)
      Surface 2:  glass1 (n1) → air (1)
      Surface 3:  air (1)  → glass2 (n2)
      Surface 4:  glass2 (n2) → air (1)

    Marginal-ray slope products (thin-lens, object at ∞, Φ = 1):
      u1 = 0
      u2 = (n1−1)·ρ1·φ1 / n1      (after surface 1)
      u3 = φ1                       (after surface 2)
      u4 = (φ1 + (n2−1)·ρ3·φ2)/n2  (after surface 3)
      u5 = φ1 + φ2 = 1             (after surface 4)
    """
    rho1 = n1 / (n1 - 1.0) + Q1
    rho3 = n2 / (n2 - 1.0) + Q2

    u1 = 0.0
    u2 = (n1 - 1.0) * rho1 * phi1 / n1
    u3 = phi1
    u4 = (phi1 + (n2 - 1.0) * rho3 * phi2) / n2
    u5 = phi1 + phi2  # = 1.0 under Phi normalization

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

    1. Solve the thin-lens Seidel equations for the two Q1 values that
       minimise spherical aberration on the coma-zero line
       (``_solve_Q_pairs``).  These are the exact analytical solutions.

    2. Apply geometric constraints (min-radius, inner-surface non-overlap).

    3. Rank by PE = Σ|Pᵢ|/4  (mean per-surface Seidel contribution).
       SA_residual is identically 0 for algebraic roots, so PE is the
       primary ranking metric; it also reflects higher-order aberration
       susceptibility and manufacturing difficulty.
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

            try:
                pairs = _solve_Q_pairs(inputs, c)
            except Exception:
                pairs = []

            # Track valid candidates for deferred thermal computation.
            valid_entries: list[tuple] = []  # (PE, Q1, Q2, Ps, radii)

            for Q1, Q2 in pairs:
                try:
                    R1, R2, R3, R4 = radii_from_Q_pair(inputs, n1, n2, phi1, phi2, Q1, Q2)
                    if not check_min_radius((R1, R2, R3, R4), inputs.D):
                        continue

                    # Non-overlap: air_gap + sag(R3, a) − sag(R2, a) ≥ 0
                    if inputs.air_gap + _sag(R3, semi_ap) - _sag(R2, semi_ap) < 0.0:
                        continue

                    Ps, PE = _Ps_and_PE(n1, n2, phi1, phi2, Q1, Q2)
                    if PE > inputs.max_PE:
                        continue

                    valid_entries.append((PE, Q1, Q2, Ps, [R1, R2, R3, R4]))
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

            # Build Candidate objects (including thermal) only when
            # at least one valid solution exists.
            if valid_entries:
                thermal = compute_thermal_metrics(
                    g1,
                    g2,
                    n1,
                    n2,
                    phi1,
                    phi2,
                    wavelength_um=inputs.lam0,
                )
                for PE, Q1, Q2, Ps, radii in valid_entries:
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
                            radii=radii,
                            PE=PE,
                            cost1=g1.relative_cost,
                            cost2=g2.relative_cost,
                            formula_id1=g1.formula_id,
                            cd1=list(g1.cd),
                            formula_id2=g2.formula_id,
                            cd2=list(g2.cd),
                            trans1=list(g1.transmission),
                            trans2=list(g2.transmission),
                            notes={"C_const_used": c["C"]},
                            thermal=thermal,
                        )
                    )

    # Sort by PE = Σ|Pᵢ|/4.  All candidates are algebraic solutions with
    # SA_residual ≡ 0 and coma ≡ 0, so PE (mean per-surface Seidel
    # contribution) is the sole distinguishing aberration / manufacturing
    # metric.
    out.sort(key=lambda c: c.PE if c.PE is not None else float("inf"))
    return out
