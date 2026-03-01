from __future__ import annotations

import logging
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


def _coeffs(n1: float, n2: float, phi1: float, phi2: float) -> dict:
    # Your spec coefficients
    A1 = (phi1**3) * (1.0 + 2.0 / n1)
    B1 = (phi1**3) * (3.0 / (n1 - 1.0))

    A2 = (phi2**3) * (1.0 + 2.0 / n2)
    B2 = (phi2**3) * (3.0 / (n2 - 1.0)) - 4.0 * phi1 * (phi2**2) * (1.0 + 1.0 / n2)

    K1 = (phi1**2) * (1.0 + 1.0 / n1)
    K2 = (phi2**2) * (1.0 + 1.0 / n2)

    L = (
        (phi1**2) * (1.0 / (n1 - 1.0))
        + (phi2**2) * (1.0 / (n2 - 1.0))
        + phi1 * phi2 * (2.0 + 1.0 / n2)
    )

    return dict(A1=A1, B1=B1, A2=A2, B2=B2, K1=K1, K2=K2, L=L)


def _C_const(n1: float, n2: float, phi1: float, phi2: float) -> float:
    """
    User-provided spaced constant term:

    C = φ1^3 * n1/(n1-1)^2
      + φ2^3 * n2/(n2-1)^2
      - φ1 φ2^2 * ( 4/(n2-1) + 1 )
      + φ1^2 φ2 * ( 3 + 2/n2 )
    """
    return (
        (phi1**3) * (n1 / ((n1 - 1.0) ** 2))
        + (phi2**3) * (n2 / ((n2 - 1.0) ** 2))
        - (phi1 * (phi2**2)) * ((4.0 / (n2 - 1.0)) + 1.0)
        + ((phi1**2) * phi2) * (3.0 + (2.0 / n2))
    )


def _solve_Q_pairs(
    inputs: Inputs, c: dict, C_const: float
) -> List[Tuple[float, float]]:
    """
    System:
      A1 Q1^2 + B1 Q1 + A2 Q2^2 + B2 Q2 + (C - P0) = 0
      K1 Q1 + K2 Q2 + (L - W0) = 0

    Eliminate Q2 using the linear equation -> quadratic in Q1
    (given your written form: quadratic in Q1 and Q2, and second equation linear).
    """
    A1, B1, A2, B2, K1, K2, L = (
        c["A1"],
        c["B1"],
        c["A2"],
        c["B2"],
        c["K1"],
        c["K2"],
        c["L"],
    )

    if abs(K2) < inputs.eps:
        return []

    # Q2 = a*Q1 + b
    a = -K1 / K2
    b = (inputs.W0 - L) / K2

    # Substitute into first equation:
    # A1 Q1^2 + B1 Q1 + A2(aQ1+b)^2 + B2(aQ1+b) + (C_const - P0)=0
    q2 = A1 + A2 * (a * a)
    q1 = B1 + A2 * (2.0 * a * b) + B2 * a
    q0 = A2 * (b * b) + B2 * b + (C_const - inputs.P0)

    roots = np.roots(np.array([q2, q1, q0], dtype=float))

    pairs: list[Tuple[float, float]] = []
    for r in roots:
        if abs(r.imag) > inputs.root_imag_tol:
            continue
        Q1 = float(r.real)
        Q2 = float(a * Q1 + b)
        pairs.append((Q1, Q2))

    # Your selection policy: smallest 2-norm first
    pairs.sort(key=lambda t: t[0] * t[0] + t[1] * t[1])
    return pairs


def _radii(
    inputs: Inputs, n1: float, n2: float, phi1: float, phi2: float, Q1: float, Q2: float
) -> Tuple[float, float, float, float]:
    f = inputs.fprime
    if abs(phi1) < inputs.eps or abs(phi2) < inputs.eps:
        raise ZeroDivisionError("phi1 or phi2 near zero in spaced _radii")
    f1 = f / phi1
    f2 = f / phi2

    R1 = f1 / (n1 / (n1 - 1.0) + Q1)
    R2 = f1 / (1.0 + Q1)
    R3 = f2 / (n2 / (n2 - 1.0) + Q2)
    R4 = f2 / (1.0 + Q2)
    return R1, R2, R3, R4


def _Ps_and_PE(
    n1: float, n2: float, phi1: float, Q1: float, Q2: float
) -> Tuple[List[float], float]:
    # Your spec u-chain (encapsulated; easy to fix later)
    u1 = 0.0
    u2 = Q1 * (1.0 - 1.0 / n1)  # u'_1 = u2
    u3 = phi1  # u'_2 = u3
    u4 = phi1 + n2 + Q2 * (n2 - 1.0)  # as you wrote
    u5 = 1.0

    P1 = ((u2 - u1) / (1.0 / n1 - 1.0)) ** 2 * (u2 / n1 - u1)
    P2 = ((u3 - u2) / (1.0 - 1.0 / n1)) ** 2 * (u3 - u2 / n1)
    P3 = ((u4 - u3) / (1.0 / n2 - 1.0)) ** 2 * (u4 / n2 - u3)
    P4 = ((u5 - u4) / (1.0 - 1.0 / n2)) ** 2 * (u5 - u4 / n2)

    Ps = [P1, P2, P3, P4]
    PE = sum(abs(x) for x in Ps) / 4.0
    return Ps, PE


def run_spaced(inputs: Inputs, glasses: list[Glass]) -> list[Candidate]:
    """
    Spaced one-file pipeline:
      enumerate pairs -> achromat -> solve (Q1,Q2) -> radii constraint -> PE -> rank by PE
    """
    gdata = prepare_glass_data(glasses, inputs.lam0, inputs.lam1, inputs.lam2)

    out: list[Candidate] = []

    for g1, n1, nu1 in gdata:
        for g2, n2, nu2 in gdata:
            if (g1.name == g2.name) and (g1.catalog == g2.catalog):
                continue

            # Abbe difference filter (same as cemented)
            if abs(nu1 - nu2) < inputs.min_delta_nu:
                continue

            try:
                phi1, phi2 = achromat_power(nu1, nu2, inputs.C0)
                c = _coeffs(n1, n2, phi1, phi2)
                Cc = _C_const(n1, n2, phi1, phi2)
                pairs = _solve_Q_pairs(inputs, c, Cc)
            except (ValueError, ZeroDivisionError):
                continue

            for Q1, Q2 in pairs:
                try:
                    R1, R2, R3, R4 = _radii(inputs, n1, n2, phi1, phi2, Q1, Q2)
                    if not check_min_radius((R1, R2, R3, R4), inputs.D):
                        continue

                    Ps, PE = _Ps_and_PE(n1, n2, phi1, Q1, Q2)
                    if PE > inputs.max_PE:
                        continue

                    thermal = compute_thermal_metrics(
                        g1, g2, n1, n2, phi1, phi2,
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
                            notes={"C_const_used": Cc},
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

    out.sort(key=lambda c: c.PE if c.PE is not None else float("inf"))
    return out
