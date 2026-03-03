from __future__ import annotations

import logging
import math
from typing import List, Tuple

from .glass_reader import Glass
from .models import Inputs, Candidate
from .optics import (
    achromat_power,
    check_min_radius,
    prepare_glass_data,
)
from .thermal import compute_thermal_metrics

logger = logging.getLogger(__name__)


def _solve_Q_roots(
    n1: float, n2: float, phi1: float, phi2: float, P0: float
) -> List[float]:
    # Quadratic: A Q^2 + B Q + (C - P0) = 0
    A = (1.0 + 2.0 / n1) * phi1 + (1.0 + 2.0 / n2) * phi2
    B = (3.0 / (n1 - 1.0)) * (phi1**2) - (3.0 / (n2 - 1.0)) * (phi2**2) - 2.0 * phi2
    C = (
        (n1 / ((n1 - 1.0) ** 2)) * (phi1**3)
        + (n2 / ((n2 - 1.0) ** 2)) * (phi2**3)
        + ((1.0 / (n2 - 1.0) + 1.0) * (phi2**2))
    )
    a, b, c = A, B, (C - P0)

    if abs(a) < 1e-15:
        if abs(b) < 1e-15:
            return []
        return [(-c) / b]

    disc = b * b - 4.0 * a * c
    if disc < 0:
        return []
    sd = math.sqrt(disc)
    return [(-b + sd) / (2.0 * a), (-b - sd) / (2.0 * a)]


def _coma_W(n1: float, n2: float, phi1: float, phi2: float, Q: float) -> float:
    # W = KQ + L ; K=(A+1)/2 ; L=(B-phi2)/3
    A = (1.0 + 2.0 / n1) * phi1 + (1.0 + 2.0 / n2) * phi2
    B = (3.0 / (n1 - 1.0)) * (phi1**2) - (3.0 / (n2 - 1.0)) * (phi2**2) - 2.0 * phi2
    K = (A + 1.0) / 2.0
    L = (B - phi2) / 3.0
    return K * Q + L


def _radii(
    inputs: Inputs, n1: float, n2: float, phi1: float, Q: float
) -> Tuple[float, float, float]:
    # Your spec:
    # R2 = f'/(phi1+Q)
    # R1 = f'/((n1/(n1-1))*phi1+Q)
    # R3 = f'/((n2/(n2-1))*phi1 + Q - 1/(n2-1))
    f = inputs.fprime
    d2 = phi1 + Q
    d1 = (n1 / (n1 - 1.0)) * phi1 + Q
    d3 = (n2 / (n2 - 1.0)) * phi1 + Q - (1.0 / (n2 - 1.0))
    if abs(d1) < inputs.eps or abs(d2) < inputs.eps or abs(d3) < inputs.eps:
        raise ZeroDivisionError("Radii denominator near zero in cemented _radii")
    R1 = f / d1
    R2 = f / d2
    R3 = f / d3
    return R1, R2, R3


def _P2_and_PE(
    inputs: Inputs,
    n1: float,
    n2: float,
    phi1: float,
    Q: float,
    W: float,
    R2: float,
    W0: float,
) -> Tuple[float, float]:
    # Your spec (encapsulated, easy to adjust later)
    u2 = Q * (1.0 - 1.0 / n1) + phi1
    u3 = Q * (1.0 - 1.0 / n2) + phi1

    denom = 1.0 / n2 - 1.0 / n1
    if abs(denom) < inputs.eps:
        raise ZeroDivisionError("P2 denom too small")

    P2 = ((u3 - u2) / denom) ** 2 * (u3 / n2 - u2 / n1)
    PE = (abs(P2) * (3.0 ** abs(W - W0))) / (R2**2)
    return P2, PE


def run_cemented(inputs: Inputs, glasses: list[Glass]) -> list[Candidate]:
    """
    Cemented one-file pipeline:
      enumerate pairs -> Abbe filter -> achromat power -> solve Q roots
      -> compute W -> radii constraint -> compute PE -> Top-N by |W-W0|
    """
    gdata = prepare_glass_data(glasses, inputs.lam0, inputs.lam1, inputs.lam2)

    kept: list[Candidate] = []

    for g1, n1, nu1 in gdata:
        for g2, n2, nu2 in gdata:
            if (g1.name == g2.name) and (g1.catalog == g2.catalog):
                continue

            # Abbe difference filter
            if abs(nu1 - nu2) < inputs.min_delta_nu:
                continue

            try:
                phi1, phi2 = achromat_power(nu1, nu2, inputs.C0)
            except (ValueError, ZeroDivisionError):
                continue

            for Q in _solve_Q_roots(n1, n2, phi1, phi2, inputs.P0):
                try:
                    W = _coma_W(n1, n2, phi1, phi2, Q)
                    R1, R2, R3 = _radii(inputs, n1, n2, phi1, Q)

                    if not check_min_radius((R1, R2, R3), inputs.D):
                        continue

                    P2, PE = _P2_and_PE(inputs, n1, n2, phi1, Q, W, R2, inputs.W0)
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

                    cand = Candidate(
                        system_type="cemented",
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
                        Q=Q,
                        W=W,
                        P2=P2,
                        radii=[R1, R2, R3],
                        PE=PE,
                        cost1=g1.relative_cost,
                        cost2=g2.relative_cost,
                        formula_id1=g1.formula_id,
                        cd1=list(g1.cd),
                        formula_id2=g2.formula_id,
                        cd2=list(g2.cd),
                        thermal=thermal,
                    )
                except (ValueError, ZeroDivisionError) as exc:
                    logger.debug(
                        "Skipping candidate %s+%s Q=%.4g – %s",
                        g1.name,
                        g2.name,
                        Q,
                        exc,
                    )
                    continue

                # Top-N selection: keep smallest |W-W0|
                if inputs.N <= 0:
                    kept.append(cand)
                elif len(kept) < inputs.N:
                    kept.append(cand)
                else:
                    worst_idx = max(
                        range(len(kept)),
                        key=lambda k: abs((kept[k].W or 0.0) - inputs.W0),
                    )
                    worst_val = abs((kept[worst_idx].W or 0.0) - inputs.W0)
                    cur_val = abs((cand.W or 0.0) - inputs.W0)
                    if cur_val < worst_val:
                        kept[worst_idx] = cand

    kept.sort(
        key=lambda c: (
            abs((c.W or 0.0) - inputs.W0),
            c.PE if c.PE is not None else float("inf"),
        )
    )
    return kept
