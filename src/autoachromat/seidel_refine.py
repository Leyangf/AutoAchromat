"""Thickness-aware Seidel refinement.

After ABCD thickening assigns physical thicknesses and corrects EFL,
the Seidel aberration balance computed under the thin-lens assumption
(h = 1 at all surfaces) is no longer accurate.  This module refines Q
by tracing the actual marginal ray through the thick system, recomputing
the Seidel SA with the real ray heights, and iterating until the thick-lens
SA matches the target P0.
"""

from __future__ import annotations

import dataclasses
import math
from dataclasses import dataclass

from .models import Candidate, Inputs, ThickPrescription

# ---------------------------------------------------------------------------
# Newton solver parameters
# ---------------------------------------------------------------------------

_REFINE_TOL = 1e-6
_REFINE_MAX_ITER = 5
_NEWTON_DQ = 1e-7
_MAX_STEP_FACTOR = 0.5  # |ΔQ| clamped to this fraction of |Q|


# ---------------------------------------------------------------------------
# Paraxial marginal-ray trace
# ---------------------------------------------------------------------------


@dataclass
class TraceState:
    """Paraxial ray state at each surface."""

    h: list[float]
    u_before: list[float]
    u_after: list[float]


def _refract(h: float, u_before: float, n_before: float, n_after: float, R: float):
    """Paraxial refraction: returns u_after."""
    nu_after = n_before * u_before - h * (n_after - n_before) / R
    return nu_after / n_after


def _transfer(h: float, u: float, t: float) -> float:
    """Paraxial transfer: returns h at next surface."""
    return h + t * u


def trace_marginal_cemented(
    R1: float,
    R2: float,
    R3: float,
    t1: float,
    t2: float,
    n1: float,
    n2: float,
) -> TraceState:
    """Trace normalised marginal ray (h=1, u=0) through thick cemented doublet."""
    h1 = 1.0
    u1_before = 0.0
    u1_after = _refract(h1, u1_before, 1.0, n1, R1)

    h2 = _transfer(h1, u1_after, t1)
    u2_before = u1_after
    u2_after = _refract(h2, u2_before, n1, n2, R2)

    h3 = _transfer(h2, u2_after, t2)
    u3_before = u2_after
    u3_after = _refract(h3, u3_before, n2, 1.0, R3)

    return TraceState(
        h=[h1, h2, h3],
        u_before=[u1_before, u2_before, u3_before],
        u_after=[u1_after, u2_after, u3_after],
    )


def trace_marginal_spaced(
    R1: float,
    R2: float,
    R3: float,
    R4: float,
    t1: float,
    t2: float,
    air_gap: float,
    n1: float,
    n2: float,
) -> TraceState:
    """Trace normalised marginal ray through thick spaced doublet (4 surfaces)."""
    # Surface 1: air(1) → glass(n1)
    h1 = 1.0
    u1_before = 0.0
    u1_after = _refract(h1, u1_before, 1.0, n1, R1)

    # Transfer through glass 1
    h2 = _transfer(h1, u1_after, t1)
    u2_before = u1_after
    u2_after = _refract(h2, u2_before, n1, 1.0, R2)

    # Transfer through air gap
    h3 = _transfer(h2, u2_after, air_gap)
    u3_before = u2_after
    u3_after = _refract(h3, u3_before, 1.0, n2, R3)

    # Transfer through glass 2
    h4 = _transfer(h3, u3_after, t2)
    u4_before = u3_after
    u4_after = _refract(h4, u4_before, n2, 1.0, R4)

    return TraceState(
        h=[h1, h2, h3, h4],
        u_before=[u1_before, u2_before, u3_before, u4_before],
        u_after=[u1_after, u2_after, u3_after, u4_after],
    )


# ---------------------------------------------------------------------------
# Thick-lens Seidel spherical aberration
# ---------------------------------------------------------------------------


def _thick_seidel_P(
    trace: TraceState, n_before: list[float], n_after: list[float]
) -> float:
    """Total third-order SA from traced ray data.

    Includes the h_k ray-height factor at each surface:
        S_I,k = h_k · A_k² · Δ(u/n)_k
    where A_k = Δu_k / Δ(1/n_k) is the Abbe refraction invariant.
    """
    P_total = 0.0
    for k in range(len(trace.h)):
        delta_u = trace.u_after[k] - trace.u_before[k]
        delta_inv_n = 1.0 / n_after[k] - 1.0 / n_before[k]
        delta_u_over_n = (
            trace.u_after[k] / n_after[k] - trace.u_before[k] / n_before[k]
        )

        if abs(delta_inv_n) < 1e-15:
            continue

        A_k_sq = (delta_u / delta_inv_n) ** 2
        P_total += trace.h[k] * A_k_sq * delta_u_over_n

    return P_total


def thick_seidel_P_cemented(
    trace: TraceState, n1: float, n2: float
) -> float:
    """Thick-lens Seidel SA for a cemented doublet (3 surfaces)."""
    return _thick_seidel_P(trace, [1.0, n1, n2], [n1, n2, 1.0])


def thick_seidel_P_spaced(
    trace: TraceState, n1: float, n2: float
) -> float:
    """Thick-lens Seidel SA for a spaced doublet (4 surfaces)."""
    return _thick_seidel_P(trace, [1.0, n1, 1.0, n2], [n1, 1.0, n2, 1.0])


# ---------------------------------------------------------------------------
# Q → P composite functions (for Newton solver)
# ---------------------------------------------------------------------------


def _P_from_Q_cemented(
    Q: float,
    phi1: float,
    n1: float,
    n2: float,
    t1: float,
    t2: float,
    f: float,
) -> float | None:
    """Compute thick-lens Seidel SA for a given Q (cemented doublet)."""
    d1 = n1 / (n1 - 1.0) * phi1 + Q
    d2 = phi1 + Q
    d3 = n2 / (n2 - 1.0) * phi1 + Q - 1.0 / (n2 - 1.0)

    if abs(d1) < 1e-12 or abs(d2) < 1e-12 or abs(d3) < 1e-12:
        return None

    R1, R2, R3 = f / d1, f / d2, f / d3

    try:
        trace = trace_marginal_cemented(R1, R2, R3, t1, t2, n1, n2)
    except (ValueError, ZeroDivisionError):
        return None

    P = thick_seidel_P_cemented(trace, n1, n2)
    return P if math.isfinite(P) else None


def _P_from_Q1_spaced(
    Q1: float,
    a: float,
    b: float,
    phi1: float,
    phi2: float,
    n1: float,
    n2: float,
    t1: float,
    t2: float,
    air_gap: float,
    f: float,
) -> float | None:
    """Compute thick-lens Seidel SA for a given Q1 on the coma line (spaced)."""
    Q2 = a * Q1 + b

    if abs(phi1) < 1e-12 or abs(phi2) < 1e-12:
        return None

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
        return None

    R1 = f1 / denom_R1
    R2 = f1 / denom_R2
    R3 = f2 / denom_R3
    R4 = f2 / denom_R4

    try:
        trace = trace_marginal_spaced(R1, R2, R3, R4, t1, t2, air_gap, n1, n2)
    except (ValueError, ZeroDivisionError):
        return None

    P = thick_seidel_P_spaced(trace, n1, n2)
    return P if math.isfinite(P) else None


# ---------------------------------------------------------------------------
# Newton iteration
# ---------------------------------------------------------------------------


def _newton_1d(eval_fn, x0: float, target: float) -> float | None:
    """Generic 1D Newton solver: find x such that eval_fn(x) = target.

    Returns refined x, or None on failure.  Falls back to current best
    if the derivative is near-zero or evaluation fails.
    """
    x = x0

    for _ in range(_REFINE_MAX_ITER):
        P = eval_fn(x)
        if P is None:
            return None

        residual = P - target
        if abs(residual) < _REFINE_TOL:
            return x

        P_plus = eval_fn(x + _NEWTON_DQ)
        if P_plus is None:
            return x

        dPdx = (P_plus - P) / _NEWTON_DQ
        if abs(dPdx) < 1e-15:
            return x

        step = residual / dPdx

        # Clamp step size to prevent wild jumps
        max_step = _MAX_STEP_FACTOR * max(abs(x), 1.0)
        if abs(step) > max_step:
            step = math.copysign(max_step, step)

        x -= step

    return x


# ---------------------------------------------------------------------------
# Cemented refinement
# ---------------------------------------------------------------------------


def _refine_cemented(
    rx: ThickPrescription,
    cand: Candidate,
    inp: Inputs,
) -> ThickPrescription:
    """Refine Q for a cemented doublet, recompute radii, re-thicken."""
    Q_old = cand.Q
    if Q_old is None:
        return rx

    e1, e2 = rx.elements[0], rx.elements[1]
    t1, t2 = e1.t_center, e2.t_center

    def eval_fn(Q):
        return _P_from_Q_cemented(Q, cand.phi1, cand.n1, cand.n2, t1, t2, inp.fprime)

    Q_new = _newton_1d(eval_fn, Q_old, inp.P0)

    if Q_new is None or abs(Q_new - Q_old) < 1e-12:
        return rx

    from .cemented import radii_from_Q

    try:
        R1, R2, R3 = radii_from_Q(inp, cand.n1, cand.n2, cand.phi1, Q_new)
    except (ValueError, ZeroDivisionError):
        return rx

    cand_refined = dataclasses.replace(cand, Q=Q_new, radii=[R1, R2, R3])

    from .thickening import thicken_cemented

    rx_new = thicken_cemented(cand_refined, inp)
    return rx_new if rx_new is not None else rx


# ---------------------------------------------------------------------------
# Spaced refinement
# ---------------------------------------------------------------------------


def _refine_spaced(
    rx: ThickPrescription,
    cand: Candidate,
    inp: Inputs,
) -> ThickPrescription:
    """Refine Q1 along the coma-zero line for a spaced doublet."""
    Q1_old = cand.Q1
    Q2_old = cand.Q2
    if Q1_old is None or Q2_old is None:
        return rx

    from .spaced import _coeffs

    c = _coeffs(cand.n1, cand.n2, cand.phi1, cand.phi2)
    K1, K2, L = c["K1"], c["K2"], c["L"]
    if abs(K2) < 1e-12:
        return rx

    a = -K1 / K2
    b = (inp.W0 - L) / K2

    e1, e2 = rx.elements[0], rx.elements[1]
    t1, t2 = e1.t_center, e2.t_center
    air_gap = rx.air_gap if rx.air_gap is not None else inp.air_gap

    def eval_fn(Q1):
        return _P_from_Q1_spaced(
            Q1, a, b, cand.phi1, cand.phi2, cand.n1, cand.n2,
            t1, t2, air_gap, inp.fprime,
        )

    Q1_new = _newton_1d(eval_fn, Q1_old, inp.P0)

    if Q1_new is None or abs(Q1_new - Q1_old) < 1e-12:
        return rx

    Q2_new = a * Q1_new + b

    from .spaced import radii_from_Q_pair

    try:
        R1, R2, R3, R4 = radii_from_Q_pair(
            inp, cand.n1, cand.n2, cand.phi1, cand.phi2, Q1_new, Q2_new,
        )
    except (ValueError, ZeroDivisionError):
        return rx

    cand_refined = dataclasses.replace(
        cand, Q1=Q1_new, Q2=Q2_new, radii=[R1, R2, R3, R4],
    )

    from .thickening import thicken_spaced

    rx_new = thicken_spaced(cand_refined, inp)
    return rx_new if rx_new is not None else rx


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------


def refine_seidel(
    rx: ThickPrescription,
    cand: Candidate,
    inp: Inputs,
) -> ThickPrescription:
    """Refine a thick prescription so the Seidel SA accounts for finite thickness.

    Dispatches to cemented or spaced refinement.
    Returns the original *rx* unchanged if refinement is unnecessary or fails.
    """
    if not isinstance(rx, ThickPrescription):
        return rx

    try:
        if cand.system_type == "cemented":
            return _refine_cemented(rx, cand, inp)
        elif cand.system_type == "spaced":
            return _refine_spaced(rx, cand, inp)
    except Exception:
        return rx

    return rx
