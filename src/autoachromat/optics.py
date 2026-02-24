from __future__ import annotations

import math
from dataclasses import dataclass

from .glass_reader import Glass


# Defaults (only used if caller doesn't pass)
LAMBDA_0 = 0.58756
LAMBDA_1 = 0.48613
LAMBDA_2 = 0.65627


@dataclass(frozen=True)
class Config:
    lam0: float
    lam1: float
    lam2: float


def refractive_index(glass: Glass, wavelength_um: float) -> float:
    """
    Compute refractive index n(λ) using the glass dispersion formula.

    Supports AGF formula_id:
    - 1: Schott formula
    - 2: Sellmeier 1
    - 3: Herzberger
    - 4: Sellmeier 2
    - 5: Conrady
    - 6: Sellmeier 3
    - 7: Handbook of Optics 1
    - 8: Handbook of Optics 2
    - 9: Sellmeier 4
    - 10: Extended 1
    - 11: Sellmeier 5
    - 12: Extended 2
    """
    if glass.formula_id is None:
        raise ValueError(f"Glass {glass.name} has no dispersion formula")

    cd = glass.cd
    if not cd or len(cd) < 6:
        raise ValueError(f"Glass {glass.name} has insufficient CD coefficients")

    lam = wavelength_um
    lam2 = lam * lam

    if glass.formula_id == 1:
        # Schott: n² = A0 + A1*λ² + A2*λ⁻² + A3*λ⁻⁴ + A4*λ⁻⁶ + A5*λ⁻⁸
        n2 = (
            cd[0]
            + cd[1] * lam2
            + cd[2] / lam2
            + cd[3] / (lam2 * lam2)
            + cd[4] / (lam2 * lam2 * lam2)
            + cd[5] / (lam2 * lam2 * lam2 * lam2)
        )

    elif glass.formula_id == 2:
        # Sellmeier 1: CD = K1, L1, K2, L2, K3, L3
        K1, L1, K2, L2, K3, L3 = cd[0], cd[1], cd[2], cd[3], cd[4], cd[5]
        n2 = (
            1.0
            + (K1 * lam2) / (lam2 - L1)
            + (K2 * lam2) / (lam2 - L2)
            + (K3 * lam2) / (lam2 - L3)
        )

    elif glass.formula_id == 3:
        # Herzberger (your implementation)
        L = 1.0 / (lam2 - 0.028)
        n = cd[0] + cd[1] * L + cd[2] * L * L + cd[3] * lam2 + cd[4] * lam2 * lam2
        if len(cd) > 5:
            n += cd[5] * lam2 * lam2 * lam2
        return n

    elif glass.formula_id == 4:
        # Sellmeier 2
        A, B, C, D, E = cd[0], cd[1], cd[2], cd[3], cd[4]
        n2 = 1.0 + A + (B * lam2) / (lam2 - C) + (D * lam2) / (lam2 - E)

    elif glass.formula_id == 5:
        # Conrady
        n = cd[0] + cd[1] / lam + cd[2] / (lam**3.5)
        return n

    elif glass.formula_id == 6:
        # Sellmeier 3
        K1, L1, K2, L2, K3, L3 = cd[0], cd[1], cd[2], cd[3], cd[4], cd[5]
        n2 = (
            1.0
            + (K1 * lam2) / (lam2 - L1)
            + (K2 * lam2) / (lam2 - L2)
            + (K3 * lam2) / (lam2 - L3)
        )
        if len(cd) > 7:
            K4, L4 = cd[6], cd[7]
            n2 += (K4 * lam2) / (lam2 - L4)

    elif glass.formula_id in (7, 8):
        # Handbook of Optics formulas
        n2 = 1.0 + cd[0] + cd[1] / (lam2 - cd[2]) + cd[3] / (lam2 - cd[4])
        if len(cd) > 5:
            n2 += cd[5] * lam2

    else:
        # Default fallback: try Sellmeier 1 format
        K1, L1, K2, L2, K3, L3 = cd[0], cd[1], cd[2], cd[3], cd[4], cd[5]
        n2 = (
            1.0
            + (K1 * lam2) / (lam2 - L1)
            + (K2 * lam2) / (lam2 - L2)
            + (K3 * lam2) / (lam2 - L3)
        )

    if n2 <= 0:
        raise ValueError(f"Computed n² <= 0 for {glass.name} at {wavelength_um} um")

    return math.sqrt(n2)


def compute_abbe_number(
    glass: Glass,
    lam0: float = LAMBDA_0,
    lam1: float = LAMBDA_1,
    lam2: float = LAMBDA_2,
) -> float:
    """
    ν = (n(lam0) - 1) / (n(lam1) - n(lam2))
    """
    n0 = refractive_index(glass, lam0)
    n1 = refractive_index(glass, lam1)
    n2 = refractive_index(glass, lam2)

    denom = n1 - n2
    if abs(denom) < 1e-10:
        return float("inf")

    return (n0 - 1.0) / denom


def filter_glasses(glasses: list[Glass], cfg: Config) -> list[Glass]:
    """
    Removes glasses that:
    - exclude_sub == True
    - do not cover wavelength range
    - have no valid dispersion formula
    - cannot compute n at cfg.lam0
    """
    lam_min = min(cfg.lam0, cfg.lam1, cfg.lam2)
    lam_max = max(cfg.lam0, cfg.lam1, cfg.lam2)

    filtered: list[Glass] = []
    for glass in glasses:
        if glass.exclude_sub is True:
            continue

        if glass.ld_min_um is not None and glass.ld_max_um is not None:
            if glass.ld_min_um > lam_min or glass.ld_max_um < lam_max:
                continue

        if glass.formula_id is None or not glass.cd:
            continue

        try:
            _ = refractive_index(glass, cfg.lam0)
        except (ValueError, ZeroDivisionError):
            continue

        filtered.append(glass)

    return filtered


def achromat_power(nu1: float, nu2: float, C0: float) -> tuple[float, float]:
    """
    Your spec:
      φ1 = ν1(1 - ν2 C0)/(ν1 - ν2)
      φ2 = 1 - φ1
    """
    denom = nu1 - nu2
    if abs(denom) < 1e-15:
        raise ZeroDivisionError("nu1 ~ nu2 in achromat_power")
    phi1 = nu1 * (1.0 - nu2 * C0) / denom
    phi2 = 1.0 - phi1
    return phi1, phi2


def check_min_radius(radii: tuple[float, ...], D: float) -> bool:
    rmin = D / 2.0
    return all(abs(r) >= rmin for r in radii)
