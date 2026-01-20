"""
Analytical synthesis pipeline for achromat design (2022 paper).

This module implements the main control flow for glass pair enumeration,
constraint solving, and candidate ranking based on:
- Nguyen 2022: "Automation of Synthesis and Ranking of Cemented and Air-Spaced Doublets"
- Ivanova et al. 2017: "Computer tool for achromatic and aplanatic cemented doublet design"

References:
- Thin-lens 3rd-order aberration theory
- Seidel aberration coefficients
- Standard Fraunhofer wavelengths: F=486.13nm, d=587.56nm, C=656.27nm
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

# from pathlib import Path
from typing import Iterator

from src.glass_reader import Glass, load_catalog
from src.score import Candidate, TopNContainer, compute_pe, rank_key
from src import io


# =============================================================================
# Constants - Default Wavelengths (micrometers)
# =============================================================================

LAMBDA_0 = 0.58756  # Primary wavelength (default: yellow)
LAMBDA_1 = 0.48613  # Short wavelength (default: blue)
LAMBDA_2 = 0.65627  # Long wavelength (default: red)


# =============================================================================
# Configuration
# =============================================================================


@dataclass
class Config:
    """Configuration for the 2022 synthesis pipeline."""

    # Glass catalog paths
    agf_paths: list[str] = field(default_factory=list)

    # System type: "cemented" or "air_spaced"
    system_type: str = "cemented"

    # Optical parameters
    f: float = 100.0  # Focal length (mm)
    D: float = 20.0  # Aperture diameter (mm)
    P0: float = 0.0  # Target spherical aberration (Seidel S_I)
    W0: float = 0.0  # Target coma (Seidel S_II)
    C0: float = 0.0  # Target chromatic aberration

    # Wavelengths (micrometers)
    lam0: float = LAMBDA_0  # Primary wavelength
    lam1: float = LAMBDA_1  # Short wavelength
    lam2: float = LAMBDA_2  # Long wavelength

    # Filtering and selection parameters
    min_delta_nu: float = 10.0  # Minimum Abbe number difference
    max_PE: float = 0.01  # Maximum PE threshold
    N_keep: int = 30  # Number of top candidates to keep
    allow_repeat: bool = False  # Allow same glass for both elements

    # Air-spaced doublet parameters
    d_air: float = 2.0  # Air gap thickness (mm) for air-spaced doublets

    # Real-ray aberration parameters (post-ranking)
    crown_lens_thickness_mm: float = 5.0  # Thickness of crown lens (higher Abbe number)
    flint_lens_thickness_mm: float = 3.0  # Thickness of flint lens (lower Abbe number)
    field_for_aberration_deg: float = 0.0  # Field angle for coma calculation (degrees)

    # Output directory
    out_dir: str = "output"


def default_cfg() -> Config:
    """
    Return a default, runnable configuration.

    Returns:
        A Config with reasonable default values for testing.
    """
    return Config(
        agf_paths=["glass_database/SCHOTT.AGF"],
        system_type="cemented",
        f=100.0,
        D=20.0,
        P0=0.0,
        W0=0.0,
        C0=0.0,
        lam0=LAMBDA_0,
        lam1=LAMBDA_1,
        lam2=LAMBDA_2,
        min_delta_nu=10.0,
        max_PE=0.01,
        N_keep=30,
        allow_repeat=False,
        d_air=2.0,
        crown_lens_thickness_mm=5.0,
        flint_lens_thickness_mm=3.0,
        field_for_aberration_deg=0.0,
        out_dir="output",
    )


# =============================================================================
# Surface Invariant Helpers (Alpha-based P and W)
# =============================================================================


def compute_alphas_from_radii(
    n_list: list[float], R_list: list[float], fprime: float, alpha_last: float = 1.0
) -> list[float]:
    """
    Compute alpha parameters from radii using backward recursion.

    Given m surfaces with radii R_1..R_m and media indices n_1..n_{m+1},
    set alpha_{m+1} = alpha_last (normalization).
    Then for i = m ... 1:
        n_{i+1}*alpha_{i+1} - n_i*alpha_i = (n_{i+1}*fprime)/R_i
        => alpha_i = (n_{i+1}*alpha_{i+1} - (n_{i+1}*fprime)/R_i) / n_i

    Args:
        n_list: List of refractive indices [n_1, n_2, ..., n_{m+1}].
                Length = number of surfaces + 1.
        R_list: List of surface radii [R_1, R_2, ..., R_m].
                Length = number of surfaces.
        fprime: Focal length (mm), used for scaling.
        alpha_last: Normalization value for alpha_{m+1}, default 1.0.

    Returns:
        List of alpha values [alpha_1, alpha_2, ..., alpha_{m+1}].
        Length = len(n_list) = len(R_list) + 1.
    """
    m = len(R_list)  # Number of surfaces
    if len(n_list) != m + 1:
        raise ValueError(
            f"n_list length ({len(n_list)}) must be R_list length + 1 ({m + 1})"
        )

    # Initialize alphas with final value
    alphas = [0.0] * (m + 1)
    alphas[m] = alpha_last  # alpha_{m+1}

    # Backward recursion: i = m, m-1, ..., 1 (but 0-indexed: i = m-1, ..., 0)
    for i in range(m - 1, -1, -1):
        # Surface i+1 (1-indexed) corresponds to R_list[i], between n_list[i] and n_list[i+1]
        n_i = n_list[i]  # n_i (medium before surface)
        n_i1 = n_list[i + 1]  # n_{i+1} (medium after surface)
        alpha_i1 = alphas[i + 1]  # alpha_{i+1}
        R_i = R_list[i]  # R_i

        # Handle large radius (nearly flat surface)
        if abs(R_i) > 1e9:
            # fprime / R_i ~ 0, so alpha_i = n_{i+1} * alpha_{i+1} / n_i
            if abs(n_i) < 1e-15:
                alphas[i] = 0.0
            else:
                alphas[i] = n_i1 * alpha_i1 / n_i
        else:
            # Normal case
            if abs(n_i) < 1e-15:
                alphas[i] = 0.0
            else:
                alphas[i] = (n_i1 * alpha_i1 - (n_i1 * fprime) / R_i) / n_i

    return alphas


def surface_PW(
    nk: float, nk1: float, alpha_k: float, alpha_k1: float
) -> tuple[float, float]:
    """
    Compute spherical aberration P_k and coma W_k for surface k.

    Using alpha-parameter formulation:
        A_k = (alpha_{k+1} - alpha_k) / ((1/n_{k+1}) - (1/n_k))
        B_k = (alpha_{k+1}/n_{k+1}) - (alpha_k/n_k)
        W_k = A_k * B_k
        P_k = A_k^2 * B_k

    Args:
        nk: Refractive index before surface (n_k).
        nk1: Refractive index after surface (n_{k+1}).
        alpha_k: Alpha parameter before surface.
        alpha_k1: Alpha parameter after surface.

    Returns:
        Tuple (P_k, W_k) of surface invariants.
        Returns (0.0, 0.0) if denominator is near zero.
    """
    # Compute 1/n terms
    if abs(nk) < 1e-15 or abs(nk1) < 1e-15:
        return (0.0, 0.0)

    inv_nk = 1.0 / nk
    inv_nk1 = 1.0 / nk1

    # Denominator for A_k
    denom = inv_nk1 - inv_nk
    if abs(denom) < 1e-15:
        # Same refractive index on both sides (no refraction)
        return (0.0, 0.0)

    # A_k = (alpha_{k+1} - alpha_k) / ((1/n_{k+1}) - (1/n_k))
    A_k = (alpha_k1 - alpha_k) / denom

    # B_k = (alpha_{k+1}/n_{k+1}) - (alpha_k/n_k)
    B_k = alpha_k1 * inv_nk1 - alpha_k * inv_nk

    # W_k = A_k * B_k
    W_k = A_k * B_k

    # P_k = A_k^2 * B_k
    P_k = A_k * A_k * B_k

    return (P_k, W_k)


# =============================================================================
# Glass Dispersion Functions
# =============================================================================


def refractive_index(glass: Glass, wavelength_um: float) -> float:
    """
    Compute refractive index n(λ) using the glass dispersion formula.

    Supports AGF formula_id:
    - 1: Schott formula
    - 2: Sellmeier 1 (most common for SCHOTT, OHARA, CDGM)
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

    Args:
        glass: Glass object with formula_id and cd coefficients.
        wavelength_um: Wavelength in micrometers.

    Returns:
        Refractive index at the given wavelength.

    Raises:
        ValueError: If formula not supported or coefficients missing.
    """
    if glass.formula_id is None:
        raise ValueError(f"Glass {glass.name} has no dispersion formula")

    cd = glass.cd
    if not cd or len(cd) < 6:
        raise ValueError(f"Glass {glass.name} has insufficient CD coefficients")

    lam = wavelength_um
    lam2 = lam * lam

    if glass.formula_id == 1:
        # Schott formula: n² = A0 + A1*λ² + A2*λ⁻² + A3*λ⁻⁴ + A4*λ⁻⁶ + A5*λ⁻⁸
        n2 = (
            cd[0]
            + cd[1] * lam2
            + cd[2] / lam2
            + cd[3] / (lam2 * lam2)
            + cd[4] / (lam2 * lam2 * lam2)
            + cd[5] / (lam2 * lam2 * lam2 * lam2)
        )

    elif glass.formula_id == 2:
        # Sellmeier 1: n² - 1 = K1*λ²/(λ²-L1) + K2*λ²/(λ²-L2) + K3*λ²/(λ²-L3)
        # CD coefficients: K1, L1, K2, L2, K3, L3
        K1, L1, K2, L2, K3, L3 = cd[0], cd[1], cd[2], cd[3], cd[4], cd[5]
        n2 = (
            1.0
            + (K1 * lam2) / (lam2 - L1)
            + (K2 * lam2) / (lam2 - L2)
            + (K3 * lam2) / (lam2 - L3)
        )

    elif glass.formula_id == 3:
        # Herzberger: n = A + B*λ² + C*λ⁴ + D*λ² + E*λ⁴ + F/(λ²-0.028)
        # Simplified: n = A + B*L + C*L² + D/(λ²-0.028) + E/(λ²-0.028)²
        # where L = 1/(λ²-0.028)
        L = 1.0 / (lam2 - 0.028)
        n = cd[0] + cd[1] * L + cd[2] * L * L + cd[3] * lam2 + cd[4] * lam2 * lam2
        if len(cd) > 5:
            n += cd[5] * lam2 * lam2 * lam2
        return n

    elif glass.formula_id == 4:
        # Sellmeier 2: n² - 1 = A + B*λ²/(λ²-C) + D*λ²/(λ²-E)
        A, B, C, D, E = cd[0], cd[1], cd[2], cd[3], cd[4]
        n2 = 1.0 + A + (B * lam2) / (lam2 - C) + (D * lam2) / (lam2 - E)

    elif glass.formula_id == 5:
        # Conrady: n = A + B/λ + C/λ^3.5
        n = cd[0] + cd[1] / lam + cd[2] / (lam**3.5)
        return n

    elif glass.formula_id == 6:
        # Sellmeier 3: n² - 1 = K1*λ²/(λ²-L1) + K2*λ²/(λ²-L2) + K3*λ²/(λ²-L3) + K4*λ²/(λ²-L4)
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
        # n² - 1 = A + B/(λ²-C) + D/(λ²-E)...
        n2 = 1.0 + cd[0] + cd[1] / (lam2 - cd[2]) + cd[3] / (lam2 - cd[4])
        if len(cd) > 5:
            n2 += cd[5] * lam2

    else:
        # Default: try Sellmeier 1 format (most common)
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


def compute_abbe_number(glass: Glass) -> float:
    """
    Compute Abbe number νd = (nd - 1) / (nF - nC).

    Uses wavelengths: lambda_0 (primary), lambda_1 (short), lambda_2 (long).

    Args:
        glass: Glass object.

    Returns:
        Abbe number (dispersion).
    """
    n0 = refractive_index(glass, LAMBDA_0)
    n1 = refractive_index(glass, LAMBDA_1)
    n2 = refractive_index(glass, LAMBDA_2)

    denom = n1 - n2
    if abs(denom) < 1e-10:
        return float("inf")

    return (n0 - 1.0) / denom


def compute_dispersion(glass: Glass) -> float:
    """
    Compute mean dispersion V = n1 - n2.

    Args:
        glass: Glass object.

    Returns:
        Mean dispersion.
    """
    n1 = refractive_index(glass, LAMBDA_1)
    n2 = refractive_index(glass, LAMBDA_2)
    return n1 - n2


# =============================================================================
# Glass Filtering
# =============================================================================


def filter_glasses(glasses: list[Glass], cfg: Config) -> list[Glass]:
    """
    Filter glasses based on configuration constraints.

    Removes glasses that:
    - Have exclude_sub == True
    - Do not cover the required wavelength range
    - Have no valid dispersion formula

    Args:
        glasses: List of Glass objects to filter.
        cfg: Configuration containing wavelength parameters.

    Returns:
        Filtered list of Glass objects.
    """
    # Determine required wavelength range
    lam_min = min(cfg.lam0, cfg.lam1, cfg.lam2)
    lam_max = max(cfg.lam0, cfg.lam1, cfg.lam2)

    filtered: list[Glass] = []

    for glass in glasses:
        # Skip if exclude_sub is True
        if glass.exclude_sub is True:
            continue

        # Check wavelength coverage
        if glass.ld_min_um is not None and glass.ld_max_um is not None:
            if glass.ld_min_um > lam_min or glass.ld_max_um < lam_max:
                continue

        # Check that dispersion formula is available
        if glass.formula_id is None or not glass.cd:
            continue

        # Verify we can compute refractive index
        try:
            _ = refractive_index(glass, LAMBDA_0)
        except (ValueError, ZeroDivisionError):
            continue

        filtered.append(glass)

    return filtered


# =============================================================================
# Glass Pair Enumeration
# =============================================================================


def enumerate_pairs(
    glasses: list[Glass], allow_repeat: bool
) -> Iterator[tuple[Glass, Glass]]:
    """
    Enumerate all valid glass pairs.

    Args:
        glasses: List of Glass objects.
        allow_repeat: If False, skip pairs where both glasses are the same.

    Yields:
        Tuples of (g1, g2) glass pairs.
    """
    n = len(glasses)
    for i in range(n):
        for j in range(n):
            if not allow_repeat and i == j:
                continue
            yield (glasses[i], glasses[j])


# =============================================================================
# Thin-Lens Aberration Coefficients
# =============================================================================
# Reference: Nguyen 2022, Ivanova 2017
# For a thin lens with refractive index n and shape factor X (bending):
#   X = (R2 + R1) / (R2 - R1)  (R1=front, R2=back radius)
# Position factor Y:
#   Y = (u' + u) / (u' - u)  for object/image conjugates
# For object at infinity: Y = 1
#
# Seidel spherical aberration coefficient P (thin lens):
#   P = a*X² + b*X + c
# where a, b, c depend on n and Y.
#
# For achromatic doublet: φ = φ1 + φ2 (total power)
# Achromatization: φ1/ν1 + φ2/ν2 = 0
# =============================================================================


def thin_lens_P_coefficients(n: float, Y: float = 1.0) -> tuple[float, float, float]:
    """
    Compute coefficients a, b, c for thin-lens spherical aberration.

    P = a*X² + b*X + c where X is the shape factor.

    Reference: Ivanova 2017, Eq. (3)-(5), Nguyen 2022 Section 2.

    For object at infinity (Y=1):
    a = (n+2) / (4*n*(n-1)²)
    b = (n+1) / (n*(n-1))
    c = (3*n+2)*(n-1) / (4*n)

    Args:
        n: Refractive index at design wavelength.
        Y: Position factor (default 1 for object at infinity).

    Returns:
        Tuple (a, b, c) of quadratic coefficients.
    """
    # Coefficients for spherical aberration P = a*X² + b*X + c
    # Reference: Kingslake "Lens Design Fundamentals", Ivanova 2017 Eq. (3-5)

    n2 = n * n
    nm1 = n - 1.0
    nm1_sq = nm1 * nm1

    # For object at infinity (Y = 1)
    a = (n + 2.0) / (4.0 * n * nm1_sq)
    b = (n + 1.0) / (n * nm1)
    c = (3.0 * n + 2.0) * nm1 / (4.0 * n)

    # General case with position factor Y
    if abs(Y - 1.0) > 1e-10:
        # Adjust for non-infinite conjugate
        # a = (n+2)/(4n(n-1)²) * Y²  (simplified)
        # These are the full expressions including Y dependence
        Y2 = Y * Y
        a = ((n + 2.0) * Y2 + 2.0 * (n2 - 1.0) * Y + (3.0 * n + 2.0) * nm1_sq) / (
            4.0 * n * nm1_sq
        )
        b = 2.0 * (n + 1.0) * Y / (n * nm1)
        c = (3.0 * n + 2.0) * nm1 / (4.0 * n)

    return a, b, c


def thin_lens_W_coefficient(n: float, Y: float = 1.0) -> float:
    """
    Compute coefficient for thin-lens coma.

    W = w * X + w0 where X is the shape factor.

    For object at infinity (Y=1):
    W depends linearly on X: W = ((n+1)/(2n(n-1))) * X + 1

    Args:
        n: Refractive index at design wavelength.
        Y: Position factor (default 1 for object at infinity).

    Returns:
        Linear coefficient w for coma.
    """
    # Coma coefficient: W = w*X + w0
    # Reference: Kingslake, Ivanova 2017
    w = (n + 1.0) / (2.0 * n * (n - 1.0))
    return w


# =============================================================================
# Achromatization Solver (solve_C)
# =============================================================================


def solve_C(g1: Glass, g2: Glass, cfg: Config) -> list[Candidate]:
    """
    Solve the achromatization constraint to determine lens powers.

    For a doublet with total power φ = 1/f:
    - Achromatization: φ1/ν1 + φ2/ν2 = 0  (for C=0)
    - Power sum: φ1 + φ2 = φ

    Solution:
    φ1 = φ * ν1 / (ν1 - ν2)
    φ2 = -φ * ν2 / (ν1 - ν2)

    Reference: Nguyen 2022 Section 2, Ivanova 2017 Eq. (1-2)

    Args:
        g1: First glass (crown).
        g2: Second glass (flint).
        cfg: Configuration.

    Returns:
        List containing one Candidate if solution exists, empty list otherwise.
    """
    try:
        # Get refractive indices at primary wavelength
        n1 = refractive_index(g1, LAMBDA_0)
        n2 = refractive_index(g2, LAMBDA_0)

        # Compute Abbe numbers
        nu1 = compute_abbe_number(g1)
        nu2 = compute_abbe_number(g2)

        # Check for valid Abbe numbers
        if not math.isfinite(nu1) or not math.isfinite(nu2):
            return []

        # Abbe number difference
        delta_nu = nu1 - nu2

        # Need sufficient difference for achromatization
        if abs(delta_nu) < cfg.min_delta_nu:
            return []

        # Total power
        phi_total = 1.0 / cfg.f

        # Solve achromatization equations
        # Reference: Nguyen 2022 Eq. (2), Ivanova 2017 Eq. (1)
        # φ1/ν1 + φ2/ν2 = C0  (C0=0 for perfect achromatization)
        # φ1 + φ2 = φ_total

        # For C0 = 0 (zero chromatic aberration at primary wavelength):
        phi1 = phi_total * nu1 / delta_nu
        phi2 = -phi_total * nu2 / delta_nu

        # Create candidate with computed optical parameters
        candidate = Candidate(
            system_type=cfg.system_type,
            P2=None,
            W=None,
            R2=None,
            Pi=[],
            PE=None,
        )

        # Store glass references and computed values in meta dict
        candidate.g1 = g1  # type: ignore
        candidate.g2 = g2  # type: ignore
        candidate.meta = {  # type: ignore
            "n1": n1,
            "n2": n2,
            "nu1": nu1,
            "nu2": nu2,
            "delta_nu": delta_nu,
            "phi1": phi1,
            "phi2": phi2,
            "phi_total": phi_total,
        }

        return [candidate]

    except (ValueError, ZeroDivisionError):
        return []


# =============================================================================
# Spherical Aberration & Coma Solver (solve_PW)
# =============================================================================


def solve_PW_cemented(candidate: Candidate, cfg: Config) -> bool:
    """
    Solve for shape factors in cemented doublet to achieve P=P0.

    For cemented doublet, there is one free parameter (Q = X1 - X2 or similar).
    The spherical aberration P is quadratic in Q:
        P(Q) = A*Q² + B*Q + C = P0

    For real solutions, discriminant must be >= 0.

    Reference: Nguyen 2022 Section 2.1, Ivanova 2017 Eq. (6-10)

    Args:
        candidate: Candidate with meta dict containing phi1, phi2, n1, n2.
        cfg: Configuration with P0, W0.

    Returns:
        True if valid solution exists, False otherwise.
    """
    meta = getattr(candidate, "meta", None)
    if meta is None:
        return False

    n1 = meta["n1"]
    n2 = meta["n2"]
    phi1 = meta["phi1"]
    phi2 = meta["phi2"]

    # Get aberration coefficients for each lens
    a1, b1, c1 = thin_lens_P_coefficients(n1)
    a2, b2, c2 = thin_lens_P_coefficients(n2)

    # For cemented doublet with shape continuity at cemented surface:
    # The system spherical aberration is:
    # P = p1*φ1³*P1(X1) + p2*φ2³*P2(X2)
    # where P1, P2 are normalized single-lens aberrations
    # and X2 is related to X1 through the cemented surface condition.

    # Reference: Ivanova 2017 Eq. (7)
    # For cemented doublet: X2 = (n2-1)/(n1-1) * X1 + (n2-n1)/((n1-1)*(n2-1))

    # Scaling factors for aberration contributions
    # P_total = h⁴ * (φ1³*P1 + φ2³*P2) where h = y/f (normalized height)
    phi1_3 = phi1 * phi1 * phi1
    phi2_3 = phi2 * phi2 * phi2

    # Combined quadratic in X1:
    # Let K = (n2-1)/(n1-1), M = (n2-n1)/((n1-1)*(n2-1))
    # X2 = K*X1 + M
    # P2(X2) = a2*(K*X1+M)² + b2*(K*X1+M) + c2

    K = (n2 - 1.0) / (n1 - 1.0)
    M = (n2 - n1) / ((n1 - 1.0) * (n2 - 1.0))

    # Expand P2(K*X1 + M):
    # = a2*K²*X1² + 2*a2*K*M*X1 + a2*M² + b2*K*X1 + b2*M + c2
    # = a2*K²*X1² + (2*a2*K*M + b2*K)*X1 + (a2*M² + b2*M + c2)

    # Total P = φ1³*(a1*X1² + b1*X1 + c1) + φ2³*(a2*K²*X1² + ...)
    # = (φ1³*a1 + φ2³*a2*K²)*X1² + (φ1³*b1 + φ2³*(2*a2*K*M + b2*K))*X1 + ...

    A_coeff = phi1_3 * a1 + phi2_3 * a2 * K * K
    B_coeff = phi1_3 * b1 + phi2_3 * (2.0 * a2 * K * M + b2 * K)
    C_coeff = phi1_3 * c1 + phi2_3 * (a2 * M * M + b2 * M + c2) - cfg.P0

    # Solve A*X1² + B*X1 + C = 0
    discriminant = B_coeff * B_coeff - 4.0 * A_coeff * C_coeff

    if discriminant < 0:
        return False

    # Store coefficients for later use
    meta["A_coeff"] = A_coeff
    meta["B_coeff"] = B_coeff
    meta["C_coeff"] = C_coeff
    meta["discriminant"] = discriminant
    meta["K"] = K
    meta["M"] = M

    # Compute solutions
    sqrt_disc = math.sqrt(discriminant)
    if abs(A_coeff) > 1e-15:
        X1_solutions = [
            (-B_coeff + sqrt_disc) / (2.0 * A_coeff),
            (-B_coeff - sqrt_disc) / (2.0 * A_coeff),
        ]
    else:
        # Linear case
        if abs(B_coeff) > 1e-15:
            X1_solutions = [-C_coeff / B_coeff]
        else:
            return False

    # Store solutions
    meta["X1_solutions"] = X1_solutions

    # For each X1, compute X2 and check coma
    valid_solutions = []
    for X1 in X1_solutions:
        X2 = K * X1 + M

        # Compute coma for this solution
        # W = w1*φ1²*X1 + w2*φ2²*X2 + constant terms
        w1 = thin_lens_W_coefficient(n1)
        w2 = thin_lens_W_coefficient(n2)

        phi1_2 = phi1 * phi1
        phi2_2 = phi2 * phi2

        # Simplified coma calculation
        W_val = phi1_2 * (w1 * X1 + 1.0) + phi2_2 * (w2 * X2 + 1.0)

        valid_solutions.append({"X1": X1, "X2": X2, "W": W_val})

    if not valid_solutions:
        return False

    # Select solution with minimum |W - W0|
    best = min(valid_solutions, key=lambda s: abs(s["W"] - cfg.W0))
    meta["X1"] = best["X1"]
    meta["X2"] = best["X2"]
    meta["W_computed"] = best["W"]
    meta["all_solutions"] = valid_solutions

    return True


def solve_PW_air_spaced(candidate: Candidate, cfg: Config) -> bool:
    """
    Solve for shape factors in air-spaced doublet.

    For air-spaced doublet, both P=P0 and W=W0 can be solved simultaneously
    since we have two independent shape factors X1 and X2.

    Reference: Nguyen 2022 Section 2.2

    Args:
        candidate: Candidate with meta dict containing phi1, phi2, n1, n2.
        cfg: Configuration with P0, W0.

    Returns:
        True if valid solution exists, False otherwise.
    """
    meta = getattr(candidate, "meta", None)
    if meta is None:
        return False

    n1 = meta["n1"]
    n2 = meta["n2"]
    phi1 = meta["phi1"]
    phi2 = meta["phi2"]

    # Get aberration coefficients
    a1, b1, c1 = thin_lens_P_coefficients(n1)
    a2, b2, c2 = thin_lens_P_coefficients(n2)
    w1 = thin_lens_W_coefficient(n1)
    w2 = thin_lens_W_coefficient(n2)

    phi1_2 = phi1 * phi1
    phi2_2 = phi2 * phi2
    phi1_3 = phi1 * phi1 * phi1
    phi2_3 = phi2 * phi2 * phi2

    # For air-spaced doublet, X1 and X2 are independent
    # We need to solve the system:
    # P(X1, X2) = φ1³*(a1*X1² + b1*X1 + c1) + φ2³*(a2*X2² + b2*X2 + c2) = P0
    # W(X1, X2) = φ1²*(w1*X1 + 1) + φ2²*(w2*X2 + 1) = W0

    # From coma equation, solve for X2 in terms of X1:
    # φ2²*w2*X2 = W0 - φ1²*(w1*X1 + 1) - φ2²
    # X2 = (W0 - φ1² - φ2² - φ1²*w1*X1) / (φ2²*w2)

    if abs(phi2_2 * w2) < 1e-15:
        return False

    # X2 = A_x2 * X1 + B_x2
    A_x2 = -phi1_2 * w1 / (phi2_2 * w2)
    B_x2 = (cfg.W0 - phi1_2 - phi2_2) / (phi2_2 * w2)

    # Substitute into P equation:
    # φ1³*(a1*X1² + b1*X1 + c1) + φ2³*(a2*(A_x2*X1+B_x2)² + b2*(A_x2*X1+B_x2) + c2) = P0

    # Expand:
    # φ1³*a1*X1² + φ2³*a2*A_x2²*X1² + ... = P0

    A_coeff = phi1_3 * a1 + phi2_3 * a2 * A_x2 * A_x2
    B_coeff = phi1_3 * b1 + phi2_3 * (2.0 * a2 * A_x2 * B_x2 + b2 * A_x2)
    C_coeff = phi1_3 * c1 + phi2_3 * (a2 * B_x2 * B_x2 + b2 * B_x2 + c2) - cfg.P0

    discriminant = B_coeff * B_coeff - 4.0 * A_coeff * C_coeff

    if discriminant < 0:
        return False

    meta["A_x2"] = A_x2
    meta["B_x2"] = B_x2
    meta["discriminant"] = discriminant

    # Solve for X1
    sqrt_disc = math.sqrt(discriminant)
    if abs(A_coeff) > 1e-15:
        X1_solutions = [
            (-B_coeff + sqrt_disc) / (2.0 * A_coeff),
            (-B_coeff - sqrt_disc) / (2.0 * A_coeff),
        ]
    else:
        if abs(B_coeff) > 1e-15:
            X1_solutions = [-C_coeff / B_coeff]
        else:
            return False

    # Compute X2 for each X1 and store
    valid_solutions = []
    for X1 in X1_solutions:
        X2 = A_x2 * X1 + B_x2
        valid_solutions.append({"X1": X1, "X2": X2})

    if not valid_solutions:
        return False

    # Use first valid solution (both achieve P=P0, W=W0 by construction)
    meta["X1"] = valid_solutions[0]["X1"]
    meta["X2"] = valid_solutions[0]["X2"]
    meta["all_solutions"] = valid_solutions

    return True


def solve_PW(candidate: Candidate, cfg: Config) -> bool:
    """
    Solve the spherical aberration and coma constraints.

    Dispatches to appropriate solver based on system type.

    Args:
        candidate: The Candidate to solve.
        cfg: Configuration.

    Returns:
        True if constraints can be satisfied, False otherwise.
    """
    if cfg.system_type == "cemented":
        return solve_PW_cemented(candidate, cfg)
    else:
        return solve_PW_air_spaced(candidate, cfg)


# =============================================================================
# Fill Candidate with Optical Values (fill_candidate)
# =============================================================================


def shape_factor_to_radii(X: float, phi: float, n: float) -> tuple[float, float]:
    """
    Convert shape factor X and power φ to surface radii R1, R2.

    Shape factor: X = (R2 + R1) / (R2 - R1)
    Power: φ = (n-1) * (1/R1 - 1/R2)

    Solving:
    Let c1 = 1/R1, c2 = 1/R2 (curvatures)
    φ = (n-1)*(c1 - c2)
    X = (c1 + c2) / (c1 - c2)  (note: X = -1/c_sum * c_diff relation)

    Actually: X = (R2+R1)/(R2-R1) = (1/c2 + 1/c1)/(1/c2 - 1/c1) = (c1+c2)/(c1-c2)

    Let S = c1 - c2 = φ/(n-1)
    Let D = c1 + c2 = X * S

    c1 = (S + D)/2 = (1 + X)*φ/(2*(n-1))
    c2 = (D - S)/2 = (X - 1)*φ/(2*(n-1))

    Args:
        X: Shape factor.
        phi: Lens power.
        n: Refractive index.

    Returns:
        Tuple (R1, R2) of surface radii. Infinite radii returned as ±1e10.
    """
    nm1 = n - 1.0

    if abs(nm1) < 1e-15:
        return (1e10, 1e10)

    c1 = (1.0 + X) * phi / (2.0 * nm1)
    c2 = (X - 1.0) * phi / (2.0 * nm1)

    # Convert curvatures to radii
    R1 = 1.0 / c1 if abs(c1) > 1e-15 else 1e10 * (1 if c1 >= 0 else -1)
    R2 = 1.0 / c2 if abs(c2) > 1e-15 else 1e10 * (1 if c2 >= 0 else -1)

    return R1, R2


def fill_candidate_cemented(candidate: Candidate, cfg: Config) -> Candidate:
    """
    Fill cemented doublet candidate with computed optical values.

    Uses alpha-parameter surface invariants:
    - P2: Spherical aberration surface invariant of the CEMENTED interface (surface 2)
    - W: Total coma (sum of W1 + W2 + W3)
    - R2: Radius of the cemented surface

    Surface layout for cemented doublet (3 surfaces):
        Surface 1: air -> glass1  (front of lens 1)
        Surface 2: glass1 -> glass2  (cemented interface)
        Surface 3: glass2 -> air  (back of lens 2)

    Reference: Nguyen 2022 Section 2.1, PE definition

    Args:
        candidate: Candidate with meta dict containing solution.
        cfg: Configuration.

    Returns:
        Updated Candidate.
    """
    meta = getattr(candidate, "meta", {})

    n1 = meta.get("n1", 1.5)
    n2 = meta.get("n2", 1.6)
    phi1 = meta.get("phi1", 0.01)
    phi2 = meta.get("phi2", -0.01)
    X1 = meta.get("X1", 0.0)
    X2 = meta.get("X2", 0.0)

    # Compute radii from shape factors
    R1_1, R1_2 = shape_factor_to_radii(X1, phi1, n1)  # First lens: front, back
    R2_1, R2_2 = shape_factor_to_radii(X2, phi2, n2)  # Second lens: front, back

    # For cemented doublet:
    # - R1 = front surface of lens 1
    # - R2 = cemented surface (back of lens 1 = front of lens 2)
    # - R3 = back surface of lens 2
    # Note: For a cemented doublet, R1_2 should equal R2_1 by construction.
    # We use R2_1 as the cemented surface radius.
    R1 = R1_1  # Front surface
    R2 = R2_1  # Cemented surface (use lens2 front = cemented interface)
    R3 = R2_2  # Back surface

    # Media indices for cemented doublet: [air, glass1, glass2, air]
    n_list = [1.0, n1, n2, 1.0]
    R_list = [R1, R2, R3]

    # Compute alphas using backward recursion
    fprime = cfg.f
    alphas = compute_alphas_from_radii(n_list, R_list, fprime, alpha_last=1.0)

    # Compute surface invariants (P, W) for each surface
    # Surface 1: between n_list[0]=1.0 and n_list[1]=n1
    P1, W1 = surface_PW(n_list[0], n_list[1], alphas[0], alphas[1])
    # Surface 2: between n_list[1]=n1 and n_list[2]=n2 (CEMENTED)
    P2_surf, W2 = surface_PW(n_list[1], n_list[2], alphas[1], alphas[2])
    # Surface 3: between n_list[2]=n2 and n_list[3]=1.0
    P3, W3 = surface_PW(n_list[2], n_list[3], alphas[2], alphas[3])

    # Store computed values in meta
    meta["R1_front"] = R1
    meta["R_cemented"] = R2
    meta["R2_back"] = R3
    meta["radii"] = [R1, R2, R3]
    meta["n_list"] = n_list
    meta["alphas"] = alphas
    meta["P_surfaces"] = [P1, P2_surf, P3]
    meta["W_surfaces"] = [W1, W2, W3]

    # Fill candidate fields for PE calculation
    # P2 = spherical aberration of cemented interface (surface 2)
    candidate.P2 = P2_surf

    # W = total system coma (sum over all surfaces)
    candidate.W = W1 + W2 + W3

    # R2 = radius of cemented surface
    candidate.R2 = R2 if abs(R2) < 1e9 else 1e3  # Cap for numerical stability

    # Ensure R2 is not zero (avoid division by zero in PE)
    if abs(candidate.R2) < 1e-6:
        candidate.R2 = 1e-3

    candidate.Pi = []

    return candidate


def fill_candidate_air_spaced(candidate: Candidate, cfg: Config) -> Candidate:
    """
    Fill air-spaced doublet candidate with computed optical values.

    Uses alpha-parameter surface invariants:
    - Pi: [P1, P2, P3, P4] spherical aberration surface invariants for each surface

    Surface layout for air-spaced doublet (4 surfaces):
        Surface 1: air -> glass1  (front of lens 1)
        Surface 2: glass1 -> air  (back of lens 1)
        Surface 3: air -> glass2  (front of lens 2)
        Surface 4: glass2 -> air  (back of lens 2)

    Reference: Nguyen 2022 Section 2.2, PE definition

    Args:
        candidate: Candidate with meta dict containing solution.
        cfg: Configuration.

    Returns:
        Updated Candidate.
    """
    meta = getattr(candidate, "meta", {})

    n1 = meta.get("n1", 1.5)
    n2 = meta.get("n2", 1.6)
    phi1 = meta.get("phi1", 0.01)
    phi2 = meta.get("phi2", -0.01)
    X1 = meta.get("X1", 0.0)
    X2 = meta.get("X2", 0.0)

    # Compute radii from shape factors
    R1_1, R1_2 = shape_factor_to_radii(X1, phi1, n1)  # First lens: front, back
    R2_1, R2_2 = shape_factor_to_radii(X2, phi2, n2)  # Second lens: front, back

    # Surface radii for air-spaced doublet
    R1 = R1_1  # Front of lens 1
    R2 = R1_2  # Back of lens 1
    R3 = R2_1  # Front of lens 2
    R4 = R2_2  # Back of lens 2

    # Media indices for air-spaced doublet: [air, glass1, air, glass2, air]
    n_list = [1.0, n1, 1.0, n2, 1.0]
    R_list = [R1, R2, R3, R4]

    # Compute alphas using backward recursion
    fprime = cfg.f
    alphas = compute_alphas_from_radii(n_list, R_list, fprime, alpha_last=1.0)

    # Compute surface invariants (P, W) for each surface
    # Surface 1: between n_list[0]=1.0 and n_list[1]=n1
    P1, W1 = surface_PW(n_list[0], n_list[1], alphas[0], alphas[1])
    # Surface 2: between n_list[1]=n1 and n_list[2]=1.0
    P2, W2 = surface_PW(n_list[1], n_list[2], alphas[1], alphas[2])
    # Surface 3: between n_list[2]=1.0 and n_list[3]=n2
    P3, W3 = surface_PW(n_list[2], n_list[3], alphas[2], alphas[3])
    # Surface 4: between n_list[3]=n2 and n_list[4]=1.0
    P4, W4 = surface_PW(n_list[3], n_list[4], alphas[3], alphas[4])

    # Store computed values in meta
    meta["R1_front"] = R1
    meta["R1_back"] = R2
    meta["R2_front"] = R3
    meta["R2_back"] = R4
    meta["radii"] = [R1, R2, R3, R4]
    meta["n_list"] = n_list
    meta["alphas"] = alphas
    meta["Pi_surfaces"] = [P1, P2, P3, P4]
    meta["Wi_surfaces"] = [W1, W2, W3, W4]

    # Fill candidate fields for PE calculation
    candidate.Pi = [P1, P2, P3, P4]
    candidate.P2 = None
    candidate.W = None
    candidate.R2 = None

    return candidate


def fill_candidate(candidate: Candidate, cfg: Config) -> Candidate:
    """
    Fill candidate with computed optical values needed for PE.

    Dispatches to appropriate filler based on system type.

    Args:
        candidate: Candidate with meta dict containing solution.
        cfg: Configuration.

    Returns:
        Updated Candidate.
    """
    if cfg.system_type == "cemented":
        return fill_candidate_cemented(candidate, cfg)
    else:
        return fill_candidate_air_spaced(candidate, cfg)


# =============================================================================
# Main Pipeline
# =============================================================================


def run(cfg: Config) -> tuple[list[Candidate], dict]:
    """
    Run the 2022 synthesis pipeline.

    Flow:
    1. Load glasses from AGF catalogs
    2. Filter glasses
    3. Initialize TopNContainer
    4. For each glass pair:
       - Check delta-nu constraint
       - Solve C = C0 (achromatization)
       - Solve P, W constraints
       - Fill candidate values
       - Compute PE
       - Add to top-N if PE <= max_PE
    5. Return sorted results and statistics

    Args:
        cfg: Configuration object.

    Returns:
        Tuple of (sorted_candidates, statistics_dict).
    """
    # Statistics
    stats: dict = {
        "glasses_total": 0,
        "glasses_filtered": 0,
        "pairs_tested": 0,
        "pairs_achromatized": 0,
        "pairs_PW_solved": 0,
        "candidates_generated": 0,
        "candidates_saved": 0,
        "best_PE": float("inf"),
    }

    # Step 1: Load glasses
    glasses = load_catalog(cfg.agf_paths)
    stats["glasses_total"] = len(glasses)

    # Step 2: Filter glasses
    filtered_glasses = filter_glasses(glasses, cfg)
    stats["glasses_filtered"] = len(filtered_glasses)

    # Step 3: Initialize top-N container
    top = TopNContainer(cfg.N_keep)

    # Step 4: Enumerate and process pairs
    for g1, g2 in enumerate_pairs(filtered_glasses, cfg.allow_repeat):
        stats["pairs_tested"] += 1

        # Solve C = C0 (achromatization)
        # This also checks delta-nu internally
        candidates = solve_C(g1, g2, cfg)

        if not candidates:
            continue

        stats["pairs_achromatized"] += 1

        for candidate in candidates:
            # Solve P, W constraints
            ok = solve_PW(candidate, cfg)
            if not ok:
                continue

            stats["pairs_PW_solved"] += 1

            # Fill candidate with optical values
            candidate = fill_candidate(candidate, cfg)

            # Compute PE
            pe = compute_pe(candidate)
            stats["candidates_generated"] += 1

            # Check PE threshold
            if candidate.PE is not None and candidate.PE > cfg.max_PE:
                continue

            # Add to top-N
            if top.push(candidate):
                stats["candidates_saved"] = len(top)

    # Step 5: Get sorted results
    results = top.sorted()

    # Update final statistics
    stats["candidates_saved"] = len(results)
    if results:
        best_pe = rank_key(results[0])
        if best_pe != float("inf"):
            stats["best_PE"] = best_pe

    return results, stats


# =============================================================================
# Entry Point
# =============================================================================

if __name__ == "__main__":
    # Run with default configuration
    cfg = default_cfg()

    print("Running 2022 synthesis pipeline...")
    print(f"  System type: {cfg.system_type}")
    print(f"  Catalogs: {cfg.agf_paths}")
    print(f"  N_keep: {cfg.N_keep}")
    print()

    # Run pipeline
    results, stats = run(cfg)

    # Ensure output directory exists
    io.ensure_out_dir(cfg.out_dir)

    # Save results
    csv_path = f"{cfg.out_dir}/results_2022.csv"
    io.save_candidates_csv(results, csv_path)
    print(f"Results saved to: {csv_path}")
    print()

    # Print summary
    io.print_summary(stats)

    # Sanity print: show top candidate's optical values
    if results:
        top_cand = results[0]
        g1_obj = getattr(top_cand, "g1", None)
        g2_obj = getattr(top_cand, "g2", None)
        g1_name = g1_obj.name if g1_obj is not None and hasattr(g1_obj, "name") else "?"
        g2_name = g2_obj.name if g2_obj is not None and hasattr(g2_obj, "name") else "?"

        print()
        print("=" * 50)
        print("Top Candidate Optical Values (Sanity Check)")
        print("=" * 50)
        print(f"  Glass pair: {g1_name} + {g2_name}")
        print(f"  PE: {top_cand.PE}")

        if cfg.system_type == "cemented":
            print(f"  P2 (cemented surface invariant): {top_cand.P2}")
            print(f"  W (system coma): {top_cand.W}")
            print(f"  R2 (cemented radius, mm): {top_cand.R2}")
            meta = getattr(top_cand, "meta", {})
            if "P_surfaces" in meta:
                print(f"  P_surfaces: {meta['P_surfaces']}")
            if "W_surfaces" in meta:
                print(f"  W_surfaces: {meta['W_surfaces']}")
            if "alphas" in meta:
                print(f"  alphas: {meta['alphas']}")
            if "radii" in meta:
                print(f"  radii (mm): {meta['radii']}")
        else:
            print(f"  Pi (surface invariants): {top_cand.Pi}")
            meta = getattr(top_cand, "meta", {})
            if "alphas" in meta:
                print(f"  alphas: {meta['alphas']}")
            if "radii" in meta:
                print(f"  radii (mm): {meta['radii']}")

        print("=" * 50)
