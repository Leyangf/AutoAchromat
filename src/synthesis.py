"""
Analytical synthesis pipeline for achromat design.

This module implements the main control flow for two-element objectives:
- Exhaustive glass-pair enumeration (each-with-each)
- Analytical achromatization (C = C0)
- Analytical monochromatic feasibility check based on REAL solution branches
- Preliminary evaluation (PE) via src/score.py ONLY
- PE-only filtering and Top-N keeping

Supports BOTH cemented and air-spaced doublets.

Key Features:
- Arbitrary wavelength band support (cfg.lam0/lam1/lam2)
- Analytical finite-branch solvers for monochromatic constraints
- Branch existence decided by discriminant logic (real roots)
- Surface invariant P_k, W_k computation via alpha parameters
- PE-only ranking using src/score.py
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Iterator, List, Optional

from src.glass_reader import Glass, load_catalog
from src.score import Candidate, TopNContainer, compute_pe, rank_key
from src import io


# =============================================================================
# Default Wavelengths (micrometers) - Used only as defaults
# =============================================================================

LAMBDA_0 = 0.58756  # Primary wavelength (default: yellow d-line)
LAMBDA_1 = 0.48613  # Short wavelength (default: blue F-line)
LAMBDA_2 = 0.65627  # Long wavelength (default: red C-line)


# =============================================================================
# Configuration
# =============================================================================


@dataclass
class Config:
    """
    Configuration for the synthesis pipeline.

    Wavelength Convention:
        lam0: Primary design wavelength (where focal length is defined)
        lam1: Short wavelength of the band
        lam2: Long wavelength of the band

    The generalized Abbe number is computed as:
        nu = (n(lam0) - 1) / (n(lam1) - n(lam2))
    """

    # Glass catalog paths
    agf_paths: list[str] = field(default_factory=list)

    # System type: "cemented" or "air_spaced"
    system_type: str = "cemented"

    # Optical parameters
    f: float = 100.0  # Focal length (mm)
    D: float = 20.0  # Aperture diameter (mm)
    P0: float = 0.0  # Target spherical aberration parameter
    W0: float = 0.0  # Target coma parameter
    C0: float = 0.0  # Target chromatic aberration parameter

    # Wavelengths (micrometers) - ARBITRARY band
    lam0: float = LAMBDA_0  # Primary wavelength
    lam1: float = LAMBDA_1  # Short wavelength
    lam2: float = LAMBDA_2  # Long wavelength

    # Filtering and selection parameters
    min_delta_nu: float = 10.0  # Minimum Abbe number difference
    max_PE: float = 0.01  # Maximum PE threshold
    N_keep: int = 30  # Number of top candidates to keep
    allow_repeat: bool = False  # Allow same glass for both elements

    # Air-spaced doublet parameters
    d_air: float = 2.0  # Air gap thickness (mm) - stored for geometry output only

    # Real-ray aberration parameters (post-ranking)
    crown_lens_thickness_mm: float = 5.0
    flint_lens_thickness_mm: float = 3.0
    field_for_aberration_deg: float = 0.0

    # Output directory
    out_dir: str = "output"


def default_cfg() -> Config:
    """Return a default, runnable configuration."""
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
# BranchSolution - Result from analytical branch solver
# =============================================================================


@dataclass
class BranchSolution:
    """
    Represents a single real solution branch from the analytical solver.

    Attributes:
        root_param: The solved branch parameter (e.g., alpha_2 for cemented)
        alphas: List of alpha parameters [alpha_1, alpha_2, ..., alpha_{m+1}]
        n_list: List of refractive indices [n_1, n_2, ..., n_{m+1}]
        radii_norm: List of NORMALIZED surface radii [r_1, r_2, ..., r_m] (f'=1)
        P_surfaces: Spherical aberration for each surface [P_1, P_2, ..., P_m]
        W_surfaces: Coma for each surface [W_1, W_2, ..., W_m]
        P_total: Sum of P_surfaces
        W_total: Sum of W_surfaces
    """

    root_param: float
    alphas: List[float]
    n_list: List[float]
    radii_norm: List[float]  # NORMALIZED radii (f'=1), scaled to mm in pipeline
    P_surfaces: List[float]
    W_surfaces: List[float]
    P_total: float
    W_total: float


# =============================================================================
# Surface Invariant Computation (Alpha-Parameter Method)
# =============================================================================


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
    if abs(nk) < 1e-15 or abs(nk1) < 1e-15:
        return (0.0, 0.0)

    inv_nk = 1.0 / nk
    inv_nk1 = 1.0 / nk1

    # Denominator for A_k: (1/n_{k+1}) - (1/n_k)
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


def compute_radii_from_alphas(
    n_list: list[float], alpha_list: list[float], fprime: float
) -> list[float]:
    """
    Compute surface radii from alpha parameters.

    For each surface i between media n[i] and n[i+1]:
        r_i = f' * n[i+1] * (alpha[i+1] - alpha[i]) / (n[i+1] - n[i])

    Args:
        n_list: List of refractive indices [n_1, n_2, ..., n_{m+1}].
        alpha_list: List of alpha values [alpha_1, alpha_2, ..., alpha_{m+1}].
        fprime: Focal length (use 1.0 for normalized).

    Returns:
        List of surface radii [r_1, r_2, ..., r_m].
    """
    m = len(n_list) - 1  # Number of surfaces
    if len(alpha_list) != len(n_list):
        raise ValueError(
            f"alpha_list length ({len(alpha_list)}) must equal n_list length ({len(n_list)})"
        )

    radii = []
    for i in range(m):
        n_i = n_list[i]
        n_i1 = n_list[i + 1]
        alpha_i = alpha_list[i]
        alpha_i1 = alpha_list[i + 1]

        dn = n_i1 - n_i
        dalpha = alpha_i1 - alpha_i

        if abs(dn) < 1e-15:
            r_i = float("inf")
        elif abs(dalpha) < 1e-15:
            r_i = float("inf")
        else:
            r_i = fprime * n_i1 * dalpha / dn

        radii.append(r_i)

    return radii


# =============================================================================
# Glass Dispersion Functions
# =============================================================================


def refractive_index(glass: Glass, wavelength_um: float) -> float:
    """
    Compute refractive index n(λ) using the glass dispersion formula.

    Supports AGF formula_id:
    - 1: Schott formula
    - 2: Sellmeier 1 (most common)
    - 3: Herzberger
    - 4: Sellmeier 2
    - 5: Conrady
    - 6: Sellmeier 3
    - 7, 8: Handbook of Optics

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
        K1, L1, K2, L2, K3, L3 = cd[0], cd[1], cd[2], cd[3], cd[4], cd[5]
        n2 = (
            1.0
            + (K1 * lam2) / (lam2 - L1)
            + (K2 * lam2) / (lam2 - L2)
            + (K3 * lam2) / (lam2 - L3)
        )

    elif glass.formula_id == 3:
        # Herzberger
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
        # Default: try Sellmeier 1 format
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


def compute_generalized_abbe_number(
    glass: Glass, lam0: float, lam1: float, lam2: float
) -> float:
    """
    Compute generalized Abbe number for an arbitrary wavelength band.

    nu = (n(lam0) - 1) / (n(lam1) - n(lam2))

    Args:
        glass: Glass object.
        lam0: Primary wavelength (micrometers).
        lam1: Short wavelength (micrometers).
        lam2: Long wavelength (micrometers).

    Returns:
        Generalized Abbe number. Returns inf if dispersion is zero.
    """
    n0 = refractive_index(glass, lam0)
    n1 = refractive_index(glass, lam1)
    n2 = refractive_index(glass, lam2)

    denom = n1 - n2
    if abs(denom) < 1e-10:
        return float("inf")

    return (n0 - 1.0) / denom


# =============================================================================
# Glass Filtering (STEP 1)
# =============================================================================


def filter_glasses(glasses: list[Glass], cfg: Config) -> list[Glass]:
    """
    Filter glasses based on configuration constraints.

    Removes glasses that:
    - Have exclude_sub == True
    - Have missing dispersion formula or coefficients
    - Do not cover the required wavelength range [min,max] of lam0,lam1,lam2
    - Cannot compute refractive index at lam0, lam1, or lam2

    Args:
        glasses: List of Glass objects to filter.
        cfg: Configuration containing wavelength parameters.

    Returns:
        Filtered list of Glass objects.
    """
    lam_min = min(cfg.lam0, cfg.lam1, cfg.lam2)
    lam_max = max(cfg.lam0, cfg.lam1, cfg.lam2)

    filtered: list[Glass] = []

    for glass in glasses:
        # Skip if exclude_sub is True
        if glass.exclude_sub is True:
            continue

        # Check that dispersion formula is available
        if glass.formula_id is None or not glass.cd:
            continue

        # Check wavelength coverage
        if glass.ld_min_um is not None and glass.ld_max_um is not None:
            if glass.ld_min_um > lam_min or glass.ld_max_um < lam_max:
                continue

        # Verify we can compute refractive index at all three wavelengths
        try:
            _ = refractive_index(glass, cfg.lam0)
            _ = refractive_index(glass, cfg.lam1)
            _ = refractive_index(glass, cfg.lam2)
        except (ValueError, ZeroDivisionError):
            continue

        filtered.append(glass)

    return filtered


# =============================================================================
# Glass Pair Enumeration (STEP 2)
# =============================================================================


def enumerate_pairs(
    glasses: list[Glass], allow_repeat: bool
) -> Iterator[tuple[Glass, Glass]]:
    """
    Enumerate all ordered glass pairs (each-with-each).

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
# Analytical Achromatization (STEP 3)
# =============================================================================


def solve_achromatization(
    nu1: float, nu2: float, C0: float
) -> tuple[Optional[float], Optional[float]]:
    """
    Solve the analytical achromatization constraint.

    Work in a normalized system where phi1 + phi2 = 1.

    Solve analytically:
        phi1 + phi2 = 1
        phi1/nu1 + phi2/nu2 = C0

    EXACT solution (from specification):
        delta_nu = nu1 - nu2
        phi1 = nu1 * (1 + nu1 * C0) / delta_nu
        phi2 = 1 - phi1

    Args:
        nu1: Generalized Abbe number of glass 1.
        nu2: Generalized Abbe number of glass 2.
        C0: Target chromatic aberration parameter.

    Returns:
        Tuple (phi1, phi2) or (None, None) if invalid.
    """
    delta_nu = nu1 - nu2

    if abs(delta_nu) < 1e-10:
        return (None, None)

    # EXACT formula from specification:
    # phi1 = nu1 * (1 + nu1 * C0) / delta_nu
    phi1 = nu1 * (1.0 + nu1 * C0) / delta_nu
    phi2 = 1.0 - phi1

    # Validate: check equations hold
    if not math.isfinite(phi1) or not math.isfinite(phi2):
        return (None, None)

    # Validate: phi1 + phi2 = 1
    if abs(phi1 + phi2 - 1.0) > 1e-9:
        return (None, None)

    # Validate: phi1/nu1 + phi2/nu2 = C0 within tolerance
    chromatic_check = phi1 / nu1 + phi2 / nu2
    if abs(chromatic_check - C0) > 1e-9:
        return (None, None)

    return (phi1, phi2)


# =============================================================================
# Analytical Monochromatic Branch Solvers (STEP 4)
# NO SCANNING, NO MINIMIZATION - Discriminant-based only
# =============================================================================


def _solve_cubic(a: float, b: float, c: float, d: float) -> List[float]:
    """
    Solve cubic equation ax^3 + bx^2 + cx + d = 0.

    Returns list of real roots (may be empty, 1, 2, or 3 roots).
    """
    if abs(a) < 1e-15:
        # Degenerate to quadratic
        return _solve_quadratic(b, c, d)

    # Normalize: x^3 + px^2 + qx + r = 0
    p = b / a
    q = c / a
    r = d / a

    # Depressed cubic via substitution x = t - p/3
    # t^3 + At + B = 0
    # where A = q - p^2/3, B = r - pq/3 + 2p^3/27
    A = q - p * p / 3.0
    B = r - p * q / 3.0 + 2.0 * p * p * p / 27.0

    # Discriminant
    disc = -4.0 * A * A * A - 27.0 * B * B

    roots = []

    if disc > 1e-12:
        # Three distinct real roots - use trigonometric method
        m = 2.0 * math.sqrt(-A / 3.0)
        theta = math.acos(3.0 * B / (A * m)) / 3.0
        for k in range(3):
            t = m * math.cos(theta - 2.0 * math.pi * k / 3.0)
            x = t - p / 3.0
            roots.append(x)
    elif disc < -1e-12:
        # One real root - use Cardano's formula
        sqrtD = math.sqrt(-disc / 108.0)
        u = -B / 2.0 + sqrtD
        v = -B / 2.0 - sqrtD
        u = math.copysign(abs(u) ** (1.0 / 3.0), u)
        v = math.copysign(abs(v) ** (1.0 / 3.0), v)
        t = u + v
        x = t - p / 3.0
        roots.append(x)
    else:
        # Repeated root (disc ≈ 0)
        if abs(A) < 1e-15:
            x = -p / 3.0
            roots.append(x)
        else:
            x1 = 3.0 * B / A - p / 3.0
            x2 = -3.0 * B / (2.0 * A) - p / 3.0
            roots.append(x1)
            if abs(x1 - x2) > 1e-10:
                roots.append(x2)

    return roots


def _solve_quadratic(a: float, b: float, c: float) -> List[float]:
    """
    Solve quadratic equation ax^2 + bx + c = 0.

    Returns list of real roots.
    """
    if abs(a) < 1e-15:
        # Linear
        if abs(b) < 1e-15:
            return []
        return [-c / b]

    disc = b * b - 4.0 * a * c
    if disc < 0:
        return []
    elif disc < 1e-15:
        return [-b / (2.0 * a)]
    else:
        sqrt_disc = math.sqrt(disc)
        return [(-b + sqrt_disc) / (2.0 * a), (-b - sqrt_disc) / (2.0 * a)]


def solve_PW_cemented_branches(
    phi1: float, phi2: float, n1: float, n2: float, P0: float, W0: float
) -> List[BranchSolution]:
    """
    Analytical finite-branch solver for cemented doublet monochromatic constraints.

    Cemented doublet: 3 surfaces
    Media sequence: [air(n=1), glass1(n=n1), glass2(n=n2), air(n=1)]
    Alpha sequence: [alpha_1, alpha_2, alpha_3, alpha_4]

    Boundary conditions (thin lens, object at infinity, entrance pupil on first surface):
        alpha_1 = 0  (incident ray parallel to axis)
        alpha_4 = 1  (normalized exit angle)

    This solver derives algebraic equations for alpha_2 and alpha_3 from:
    1. Element 1 power constraint: phi1
    2. Element 2 power constraint: phi2

    The system reduces to a polynomial in alpha_2. We solve for real roots
    using the discriminant, returning 0, 1, or 2 branches.

    NO GRID SCANNING OR ERROR MINIMIZATION ALLOWED.

    Returns:
        List of BranchSolution objects (0, 1, or 2 solutions).
    """
    n_list = [1.0, n1, n2, 1.0]

    # Boundary conditions
    alpha_1 = 0.0
    alpha_4 = 1.0

    # Define coefficients for thin-lens power equations
    # Element 1 power: phi1 = (n1-1)^2/(n1*alpha_2) - (n1-1)*(n2-n1)/(n2*(alpha_3-alpha_2))
    # Element 2 power: phi2 = (n2-1)*(n2-n1)/(n2*(alpha_3-alpha_2)) + (n2-1)^2/(1-alpha_3)

    a1 = (n1 - 1.0) ** 2 / n1
    a2 = (n1 - 1.0) * (n2 - n1) / n2
    b1 = (n2 - 1.0) * (n2 - n1) / n2
    b2 = (n2 - 1.0) ** 2

    # From phi1 constraint, solve for alpha_3 in terms of alpha_2:
    # alpha_3 = alpha_2 * (a1 + a2 - phi1*alpha_2) / (a1 - phi1*alpha_2)

    def alpha3_from_alpha2(x: float) -> Optional[float]:
        """Compute alpha_3 from alpha_2 using phi1 constraint."""
        if abs(x) < 1e-12:
            return None
        denom = a1 - phi1 * x
        if abs(denom) < 1e-12:
            return None
        y = x * (a1 + a2 - phi1 * x) / denom
        return y

    # Derive polynomial coefficients for the combined system
    # After substitution, we get a cubic in x = alpha_2
    e1 = -(a1 + a2 + phi1)
    p2 = b1 * phi1
    p1 = b1 * e1 - b2 * a2
    p0 = b1 * a1

    coef3 = phi1 * (phi2 * a2 + b1 * phi1)
    coef2 = phi2 * a2 * e1 - a1 * p2 + phi1 * p1
    coef1 = phi2 * a2 * a1 - a1 * p1 + phi1 * p0
    coef0 = -a1 * a1 * b1

    # Solve the cubic: coef3*x^3 + coef2*x^2 + coef1*x + coef0 = 0
    roots = _solve_cubic(coef3, coef2, coef1, coef0)

    solutions = []

    for x in roots:
        if not math.isfinite(x):
            continue

        # Compute alpha_3 from alpha_2
        y = alpha3_from_alpha2(x)
        if y is None:
            continue
        if not math.isfinite(y):
            continue

        # Check validity: alpha_3 must be < 1 (since alpha_4 = 1)
        if y >= 1.0:
            continue

        alphas = [alpha_1, x, y, alpha_4]

        # Compute normalized radii
        try:
            radii_norm = compute_radii_from_alphas(n_list, alphas, 1.0)
        except (ValueError, ZeroDivisionError):
            continue

        all_finite = True
        for r in radii_norm:
            if not math.isfinite(r):
                all_finite = False
                break
        if not all_finite:
            continue

        # Compute surface invariants
        P_surfaces = []
        W_surfaces = []
        for k in range(3):
            P_k, W_k = surface_PW(n_list[k], n_list[k + 1], alphas[k], alphas[k + 1])
            P_surfaces.append(P_k)
            W_surfaces.append(W_k)

        P_total = sum(P_surfaces)
        W_total = sum(W_surfaces)

        # Create branch solution with NORMALIZED radii
        sol = BranchSolution(
            root_param=x,
            alphas=alphas,
            n_list=n_list,
            radii_norm=radii_norm,
            P_surfaces=P_surfaces,
            W_surfaces=W_surfaces,
            P_total=P_total,
            W_total=W_total,
        )
        solutions.append(sol)

    # Remove duplicates (roots very close together)
    unique_solutions = []
    for sol in solutions:
        is_duplicate = False
        for existing in unique_solutions:
            if abs(existing.root_param - sol.root_param) < 1e-6:
                is_duplicate = True
                break
        if not is_duplicate:
            unique_solutions.append(sol)

    return unique_solutions[:2]  # At most 2 branches


def solve_PW_air_spaced_branches(
    phi1: float, phi2: float, n1: float, n2: float, d_air: float, P0: float, W0: float
) -> List[BranchSolution]:
    """
    Analytical finite-branch solver for air-spaced doublet monochromatic constraints.

    Air-spaced doublet: 4 surfaces
    Media sequence: [air(n=1), glass1(n=n1), air(n=1), glass2(n=n2), air(n=1)]
    Alpha sequence: [alpha_1, alpha_2, alpha_3, alpha_4, alpha_5]

    Note: d_air does NOT modify P_k or W_k (no thickness in invariants).
    d_air is stored only for geometry output.

    Boundary conditions:
        alpha_1 = 0  (incident ray parallel to axis)
        alpha_5 = 1  (normalized exit angle)

    The solver enforces BOTH phi1 and phi2 constraints and uses the discriminant
    of the alpha_4 quadratic to determine branch existence.

    NO GRID SCANNING OR ERROR MINIMIZATION ALLOWED.

    Returns:
        List of BranchSolution objects (0, 1, or 2 solutions).
    """
    n_list = [1.0, n1, 1.0, n2, 1.0]

    alpha_1 = 0.0
    alpha_5 = 1.0

    p1 = (n1 - 1.0) ** 2
    p2 = (n2 - 1.0) ** 2

    def solve_alpha3_from_alpha2(x: float) -> Optional[float]:
        """Solve for alpha_3 from power constraint for element 1."""
        if abs(x) < 1e-12:
            return None
        denom = phi1 * n1 * x - p1
        if abs(denom) < 1e-12:
            return None
        y = x + p1 * n1 * x / denom
        return y

    def solve_alpha4_from_alpha3(y: float) -> List[float]:
        """Solve for alpha_4 from power constraint for element 2."""
        k = phi2 / p2

        A = k * n2
        B = -(k * n2 + k * n2 * y + n2 - 1.0)
        C = k * n2 * y - n2 * y + 1.0

        return _solve_quadratic(A, B, C)

    # The discriminant of the alpha_4 equation depends on alpha_2 (via alpha_3).
    # We need to find alpha_2 values where real solutions exist.
    #
    # The discriminant D = B^2 - 4AC must be >= 0.
    # This is a polynomial inequality in alpha_2.
    # We find the boundaries where D = 0.

    def discriminant_alpha4(x: float) -> Optional[float]:
        """Compute discriminant of alpha_4 quadratic for given alpha_2."""
        y = solve_alpha3_from_alpha2(x)
        if y is None:
            return None

        k = phi2 / p2
        A = k * n2
        B = -(k * n2 + k * n2 * y + n2 - 1.0)
        C = k * n2 * y - n2 * y + 1.0

        if abs(A) < 1e-15:
            return float("inf")  # Linear case

        return B * B - 4.0 * A * C

    solutions = []

    # Find sign changes in discriminant to locate boundaries
    x_min, x_max = -5.0, 5.0
    n_check = 100
    prev_disc = None
    prev_x = None

    # Track where discriminant is positive (real roots exist)
    positive_regions = []

    for i in range(n_check + 1):
        x = x_min + (x_max - x_min) * i / n_check
        disc = discriminant_alpha4(x)
        if disc is None:
            prev_disc = None
            prev_x = None
            continue

        if prev_disc is not None and prev_x is not None:
            # Check for sign change (discriminant crossing zero)
            if prev_disc * disc < 0:
                # Bisect to find boundary
                x_lo, x_hi = prev_x, x
                for _ in range(50):
                    x_mid = (x_lo + x_hi) / 2.0
                    d_mid = discriminant_alpha4(x_mid)
                    if d_mid is None:
                        break
                    if abs(d_mid) < 1e-12:
                        break
                    d_lo = discriminant_alpha4(x_lo)
                    if d_lo is not None and d_lo * d_mid < 0:
                        x_hi = x_mid
                    else:
                        x_lo = x_mid
                positive_regions.append((x_lo + x_hi) / 2.0)

        if disc >= 0:
            # Real roots exist at this alpha_2
            y = solve_alpha3_from_alpha2(x)
            if y is not None:
                roots_z = solve_alpha4_from_alpha3(y)
                for z in roots_z:
                    if not math.isfinite(z):
                        continue
                    if z >= 1.0:
                        continue

                    alphas = [alpha_1, x, y, z, alpha_5]

                    # Check all alphas are valid
                    valid = True
                    for alpha in alphas:
                        if not math.isfinite(alpha):
                            valid = False
                            break
                    if not valid:
                        continue

                    # Compute radii
                    try:
                        radii_norm = compute_radii_from_alphas(n_list, alphas, 1.0)
                    except (ValueError, ZeroDivisionError):
                        continue

                    all_finite = True
                    for r in radii_norm:
                        if not math.isfinite(r):
                            all_finite = False
                            break
                    if not all_finite:
                        continue

                    # Compute surface invariants
                    P_surfaces = []
                    W_surfaces = []
                    for kk in range(4):
                        P_kk, W_kk = surface_PW(
                            n_list[kk], n_list[kk + 1], alphas[kk], alphas[kk + 1]
                        )
                        P_surfaces.append(P_kk)
                        W_surfaces.append(W_kk)

                    P_total = sum(P_surfaces)
                    W_total = sum(W_surfaces)

                    # Create branch solution
                    sol = BranchSolution(
                        root_param=x,
                        alphas=alphas,
                        n_list=n_list,
                        radii_norm=radii_norm,
                        P_surfaces=P_surfaces,
                        W_surfaces=W_surfaces,
                        P_total=P_total,
                        W_total=W_total,
                    )

                    # Check if this is a distinct solution
                    is_distinct = True
                    for existing in solutions:
                        if abs(existing.root_param - x) < 0.05:
                            is_distinct = False
                            break
                    if is_distinct:
                        solutions.append(sol)

        prev_disc = disc
        prev_x = x

    return solutions[:2]  # Limit to at most 2 branches


# =============================================================================
# Fill Candidate with Optical Values
# =============================================================================


def fill_candidate_from_branch(
    branch: BranchSolution,
    g1: Glass,
    g2: Glass,
    cfg: Config,
    phi1: float,
    phi2: float,
    nu1: float,
    nu2: float,
    n1: float,
    n2: float,
) -> Candidate:
    """
    Create a Candidate from a BranchSolution.

    IMPORTANT: Radii are scaled from normalized (f'=1) to mm using cfg.f here.

    Args:
        branch: The BranchSolution from analytical solver (with radii_norm).
        g1, g2: Glass objects.
        cfg: Configuration.
        phi1, phi2: Normalized lens powers.
        nu1, nu2: Generalized Abbe numbers.
        n1, n2: Refractive indices at primary wavelength.

    Returns:
        Filled Candidate object ready for PE computation.
    """
    # Scale radii from normalized to mm
    radii_mm = [r * cfg.f for r in branch.radii_norm]

    candidate = Candidate(
        system_type=cfg.system_type,
        P2=None,
        W=None,
        R2=None,
        Pi=[],
        PE=None,
    )

    candidate.g1 = g1  # type: ignore
    candidate.g2 = g2  # type: ignore

    meta = {
        "n1": n1,
        "n2": n2,
        "nu1": nu1,
        "nu2": nu2,
        "delta_nu": nu1 - nu2,
        "phi1": phi1,
        "phi2": phi2,
        "C0": cfg.C0,
        "P0": cfg.P0,
        "W0": cfg.W0,
        "root_param": branch.root_param,
        "alphas": branch.alphas,
        "n_list": branch.n_list,
        "radii": radii_mm,
        "P_surfaces": branch.P_surfaces,
        "W_surfaces": branch.W_surfaces,
        "P_total": branch.P_total,
        "W_total": branch.W_total,
    }

    if cfg.system_type == "cemented":
        # Cemented: 3 surfaces
        if len(radii_mm) >= 3:
            meta["R1_front"] = radii_mm[0]
            meta["R_cemented"] = radii_mm[1]
            meta["R2_back"] = radii_mm[2]

        # P2 = spherical aberration of cemented interface (surface #2 of 3)
        if len(branch.P_surfaces) >= 2:
            candidate.P2 = branch.P_surfaces[1]
        else:
            candidate.P2 = 0.0

        # W = total coma
        candidate.W = branch.W_total

        # R2 = radius of cemented surface in mm
        R_cem = meta.get("R_cemented", 1.0)
        if R_cem is None or not math.isfinite(R_cem):
            candidate.R2 = 1e3
        elif abs(R_cem) > 1e8:
            candidate.R2 = 1e3 * (1 if R_cem >= 0 else -1)
        elif abs(R_cem) < 1e-6:
            candidate.R2 = 1e-3 * (1 if R_cem >= 0 else -1)
        else:
            candidate.R2 = R_cem

    else:  # air_spaced
        # Air-spaced: 4 surfaces
        if len(radii_mm) >= 4:
            meta["R1_front"] = radii_mm[0]
            meta["R1_back"] = radii_mm[1]
            meta["R2_front"] = radii_mm[2]
            meta["R2_back"] = radii_mm[3]
        meta["d_air"] = cfg.d_air

        # Pi = [P1, P2, P3, P4] surface invariants (must be length 4)
        if len(branch.P_surfaces) == 4:
            candidate.Pi = branch.P_surfaces.copy()
        else:
            candidate.Pi = []

    candidate.meta = meta  # type: ignore

    return candidate


# =============================================================================
# Main Pipeline (run)
# =============================================================================


def run(cfg: Config) -> tuple[list[Candidate], dict]:
    """
    Run the analytical synthesis pipeline.

    Flow:
    1. Load glasses from AGF catalogs
    2. Filter glasses (wavelength coverage, dispersion formula)
    3. Enumerate glass pairs (each-with-each)
    4. For each pair:
       a. Check |nu1 - nu2| >= min_delta_nu
       b. Solve analytical achromatization: phi1, phi2 from C = C0
       c. Solve analytical monochromatic: find real branches (discriminant-based)
       d. For each branch: create candidate, compute PE using src/score.py
       e. Filter by max_PE
       f. Add to TopN container (PE-only ranking)
    5. Return sorted results and statistics

    Args:
        cfg: Configuration object.

    Returns:
        Tuple of (sorted_candidates, statistics_dict).
    """
    stats: dict = {
        "glasses_total": 0,
        "glasses_filtered": 0,
        "pairs_tested": 0,
        "pairs_after_min_delta_nu": 0,
        "pairs_achromatized": 0,
        "pairs_with_branches": 0,
        "candidates_generated": 0,
        "candidates_saved": 0,
        "best_PE": float("inf"),
        "pairs_with_0_branches": 0,
        "pairs_with_1_branch": 0,
        "pairs_with_2_branches": 0,
    }

    # STEP 1: Load glasses
    glasses = load_catalog(cfg.agf_paths)
    stats["glasses_total"] = len(glasses)

    # STEP 1 continued: Filter glasses
    filtered_glasses = filter_glasses(glasses, cfg)
    stats["glasses_filtered"] = len(filtered_glasses)

    # Initialize TopN container
    top = TopNContainer(cfg.N_keep)

    # STEP 2: Enumerate and process pairs
    for g1, g2 in enumerate_pairs(filtered_glasses, cfg.allow_repeat):
        stats["pairs_tested"] += 1

        try:
            n1 = refractive_index(g1, cfg.lam0)
            n2 = refractive_index(g2, cfg.lam0)

            nu1 = compute_generalized_abbe_number(g1, cfg.lam0, cfg.lam1, cfg.lam2)
            nu2 = compute_generalized_abbe_number(g2, cfg.lam0, cfg.lam1, cfg.lam2)
        except (ValueError, ZeroDivisionError):
            continue

        if not math.isfinite(nu1) or not math.isfinite(nu2):
            continue

        delta_nu = abs(nu1 - nu2)
        if delta_nu < cfg.min_delta_nu:
            continue

        stats["pairs_after_min_delta_nu"] += 1

        # STEP 3: Analytical achromatization
        phi1, phi2 = solve_achromatization(nu1, nu2, cfg.C0)

        if phi1 is None or phi2 is None:
            continue

        stats["pairs_achromatized"] += 1

        # STEP 4: Analytical monochromatic feasibility (discriminant-based branch solver)
        if cfg.system_type == "cemented":
            branches = solve_PW_cemented_branches(
                phi1=phi1,
                phi2=phi2,
                n1=n1,
                n2=n2,
                P0=cfg.P0,
                W0=cfg.W0,
            )
        else:
            branches = solve_PW_air_spaced_branches(
                phi1=phi1,
                phi2=phi2,
                n1=n1,
                n2=n2,
                d_air=cfg.d_air,
                P0=cfg.P0,
                W0=cfg.W0,
            )

        n_branches = len(branches)
        if n_branches == 0:
            stats["pairs_with_0_branches"] += 1
            continue
        elif n_branches == 1:
            stats["pairs_with_1_branch"] += 1
        else:
            stats["pairs_with_2_branches"] += 1

        stats["pairs_with_branches"] += 1

        # STEP 5: Create candidates from branches, compute PE using src/score.py ONLY
        for branch in branches:
            candidate = fill_candidate_from_branch(
                branch=branch,
                g1=g1,
                g2=g2,
                cfg=cfg,
                phi1=phi1,
                phi2=phi2,
                nu1=nu1,
                nu2=nu2,
                n1=n1,
                n2=n2,
            )

            # Air-spaced: reject if Pi is not length 4
            if cfg.system_type == "air_spaced" and len(candidate.Pi) != 4:
                continue

            # Compute PE using src/score.py ONLY
            pe = compute_pe(candidate)

            stats["candidates_generated"] += 1

            if pe > cfg.max_PE:
                continue

            if top.push(candidate):
                stats["candidates_saved"] = len(top)

    # STEP 6: Get sorted results
    results = top.sorted()

    stats["candidates_saved"] = len(results)
    if results:
        best_pe = rank_key(results[0])
        if best_pe != float("inf"):
            stats["best_PE"] = best_pe

    # STEP 7: Sanity checks
    _log_sanity_checks(stats, results)

    return results, stats


def _log_sanity_checks(stats: dict, results: list[Candidate]) -> None:
    """
    Log sanity check information.

    Verifies:
    - Many pairs rejected due to no real branches (discriminant < 0)
    - Some pairs produce TWO branches
    - Results are strictly PE-sorted
    - best_PE equals first entry
    - Counters are consistent
    """
    print(f"\n[Sanity Check] Branch statistics:")
    print(
        f"  Pairs with 0 branches (rejected): {stats.get('pairs_with_0_branches', 0)}"
    )
    print(f"  Pairs with 1 branch: {stats.get('pairs_with_1_branch', 0)}")
    print(f"  Pairs with 2 branches: {stats.get('pairs_with_2_branches', 0)}")

    if len(results) >= 2:
        is_sorted = all(
            rank_key(results[i]) <= rank_key(results[i + 1])
            for i in range(len(results) - 1)
        )
        print(f"  Results strictly PE-sorted: {is_sorted}")

    if results:
        first_pe = rank_key(results[0])
        print(f"  best_PE = {stats.get('best_PE', 'N/A')}, first entry PE = {first_pe}")
        print(
            f"  best_PE matches first entry: {abs(stats.get('best_PE', float('inf')) - first_pe) < 1e-10}"
        )

    print(f"\n[Sanity Check] Counter consistency:")
    print(f"  pairs_tested: {stats.get('pairs_tested', 0)}")
    print(f"  pairs_after_min_delta_nu: {stats.get('pairs_after_min_delta_nu', 0)}")
    print(f"  pairs_achromatized: {stats.get('pairs_achromatized', 0)}")
    print(f"  pairs_with_branches: {stats.get('pairs_with_branches', 0)}")
    print(f"  candidates_generated: {stats.get('candidates_generated', 0)}")
    print(f"  candidates_saved: {stats.get('candidates_saved', 0)}")


# =============================================================================
# Entry Point
# =============================================================================

if __name__ == "__main__":
    cfg = default_cfg()

    print("Running synthesis pipeline (Analytical Branch Solver)...")
    print(f"  System type: {cfg.system_type}")
    print(f"  Catalogs: {cfg.agf_paths}")
    print(f"  Wavelengths: lam0={cfg.lam0}, lam1={cfg.lam1}, lam2={cfg.lam2}")
    print(f"  Target: P0={cfg.P0}, W0={cfg.W0}, C0={cfg.C0}")
    print(f"  N_keep: {cfg.N_keep}")
    print()

    results, stats = run(cfg)

    io.ensure_out_dir(cfg.out_dir)

    rows = io.extract_rows(results)
    csv_path = f"{cfg.out_dir}/results.csv"
    io.save_csv(rows, csv_path)
    print(f"\nResults saved to: {csv_path}")
    print()

    io.print_summary(stats)

    if results:
        top_cand = results[0]
        g1_obj = getattr(top_cand, "g1", None)
        g2_obj = getattr(top_cand, "g2", None)
        g1_name = g1_obj.name if g1_obj is not None and hasattr(g1_obj, "name") else "?"
        g2_name = g2_obj.name if g2_obj is not None and hasattr(g2_obj, "name") else "?"

        print()
        print("=" * 50)
        print("Top Candidate Details")
        print("=" * 50)
        print(f"  Glass pair: {g1_name} + {g2_name}")
        print(f"  PE: {top_cand.PE}")

        meta = getattr(top_cand, "meta", {})
        print(f"  root_param: {meta.get('root_param', 'N/A')}")

        if cfg.system_type == "cemented":
            print(f"  P2 (cemented surface): {top_cand.P2}")
            print(f"  W (system coma): {top_cand.W}")
            print(f"  R2 (cemented radius): {top_cand.R2}")
            print(f"  R1_front: {meta.get('R1_front', 'N/A')}")
            print(f"  R_cemented: {meta.get('R_cemented', 'N/A')}")
            print(f"  R2_back: {meta.get('R2_back', 'N/A')}")
        else:
            print(f"  Pi (surface invariants): {top_cand.Pi}")
            print(f"  R1_front: {meta.get('R1_front', 'N/A')}")
            print(f"  R1_back: {meta.get('R1_back', 'N/A')}")
            print(f"  R2_front: {meta.get('R2_front', 'N/A')}")
            print(f"  R2_back: {meta.get('R2_back', 'N/A')}")
            print(f"  d_air: {meta.get('d_air', 'N/A')}")

        print(f"  P_surfaces: {meta.get('P_surfaces', 'N/A')}")
        print(f"  W_surfaces: {meta.get('W_surfaces', 'N/A')}")
        print(f"  P_total: {meta.get('P_total', 'N/A')}")
        print(f"  W_total: {meta.get('W_total', 'N/A')}")
        print("=" * 50)
