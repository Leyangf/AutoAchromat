"""
thermal.py – First-order thermal (athermal) analysis for doublets.

Implements the Schott thermal dispersion formula and the thin-lens
passive athermalization condition.  All public functions return None
when glass TD/ED data is missing rather than raising exceptions.

Physical background
-------------------
Lens power:  Φ = (n−1)(C₁−C₂)

Temperature derivative (radii expand at rate α, index changes at dn/dT):

    dΦ/dT = Φ · V      where  V = dn/dT / (n−1) − α    [1/K]

Schott first-order formula at reference temperature T_ref:

    dn/dT = (n²−1)/(2n) × [D₀ + E₀/(λ²−λ_tk²)]    [1/K]

Doublet (normalised powers phi1+phi2=1, thin-lens approximation):

    dΦ/dT_norm = V₁·phi1 + V₂·phi2    [1/K]

Passive athermalization (image tracks detector):

    df'/dT = dl/dT  →  α_h,required = −dΦ/dT_norm    [1/K]

V < 0 for most glasses  →  α_h,required > 0  (physically achievable).
"""
from __future__ import annotations

from dataclasses import dataclass

from .glass_reader import Glass


# ---------------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------------


@dataclass
class ThermalMetrics:
    """First-order thermal analysis result for a doublet glass pair.

    Internal quantities are in SI [1/K].
    Multiply by 1e6 to convert to the more readable ppm/K.
    """

    dn_dT_1: float | None = None
    """dn/dT for element 1 at design wavelength [1/K]."""

    dn_dT_2: float | None = None
    """dn/dT for element 2 at design wavelength [1/K]."""

    V1: float | None = None
    """Thermo-optical coefficient element 1: dn/dT/(n−1) − α [1/K]."""

    V2: float | None = None
    """Thermo-optical coefficient element 2 [1/K]."""

    dphi_dT_norm: float | None = None
    """Normalised system thermal power derivative V₁·phi1+V₂·phi2 [1/K].
    Typically negative (focal length increases with temperature)."""

    alpha_housing_required: float | None = None
    """Required housing CTE for passive athermalization: −dphi_dT_norm [1/K].
    Typically positive (achievable with real materials)."""

    thermal_data_available: bool = False
    """True when both glasses provide complete TD and ED thermal data."""


# ---------------------------------------------------------------------------
# Elementary functions
# ---------------------------------------------------------------------------


def dn_dT(glass: Glass, wavelength_um: float, n: float) -> float | None:
    """Schott first-order thermal dispersion at reference temperature.

    Formula: dn/dT = (n²−1)/(2n) × [D₀ + E₀/(λ²−λ_tk²)]    [1/K]

    Parameters
    ----------
    glass :
        Glass record with td_D0, td_E0, td_ltk populated from AGF TD line.
    wavelength_um :
        Evaluation wavelength [µm].
    n :
        Refractive index at *wavelength_um* (passed in to avoid recomputing).

    Returns
    -------
    float | None
        dn/dT in [1/K], or None when any required TD field is missing.
    """
    if glass.td_D0 is None or glass.td_E0 is None or glass.td_ltk is None:
        return None
    lam2 = wavelength_um * wavelength_um
    ltk2 = glass.td_ltk * glass.td_ltk
    denom = lam2 - ltk2
    if abs(denom) < 1e-15:
        return None
    factor = (n * n - 1.0) / (2.0 * n)
    return factor * (glass.td_D0 + glass.td_E0 / denom)


def thermo_optical_coeff(
    glass: Glass, wavelength_um: float, n: float
) -> float | None:
    """Thermo-optical coefficient V = dn/dT / (n−1) − α    [1/K].

    Uses CTE from −40 to +20°C (AGF ED parts[1], stored in ppm/K).
    The ppm/K value is converted to [1/K] internally (× 1e-6).

    Physical meaning
    ----------------
    V < 0 (most glasses): temperature rise → focal length increases
                          → positive-CTE housing can compensate.
    V > 0 (some high-index): temperature rise → focal length decreases.

    Returns None when dn/dT or CTE data is unavailable.
    """
    dndt = dn_dT(glass, wavelength_um, n)
    if dndt is None:
        return None
    if glass.cte_m40_20 is None:
        return None
    alpha = glass.cte_m40_20 * 1e-6   # ppm/K → 1/K
    n_minus_1 = n - 1.0
    if abs(n_minus_1) < 1e-15:
        return None
    return dndt / n_minus_1 - alpha


def system_thermal_power_derivative(
    V1: float | None,
    phi1: float,
    V2: float | None,
    phi2: float,
) -> float | None:
    """Normalised system thermal power derivative [1/K].

    dΦ/dT_norm = V₁·phi1 + V₂·phi2

    Exact for cemented doublets (d = 0).  For air-spaced doublets the
    correction terms are O(d/f') ≈ 3–6 % for d ≤ 2 mm, f' ≥ 100 mm,
    which is within the thin-lens approximation accuracy.

    Returns None when either V is None.
    """
    if V1 is None or V2 is None:
        return None
    return V1 * phi1 + V2 * phi2


def required_housing_cte(dphi_dT_norm: float | None) -> float | None:
    """Required housing CTE for passive athermalization [1/K].

    Derived from the athermalization condition df'/dT = dl/dT:

        α_h,required = −dΦ/dT_norm

    Sign behaviour
    --------------
    dΦ/dT_norm < 0 (typical)  →  α_h,required > 0  (physically achievable).
    α_h,required < 0           →  cannot be met with standard materials.

    Returns None when *dphi_dT_norm* is None.
    """
    if dphi_dT_norm is None:
        return None
    return -dphi_dT_norm


def thermal_defocus(
    f_prime: float,
    alpha_h_required: float | None,
    alpha_h_actual: float,
    delta_T: float,
) -> float | None:
    """Thermal defocus for a given housing CTE [mm].

    δf' = f' × (α_h,required − α_h_actual) × ΔT

    Sign convention
    ---------------
    Positive: focal point moves behind detector (housing expands too little).
    Negative: focal point moves in front of detector (housing over-expands;
              common for aluminium housings).

    Returns None when *alpha_h_required* is None.
    """
    if alpha_h_required is None:
        return None
    return f_prime * (alpha_h_required - alpha_h_actual) * delta_T


# ---------------------------------------------------------------------------
# Integrated entry point
# ---------------------------------------------------------------------------


def compute_thermal_metrics(
    glass1: Glass,
    glass2: Glass,
    n1: float,
    n2: float,
    phi1: float,
    phi2: float,
    wavelength_um: float,
) -> ThermalMetrics:
    """Compute all first-order thermal metrics for a doublet glass pair.

    Called from ``cemented.py`` and ``spaced.py`` during thin-lens synthesis.
    Never raises; returns a ``ThermalMetrics`` with ``None`` fields when
    glass TD/ED data is incomplete.

    Parameters
    ----------
    glass1, glass2 :
        Glass records from the AGF catalog (must have td_* and cte_* fields
        populated for meaningful results).
    n1, n2 :
        Refractive indices at *wavelength_um* (already computed by
        ``prepare_glass_data``).
    phi1, phi2 :
        Normalised element powers from ``achromat_power()``; phi1+phi2 = 1.
    wavelength_um :
        Design wavelength [µm] (typically ``inputs.lam0``).
    """
    try:
        dndt1 = dn_dT(glass1, wavelength_um, n1)
        dndt2 = dn_dT(glass2, wavelength_um, n2)
        V1 = thermo_optical_coeff(glass1, wavelength_um, n1)
        V2 = thermo_optical_coeff(glass2, wavelength_um, n2)
        dphi = system_thermal_power_derivative(V1, phi1, V2, phi2)
        alpha_req = required_housing_cte(dphi)
        available = (
            dndt1 is not None
            and dndt2 is not None
            and V1 is not None
            and V2 is not None
        )
        return ThermalMetrics(
            dn_dT_1=dndt1,
            dn_dT_2=dndt2,
            V1=V1,
            V2=V2,
            dphi_dT_norm=dphi,
            alpha_housing_required=alpha_req,
            thermal_data_available=available,
        )
    except Exception:
        return ThermalMetrics()
