"""
Optiland interface for real-ray aberration calculation.

This module builds optical systems from ranked candidates and computes
real-ray aberrations using the Optiland ray tracing library.

References:
- ITMO doublet synthesis software behavior
- Optiland ray tracing library documentation
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, TYPE_CHECKING

if TYPE_CHECKING:
    from ..synthesis import Config

# Check if Optiland is available
try:
    from optiland.optic import Optic
    from optiland.materials import IdealMaterial
    import optiland.backend as be

    OPTILAND_AVAILABLE = True
except ImportError:
    Optic = None  # type: ignore[misc, assignment]
    IdealMaterial = None  # type: ignore[misc, assignment]
    be = None  # type: ignore[assignment]
    OPTILAND_AVAILABLE = False

# Glass catalog cache
_glass_cache: dict[str, Any] = {}


def _load_glass_catalogs(agf_paths: list[str]) -> None:
    """Load glass catalogs and cache them."""
    global _glass_cache

    if not _glass_cache:
        from ..glass_reader import read_agf

        for agf_path in agf_paths:
            try:
                _, glasses = read_agf(agf_path)
                for g in glasses:
                    # Key by name and catalog
                    key = f"{g.name}|{g.catalog}".upper()
                    _glass_cache[key] = g
            except Exception:
                pass  # Silently skip missing catalogs


def _find_glass(name: str, catalog: str, agf_paths: list[str]) -> Any:
    """Find a glass by name and catalog."""
    _load_glass_catalogs(agf_paths)

    key = f"{name}|{catalog}".upper()
    return _glass_cache.get(key)


def _compute_n_at_wavelength(glass: Any, wavelength_um: float) -> float | None:
    """Compute refractive index at a specific wavelength using glass dispersion."""
    if glass is None:
        return None

    try:
        from ..synthesis import refractive_index

        return refractive_index(glass, wavelength_um)
    except Exception:
        return None


@dataclass
class RealRayResults:
    """
    Container for real-ray aberration calculation results.

    All values are in mm unless otherwise noted.

    Attributes:
        transverse_aberration: Ray height at paraxial focus plane (marginal ray)
        longitudinal_aberration: Distance from axis crossing to paraxial focus
        axial_color: z_focus(lam1) - z_focus(lam2)
        secondary_spectrum: Residual chromatic focus relative to lam0
        coma: y_marginal - y_chief at paraxial focus (off-axis)
        success: Whether calculation was successful
        error_message: Error message if calculation failed
    """

    transverse_aberration: float | None = None
    longitudinal_aberration: float | None = None
    axial_color: float | None = None
    secondary_spectrum: float | None = None
    coma: float | None = None
    success: bool = False
    error_message: str | None = None


def _check_optiland() -> str | None:
    """
    Check if Optiland is available.

    Returns:
        Error message if Optiland is not available, None otherwise.
    """
    if not OPTILAND_AVAILABLE:
        return (
            "Optiland is not installed. Please install it with:\n"
            "  pip install optiland\n"
            "Real-ray aberration calculation requires Optiland."
        )
    return None


def _assign_thicknesses(
    nu1: float,
    nu2: float,
    crown_thickness: float,
    flint_thickness: float,
) -> tuple[float, float]:
    """
    Assign lens thicknesses based on crown/flint rule.

    Crown = glass with higher Abbe number (nu)
    Flint = glass with lower Abbe number

    Args:
        nu1: Abbe number of glass 1
        nu2: Abbe number of glass 2
        crown_thickness: User-specified crown lens thickness (mm)
        flint_thickness: User-specified flint lens thickness (mm)

    Returns:
        Tuple of (lens1_thickness, lens2_thickness) in mm
    """
    if nu1 >= nu2:
        # Glass 1 is crown, glass 2 is flint
        return crown_thickness, flint_thickness
    else:
        # Glass 1 is flint, glass 2 is crown
        return flint_thickness, crown_thickness


def _safe_radius(r: float) -> float:
    """Convert very large radii to infinity."""
    if not OPTILAND_AVAILABLE or be is None:
        return float("inf")
    if r is None or abs(r) > 1e9:
        return be.inf
    return r


def _build_cemented_doublet(
    candidate_row: dict[str, Any],
    cfg: "Config",
    n1_override: float | None = None,
    n2_override: float | None = None,
) -> Any:
    """
    Build an Optiland optical system for a cemented doublet.

    Surface layout (propagation order):
        0) Object at infinity
        1) R1_front: air -> lens1 (stop)
        2) R_cemented: lens1 -> lens2
        3) R2_back: lens2 -> air
        4) Image surface

    Entrance pupil (stop) is located on the first surface.
    Object at infinity.

    Args:
        candidate_row: Row dictionary from results table
        cfg: Configuration object
        n1_override: Override n1 for wavelength-specific calculation
        n2_override: Override n2 for wavelength-specific calculation

    Returns:
        Optiland Optic object
    """
    # Extract geometry from candidate
    R1_front = candidate_row.get("R1_front")
    R_cemented = candidate_row.get("R_cemented")
    R2_back = candidate_row.get("R2_back")

    if R1_front is None or R_cemented is None or R2_back is None:
        raise ValueError("Missing radius data for cemented doublet")

    # Extract glass properties
    n1 = n1_override if n1_override is not None else candidate_row.get("n1")
    n2 = n2_override if n2_override is not None else candidate_row.get("n2")
    nu1 = candidate_row.get("nu1")
    nu2 = candidate_row.get("nu2")

    if n1 is None or n2 is None or nu1 is None or nu2 is None:
        raise ValueError("Missing glass property data")

    # Assign thicknesses based on crown/flint rule
    t1, t2 = _assign_thicknesses(
        nu1,
        nu2,
        cfg.crown_lens_thickness_mm,
        cfg.flint_lens_thickness_mm,
    )

    # Handle very large radii
    R1 = _safe_radius(R1_front)
    R2 = _safe_radius(R_cemented)
    R3 = _safe_radius(R2_back)

    # Create optical system
    if Optic is None or IdealMaterial is None or be is None:
        raise RuntimeError("Optiland is not installed")
    lens = Optic()

    # Set aperture (entrance pupil diameter)
    lens.set_aperture(aperture_type="EPD", value=cfg.D)

    # Set field type and add on-axis field
    lens.set_field_type(field_type="angle")
    lens.add_field(y=0.0, x=0.0)

    # Add wavelengths (in micrometers)
    lens.add_wavelength(value=cfg.lam0, is_primary=True, unit="um")
    lens.add_wavelength(value=cfg.lam1, is_primary=False, unit="um")
    lens.add_wavelength(value=cfg.lam2, is_primary=False, unit="um")

    # Add object surface (at infinity)
    lens.add_surface(index=0, radius=be.inf, thickness=be.inf)

    # Surface 1: Front of lens 1 (R1_front), stop surface
    lens.add_surface(
        index=1,
        radius=R1,
        thickness=t1,
        material=IdealMaterial(n=n1),
        is_stop=True,
    )

    # Surface 2: Cemented interface (R_cemented)
    lens.add_surface(
        index=2,
        radius=R2,
        thickness=t2,
        material=IdealMaterial(n=n2),
    )

    # Surface 3: Back of lens 2 (R2_back) - air after
    lens.add_surface(
        index=3,
        radius=R3,
        thickness=100.0,  # Initial guess, will be solved for focus
    )

    # Add image surface
    lens.add_surface(index=4, radius=be.inf, thickness=0.0)

    # Solve for focus: move image plane so marginal ray crosses axis
    lens.update()
    lens.image_solve()

    return lens


def _build_air_spaced_doublet(
    candidate_row: dict[str, Any],
    cfg: "Config",
    n1_override: float | None = None,
    n2_override: float | None = None,
) -> Any:
    """
    Build an Optiland optical system for an air-spaced doublet.

    Surface layout (propagation order):
        0) Object at infinity
        1) R1_front: air -> lens1 (stop)
        2) R1_back: lens1 -> air
        3) R2_front: air -> lens2
        4) R2_back: lens2 -> air
        5) Image surface

    Thicknesses:
        - lens1 thickness (assigned by crown/flint rule)
        - air gap = d_air
        - lens2 thickness (assigned by crown/flint rule)

    Entrance pupil (stop) is located on the first surface.
    Object at infinity.

    Args:
        candidate_row: Row dictionary from results table
        cfg: Configuration object
        n1_override: Override n1 for wavelength-specific calculation
        n2_override: Override n2 for wavelength-specific calculation

    Returns:
        Optiland Optic object
    """
    # Extract geometry from candidate
    R1_front = candidate_row.get("R1_front")
    R1_back = candidate_row.get("R1_back")
    R2_front = candidate_row.get("R2_front")
    R2_back = candidate_row.get("R2_back")

    if R1_front is None or R1_back is None or R2_front is None or R2_back is None:
        raise ValueError("Missing radius data for air-spaced doublet")

    # Extract glass properties
    n1 = n1_override if n1_override is not None else candidate_row.get("n1")
    n2 = n2_override if n2_override is not None else candidate_row.get("n2")
    nu1 = candidate_row.get("nu1")
    nu2 = candidate_row.get("nu2")

    if n1 is None or n2 is None or nu1 is None or nu2 is None:
        raise ValueError("Missing glass property data")

    # Assign thicknesses based on crown/flint rule
    t1, t2 = _assign_thicknesses(
        nu1,
        nu2,
        cfg.crown_lens_thickness_mm,
        cfg.flint_lens_thickness_mm,
    )

    # Air gap thickness
    d_air = cfg.d_air

    # Handle very large radii
    R1 = _safe_radius(R1_front)
    R2 = _safe_radius(R1_back)
    R3 = _safe_radius(R2_front)
    R4 = _safe_radius(R2_back)

    # Create optical system
    if Optic is None or IdealMaterial is None or be is None:
        raise RuntimeError("Optiland is not installed")
    lens = Optic()

    # Set aperture (entrance pupil diameter)
    lens.set_aperture(aperture_type="EPD", value=cfg.D)

    # Set field type and add on-axis field
    lens.set_field_type(field_type="angle")
    lens.add_field(y=0.0, x=0.0)

    # Add wavelengths (in micrometers)
    lens.add_wavelength(value=cfg.lam0, is_primary=True, unit="um")
    lens.add_wavelength(value=cfg.lam1, is_primary=False, unit="um")
    lens.add_wavelength(value=cfg.lam2, is_primary=False, unit="um")

    # Add object surface (at infinity)
    lens.add_surface(index=0, radius=be.inf, thickness=be.inf)

    # Surface 1: Front of lens 1 (R1_front), stop surface
    lens.add_surface(
        index=1,
        radius=R1,
        thickness=t1,
        material=IdealMaterial(n=n1),
        is_stop=True,
    )

    # Surface 2: Back of lens 1 (R1_back) - air gap after
    lens.add_surface(
        index=2,
        radius=R2,
        thickness=d_air,
    )

    # Surface 3: Front of lens 2 (R2_front)
    lens.add_surface(
        index=3,
        radius=R3,
        thickness=t2,
        material=IdealMaterial(n=n2),
    )

    # Surface 4: Back of lens 2 (R2_back)
    lens.add_surface(
        index=4,
        radius=R4,
        thickness=100.0,  # Initial guess, will be solved for focus
    )

    # Add image surface
    lens.add_surface(index=5, radius=be.inf, thickness=0.0)

    # Solve for focus: move image plane so marginal ray crosses axis
    lens.update()
    lens.image_solve()

    return lens


def _get_ray_at_image(
    lens: Any, wavelength: float, Hx: float, Hy: float, Px: float, Py: float
) -> tuple[float, float, float]:
    """
    Trace a ray and get its position at the image surface.

    Args:
        lens: Optiland Optic object
        wavelength: Wavelength in micrometers
        Hx, Hy: Normalized field coordinates
        Px, Py: Normalized pupil coordinates

    Returns:
        Tuple of (y_position, slope_M, slope_N) at image surface
    """
    lens.trace_generic(Hx=Hx, Hy=Hy, Px=Px, Py=Py, wavelength=wavelength)

    img = lens.image_surface

    # Get y position at image
    y = 0.0
    if hasattr(img, "y") and img.y is not None:
        y_arr = img.y
        if hasattr(y_arr, "__len__") and len(y_arr) > 0:
            y = float(y_arr[0])
        elif hasattr(y_arr, "item"):
            y = float(y_arr.item())
        else:
            y = float(y_arr)

    # Get direction cosines M (y) and N (z)
    M = 0.0
    N = 1.0
    if hasattr(img, "M") and img.M is not None:
        M_arr = img.M
        if hasattr(M_arr, "__len__") and len(M_arr) > 0:
            M = float(M_arr[0])
        elif hasattr(M_arr, "item"):
            M = float(M_arr.item())
        else:
            M = float(M_arr)

    if hasattr(img, "N") and img.N is not None:
        N_arr = img.N
        if hasattr(N_arr, "__len__") and len(N_arr) > 0:
            N = float(N_arr[0])
        elif hasattr(N_arr, "item"):
            N = float(N_arr.item())
        else:
            N = float(N_arr)

    return y, M, N


def _build_system_for_wavelength(
    candidate_row: dict[str, Any],
    cfg: "Config",
    wavelength_um: float,
) -> Any:
    """
    Build an optical system with refractive indices at a specific wavelength.

    Uses glass dispersion data to get wavelength-specific n values.

    Args:
        candidate_row: Row dictionary from results table
        cfg: Configuration object
        wavelength_um: Wavelength in micrometers

    Returns:
        Optiland Optic object with correct n for the wavelength
    """
    system_type = candidate_row.get("system_type", "cemented")

    # Try to get wavelength-specific n values from glass data
    glass1_name = candidate_row.get("glass1_name")
    glass1_catalog = candidate_row.get("glass1_catalog")
    glass2_name = candidate_row.get("glass2_name")
    glass2_catalog = candidate_row.get("glass2_catalog")

    n1_at_lam = None
    n2_at_lam = None

    if glass1_name and glass1_catalog:
        glass1 = _find_glass(glass1_name, glass1_catalog, cfg.agf_paths)
        n1_at_lam = _compute_n_at_wavelength(glass1, wavelength_um)

    if glass2_name and glass2_catalog:
        glass2 = _find_glass(glass2_name, glass2_catalog, cfg.agf_paths)
        n2_at_lam = _compute_n_at_wavelength(glass2, wavelength_um)

    # Build system
    if system_type == "cemented":
        return _build_cemented_doublet(candidate_row, cfg, n1_at_lam, n2_at_lam)
    else:
        return _build_air_spaced_doublet(candidate_row, cfg, n1_at_lam, n2_at_lam)


def compute_realray_aberrations(
    candidate_row: dict[str, Any],
    cfg: "Config",
) -> RealRayResults:
    """
    Compute real-ray aberrations for a single candidate.

    This function is called ONLY after ranking, for a user-selected candidate.
    It builds Optiland optical systems with correct wavelength-specific
    refractive indices and computes:

    - Transverse Aberration [mm]: ray height at z_ref (marginal axial ray)
    - Longitudinal Aberration [mm]: z_cross - z_ref
    - Axial Color [mm]: z_focus(lam1) - z_focus(lam2)
    - Secondary Spectrum [mm]: residual chromatic focus relative to lam0
    - Coma [mm]: y_marginal - y_chief at z_ref (off-axis)

    Args:
        candidate_row: Row dictionary from results table containing geometry
        cfg: Configuration object with wavelengths and new thickness parameters

    Returns:
        RealRayResults with computed aberrations or error message
    """
    results = RealRayResults()

    # Check Optiland availability
    error = _check_optiland()
    if error:
        results.error_message = error
        return results

    try:
        # Wavelengths (in micrometers)
        lam0 = cfg.lam0
        lam1 = cfg.lam1
        lam2 = cfg.lam2

        # Build systems for each wavelength (with correct n values)
        lens_lam0 = _build_system_for_wavelength(candidate_row, cfg, lam0)
        lens_lam1 = _build_system_for_wavelength(candidate_row, cfg, lam1)
        lens_lam2 = _build_system_for_wavelength(candidate_row, cfg, lam2)

        # ------------------------------------------------------------------
        # Transverse and Longitudinal Aberration at lam0
        # Trace marginal ray: Hx=0, Hy=0, Px=0, Py=1.0
        # ------------------------------------------------------------------
        y_marginal, M, N = _get_ray_at_image(
            lens_lam0, lam0, Hx=0.0, Hy=0.0, Px=0.0, Py=1.0
        )

        # Transverse aberration = ray height at image plane
        results.transverse_aberration = y_marginal

        # Longitudinal aberration: estimate from ray slope
        # If M is the y-direction cosine, slope = M/N
        # Longitudinal = -y / slope (where ray crosses axis)
        if abs(M) > 1e-10 and abs(N) > 1e-10:
            slope = M / N
            results.longitudinal_aberration = (
                -y_marginal / slope if abs(slope) > 1e-10 else 0.0
            )
        else:
            results.longitudinal_aberration = 0.0

        # ------------------------------------------------------------------
        # Axial Color: difference in focus position between wavelengths
        # Build separate systems for each wavelength with correct n values
        # ------------------------------------------------------------------
        y_marg_lam1, M1, N1 = _get_ray_at_image(
            lens_lam1, lam1, Hx=0.0, Hy=0.0, Px=0.0, Py=1.0
        )
        y_marg_lam2, M2, N2 = _get_ray_at_image(
            lens_lam2, lam2, Hx=0.0, Hy=0.0, Px=0.0, Py=1.0
        )

        # Compute focus shift from marginal ray intercept
        def focus_shift(y, M, N):
            if abs(M) > 1e-10 and abs(N) > 1e-10:
                return -y * N / M
            return 0.0

        delta_z_lam0 = focus_shift(y_marginal, M, N)
        delta_z_lam1 = focus_shift(y_marg_lam1, M1, N1)
        delta_z_lam2 = focus_shift(y_marg_lam2, M2, N2)

        results.axial_color = delta_z_lam1 - delta_z_lam2

        # ------------------------------------------------------------------
        # Secondary Spectrum: residual chromatic focus relative to lam0
        # ------------------------------------------------------------------
        mean_focus = (delta_z_lam1 + delta_z_lam2) / 2.0
        results.secondary_spectrum = mean_focus - delta_z_lam0

        # ------------------------------------------------------------------
        # Coma: at field_for_aberration_deg and lam0
        # coma = y_marginal(z_ref) - y_chief(z_ref)
        # ------------------------------------------------------------------
        field_deg = cfg.field_for_aberration_deg

        if abs(field_deg) > 1e-6:
            # Add off-axis field to system
            lens_lam0.add_field(y=field_deg, x=0.0)
            lens_lam0.update()

            # Trace marginal ray at field (Hy=1.0 for full field, Py=1.0 for edge of pupil)
            y_marg_field, _, _ = _get_ray_at_image(
                lens_lam0, lam0, Hx=0.0, Hy=1.0, Px=0.0, Py=1.0
            )

            # Trace chief ray at field (Hy=1.0 for full field, Py=0.0 for center of pupil)
            y_chief_field, _, _ = _get_ray_at_image(
                lens_lam0, lam0, Hx=0.0, Hy=1.0, Px=0.0, Py=0.0
            )

            results.coma = y_marg_field - y_chief_field
        else:
            # No field specified, coma is zero on-axis
            results.coma = 0.0

        results.success = True

    except Exception as exc:
        results.error_message = f"Real-ray calculation failed: {exc}"
        results.success = False

    return results


def compute_realray_for_candidate(
    candidate_row: dict[str, Any],
    cfg: "Config",
) -> dict[str, Any]:
    """
    Compute real-ray aberrations and return as dictionary for table integration.

    This is the main entry point called by the CLI and GUI after ranking.

    Args:
        candidate_row: Row dictionary from results table
        cfg: Configuration object

    Returns:
        Dictionary with aberration column keys and values (mm)
    """
    results = compute_realray_aberrations(candidate_row, cfg)

    return {
        "Transverse Aberration [mm]": results.transverse_aberration,
        "Longitudinal Aberration [mm]": results.longitudinal_aberration,
        "Axial Color [mm]": results.axial_color,
        "Secondary Spectrum [mm]": results.secondary_spectrum,
        "Coma [mm]": results.coma,
        "_realray_success": results.success,
        "_realray_error": results.error_message,
    }
