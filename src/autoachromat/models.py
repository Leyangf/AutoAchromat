from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Optional, Literal, TYPE_CHECKING

if TYPE_CHECKING:
    from .thermal import ThermalMetrics

SystemType = Literal["cemented", "spaced"]


@dataclass(frozen=True)
class Inputs:
    # Spectral (um)
    lam0: float
    lam1: float
    lam2: float

    # Geometry (same unit for D and fprime, e.g. mm)
    D: float
    fprime: float

    # Targets
    C0: float
    P0: float
    W0: float

    # Screening constraints
    min_delta_nu: float
    max_PE: float
    N: int  # cemented TopN; spaced can ignore

    # Type
    system_type: SystemType

    # Spaced doublet air gap (mm); ignored for cemented
    air_gap: float = 1.0

    # Half field angle [degrees] for off-axis evaluation in Stage B
    half_field_angle: float = 1.0

    # Manufacturing minimum thickness [mm] (0 = auto from Table 10-3)
    # te_min: min edge thickness for positive lenses
    # tc_min: min centre thickness for negative lenses
    # Auto values depend on D; GUI auto-fills from lookup tables.
    te_min: float = 0.0
    tc_min: float = 0.0

    # Stage B centre-thickness bounds [mm] (0 = auto)
    # Auto: t_min = manufacturing minimum, t_max = max(5×t_min, 20)
    # User values are clamped to respect manufacturing limits.
    t1_min: float = 0.0
    t1_max: float = 0.0
    t2_min: float = 0.0
    t2_max: float = 0.0

    # Stage B air-gap bounds [mm] (0 = auto); spaced only
    gap_min: float = 0.0
    gap_max: float = 0.0

    # Numerical
    eps: float = 1e-12
    root_imag_tol: float = 1e-9

    # Thermal analysis
    T_ref: float = 20.0  # reference temperature [°C]
    T_delta: float = 20.0  # ΔT for thermal defocus estimate [K]
    alpha_housing: Optional[float] = (
        None  # actual housing CTE [1/K]; None = report only
    )

    @classmethod
    def with_defaults(cls, **overrides) -> "Inputs":
        """Create an ``Inputs`` with sensible defaults, overriding any field.

        Useful for tests and quick prototyping where only a few
        parameters differ from the standard d-line configuration.

        >>> inp = Inputs.with_defaults(fprime=100.0, system_type="spaced")
        """
        defaults: dict[str, Any] = dict(
            lam0=0.58756,
            lam1=0.48613,
            lam2=0.65627,
            D=50.0,
            fprime=200.0,
            C0=0.0,
            P0=0.0,
            W0=0.0,
            min_delta_nu=10.0,
            max_PE=100.0,
            N=20,
            system_type="cemented",
            air_gap=1.0,
            half_field_angle=1.0,
        )
        defaults.update(overrides)
        return cls(**defaults)


@dataclass
class Candidate:
    system_type: SystemType

    glass1: str
    glass2: str
    catalog1: str
    catalog2: str

    # n at lam0, nu computed from (lam0, lam1, lam2)
    n1: float
    n2: float
    nu1: float
    nu2: float

    phi1: float
    phi2: float

    # cemented
    Q: Optional[float] = None
    W: Optional[float] = None
    P2: Optional[float] = None

    # spaced
    Q1: Optional[float] = None
    Q2: Optional[float] = None
    P_surfaces: list[float] = field(default_factory=list)

    radii: list[float] = field(default_factory=list)
    PE: Optional[float] = None

    # Glass relative cost (from AGF OD line; None = unknown)
    cost1: Optional[float] = None
    cost2: Optional[float] = None

    # Glass dispersion data (for optiland builder)
    formula_id1: Optional[int] = None
    cd1: list[float] = field(default_factory=list)
    formula_id2: Optional[int] = None
    cd2: list[float] = field(default_factory=list)

    notes: dict = field(default_factory=dict)

    # Thermal analysis result (None when glass lacks TD/ED data)
    thermal: Optional["ThermalMetrics"] = None


# ---------------------------------------------------------------------------
# Thick-lens prescription  (output of thickening, input to builder)
# ---------------------------------------------------------------------------


@dataclass
class ElementRx:
    """Thick-lens prescription for a single element."""

    R_front: float  # front radius [mm] (unchanged from thin-lens)
    R_back: float  # back radius [mm] (unchanged from thin-lens)
    t_center: float  # center thickness [mm]
    t_edge: float  # edge thickness [mm]
    nd: float  # refractive index at design wavelength
    vd: float  # Abbe number

    # Glass dispersion data (optional, for accurate ray tracing)
    formula_id: Optional[int] = None
    cd: list[float] = field(default_factory=list)


@dataclass
class ThickPrescription:
    """Complete thick-lens prescription for a doublet.

    Produced by ``thickening.thicken()`` and consumed by
    ``optiland_bridge.builder.build_optic()``.
    """

    system_type: SystemType  # "cemented" | "spaced"
    elements: list[ElementRx]  # always 2 elements
    air_gap: float | None  # spaced doublet only (mm)
    back_focus_guess: float  # initial back focal distance for image_solve [mm]
    D: float  # entrance pupil diameter [mm]
    wavelengths: tuple[float, float, float]  # (lam1, lam0, lam2) [µm]
    half_field_angle: float = 1.0  # off-axis field angle [degrees]
    actual_efl: float | None = None  # thick-lens EFL from ABCD [mm]
    efl_deviation: float | None = None  # (actual - target) / target
