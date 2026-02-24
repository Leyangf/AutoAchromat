from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Literal

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

    # Numerical
    eps: float = 1e-12
    root_imag_tol: float = 1e-9


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

    notes: dict = field(default_factory=dict)
