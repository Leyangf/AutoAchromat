"""
evaluator.py – Extract optical performance metrics from an optiland ``Optic``.

All metrics are returned as a plain ``dict`` so they can be serialised to
JSON, sorted, or fed into a ranking function without depending on optiland
at consumption time.

The same ``evaluate()`` function works for both Stage A and Stage B optics
(they are both ``optiland.optic.Optic`` instances).
"""

from __future__ import annotations

import logging
import time
import types
from contextlib import contextmanager
from dataclasses import dataclass, field, asdict
from typing import Optional, List, Dict, Any

import numpy as np
from optiland import optic
from optiland.analysis import SpotDiagram

from ..models import Candidate, Inputs
from .builder import build_optic

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Metrics data container
# ---------------------------------------------------------------------------


@dataclass
class OpticMetrics:
    """Optical performance metrics extracted from an optiland ``Optic``."""

    # Identification (carried over from the candidate)
    glass1: str = ""
    glass2: str = ""
    catalog1: str = ""
    catalog2: str = ""
    system_type: str = ""

    # Paraxial
    efl: Optional[float] = None  # effective focal length (f2)
    fno: Optional[float] = None  # image-space F/#
    bfd: Optional[float] = None  # back focal distance (last surface → image)

    # Spot diagram (on-axis, across all wavelengths)
    rms_spot_radius: Optional[float] = None  # mean RMS across wavelengths
    rms_spot_per_wl: List[float] = field(default_factory=list)
    geo_spot_radius: Optional[float] = None  # mean geometric radius
    geo_spot_per_wl: List[float] = field(default_factory=list)

    # Seidel aberration coefficients (sum over all surfaces)
    SA: Optional[float] = None  # spherical aberration  (SC)
    CC: Optional[float] = None  # coma                  (CC)
    AC: Optional[float] = None  # astigmatism            (AC)
    PC: Optional[float] = None  # Petzval curvature      (PC)
    DC: Optional[float] = None  # distortion             (DC)

    # Transverse Seidel coefficients
    TSC: Optional[float] = None  # transverse spherical
    TCC: Optional[float] = None  # transverse coma
    TAC: Optional[float] = None  # transverse astigmatism
    TPC: Optional[float] = None  # transverse Petzval
    TDC: Optional[float] = None  # transverse distortion (not standard name)

    # Chromatic
    LchC: Optional[float] = None  # longitudinal chromatic aberration
    TchC: Optional[float] = None  # transverse chromatic aberration
    secondary_spectrum: Optional[float] = None  # BFD(λ₀) − avg(BFD(λ₁),BFD(λ₂)) [mm]

    # Thin-lens synthesis values (for comparison / ranking)
    PE: Optional[float] = None
    Q: Optional[float] = None
    W: Optional[float] = None
    Q1: Optional[float] = None
    Q2: Optional[float] = None
    cost1: Optional[float] = None
    cost2: Optional[float] = None

    # Build / trace status
    success: bool = False
    error_msg: str = ""
    build_time_ms: float = 0.0
    eval_time_ms: float = 0.0

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


# ---------------------------------------------------------------------------
# Monkey-patch helper for optiland chromatic dispersion
# ---------------------------------------------------------------------------


@contextmanager
def _patched_dispersion(ab, op, inputs: Inputs):
    """Temporarily patch optiland's _precalculations to use design wavelengths.

    optiland hardcodes F/C Fraunhofer lines for _dn computation.  For
    non-visible designs this is wrong.  This context manager installs a
    patched version and guarantees cleanup via ``finally``.
    """
    _orig_precalc = type(ab)._precalculations

    def _patched_precalc(self_ab):
        _orig_precalc(self_ab)
        self_ab._dn = op.n(inputs.lam1) - op.n(inputs.lam2)

    ab._precalculations = types.MethodType(_patched_precalc, ab)
    try:
        yield
    finally:
        # Restore original behaviour by removing instance override
        if hasattr(ab, "_precalculations"):
            del ab._precalculations


# ---------------------------------------------------------------------------
# Secondary spectrum via ABCD
# ---------------------------------------------------------------------------


def _bfd_at_wavelength(
    rx, candidate: Candidate, lam_um: float,
) -> float | None:
    """Back focal distance at a single wavelength via ABCD matrix."""
    from ..optics import refractive_index
    from ..glass_reader import Glass
    from ..thickening import (
        _system_matrix_cemented, _system_matrix_spaced, _bfd_from_matrix,
    )

    if not candidate.cd1 or not candidate.cd2:
        return None
    if candidate.formula_id1 is None or candidate.formula_id2 is None:
        return None

    g1 = Glass(
        name=candidate.glass1, catalog=candidate.catalog1,
        formula_id=candidate.formula_id1, cd=candidate.cd1,
    )
    g2 = Glass(
        name=candidate.glass2, catalog=candidate.catalog2,
        formula_id=candidate.formula_id2, cd=candidate.cd2,
    )

    try:
        n1 = refractive_index(g1, lam_um)
        n2 = refractive_index(g2, lam_um)
    except (ValueError, ZeroDivisionError):
        return None

    e1, e2 = rx.elements

    if rx.system_type == "cemented":
        M = _system_matrix_cemented(
            e1.R_front, e1.R_back, e2.R_back,
            e1.t_center, e2.t_center, n1, n2,
        )
    else:
        air_gap = rx.air_gap if rx.air_gap is not None else 1.0
        M = _system_matrix_spaced(
            e1.R_front, e1.R_back, e2.R_front, e2.R_back,
            e1.t_center, e2.t_center, air_gap, n1, n2,
        )

    return _bfd_from_matrix(M)


def _secondary_spectrum(
    rx, candidate: Candidate, inputs: Inputs,
) -> float | None:
    """BFD(λ₀) − avg(BFD(λ₁), BFD(λ₂))  [mm].

    Uses ABCD matrices with the candidate's dispersion data.
    Returns None if dispersion data is unavailable.
    """
    if rx is None:
        return None

    bfd0 = _bfd_at_wavelength(rx, candidate, inputs.lam0)
    bfd1 = _bfd_at_wavelength(rx, candidate, inputs.lam1)
    bfd2 = _bfd_at_wavelength(rx, candidate, inputs.lam2)

    if bfd0 is None or bfd1 is None or bfd2 is None:
        return None

    return bfd0 - (bfd1 + bfd2) / 2.0


def chromatic_focal_shift(
    rx,
    candidate: Candidate,
    inputs: Inputs,
    n_points: int = 20,
) -> tuple[list[float], list[float]] | None:
    """Compute BFD(λ) curve across the design waveband.

    Returns (wavelengths_um, delta_bfd_mm) where delta_bfd is relative
    to BFD at the design wavelength λ₀.  Returns None if dispersion
    data is unavailable.
    """
    if rx is None:
        return None

    bfd0 = _bfd_at_wavelength(rx, candidate, inputs.lam0)
    if bfd0 is None:
        return None

    lam_min = min(inputs.lam1, inputs.lam2)
    lam_max = max(inputs.lam1, inputs.lam2)
    margin = (lam_max - lam_min) * 0.05
    lam_lo = lam_min - margin
    lam_hi = lam_max + margin

    wavelengths: list[float] = []
    delta_bfd: list[float] = []

    for i in range(n_points):
        lam = lam_lo + (lam_hi - lam_lo) * i / (n_points - 1)
        bfd = _bfd_at_wavelength(rx, candidate, lam)
        if bfd is None:
            continue
        wavelengths.append(lam)
        delta_bfd.append(bfd - bfd0)

    if len(wavelengths) < 2:
        return None

    return wavelengths, delta_bfd


# ---------------------------------------------------------------------------
# Single-candidate evaluation
# ---------------------------------------------------------------------------


def evaluate(
    op: optic.Optic,
    candidate: Candidate,
    inputs: Inputs,
    rx=None,
) -> OpticMetrics:
    """Extract performance metrics from an already-built ``Optic``.

    Parameters
    ----------
    op : optiland.optic.Optic
        A fully configured optic (aperture, field, wavelengths set;
        ``image_solve`` already called).
    candidate : Candidate
        The thin-lens candidate this optic was built from.
    inputs : Inputs
        System specification.

    Returns
    -------
    OpticMetrics
        Populated metrics dict.  ``success`` is True if all extractions
        succeeded; partial results are still stored on failure.
    """
    m = OpticMetrics(
        glass1=candidate.glass1,
        glass2=candidate.glass2,
        catalog1=candidate.catalog1,
        catalog2=candidate.catalog2,
        system_type=candidate.system_type,
        PE=candidate.PE,
        Q=candidate.Q,
        W=candidate.W,
        Q1=candidate.Q1,
        Q2=candidate.Q2,
        cost1=candidate.cost1,
        cost2=candidate.cost2,
    )

    t0 = time.perf_counter()

    try:
        # ---- Paraxial first-order data ----
        m.efl = float(op.paraxial.f2())
        m.fno = float(op.paraxial.FNO())

        # ---- Back focal distance (after image_solve) ----
        m.bfd = float(op.surface_group.surfaces[-2].thickness)

        # ---- Spot diagram (on-axis, all wavelengths) ----
        spot = SpotDiagram(
            op, fields="all", wavelengths="all", num_rings=6, distribution="hexapolar"
        )
        rms_all = spot.rms_spot_radius()  # list[list[float]]
        geo_all = spot.geometric_spot_radius()  # list[list[float]]

        # rms_all is [[val_wl0, val_wl1, val_wl2]] for 1 field
        # Filter out NaN values (from failed marginal rays)
        rms_flat = [float(v) for v in rms_all[0]]
        geo_flat = [float(v) for v in geo_all[0]]
        rms_valid = [v for v in rms_flat if np.isfinite(v)]
        geo_valid = [v for v in geo_flat if np.isfinite(v)]

        m.rms_spot_per_wl = rms_flat
        m.geo_spot_per_wl = geo_flat
        m.rms_spot_radius = float(np.mean(rms_valid)) if rms_valid else None
        m.geo_spot_radius = float(np.mean(geo_valid)) if geo_valid else None

        # ---- Seidel / third-order aberrations ----
        ab = op.aberrations

        # Patch optiland's hardcoded chromatic dispersion (F/C lines) to use
        # the actual design wavelengths.  optiland computes
        #   _dn = n(0.4861) - n(0.6563)
        # inside _precalculations(), but for non-visible designs this is wrong.
        with _patched_dispersion(ab, op, inputs):
            m.SA = _scalar(ab.SC())
            m.CC = _scalar(ab.CC())
            m.AC = _scalar(ab.AC())
            m.PC = _scalar(ab.PC())
            m.DC = _scalar(ab.DC())

            # Transverse Seidel
            m.TSC = _scalar(ab.TSC())
            m.TCC = _scalar(ab.TCC())
            m.TAC = _scalar(ab.TAC())
            m.TPC = _scalar(ab.TPC())

            # ---- Chromatic ----
            m.LchC = _scalar(ab.LchC())
            m.TchC = _scalar(ab.TchC())

        # ---- Secondary spectrum (ABCD, independent of optiland) ----
        m.secondary_spectrum = _secondary_spectrum(rx, candidate, inputs)

        m.success = True

    except Exception as exc:
        m.error_msg = f"{type(exc).__name__}: {exc}"
        logger.debug(
            "evaluate() failed for %s:%s + %s:%s – %s",
            candidate.catalog1,
            candidate.glass1,
            candidate.catalog2,
            candidate.glass2,
            m.error_msg,
            exc_info=True,
        )

    m.eval_time_ms = (time.perf_counter() - t0) * 1000.0
    return m


# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# Batch convenience  (DEPRECATED – use pipeline.run_pipeline instead)
# ---------------------------------------------------------------------------


def batch_evaluate(
    candidates: List[Candidate],
    inputs: Inputs,
    *,
    stage: str = "A",
    max_n: Optional[int] = None,
) -> List[OpticMetrics]:
    """Build + evaluate a list of candidates in one call.

    .. deprecated::
        Use :func:`autoachromat.pipeline.run_pipeline` instead.
        This function is kept only for backward compatibility with
        ``test_bridge.py`` and will be removed in a future release.

    Parameters
    ----------
    candidates : list[Candidate]
        From ``run_cemented`` or ``run_spaced``.
    inputs : Inputs
        System specification.
    stage : str
        ``"A"`` or (future) ``"B"``.
    max_n : int | None
        If given, only process the first *max_n* candidates.

    Returns
    -------
    list[OpticMetrics]
        One entry per candidate (in the same order).  Candidates that
        failed to build have ``success=False`` and ``error_msg`` set.
    """
    import warnings

    warnings.warn(
        "batch_evaluate() is deprecated; use pipeline.run_pipeline() instead.",
        DeprecationWarning,
        stacklevel=2,
    )

    todo = candidates[:max_n] if max_n else candidates
    results: List[OpticMetrics] = []
    n_ok = 0

    for i, cand in enumerate(todo):
        t_build = time.perf_counter()

        op = build_optic(cand, inputs, stage=stage)

        build_ms = (time.perf_counter() - t_build) * 1000.0

        if op is None:
            m = OpticMetrics(
                glass1=cand.glass1,
                glass2=cand.glass2,
                catalog1=cand.catalog1,
                catalog2=cand.catalog2,
                system_type=cand.system_type,
                PE=cand.PE,
                Q=cand.Q,
                W=cand.W,
                Q1=cand.Q1,
                Q2=cand.Q2,
                cost1=cand.cost1,
                cost2=cand.cost2,
                success=False,
                error_msg="build_optic returned None",
                build_time_ms=build_ms,
            )
        else:
            m = evaluate(op, cand, inputs)
            m.build_time_ms = build_ms
            if m.success:
                n_ok += 1

        results.append(m)

    logger.info(
        "batch_evaluate: %d / %d succeeded  (total %.1f ms)",
        n_ok,
        len(todo),
        sum(r.build_time_ms + r.eval_time_ms for r in results),
    )
    return results


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _scalar(arr) -> Optional[float]:
    """Safely convert a numpy array / scalar to a Python float."""
    try:
        if hasattr(arr, "__len__"):
            return float(np.sum(arr))  # sum over surfaces
        return float(arr)
    except (TypeError, ValueError):
        return None
