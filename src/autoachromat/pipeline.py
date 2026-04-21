"""
pipeline.py – Shared thicken → build → evaluate pipeline.

Extracts the common processing loop that was duplicated across
``cli.py``, ``gui.py``, and ``evaluator.batch_evaluate``.
"""

from __future__ import annotations

import logging
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, asdict
from typing import Optional, Callable, Any

from .glass_reader import load_catalog
from .models import Candidate, Inputs, ThickPrescription
from .cemented import run_cemented
from .spaced import run_spaced
from .thickening import thicken
from .seidel_refine import refine_seidel
from .optiland_bridge.builder import (
    build_optic_from_prescription,
    rx_from_optic,
)
from .optiland_bridge.builder import _build_from_prescription
from .optiland_bridge.optimizer import optimize_optic
from .optiland_bridge.evaluator import evaluate, OpticMetrics

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Pipeline result container
# ---------------------------------------------------------------------------


@dataclass
class PipelineResult:
    """One candidate fully processed through thicken → build → evaluate."""

    candidate: Candidate
    rx: Optional[ThickPrescription]
    metrics: OpticMetrics

    # -- Serialisation ---------------------------------------------------

    def to_dict(self) -> dict[str, Any]:
        """JSON-serialisable flat dict merging candidate, rx, and metrics."""
        d = asdict(self.metrics)

        # Thin-lens radii
        d["thin_radii"] = self.candidate.radii

        # Thick prescription details
        if self.rx is not None:
            d["thick_radii"] = [[e.R_front, e.R_back] for e in self.rx.elements]
            d["t_center"] = [e.t_center for e in self.rx.elements]
            d["t_edge"] = [e.t_edge for e in self.rx.elements]
            d["air_gap"] = self.rx.air_gap
            d["actual_efl"] = self.rx.actual_efl
            d["efl_deviation"] = self.rx.efl_deviation
        else:
            d["thick_radii"] = None
            d["t_center"] = None
            d["t_edge"] = None
            d["air_gap"] = None
            d["actual_efl"] = None
            d["efl_deviation"] = None

        # Thermal metrics (computed at thin-lens synthesis stage)
        th = self.candidate.thermal
        if th is not None:
            d["thermal_available"] = th.thermal_data_available
            d["dn_dT_1"] = th.dn_dT_1
            d["dn_dT_2"] = th.dn_dT_2
            d["V1_ppm_K"] = th.V1 * 1e6 if th.V1 is not None else None
            d["V2_ppm_K"] = th.V2 * 1e6 if th.V2 is not None else None
            d["dphi_dT_norm"] = th.dphi_dT_norm
            d["alpha_housing_required_ppm_K"] = (
                th.alpha_housing_required * 1e6
                if th.alpha_housing_required is not None else None
            )
        else:
            for k in [
                "thermal_available", "dn_dT_1", "dn_dT_2",
                "V1_ppm_K", "V2_ppm_K", "dphi_dT_norm",
                "alpha_housing_required_ppm_K",
            ]:
                d[k] = None

        return d


# ---------------------------------------------------------------------------
# Single-candidate processing
# ---------------------------------------------------------------------------


def process_candidate(
    cand: Candidate,
    inputs: Inputs,
) -> PipelineResult:
    """Thicken → build → evaluate a single candidate.

    Always returns a ``PipelineResult``; on failure the ``metrics.success``
    flag is ``False`` and ``metrics.error_msg`` describes the problem.
    """
    t0 = time.perf_counter()

    try:
        # -- Step 1: Thicken --
        rx = thicken(cand, inputs)
    except Exception as exc:
        m = OpticMetrics(
            glass1=cand.glass1,
            glass2=cand.glass2,
            catalog1=cand.catalog1,
            catalog2=cand.catalog2,
            system_type=cand.system_type,
            PE=cand.PE,
            success=False,
            error_msg=f"thickening error: {type(exc).__name__}: {exc}",
            build_time_ms=(time.perf_counter() - t0) * 1000.0,
        )
        return PipelineResult(cand, None, m)

    if rx is None:
        m = OpticMetrics(
            glass1=cand.glass1,
            glass2=cand.glass2,
            catalog1=cand.catalog1,
            catalog2=cand.catalog2,
            system_type=cand.system_type,
            PE=cand.PE,
            success=False,
            error_msg="thickening failed (sag / geometry rejection)",
            build_time_ms=(time.perf_counter() - t0) * 1000.0,
        )
        return PipelineResult(cand, None, m)

    # -- Step 1.5: Thickness-aware Seidel refinement --
    rx = refine_seidel(rx, cand, inputs)

    # -- Step 2: Build optiland Optic --
    op = build_optic_from_prescription(rx)
    build_ms = (time.perf_counter() - t0) * 1000.0

    if op is None:
        m = OpticMetrics(
            glass1=cand.glass1,
            glass2=cand.glass2,
            catalog1=cand.catalog1,
            catalog2=cand.catalog2,
            system_type=cand.system_type,
            PE=cand.PE,
            success=False,
            error_msg="build_optic returned None",
            build_time_ms=build_ms,
        )
        return PipelineResult(cand, rx, m)

    # -- Step 3: Evaluate (ray tracing + aberrations) --
    m = evaluate(op, cand, inputs, rx=rx)
    m.build_time_ms = build_ms
    return PipelineResult(cand, rx, m)


# ---------------------------------------------------------------------------
# Batch processing
# ---------------------------------------------------------------------------


def run_pipeline(
    candidates: list[Candidate],
    inputs: Inputs,
    *,
    max_n: Optional[int] = None,
    on_progress: Optional[Callable[[int, int, PipelineResult], None]] = None,
) -> list[PipelineResult]:
    """Process a batch of candidates with optional progress callback.

    Parameters
    ----------
    candidates :
        Thin-lens candidates from ``run_cemented`` / ``run_spaced``.
    inputs :
        System specification.
    max_n :
        If given, only process the first *max_n* candidates.
    on_progress :
        ``(current_1based, total, result)`` called after each candidate.

    Returns
    -------
    list[PipelineResult]
        One entry per processed candidate, in order.
    """
    todo = candidates[:max_n] if max_n else candidates
    results: list[PipelineResult] = []

    for i, cand in enumerate(todo):
        result = process_candidate(cand, inputs)
        results.append(result)
        if on_progress is not None:
            on_progress(i + 1, len(todo), result)

    return results


# ---------------------------------------------------------------------------
# Unified entry point: catalog → synthesis → pipeline
# ---------------------------------------------------------------------------


@dataclass
class DesignResult:
    """Bundles synthesis output with pipeline results."""

    candidates: list[Candidate]
    results: list[PipelineResult]
    n_glasses: int
    synth_time_ms: float


# ---------------------------------------------------------------------------
# Stage B: optical-equivalence grouping + thick-lens optimisation
# ---------------------------------------------------------------------------


def _optical_fingerprint(c: Candidate) -> tuple:
    """Group key: (n1, n2, ν1, ν2) with n rounded to 3 decimals, ν to 1.

    Bucket grouping — not a continuous tolerance test. Candidates sharing
    a fingerprint differ in glass properties only at the rounding scale
    (n to 0.001, ν to 0.1) and therefore produce the same Stage B geometry.
    """
    return (round(c.n1, 3), round(c.n2, 3), round(c.nu1, 1), round(c.nu2, 1))


def deduplicate_candidates(
    candidates: list[Candidate],
) -> dict[tuple, list[Candidate]]:
    """Group candidates by optical equivalence.

    Returns a dict mapping fingerprint → list of candidates.
    Stage B optimisation runs once per key (shared geometry); thermal and
    full-spectrum chromatic evaluation uses each candidate's own data.
    """
    groups: dict[tuple, list[Candidate]] = {}
    for c in candidates:
        key = _optical_fingerprint(c)
        groups.setdefault(key, []).append(c)
    return groups


def process_candidate_stage_b(
    cand: Candidate,
    inputs: Inputs,
    rx_a: Optional[ThickPrescription] = None,
) -> PipelineResult:
    """Stage B: build → optimise → evaluate a single candidate.

    Runs ``optimize_optic`` on the Stage-A starting design, then reads back
    the optimised geometry via ``rx_from_optic`` so that ``PipelineResult.rx``
    reflects the actual optimised surface parameters.  Falls back to the
    Stage-A optic (and Stage-A rx) if optimisation fails.

    If *rx_a* is provided (from Stage A), thickening is skipped.
    """
    t0 = time.perf_counter()

    # -- Step 1: Reuse Stage-A prescription or thicken from scratch --
    if rx_a is None:
        try:
            rx_a = thicken(cand, inputs)
        except Exception as exc:
            m = OpticMetrics(
                glass1=cand.glass1, glass2=cand.glass2,
                catalog1=cand.catalog1, catalog2=cand.catalog2,
                system_type=cand.system_type, PE=cand.PE,
                success=False,
                error_msg=f"thickening error: {type(exc).__name__}: {exc}",
                build_time_ms=(time.perf_counter() - t0) * 1000.0,
            )
            return PipelineResult(cand, None, m)

        if rx_a is None:
            m = OpticMetrics(
                glass1=cand.glass1, glass2=cand.glass2,
                catalog1=cand.catalog1, catalog2=cand.catalog2,
                system_type=cand.system_type, PE=cand.PE,
                success=False,
                error_msg="thickening failed (geometry rejection)",
                build_time_ms=(time.perf_counter() - t0) * 1000.0,
            )
            return PipelineResult(cand, None, m)

    build_ms_start = time.perf_counter()

    # -- Step 2: Build Stage-A optic --
    op_a = _build_from_prescription(rx_a)
    if op_a is None:
        build_ms = (time.perf_counter() - build_ms_start) * 1000.0
        m = OpticMetrics(
            glass1=cand.glass1, glass2=cand.glass2,
            catalog1=cand.catalog1, catalog2=cand.catalog2,
            system_type=cand.system_type, PE=cand.PE,
            success=False,
            error_msg="build_optic (Stage A) returned None",
            build_time_ms=build_ms,
        )
        return PipelineResult(cand, rx_a, m)

    # -- Step 3: Optimise (Stage B) — fall back to Stage-A optic on failure --
    op_b = optimize_optic(op_a, rx_a, inputs)
    optimised = op_b is not None
    if not optimised:
        op_b = op_a  # fall back to Stage-A optic

    build_ms = (time.perf_counter() - build_ms_start) * 1000.0

    # -- Step 4: Read back optimised geometry → update prescription --
    try:
        rx_b = rx_from_optic(op_b, rx_a) if optimised else rx_a
    except Exception:
        rx_b = rx_a  # read-back failed; keep Stage-A prescription

    # -- Step 5: Evaluate --
    m = evaluate(op_b, cand, inputs, rx=rx_b)
    m.build_time_ms = build_ms
    return PipelineResult(cand, rx_b, m)


def run_stage_b(
    stage_a_results: list[PipelineResult],
    inputs: Inputs,
    *,
    top_n: int = 10,
    on_progress: Optional[Callable[[int, int, PipelineResult], None]] = None,
) -> list[PipelineResult]:
    """Run Stage B optimisation on the top-N Stage A results.

    Groups successful Stage A results by optical fingerprint
    ``(n1, n2, ν1, ν2)``, runs the optiland optimiser once per unique
    optical group, then returns the optimised results sorted by RMS spot
    radius.

    Parameters
    ----------
    stage_a_results :
        Stage A results in ranking order (e.g. from ``run_pipeline``).
    inputs :
        System specification.
    top_n :
        How many top Stage A candidates to promote to Stage B.
    on_progress :
        ``(current_1based, total, result)`` called after each optimisation.

    Returns
    -------
    list[PipelineResult]
        One entry per unique optical group (≤ top_n), sorted by
        ``rms_spot_radius`` ascending.
    """
    successful = [r for r in stage_a_results if r.metrics.success]
    todo = successful[:top_n]

    if not todo:
        logger.warning("run_stage_b: no successful Stage A results to optimise")
        return []

    # Deduplicate: keep only the best representative per optical group
    groups: dict[tuple, PipelineResult] = {}
    for r in todo:
        key = _optical_fingerprint(r.candidate)
        if key not in groups:
            groups[key] = r  # first occurrence = highest Stage A rank

    group_list = list(groups.values())
    logger.info(
        "Stage B: %d candidates → %d unique optical groups",
        len(todo), len(group_list),
    )

    total = len(group_list)
    results_b: list[PipelineResult] = []

    # Parallel processing — each candidate builds its own optiland Optic
    # internally, so there is no shared state between workers.
    import os
    max_workers = min(total, max(1, os.cpu_count() or 1))

    if max_workers == 1:
        # Single candidate — no overhead from multiprocessing
        for i, stage_a_r in enumerate(group_list):
            r_b = process_candidate_stage_b(
                stage_a_r.candidate, inputs, rx_a=stage_a_r.rx,
            )
            results_b.append(r_b)
            if on_progress is not None:
                on_progress(i + 1, total, r_b)
    else:
        with ProcessPoolExecutor(max_workers=max_workers) as pool:
            future_map = {
                pool.submit(
                    process_candidate_stage_b,
                    stage_a_r.candidate,
                    inputs,
                    stage_a_r.rx,
                ): stage_a_r
                for stage_a_r in group_list
            }
            done_count = 0
            for future in as_completed(future_map):
                done_count += 1
                r_b = future.result()
                results_b.append(r_b)
                if on_progress is not None:
                    on_progress(done_count, total, r_b)

    results_b.sort(
        key=lambda r: r.metrics.rms_spot_radius
        if r.metrics.rms_spot_radius is not None
        else float("inf")
    )

    n_ok = sum(1 for r in results_b if r.metrics.success)
    logger.info("Stage B: %d/%d succeeded", n_ok, len(results_b))
    return results_b


def run_design(
    inputs: Inputs,
    catalog_paths: list[str],
    *,
    on_progress: Optional[Callable[[int, int, PipelineResult], None]] = None,
) -> DesignResult:
    """Load catalogs → thin-lens synthesis → thicken/build/evaluate.

    This is the single high-level function that both CLI and GUI call.

    Parameters
    ----------
    inputs :
        System specification (wavelengths, geometry, constraints, type).
    catalog_paths :
        Paths to ``.agf`` glass catalog files.
    on_progress :
        ``(current_1based, total, result)`` called after each candidate
        has been processed through the pipeline.

    Returns
    -------
    DesignResult
        Contains the raw candidates, pipeline results, glass count
        and synthesis timing.
    """
    # Step 1: Load catalogs
    glasses = load_catalog(catalog_paths)
    n_glasses = len(glasses)
    logger.info("Loaded %d glasses from %d catalog(s)", n_glasses, len(catalog_paths))

    # Step 2: Thin-lens synthesis
    t0 = time.perf_counter()
    if inputs.system_type == "cemented":
        cands = run_cemented(inputs, glasses)
    elif inputs.system_type == "spaced":
        cands = run_spaced(inputs, glasses)
    else:
        raise ValueError(f"Unknown system_type={inputs.system_type!r}")
    synth_ms = (time.perf_counter() - t0) * 1000.0
    logger.info("Synthesis: %d candidates in %.0f ms", len(cands), synth_ms)

    # Step 3-5: Thicken → Build → Evaluate (top-N)
    top_n = inputs.N
    results = run_pipeline(
        cands,
        inputs,
        max_n=top_n,
        on_progress=on_progress,
    )

    n_ok = sum(1 for r in results if r.metrics.success)
    logger.info("Pipeline: %d/%d succeeded", n_ok, len(results))

    return DesignResult(
        candidates=cands,
        results=results,
        n_glasses=n_glasses,
        synth_time_ms=synth_ms,
    )
