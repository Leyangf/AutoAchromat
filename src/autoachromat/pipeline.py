"""
pipeline.py – Shared thicken → build → evaluate pipeline.

Extracts the common processing loop that was duplicated across
``cli.py``, ``gui.py``, and ``evaluator.batch_evaluate``.
"""

from __future__ import annotations

import logging
import time
from dataclasses import dataclass, asdict
from typing import Optional, Callable, Any

from .glass_reader import load_catalog
from .models import Candidate, Inputs, ThickPrescription
from .cemented import run_cemented
from .spaced import run_spaced
from .thickening import thicken
from .optiland_bridge.builder import build_optic_from_prescription
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
    m = evaluate(op, cand, inputs)
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
    n_glasses = (
        sum(len(v) for v in glasses.values())
        if isinstance(glasses, dict)
        else len(glasses)
    )
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
