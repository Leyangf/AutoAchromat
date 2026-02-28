"""
smoke_test_bridge.py – End-to-end smoke test for the optiland bridge.

Usage:
    python scripts/smoke_test_bridge.py --config config_example.json
    python scripts/smoke_test_bridge.py --config config_example.json --max-n 10

Runs the thin-lens synthesis, then builds + evaluates every candidate
via the shared pipeline and prints a summary table.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

# Ensure the src directory is importable when running as a script
_src = str((Path(__file__).resolve().parent.parent / "src"))
if _src not in sys.path:
    sys.path.insert(0, _src)

from autoachromat.cli import load_inputs
from autoachromat.glass_reader import load_catalog
from autoachromat.cemented import run_cemented
from autoachromat.spaced import run_spaced
from autoachromat.pipeline import run_pipeline, PipelineResult
from autoachromat.optiland_bridge.evaluator import OpticMetrics


def _fmt(v, fmt=".4g") -> str:
    if v is None:
        return "N/A"
    return f"{v:{fmt}}"


def print_summary(results: list[OpticMetrics]) -> None:
    n_total = len(results)
    n_ok = sum(1 for r in results if r.success)
    n_fail = n_total - n_ok

    print("=" * 90)
    print(
        f"  Stage A  evaluation summary:  {n_ok}/{n_total} succeeded,  {n_fail} failed"
    )
    print("=" * 90)
    print()

    # Header
    hdr = (
        f"{'#':>3}  {'OK':>2}  {'Catalog1':>8}:{'Glass1':<12}  "
        f"{'Catalog2':>8}:{'Glass2':<12}  "
        f"{'EFL':>8}  {'FNO':>6}  {'RMS':>8}  {'GEO':>8}  "
        f"{'SA':>9}  {'LchC':>9}  {'TchC':>9}  "
        f"{'build_ms':>8}  {'eval_ms':>8}"
    )
    print(hdr)
    print("-" * len(hdr))

    for i, m in enumerate(results):
        ok_str = "OK" if m.success else "!!"
        line = (
            f"{i:3d}  {ok_str:>2}  {m.catalog1:>8}:{m.glass1:<12}  "
            f"{m.catalog2:>8}:{m.glass2:<12}  "
            f"{_fmt(m.efl):>8}  {_fmt(m.fno):>6}  "
            f"{_fmt(m.rms_spot_radius):>8}  {_fmt(m.geo_spot_radius):>8}  "
            f"{_fmt(m.SA):>9}  {_fmt(m.LchC):>9}  {_fmt(m.TchC):>9}  "
            f"{m.build_time_ms:8.1f}  {m.eval_time_ms:8.1f}"
        )
        print(line)

    if n_fail:
        print()
        print("  FAILED candidates:")
        for i, m in enumerate(results):
            if not m.success:
                print(
                    f"    [{i:3d}] {m.catalog1}:{m.glass1} + "
                    f"{m.catalog2}:{m.glass2} – {m.error_msg}"
                )

    t_total = sum(r.build_time_ms + r.eval_time_ms for r in results)
    print()
    print(
        f"  Total wall time: {t_total:.1f} ms  "
        f"({t_total / max(n_total, 1):.1f} ms/candidate)"
    )
    print()


def main() -> None:
    ap = argparse.ArgumentParser(description="Stage A optiland bridge test")
    ap.add_argument("--config", required=True, help="JSON config path")
    ap.add_argument(
        "--max-n",
        type=int,
        default=None,
        help="Max # candidates to evaluate (default: all)",
    )
    args = ap.parse_args()

    # ---- 1. Load config & glasses ----
    inputs, catalog_paths = load_inputs(args.config)
    config_dir = Path(args.config).resolve().parent
    resolved = []
    for p in catalog_paths:
        pp = Path(p)
        if not pp.is_absolute():
            pp = config_dir / pp
        resolved.append(str(pp))

    glasses = load_catalog(resolved)
    print(f"Loaded {len(glasses)} glasses from {len(resolved)} catalog(s)\n")

    # ---- 2. Run thin-lens synthesis ----
    t0 = time.perf_counter()
    if inputs.system_type == "cemented":
        cands = run_cemented(inputs, glasses)
    elif inputs.system_type == "spaced":
        cands = run_spaced(inputs, glasses)
    else:
        raise ValueError(f"Unknown system_type={inputs.system_type}")
    t_synth = (time.perf_counter() - t0) * 1000.0
    print(f"Synthesis: {len(cands)} candidates in {t_synth:.1f} ms\n")

    if not cands:
        print("No candidates found – nothing to evaluate.")
        return

    # ---- 3. Build + evaluate via pipeline ----
    pipeline_results = run_pipeline(cands, inputs, max_n=args.max_n)
    metrics = [pr.metrics for pr in pipeline_results]

    # ---- 4. Print summary ----
    print_summary(metrics)


if __name__ == "__main__":
    main()
