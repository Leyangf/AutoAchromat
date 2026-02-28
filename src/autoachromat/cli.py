from __future__ import annotations

import argparse
import json
import math
import time
from dataclasses import replace
from pathlib import Path
from typing import Optional

from .models import Inputs, Candidate
from .pipeline import run_design, PipelineResult
from .optiland_bridge.evaluator import OpticMetrics


def load_inputs(path: str) -> tuple[Inputs, list[str]]:
    """Return (Inputs, catalog_paths) from a JSON config file."""
    d = json.loads(Path(path).read_text(encoding="utf-8"))
    inputs = Inputs(
        lam0=d["lam0"],
        lam1=d["lam1"],
        lam2=d["lam2"],
        D=d["D"],
        fprime=d["fprime"],
        C0=d["C0"],
        P0=d["P0"],
        W0=d["W0"],
        min_delta_nu=d["min_delta_nu"],
        max_PE=d["max_PE"],
        N=d.get("N", 50),
        system_type=d["system_type"],
        air_gap=d.get("air_gap", 1.0),
        eps=d.get("eps", 1e-12),
        root_imag_tol=d.get("root_imag_tol", 1e-9),
    )
    catalog_paths = d.get("catalogs")
    if not isinstance(catalog_paths, list) or not catalog_paths:
        raise ValueError(
            "Config must include a non-empty 'catalogs' list with AGF paths."
        )
    return inputs, catalog_paths


# ---------------------------------------------------------------------------
# Formatting helpers
# ---------------------------------------------------------------------------


def _val(v: Optional[float], fmt: str = ".4f") -> str:
    """Format a nullable float, returning 'N/A' for None / NaN."""
    if v is None:
        return "N/A"
    try:
        if not math.isfinite(v):
            return "N/A"
    except (TypeError, ValueError):
        return "N/A"
    return f"{v:{fmt}}"


def _print_header() -> None:
    """Print the table header for the results."""
    header = (
        f"{'#':>3}  {'OK':>2}  {'Glass1':<22} {'Glass2':<22}"
        f"  {'EFL':>8} {'FNO':>7} {'BFD':>8}"
        f"  {'t1':>6} {'t2':>6} {'gap':>5}"
        f"  {'RMS':>8} {'GEO':>8}"
        f"  {'SA':>10} {'LchC':>10} {'TchC':>10}"
        f"  {'build':>7} {'eval':>7}"
    )
    print(header)
    print("-" * len(header))


def _print_row(
    idx: int,
    cand: Candidate,
    m: OpticMetrics,
    rx_t1: Optional[float],
    rx_t2: Optional[float],
    rx_gap: Optional[float],
) -> None:
    """Print one result row."""
    ok = "OK" if m.success else "!!"
    g1 = f"{cand.catalog1}:{cand.glass1}"
    g2 = f"{cand.catalog2}:{cand.glass2}"
    print(
        f"{idx:3d}  {ok:>2}  {g1:<22} {g2:<22}"
        f"  {_val(m.efl, '.1f'):>8} {_val(m.fno, '.3f'):>7} {_val(m.bfd, '.1f'):>8}"
        f"  {_val(rx_t1, '.2f'):>6} {_val(rx_t2, '.2f'):>6} {_val(rx_gap, '.2f'):>5}"
        f"  {_val(m.rms_spot_radius, '.2f'):>8} {_val(m.geo_spot_radius, '.2f'):>8}"
        f"  {_val(m.SA, '.3f'):>10} {_val(m.LchC, '.3f'):>10} {_val(m.TchC, '.3f'):>10}"
        f"  {m.build_time_ms:7.1f} {m.eval_time_ms:7.1f}"
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    ap = argparse.ArgumentParser(
        description="AutoAchromat – thin-lens synthesis → thickening → ray tracing"
    )
    ap.add_argument("--config", required=True, help="JSON config path")
    ap.add_argument("--out", default="", help="Optional output JSON path")
    ap.add_argument(
        "--max-n",
        type=int,
        default=0,
        help="Max candidates to evaluate (0 = use config N)",
    )
    ap.add_argument(
        "--thin-only",
        action="store_true",
        help="Only run thin-lens synthesis (no thickening / ray tracing)",
    )
    args = ap.parse_args()

    try:
        inputs, catalog_paths = load_inputs(args.config)
    except ValueError as e:
        ap.error(str(e))

    # Resolve catalog paths relative to config file directory
    config_dir = Path(args.config).resolve().parent
    resolved = []
    for p in catalog_paths:
        pp = Path(p)
        if not pp.is_absolute():
            pp = config_dir / pp
        resolved.append(str(pp))

    # ---- Thin-only mode: still needs direct synthesis ----
    if args.thin_only:
        from .glass_reader import load_catalog
        from .cemented import run_cemented
        from .spaced import run_spaced

        glasses = load_catalog(resolved)
        t_synth = time.perf_counter()
        if inputs.system_type == "cemented":
            cands = run_cemented(inputs, glasses)
        elif inputs.system_type == "spaced":
            cands = run_spaced(inputs, glasses)
        else:
            raise ValueError(f"Unknown system_type={inputs.system_type}")
        synth_ms = (time.perf_counter() - t_synth) * 1000.0
        print(f"\nSynthesis: {len(cands)} candidates in {synth_ms:.1f} ms\n")

        max_n = args.max_n if args.max_n > 0 else inputs.N
        todo = cands[:max_n]

        def _cost_str(cost: Optional[float]) -> str:
            return f"{cost:.1f}" if cost is not None else "N/A"

        for i, c in enumerate(todo):
            cost_info = f"cost1={_cost_str(c.cost1)} cost2={_cost_str(c.cost2)}"
            if c.system_type == "cemented":
                print(
                    f"[{i:02d}] {c.catalog1}:{c.glass1} + {c.catalog2}:{c.glass2} "
                    f"Q={c.Q or 0.0:.6g} W={c.W or 0.0:.6g} PE={c.PE or 0.0:.6g} "
                    f"{cost_info} R={['%.3g' % r for r in c.radii]}"
                )
            else:
                print(
                    f"[{i:02d}] {c.catalog1}:{c.glass1} + {c.catalog2}:{c.glass2} "
                    f"Q1={c.Q1 or 0.0:.6g} Q2={c.Q2 or 0.0:.6g} PE={c.PE or 0.0:.6g} "
                    f"{cost_info} R={['%.3g' % r for r in c.radii]}"
                )
        return

    # ---- Full pipeline via run_design() ----
    print("=" * 90)
    _print_header()

    results_dicts: list[dict] = []
    n_ok = 0

    def _on_progress(current: int, total: int, pr: PipelineResult) -> None:
        nonlocal n_ok
        cand = pr.candidate
        m = pr.metrics
        rx = pr.rx
        if m.success:
            n_ok += 1
        rx_t1 = rx.elements[0].t_center if rx else None
        rx_t2 = rx.elements[1].t_center if rx else None
        rx_gap = rx.air_gap if rx else None
        _print_row(current - 1, cand, m, rx_t1, rx_t2, rx_gap)
        results_dicts.append(pr.to_dict())

    if args.max_n > 0:
        inputs = replace(inputs, N=args.max_n)

    dr = run_design(inputs, resolved, on_progress=_on_progress)

    print(f"\nSynthesis: {len(dr.candidates)} candidates in {dr.synth_time_ms:.1f} ms")

    # ---- Summary ----
    total_ms = sum(
        r.get("build_time_ms", 0) + r.get("eval_time_ms", 0) for r in results_dicts
    )
    print("=" * 90)
    print(
        f"\n  {n_ok}/{len(dr.results)} succeeded  |  Total: {total_ms:.1f} ms  "
        f"({total_ms / max(len(dr.results), 1):.1f} ms/candidate)\n"
    )

    # ---- JSON export ----
    if args.out:
        Path(args.out).write_text(
            json.dumps(results_dicts, indent=2, ensure_ascii=False),
            encoding="utf-8",
        )
        print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
