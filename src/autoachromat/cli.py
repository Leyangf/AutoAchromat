from __future__ import annotations

import argparse
import json
import time
from dataclasses import asdict
from pathlib import Path
from typing import Optional

from .glass_reader import load_catalog
from .models import Inputs, Candidate
from .cemented import run_cemented
from .spaced import run_spaced
from .thickening import thicken
from .optiland_bridge.builder import build_optic_from_prescription
from .optiland_bridge.evaluator import evaluate, OpticMetrics


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
        import math

        if not math.isfinite(v):
            return "N/A"
    except (TypeError, ValueError):
        return "N/A"
    return f"{v:{fmt}}"


def _print_header() -> None:
    """Print the table header for the results."""
    header = (
        f"{'#':>3}  {'OK':>2}  {'Glass1':<22} {'Glass2':<22}"
        f"  {'EFL':>8} {'FNO':>7}"
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
        f"  {_val(m.efl, '.1f'):>8} {_val(m.fno, '.3f'):>7}"
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

    # ---- Step 1: Load catalogs ----
    glasses = load_catalog(resolved)

    # ---- Step 2: Thin-lens synthesis ----
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

    # ---- Thin-only mode ----
    if args.thin_only:

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

    # ---- Step 3–5: Thicken → Build → Evaluate ----
    print("=" * 90)
    _print_header()

    results: list[dict] = []
    n_ok = 0

    for i, cand in enumerate(todo):
        # Step 3: Thicken
        t0 = time.perf_counter()
        rx = thicken(cand, inputs)
        thicken_ms = (time.perf_counter() - t0) * 1000.0

        if rx is None:
            # Thickening failed (sag rejection etc.)
            m = OpticMetrics(
                glass1=cand.glass1,
                glass2=cand.glass2,
                catalog1=cand.catalog1,
                catalog2=cand.catalog2,
                system_type=cand.system_type,
                PE=cand.PE,
                success=False,
                error_msg="thickening failed (sag / geometry rejection)",
                build_time_ms=thicken_ms,
            )
            _print_row(i, cand, m, None, None, None)
            results.append(_build_result_dict(cand, m, None))
            continue

        # Step 4: Build optiland Optic
        t1 = time.perf_counter()
        op = build_optic_from_prescription(rx)
        build_ms = (time.perf_counter() - t1) * 1000.0 + thicken_ms

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
            _print_row(
                i, cand, m, rx.elements[0].t_center, rx.elements[1].t_center, rx.air_gap
            )
            results.append(_build_result_dict(cand, m, rx))
            continue

        # Step 5: Evaluate (ray tracing + aberrations)
        m = evaluate(op, cand, inputs)
        m.build_time_ms = build_ms

        if m.success:
            n_ok += 1

        _print_row(
            i,
            cand,
            m,
            rx.elements[0].t_center,
            rx.elements[1].t_center,
            rx.air_gap,
        )
        results.append(_build_result_dict(cand, m, rx))

    # ---- Summary ----
    total_ms = sum(
        r.get("build_time_ms", 0) + r.get("eval_time_ms", 0) for r in results
    )
    print("=" * 90)
    print(
        f"\n  {n_ok}/{len(todo)} succeeded  |  Total: {total_ms:.1f} ms  "
        f"({total_ms / max(len(todo), 1):.1f} ms/candidate)\n"
    )

    # ---- JSON export ----
    if args.out:
        Path(args.out).write_text(
            json.dumps(results, indent=2, ensure_ascii=False),
            encoding="utf-8",
        )
        print(f"Wrote {args.out}")


def _build_result_dict(
    cand: Candidate,
    m: OpticMetrics,
    rx,
) -> dict:
    """Merge candidate info, metrics, and prescription into one flat dict."""
    d = asdict(m)

    # Thin-lens radii
    d["thin_radii"] = cand.radii

    # Thick prescription details
    if rx is not None:
        elems = rx.elements
        d["thick_radii"] = [
            [elems[0].R_front, elems[0].R_back],
            [elems[1].R_front, elems[1].R_back],
        ]
        d["t_center"] = [elems[0].t_center, elems[1].t_center]
        d["t_edge"] = [elems[0].t_edge, elems[1].t_edge]
        d["air_gap"] = rx.air_gap
        d["back_focus_guess"] = rx.back_focus_guess
    else:
        d["thick_radii"] = None
        d["t_center"] = None
        d["t_edge"] = None
        d["air_gap"] = None
        d["back_focus_guess"] = None

    return d


if __name__ == "__main__":
    main()
