from __future__ import annotations

import argparse
import json
from pathlib import Path

from .glass_reader import load_catalog
from .models import Inputs
from .cemented import run_cemented
from .spaced import run_spaced


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
    catalog_paths: list[str] = d.get("catalogs", [])
    return inputs, catalog_paths


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True, help="JSON config path")
    ap.add_argument(
        "--catalog",
        nargs="+",
        default=None,
        help="AGF catalog paths (overrides 'catalogs' in JSON)",
    )
    ap.add_argument("--out", default="", help="Optional output json path")
    args = ap.parse_args()

    inputs, catalog_paths = load_inputs(args.config)

    # --catalog CLI flag takes priority; otherwise use JSON 'catalogs'
    if args.catalog:
        catalog_paths = args.catalog
    if not catalog_paths:
        ap.error(
            "No catalog paths provided. Set 'catalogs' in the JSON config "
            "or pass --catalog on the command line."
        )

    # Resolve paths relative to the config file directory
    config_dir = Path(args.config).resolve().parent
    resolved = []
    for p in catalog_paths:
        pp = Path(p)
        if not pp.is_absolute():
            pp = config_dir / pp
        resolved.append(str(pp))

    glasses = load_catalog(resolved)

    if inputs.system_type == "cemented":
        cands = run_cemented(inputs, glasses)
    elif inputs.system_type == "spaced":
        cands = run_spaced(inputs, glasses)
    else:
        raise ValueError(f"Unknown system_type={inputs.system_type}")

    print(f"Accepted: {len(cands)}")

    def _cost_str(cost) -> str:
        return f"{cost:.1f}" if cost is not None else "N/A"

    for i, c in enumerate(cands[:20]):
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

    if args.out:
        out = [c.__dict__ for c in cands]
        Path(args.out).write_text(json.dumps(out, indent=2), encoding="utf-8")
        print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
