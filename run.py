#!/usr/bin/env python3
"""
AutoAchromat - Automated achromat doublet lens design.

Main entry point for running the synthesis pipeline.

Usage examples:
    python run.py --agf glass_database/SCHOTT.AGF --system cemented --N 30
    python run.py --agf glass_database/SCHOTT.AGF --agf glass_database/OHARA.agf --system air_spaced --N 20
    python run.py --config configs/example_config.json
    python run.py --config configs/example_config.json --N 50  # CLI overrides config
    python run.py --realray-index 1  # Compute real-ray aberrations for top candidate
    python run.py --help
"""

from __future__ import annotations

import sys
from pathlib import Path

# Add src to path for direct script execution
sys.path.insert(0, str(Path(__file__).parent))

from src import cli, io
from src.synthesis import run


def main() -> int:
    """Main entry point."""
    # Parse CLI arguments
    args = cli.parse_args()

    # Build configuration (defaults → config file → CLI overrides)
    result = cli.build_config(args)
    if result[0] is None:
        print(result[1])  # Error message
        return 1

    cfg, display_count = result

    # Print configuration summary
    cli.print_config_summary(cfg, display_count, args.config)

    # Run synthesis pipeline
    print("Running synthesis pipeline...")
    results, stats = run(cfg)
    print()

    # Print summary
    io.print_summary(stats)
    print()

    # Extract rows from candidates
    rows = io.extract_rows(results)

    # Ensure output directory exists
    io.ensure_out_dir(cfg.out_dir)

    # Save outputs
    io.save_csv(rows, f"{cfg.out_dir}/results_table.csv")
    io.save_report_md(rows, stats, cfg, display_count, f"{cfg.out_dir}/report.md")
    io.save_resolved_config(cfg, display_count, f"{cfg.out_dir}/resolved_config.json")
    print()

    # Print table to console
    print(f"Top {min(display_count, len(rows))} candidates:")
    print()
    io.print_table(rows, max_rows=display_count)

    if not rows:
        print()
        print("No candidates found matching the criteria.")
        return 0

    # Handle real-ray calculation if requested
    realray_index = getattr(args, "realray_index", None)
    if realray_index is not None:
        print()
        _compute_realray_for_index(realray_index, rows, cfg)

    return 0


def _compute_realray_for_index(
    index: int,
    rows: list[dict],
    cfg,
) -> None:
    """
    Compute real-ray aberrations for a specific candidate index.

    Args:
        index: 1-based candidate rank
        rows: List of row dictionaries
        cfg: Configuration object
    """
    # Import raytrace module (lazy import to handle missing Optiland gracefully)
    try:
        from src.raytrace import OPTILAND_AVAILABLE
        from src.raytrace.optiland_interface import compute_realray_for_candidate
    except ImportError as e:
        print(f"ERROR: Failed to import raytrace module: {e}")
        return

    # Check Optiland availability
    if not OPTILAND_AVAILABLE:
        print()
        print("=" * 60)
        print("ERROR: Optiland is not installed")
        print("=" * 60)
        print("Real-ray aberration calculation requires Optiland.")
        print("Install it with: pip install optiland")
        print("=" * 60)
        return

    # Validate index
    if index < 1 or index > len(rows):
        print()
        print(f"ERROR: Invalid --realray-index {index}")
        print(f"       Valid range: 1 to {len(rows)}")
        return

    row_index = index - 1  # Convert to 0-based
    row = rows[row_index]

    print(f"Computing real-ray aberrations for candidate #{index}...")

    # Compute real-ray aberrations
    realray_results = compute_realray_for_candidate(row, cfg)

    # Print results
    io.print_realray_results(index, row, realray_results)

    # Update CSV and report
    csv_path = f"{cfg.out_dir}/results_table.csv"
    report_path = f"{cfg.out_dir}/report.md"

    io.update_csv_row(csv_path, row_index, realray_results)
    io.append_realray_to_report(report_path, index, row, realray_results)


if __name__ == "__main__":
    sys.exit(main())
