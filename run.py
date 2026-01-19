#!/usr/bin/env python3
"""
AutoAchromat - Automated achromat doublet lens design.

Main entry point for running the synthesis pipeline.

Usage examples:
    python run.py --agf glass_database/SCHOTT.AGF --system cemented --N 30
    python run.py --agf glass_database/SCHOTT.AGF --agf glass_database/OHARA.agf --system air_spaced --N 20
    python run.py --config configs/example_config.json
    python run.py --config configs/example_config.json --N 50  # CLI overrides config
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


if __name__ == "__main__":
    sys.exit(main())
