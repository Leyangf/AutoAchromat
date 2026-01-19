"""
I/O utilities for AutoAchromat.

Provides functions for:
- Directory management
- Candidate extraction and formatting
- Console table printing
- CSV export
- Markdown report export
- Configuration export
"""

from __future__ import annotations

import csv
import dataclasses
import json
import math
from datetime import datetime
from pathlib import Path
from typing import Any, TYPE_CHECKING

if TYPE_CHECKING:
    from .synthesis import Config


# =============================================================================
# Column definitions for output tables
# =============================================================================

# Column definitions: (key, header, width, is_float)
COLUMNS = [
    ("rank", "Rank", 5, False),
    ("system_type", "Type", 10, False),
    ("glass1_name", "Glass1", 14, False),
    ("glass1_catalog", "Cat1", 8, False),
    ("glass1_cost", "Cost1", 7, True),
    ("glass2_name", "Glass2", 14, False),
    ("glass2_catalog", "Cat2", 8, False),
    ("glass2_cost", "Cost2", 7, True),
    ("total_cost", "TotCost", 8, True),
    ("phi1", "phi1", 10, True),
    ("phi2", "phi2", 10, True),
    ("n1", "n1", 8, True),
    ("n2", "n2", 8, True),
    ("nu1", "nu1", 8, True),
    ("nu2", "nu2", 8, True),
    ("delta_nu", "d_nu", 8, True),
    ("X1", "X1", 10, True),
    ("X2", "X2", 10, True),
    ("R1_front", "R1f", 10, True),
    ("R1_back", "R1b", 10, True),
    ("R2_front", "R2f", 10, True),
    ("R2_back", "R2b", 10, True),
    ("R_cemented", "R_cem", 10, True),
    ("P2", "P2", 12, True),
    ("W", "W", 10, True),
    ("R2", "R2", 10, True),
    ("Pi_len", "Pi#", 4, False),
    ("PE", "PE", 14, True),
]


# =============================================================================
# Directory management
# =============================================================================


def ensure_out_dir(path: str) -> str:
    """
    Create the output directory if it doesn't exist.

    Args:
        path: Path to the directory.

    Returns:
        The normalized absolute path as a string.
    """
    out_path = Path(path)
    out_path.mkdir(parents=True, exist_ok=True)
    return str(out_path.resolve())


# =============================================================================
# Formatting helpers
# =============================================================================


def fmt_float(value: Any, precision: int = 6) -> str:
    """
    Format a float value for display/CSV.

    Args:
        value: Value to format (may be None).
        precision: Number of significant digits.

    Returns:
        Formatted string.
    """
    if value is None:
        return ""
    try:
        f = float(value)
        if math.isnan(f) or math.isinf(f):
            return str(f)
        return f"{f:.{precision}g}"
    except (TypeError, ValueError):
        return str(value)


def _fmt_str(value: Any, max_len: int = 20) -> str:
    """
    Format a string value, truncating if needed.

    Args:
        value: Value to format.
        max_len: Maximum length before truncation.

    Returns:
        Formatted string.
    """
    if value is None:
        return ""
    s = str(value)
    if len(s) > max_len:
        return s[: max_len - 2] + ".."
    return s


def _get_glass_attr(candidate: Any, glass_num: int, attr: str) -> Any:
    """
    Get attribute from glass1 or glass2.

    Args:
        candidate: Candidate object.
        glass_num: 1 or 2.
        attr: Attribute name to retrieve.

    Returns:
        Attribute value or None.
    """
    glass_attr = f"g{glass_num}"
    glass = getattr(candidate, glass_attr, None)
    if glass is None:
        return None
    return getattr(glass, attr, None)


# =============================================================================
# Candidate extraction
# =============================================================================


def extract_candidate_row(candidate: Any, rank: int) -> dict[str, Any]:
    """
    Extract all relevant fields from a candidate into a dict.

    Args:
        candidate: Candidate object from synthesis pipeline.
        rank: Ranking position (1-based).

    Returns:
        Dictionary with all candidate fields.
    """
    meta = getattr(candidate, "meta", {}) or {}

    # Glass 1 info
    g1_name = _get_glass_attr(candidate, 1, "name")
    g1_catalog = _get_glass_attr(candidate, 1, "catalog")
    g1_cost = _get_glass_attr(candidate, 1, "relative_cost")

    # Glass 2 info
    g2_name = _get_glass_attr(candidate, 2, "name")
    g2_catalog = _get_glass_attr(candidate, 2, "catalog")
    g2_cost = _get_glass_attr(candidate, 2, "relative_cost")

    # Total cost
    total_cost = None
    if g1_cost is not None and g2_cost is not None:
        try:
            total_cost = float(g1_cost) + float(g2_cost)
        except (TypeError, ValueError):
            pass

    # Meta fields
    phi1 = meta.get("phi1")
    phi2 = meta.get("phi2")
    phi_total = meta.get("phi_total")
    n1 = meta.get("n1")
    n2 = meta.get("n2")
    nu1 = meta.get("nu1")
    nu2 = meta.get("nu2")
    delta_nu = meta.get("delta_nu")
    X1 = meta.get("X1")
    X2 = meta.get("X2")

    # Radii from meta
    R1_front = meta.get("R1_front")
    R1_back = meta.get("R1_back")
    R2_front = meta.get("R2_front")
    R2_back = meta.get("R2_back")
    R_cemented = meta.get("R_cemented")

    # Candidate direct fields
    P2 = candidate.P2
    W = candidate.W
    R2 = candidate.R2
    Pi = candidate.Pi
    Pi_len = len(Pi) if Pi else 0
    PE = candidate.PE

    return {
        "rank": rank,
        "system_type": candidate.system_type,
        "glass1_name": g1_name,
        "glass1_catalog": g1_catalog,
        "glass1_cost": g1_cost,
        "glass2_name": g2_name,
        "glass2_catalog": g2_catalog,
        "glass2_cost": g2_cost,
        "total_cost": total_cost,
        "phi1": phi1,
        "phi2": phi2,
        "phi_total": phi_total,
        "n1": n1,
        "n2": n2,
        "nu1": nu1,
        "nu2": nu2,
        "delta_nu": delta_nu,
        "X1": X1,
        "X2": X2,
        "R1_front": R1_front,
        "R1_back": R1_back,
        "R2_front": R2_front,
        "R2_back": R2_back,
        "R_cemented": R_cemented,
        "P2": P2,
        "W": W,
        "R2": R2,
        "Pi_len": Pi_len,
        "PE": PE,
    }


def extract_rows(candidates: list[Any]) -> list[dict[str, Any]]:
    """
    Extract all candidates into row dictionaries.

    Args:
        candidates: List of Candidate objects.

    Returns:
        List of row dictionaries.
    """
    return [extract_candidate_row(c, i) for i, c in enumerate(candidates, start=1)]


# =============================================================================
# Console output
# =============================================================================


def print_table(rows: list[dict[str, Any]], max_rows: int | None = None) -> None:
    """
    Print a formatted table to console.

    Args:
        rows: List of row dictionaries.
        max_rows: Maximum rows to display (None for all).
    """
    if not rows:
        print("No candidates found.")
        return

    display_rows = rows[:max_rows] if max_rows else rows

    # Print header
    header_parts = []
    for key, header, width, _ in COLUMNS:
        header_parts.append(f"{header:>{width}}")
    print(" ".join(header_parts))
    print("-" * (sum(c[2] for c in COLUMNS) + len(COLUMNS) - 1))

    # Print rows
    for row in display_rows:
        row_parts = []
        for key, _, width, is_float in COLUMNS:
            val = row.get(key)
            if is_float:
                s = fmt_float(val, 4)
            else:
                s = _fmt_str(val, width)
            # Right-align, truncate if needed
            if len(s) > width:
                s = s[: width - 1] + "~"
            row_parts.append(f"{s:>{width}}")
        print(" ".join(row_parts))

    if max_rows and len(rows) > max_rows:
        print(f"... ({len(rows) - max_rows} more rows not shown)")


def print_summary(stats: dict) -> None:
    """
    Print a compact summary to console.

    Args:
        stats: Dictionary with statistics.
    """
    print("=" * 50)
    print("AUTOACHROMAT Summary")
    print("=" * 50)

    # Define display order and labels
    key_labels = [
        ("glasses_total", "Total glasses loaded"),
        ("glasses_filtered", "Glasses after filtering"),
        ("pairs_tested", "Glass pairs tested"),
        ("candidates_generated", "Candidates generated"),
        ("candidates_saved", "Candidates saved"),
        ("best_PE", "Best PE score"),
    ]

    # Print known keys in order
    printed_keys = set()
    for key, label in key_labels:
        if key in stats:
            value = stats[key]
            if isinstance(value, float):
                print(f"  {label}: {value:.6g}")
            else:
                print(f"  {label}: {value}")
            printed_keys.add(key)

    # Print any remaining keys not in the predefined list
    for key, value in stats.items():
        if key not in printed_keys:
            if isinstance(value, float):
                print(f"  {key}: {value:.6g}")
            else:
                print(f"  {key}: {value}")

    print("=" * 50)


# =============================================================================
# CSV export
# =============================================================================


def save_csv(rows: list[dict[str, Any]], csv_path: str) -> None:
    """
    Save rows to CSV file.

    Args:
        rows: List of row dictionaries.
        csv_path: Output file path.
    """
    # Ensure parent directory exists
    Path(csv_path).parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [c[0] for c in COLUMNS]

    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for row in rows:
            # Format values for CSV
            csv_row = {}
            for key, _, _, is_float in COLUMNS:
                val = row.get(key)
                if val is None:
                    csv_row[key] = ""
                elif is_float:
                    csv_row[key] = fmt_float(val, 6)
                else:
                    csv_row[key] = str(val)
            writer.writerow(csv_row)

    print(f"Results saved to: {csv_path}")


def save_candidates_csv(candidates: list[Any], csv_path: str) -> None:
    """
    Save candidates to CSV file (backward compatibility wrapper).

    Extracts candidates to rows and saves using save_csv().

    Args:
        candidates: List of Candidate objects.
        csv_path: Output file path.
    """
    rows = extract_rows(candidates)
    save_csv(rows, csv_path)


# =============================================================================
# Markdown report export
# =============================================================================


def save_report_md(
    rows: list[dict[str, Any]],
    stats: dict[str, Any],
    cfg: "Config",
    display_count: int,
    report_path: str,
) -> None:
    """
    Save a human-readable Markdown report.

    Args:
        rows: List of candidate row dictionaries.
        stats: Statistics dictionary from the pipeline run.
        cfg: The resolved Config object.
        display_count: Number of rows to display.
        report_path: Path to save the report file.
    """
    # Ensure parent directory exists
    Path(report_path).parent.mkdir(parents=True, exist_ok=True)

    lines: list[str] = []

    # -------------------------------------------------------------------------
    # 1) Title + timestamp
    # -------------------------------------------------------------------------
    lines.append("# AutoAchromat Report")
    lines.append("")
    timestamp = datetime.now().isoformat(timespec="seconds")
    lines.append(f"Generated: {timestamp}")
    lines.append("")

    # -------------------------------------------------------------------------
    # 2) Configuration summary
    # -------------------------------------------------------------------------
    lines.append("## Configuration")
    lines.append("")
    lines.append(f"- **system_type**: {cfg.system_type}")
    lines.append(f"- **N_keep**: {cfg.N_keep}")
    lines.append(f"- **f**: {cfg.f}")
    lines.append(f"- **D**: {cfg.D}")
    lines.append(f"- **P0**: {cfg.P0}")
    lines.append(f"- **W0**: {cfg.W0}")
    lines.append(f"- **C0**: {cfg.C0}")
    lines.append(f"- **lam0**: {cfg.lam0}")
    lines.append(f"- **lam1**: {cfg.lam1}")
    lines.append(f"- **lam2**: {cfg.lam2}")
    lines.append(f"- **min_delta_nu**: {cfg.min_delta_nu}")
    lines.append(f"- **max_PE**: {cfg.max_PE}")
    lines.append(f"- **allow_repeat**: {cfg.allow_repeat}")
    if cfg.system_type == "air_spaced":
        lines.append(f"- **d_air**: {cfg.d_air}")
    lines.append("- **agf_paths**:")
    for agf_path in cfg.agf_paths:
        lines.append(f"  - {agf_path}")
    lines.append(f"- **out_dir**: {cfg.out_dir}")
    lines.append("")

    # -------------------------------------------------------------------------
    # 3) Statistics summary
    # -------------------------------------------------------------------------
    lines.append("## Statistics")
    lines.append("")
    lines.append("| Key | Value |")
    lines.append("|-----|-------|")
    for key, value in stats.items():
        # Format value nicely
        if isinstance(value, float):
            formatted_value = fmt_float(value, 6)
        else:
            formatted_value = str(value)
        lines.append(f"| {key} | {formatted_value} |")
    lines.append("")

    # -------------------------------------------------------------------------
    # 4) Results table
    # -------------------------------------------------------------------------
    lines.append("## Results")
    lines.append("")

    display_rows = rows[:display_count] if display_count else rows

    if not display_rows:
        lines.append("No candidates found.")
    else:
        # Build header row
        headers = [header for _, header, _, _ in COLUMNS]
        lines.append("| " + " | ".join(headers) + " |")
        lines.append("| " + " | ".join(["---"] * len(COLUMNS)) + " |")

        # Build data rows
        for row in display_rows:
            row_cells = []
            for key, _, _, is_float in COLUMNS:
                val = row.get(key)
                if val is None:
                    cell = ""
                elif is_float:
                    cell = fmt_float(val, 6)
                else:
                    cell = str(val)
                row_cells.append(cell)
            lines.append("| " + " | ".join(row_cells) + " |")

        if display_count and len(rows) > display_count:
            lines.append("")
            lines.append(f"*({len(rows) - display_count} more rows not shown)*")

    lines.append("")

    # Write to file
    with open(report_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))

    print(f"Report saved to: {report_path}")


# =============================================================================
# Configuration export
# =============================================================================


def save_resolved_config(cfg: "Config", display: int, out_path: str) -> None:
    """
    Save the final resolved configuration to a JSON file.

    Args:
        cfg: The resolved Config object.
        display: The display row count.
        out_path: Path to save the JSON file.
    """
    resolved = {
        "agf_paths": cfg.agf_paths,
        "system_type": cfg.system_type,
        "N_keep": cfg.N_keep,
        "f": cfg.f,
        "D": cfg.D,
        "P0": cfg.P0,
        "W0": cfg.W0,
        "C0": cfg.C0,
        "lam0": cfg.lam0,
        "lam1": cfg.lam1,
        "lam2": cfg.lam2,
        "min_delta_nu": cfg.min_delta_nu,
        "max_PE": cfg.max_PE,
        "allow_repeat": cfg.allow_repeat,
        "d_air": cfg.d_air,
        "out_dir": cfg.out_dir,
        "display": display,
    }

    # Ensure parent directory exists
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)

    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(resolved, f, indent=2)

    print(f"Resolved config saved to: {out_path}")


def dump_config(cfg: Any, path: str) -> None:
    """
    Write configuration to a JSON file (best effort).

    Args:
        cfg: Configuration object (dataclass or dict).
        path: Path to the output JSON file.
    """
    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Convert to dict if dataclass
    if dataclasses.is_dataclass(cfg) and not isinstance(cfg, type):
        data = dataclasses.asdict(cfg)
    elif isinstance(cfg, dict):
        data = cfg
    else:
        # Best effort: use __dict__ if available
        data = getattr(cfg, "__dict__", {"value": str(cfg)})

    # Custom JSON encoder for non-serializable types
    def default_encoder(obj: Any) -> Any:
        if dataclasses.is_dataclass(obj) and not isinstance(obj, type):
            return dataclasses.asdict(obj)
        if hasattr(obj, "__dict__"):
            return obj.__dict__
        return str(obj)

    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, default=default_encoder)
