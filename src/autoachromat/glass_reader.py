"""
AGF (ANSI Glass Format) reader for Zemax glass catalogs.

This module provides functionality to parse AGF files containing glass data
from manufacturers like SCHOTT, OHARA, and CDGM.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class Glass:
    """Representation of a single glass material from an AGF catalog."""

    # Required fields
    name: str
    catalog: str

    # Dispersion formula
    formula_id: int | None = None
    cd: list[float] = field(default_factory=list)

    # Wavelength limits (micrometers)
    ld_min_um: float | None = None
    ld_max_um: float | None = None

    # NM line fields
    exclude_sub: bool | None = None
    status: int | None = None
    melt_freq: int | None = None

    # ED line field
    dpgf: float | None = None

    # Thermal expansion (ED line, stored as ppm/K in AGF file)
    cte_m40_20: float | None = None  # CTE -40 to +20°C [ppm/K]
    cte_20_300: float | None = None  # CTE +20 to +300°C [ppm/K]

    # Schott thermal dispersion constants (TD line)
    td_D0: float | None = None  # [1/K]
    td_D1: float | None = None  # [1/K²]
    td_D2: float | None = None  # [1/K³]
    td_E0: float | None = None  # [µm²/K]
    td_E1: float | None = None  # [µm²/K²]
    td_ltk: float | None = None  # thermal knee wavelength [µm]
    td_T_ref: float | None = None  # reference temperature [°C]

    # IT lines: internal transmittance data [(wavelength_um, transmittance, thickness_mm)]
    transmission: list[tuple[float, float, float]] = field(default_factory=list)

    # OD line field
    relative_cost: float | None = None

    # Reference values (not used in calculations)
    nd_ref: float | None = None
    vd_ref: float | None = None


def _parse_float(value: str) -> float | None:
    """Parse a float value, returning None for missing values ('_')."""
    value = value.strip()
    if value in {"_", ""}:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def _parse_int(value: str) -> int | None:
    """Parse an integer value, returning None for missing values ('_')."""
    value = value.strip()
    if value in {"_", ""}:
        return None
    try:
        return int(value)
    except ValueError:
        # Try parsing as float first, then convert to int
        try:
            return int(float(value))
        except ValueError:
            return None


def _parse_bool(value: str) -> bool | None:
    """Parse a boolean value from '0' or '1', returning None for missing."""
    value = value.strip()
    if value in {"_", ""}:
        return None
    try:
        val = int(float(value))
        return val == 1
    except ValueError:
        return None


def read_agf(path: str) -> tuple[str, list[Glass]]:
    """
    Read an AGF (ANSI Glass Format) file and return parsed glass data.

    Args:
        path: Path to the AGF file.

    Returns:
        A tuple of (catalog_comment, list_of_glass_records).
        catalog_comment is from the CC line (may be empty).
    """
    glasses: list[Glass] = []
    catalog_comment = ""

    # Derive catalog name from filename
    catalog_name = Path(path).stem.upper()

    current_glass: Glass | None = None

    # Try different encodings commonly used in AGF files
    encodings = ["utf-8", "latin-1", "cp1252"]
    content: str | None = None

    for encoding in encodings:
        try:
            with open(path, "r", encoding=encoding) as f:
                content = f.read()
            break
        except UnicodeDecodeError:
            continue

    if content is None:
        raise ValueError(f"Could not decode file {path} with any known encoding")

    lines = content.splitlines()

    for line in lines:
        line = line.strip()
        if not line:
            continue

        # Split on whitespace, handling multiple spaces
        parts = line.split()
        if not parts:
            continue

        tag = parts[0].upper()

        if tag == "CC":
            # Catalog comment - join remaining parts
            catalog_comment = " ".join(parts[1:]) if len(parts) > 1 else ""

        elif tag == "NM":
            # Start of new glass record
            # NM <name> <formula#> <MIL#> <N(d)> <V(d)> <ExcludeSub> <status> <melt_freq>
            if current_glass is not None:
                glasses.append(current_glass)

            name = parts[1] if len(parts) > 1 else ""
            formula_id = _parse_int(parts[2]) if len(parts) > 2 else None
            # parts[3] is MIL# - we skip it
            nd_ref = _parse_float(parts[4]) if len(parts) > 4 else None
            vd_ref = _parse_float(parts[5]) if len(parts) > 5 else None
            exclude_sub = _parse_bool(parts[6]) if len(parts) > 6 else None
            status = _parse_int(parts[7]) if len(parts) > 7 else None
            melt_freq = _parse_int(parts[8]) if len(parts) > 8 else None

            current_glass = Glass(
                name=name,
                catalog=catalog_name,
                formula_id=formula_id,
                nd_ref=nd_ref,
                vd_ref=vd_ref,
                exclude_sub=exclude_sub,
                status=status,
                melt_freq=melt_freq,
            )

        elif tag == "ED" and current_glass is not None:
            # ED line: TCE(-40 to 20), TCE(20 to 300), density, dPgF, ignore_thermal_exp, ...
            if len(parts) > 1:
                current_glass.cte_m40_20 = _parse_float(parts[1])  # ppm/K
            if len(parts) > 2:
                current_glass.cte_20_300 = _parse_float(parts[2])  # ppm/K
            # parts[3] = density (unused)
            if len(parts) > 4:
                current_glass.dpgf = _parse_float(parts[4])

        elif tag == "CD" and current_glass is not None:
            # CD line: dispersion coefficients.
            # IMPORTANT: positions are meaningful (K1,L1,K2,L2,... for
            # Sellmeier).  A '_' placeholder means "unknown" — we store
            # 0.0 to preserve positional alignment.  Glasses with missing
            # coefficients will later produce wrong n(λ) and be rejected
            # by the n-range sanity check in filter_glasses().
            coeffs: list[float] = []
            for val in parts[1:]:
                parsed = _parse_float(val)
                coeffs.append(parsed if parsed is not None else 0.0)
            current_glass.cd = coeffs

        elif tag == "LD" and current_glass is not None:
            # LD line: wavelength limits (min, max) in micrometers
            if len(parts) > 1:
                current_glass.ld_min_um = _parse_float(parts[1])
            if len(parts) > 2:
                current_glass.ld_max_um = _parse_float(parts[2])

        elif tag == "OD" and current_glass is not None:
            # OD line: rel_cost, CR, FR, SR, AR, PR
            # First value is relative cost; -1 or '_' means missing
            if len(parts) > 1:
                cost = _parse_float(parts[1])
                if cost is not None and cost == -1:
                    cost = None
                current_glass.relative_cost = cost

        elif tag == "TD" and current_glass is not None:
            # TD line: D0 D1 D2 E0 E1 lambda_tk T_ref
            if len(parts) > 1:
                current_glass.td_D0 = _parse_float(parts[1])
            if len(parts) > 2:
                current_glass.td_D1 = _parse_float(parts[2])
            if len(parts) > 3:
                current_glass.td_D2 = _parse_float(parts[3])
            if len(parts) > 4:
                current_glass.td_E0 = _parse_float(parts[4])
            if len(parts) > 5:
                current_glass.td_E1 = _parse_float(parts[5])
            if len(parts) > 6:
                current_glass.td_ltk = _parse_float(parts[6])
            if len(parts) > 7:
                current_glass.td_T_ref = _parse_float(parts[7])

        elif tag == "IT" and current_glass is not None:
            # IT line: wavelength(µm) transmittance thickness(mm)
            if len(parts) >= 4:
                lam = _parse_float(parts[1])
                trans = _parse_float(parts[2])
                thick = _parse_float(parts[3])
                if lam is not None and trans is not None and thick is not None:
                    current_glass.transmission.append((lam, trans, thick))

        # Ignore other tags: GC, BD, MD, etc.

    # Don't forget the last glass
    if current_glass is not None:
        glasses.append(current_glass)

    return catalog_comment, glasses


def load_catalog(paths: list[str]) -> list[Glass]:
    """
    Load multiple AGF files and return a merged list of glasses.

    The catalog name for each glass is derived from the file stem
    (e.g., "SCHOTT" from "SCHOTT.AGF").

    Duplicates (same name across catalogs) are kept as separate entries.

    Args:
        paths: List of paths to AGF files.

    Returns:
        A merged list of all Glass records from all catalogs.
    """
    all_glasses: list[Glass] = []

    for path in paths:
        _, glasses = read_agf(path)
        all_glasses.extend(glasses)

    return all_glasses


if __name__ == "__main__":
    # Self-test: load catalogs from ../glass_database relative to this file
    script_dir = Path(__file__).parent
    glass_db_dir = script_dir.parent / "glass_database"

    catalog_files = [
        glass_db_dir / "SCHOTT.AGF",
        glass_db_dir / "OHARA.agf",
        glass_db_dir / "CDGM.AGF",
    ]

    # Check files exist
    for f in catalog_files:
        if not f.exists():
            print(f"Warning: {f} not found")

    # Load all catalogs
    all_glasses = load_catalog([str(f) for f in catalog_files if f.exists()])

    # Count glasses per catalog
    catalog_counts: dict[str, int] = {}
    catalog_glasses: dict[str, list[Glass]] = {}

    for glass in all_glasses:
        catalog_counts[glass.catalog] = catalog_counts.get(glass.catalog, 0) + 1
        if glass.catalog not in catalog_glasses:
            catalog_glasses[glass.catalog] = []
        catalog_glasses[glass.catalog].append(glass)

    print("=" * 60)
    print("AGF Reader Self-Test")
    print("=" * 60)
    print()
    print("Number of glasses per catalog:")
    for catalog, count in sorted(catalog_counts.items()):
        print(f"  {catalog}: {count} glasses")

    print()
    print("First 3 glasses from each catalog:")
    for catalog in sorted(catalog_glasses.keys()):
        print(f"\n  {catalog}:")
        for glass in catalog_glasses[catalog][:3]:
            print(
                f"    {glass.name}: "
                f"(formula_id={glass.formula_id}, "
                f"ld_min_um={glass.ld_min_um}, "
                f"ld_max_um={glass.ld_max_um}, "
                f"relative_cost={glass.relative_cost})"
            )

    print()
    print("=" * 60)
