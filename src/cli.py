"""
CLI and configuration handling for AutoAchromat.

This module provides:
- Command-line argument parsing
- JSON config file loading
- Configuration resolution (defaults → config → CLI overrides)
- AGF path validation and fallback
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, TYPE_CHECKING

if TYPE_CHECKING:
    from .synthesis import Config


def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns:
        Parsed arguments namespace.
    """
    parser = argparse.ArgumentParser(
        description="Run the AutoAchromat synthesis pipeline with customizable parameters.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python run.py --agf glass_database/SCHOTT.AGF --system cemented --N 30
  python run.py --agf glass_database/SCHOTT.AGF --agf glass_database/OHARA.agf --system air_spaced --N 20
  python run.py --config configs/example_config.json
  python run.py --config configs/example_config.json --N 50  # CLI overrides config
  python run.py --f 150 --D 30 --system cemented
""",
    )

    # Config file (loaded first, then CLI overrides)
    parser.add_argument(
        "--config",
        type=str,
        default=None,
        metavar="PATH",
        help="Path to JSON config file (CLI args override config values)",
    )

    # AGF catalog paths (can specify multiple)
    parser.add_argument(
        "--agf",
        action="append",
        dest="agf_paths",
        metavar="PATH",
        help="Path to AGF catalog file (can specify multiple times)",
    )

    # System type (default None so we can detect if user provided it)
    parser.add_argument(
        "--system",
        choices=["cemented", "air_spaced"],
        default=None,
        help="System type: cemented or air_spaced",
    )

    # Number of candidates (default None)
    parser.add_argument(
        "--N",
        type=int,
        default=None,
        help="Number of top candidates to keep/display (default: 30)",
    )

    # Optical parameters (all default None)
    parser.add_argument(
        "--f", type=float, default=None, help="Focal length in mm (default: 100.0)"
    )
    parser.add_argument(
        "--D", type=float, default=None, help="Aperture diameter in mm (default: 25.0)"
    )
    parser.add_argument(
        "--P0",
        type=float,
        default=None,
        help="Target spherical aberration (default: 0.0)",
    )
    parser.add_argument(
        "--W0", type=float, default=None, help="Target coma (default: 0.0)"
    )
    parser.add_argument(
        "--C0",
        type=float,
        default=None,
        help="Target chromatic aberration (default: 0.0)",
    )

    # Wavelengths (all default None)
    parser.add_argument(
        "--lam0",
        type=float,
        default=None,
        help="Primary wavelength in um (default: 0.58756, d-line)",
    )
    parser.add_argument(
        "--lam1",
        type=float,
        default=None,
        help="Blue wavelength in um (default: 0.48613, F-line)",
    )
    parser.add_argument(
        "--lam2",
        type=float,
        default=None,
        help="Red wavelength in um (default: 0.65627, C-line)",
    )

    # Filtering parameters (default None)
    parser.add_argument(
        "--min_delta_nu",
        type=float,
        default=None,
        help="Minimum Abbe number difference (default: 10.0)",
    )
    parser.add_argument(
        "--max_PE",
        type=float,
        default=None,
        help="Maximum PE threshold (default: 1e10)",
    )

    # allow_repeat: use mutually exclusive group for --allow_repeat / --no_allow_repeat
    repeat_group = parser.add_mutually_exclusive_group()
    repeat_group.add_argument(
        "--allow_repeat",
        action="store_true",
        dest="allow_repeat",
        default=None,
        help="Allow same glass for both elements",
    )
    repeat_group.add_argument(
        "--no_allow_repeat",
        action="store_false",
        dest="allow_repeat",
        help="Disallow same glass for both elements (default)",
    )

    # Air-spaced specific (default None)
    parser.add_argument(
        "--d_air",
        type=float,
        default=None,
        help="Air gap thickness in mm for air-spaced (default: 5.0)",
    )

    # Output (default None)
    parser.add_argument(
        "--out_dir",
        type=str,
        default=None,
        help="Output directory (default: output)",
    )

    # Display options (default None)
    parser.add_argument(
        "--display",
        type=int,
        default=None,
        help="Number of rows to display in console (default: N_keep)",
    )

    return parser.parse_args()


def load_config_file(config_path: str) -> dict[str, Any]:
    """
    Load configuration from a JSON file.

    Args:
        config_path: Path to the JSON config file.

    Returns:
        Dictionary of configuration values.

    Raises:
        FileNotFoundError: If config file doesn't exist.
        json.JSONDecodeError: If config file is invalid JSON.
    """
    path = Path(config_path)
    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")

    with open(path, "r", encoding="utf-8") as f:
        config = json.load(f)

    return config


def _apply_config_to_cfg(cfg: "Config", config: dict[str, Any]) -> int | None:
    """
    Apply configuration dictionary values to Config object.

    Args:
        cfg: The Config object to modify.
        config: Dictionary of configuration values from JSON.

    Returns:
        display value if present in config, else None.
    """
    display = None

    # Map config keys to cfg attributes
    if "agf_paths" in config and config["agf_paths"]:
        cfg.agf_paths = config["agf_paths"]
    if "system_type" in config and config["system_type"]:
        cfg.system_type = config["system_type"]
    if "N_keep" in config and config["N_keep"] is not None:
        cfg.N_keep = int(config["N_keep"])
    if "f" in config and config["f"] is not None:
        cfg.f = float(config["f"])
    if "D" in config and config["D"] is not None:
        cfg.D = float(config["D"])
    if "P0" in config and config["P0"] is not None:
        cfg.P0 = float(config["P0"])
    if "W0" in config and config["W0"] is not None:
        cfg.W0 = float(config["W0"])
    if "C0" in config and config["C0"] is not None:
        cfg.C0 = float(config["C0"])
    if "lam0" in config and config["lam0"] is not None:
        cfg.lam0 = float(config["lam0"])
    if "lam1" in config and config["lam1"] is not None:
        cfg.lam1 = float(config["lam1"])
    if "lam2" in config and config["lam2"] is not None:
        cfg.lam2 = float(config["lam2"])
    if "min_delta_nu" in config and config["min_delta_nu"] is not None:
        cfg.min_delta_nu = float(config["min_delta_nu"])
    if "max_PE" in config and config["max_PE"] is not None:
        cfg.max_PE = float(config["max_PE"])
    if "allow_repeat" in config and config["allow_repeat"] is not None:
        cfg.allow_repeat = bool(config["allow_repeat"])
    if "d_air" in config and config["d_air"] is not None:
        cfg.d_air = float(config["d_air"])
    if "out_dir" in config and config["out_dir"]:
        cfg.out_dir = config["out_dir"]

    # display is harness-specific, not part of Config
    if "display" in config and config["display"] is not None:
        display = int(config["display"])

    return display


def _apply_cli_overrides(cfg: "Config", args: argparse.Namespace) -> None:
    """
    Apply CLI argument overrides to Config object.

    Only applies values that are not None (i.e., explicitly provided by user).

    Args:
        cfg: The Config object to modify.
        args: Parsed CLI arguments.
    """
    if args.agf_paths:
        cfg.agf_paths = args.agf_paths
    if args.system is not None:
        cfg.system_type = args.system
    if args.N is not None:
        cfg.N_keep = args.N
    if args.f is not None:
        cfg.f = args.f
    if args.D is not None:
        cfg.D = args.D
    if args.P0 is not None:
        cfg.P0 = args.P0
    if args.W0 is not None:
        cfg.W0 = args.W0
    if args.C0 is not None:
        cfg.C0 = args.C0
    if args.lam0 is not None:
        cfg.lam0 = args.lam0
    if args.lam1 is not None:
        cfg.lam1 = args.lam1
    if args.lam2 is not None:
        cfg.lam2 = args.lam2
    if args.min_delta_nu is not None:
        cfg.min_delta_nu = args.min_delta_nu
    if args.max_PE is not None:
        cfg.max_PE = args.max_PE
    if args.allow_repeat is not None:
        cfg.allow_repeat = args.allow_repeat
    if args.d_air is not None:
        cfg.d_air = args.d_air
    if args.out_dir is not None:
        cfg.out_dir = args.out_dir


def _validate_agf_paths(cfg: "Config") -> str | None:
    """
    Validate and resolve AGF paths with fallback to defaults.

    Modifies cfg.agf_paths in place.

    Args:
        cfg: The Config object to validate/modify.

    Returns:
        Error message string if validation fails, None if successful.
    """
    # Default catalogs to try if none specified
    default_catalogs = [
        "glass_database/SCHOTT.AGF",
        "glass_database/OHARA.agf",
        "glass_database/CDGM.AGF",
    ]

    # If no paths specified, try defaults
    if not cfg.agf_paths:
        cfg.agf_paths = [p for p in default_catalogs if Path(p).exists()]
        if not cfg.agf_paths:
            return (
                "ERROR: No AGF catalogs found. Please specify with --agf <path> or in config file.\n"
                "       Default locations checked:\n"
                + "\n".join(f"         - {p}" for p in default_catalogs)
            )
        print(f"Using default catalogs: {cfg.agf_paths}")

    # Validate specified paths - warn for missing, continue if at least one exists
    valid_paths = []
    for p in cfg.agf_paths:
        if Path(p).exists():
            valid_paths.append(p)
        else:
            print(f"WARNING: AGF file not found: {p}")

    if not valid_paths:
        return "ERROR: No valid AGF files found."

    cfg.agf_paths = valid_paths
    return None


def build_config(args: argparse.Namespace) -> tuple["Config", int] | tuple[None, str]:
    """
    Build final configuration from defaults, config file, and CLI overrides.

    Resolution order (lowest to highest priority):
    1. default_cfg() - hardcoded defaults
    2. JSON config file values (if --config provided)
    3. CLI argument overrides

    Args:
        args: Parsed CLI arguments.

    Returns:
        Tuple of (Config, display_count) on success.
        Tuple of (None, error_message) on failure.
    """
    # Import here to avoid circular imports
    from .synthesis import default_cfg

    # Step 1: Start with defaults
    cfg = default_cfg()
    config_display: int | None = None

    # Step 2: Load config file if provided
    if args.config:
        try:
            print(f"Loading config from: {args.config}")
            config_data = load_config_file(args.config)
            config_display = _apply_config_to_cfg(cfg, config_data)
        except FileNotFoundError as e:
            return None, f"ERROR: {e}"
        except json.JSONDecodeError as e:
            return None, f"ERROR: Invalid JSON in config file: {e}"

    # Step 3: Apply CLI overrides
    _apply_cli_overrides(cfg, args)

    # Step 4: Validate AGF paths
    error = _validate_agf_paths(cfg)
    if error:
        return None, error

    # Step 5: Determine display count (CLI > config > N_keep)
    if args.display is not None:
        display_count = args.display
    elif config_display is not None:
        display_count = config_display
    else:
        display_count = cfg.N_keep

    return cfg, display_count


def print_config_summary(
    cfg: "Config", display_count: int, config_path: str | None = None
) -> None:
    """
    Print configuration summary to console.

    Args:
        cfg: The resolved Config object.
        display_count: Number of rows to display.
        config_path: Path to config file if used, None otherwise.
    """
    print()
    print("=" * 60)
    print("AutoAchromat Pipeline")
    print("=" * 60)
    if config_path:
        print(f"  Config file:    {config_path}")
    print(f"  System type:    {cfg.system_type}")
    print(f"  Catalogs:       {cfg.agf_paths}")
    print(f"  N_keep:         {cfg.N_keep}")
    print(f"  Focal length:   {cfg.f} mm")
    print(f"  Aperture:       {cfg.D} mm")
    print(f"  P0, W0, C0:     {cfg.P0}, {cfg.W0}, {cfg.C0}")
    print(f"  Wavelengths:    {cfg.lam0}, {cfg.lam1}, {cfg.lam2} um")
    print(f"  min_delta_nu:   {cfg.min_delta_nu}")
    print(f"  max_PE:         {cfg.max_PE}")
    print(f"  allow_repeat:   {cfg.allow_repeat}")
    if cfg.system_type == "air_spaced":
        print(f"  d_air:          {cfg.d_air} mm")
    print(f"  Output dir:     {cfg.out_dir}")
    print(f"  Display rows:   {display_count}")
    print("=" * 60)
    print()
