# AutoAchromat

**Automated achromat doublet lens synthesis engine** — Design and rank optical doublet lenses using analytical thin-lens theory with alpha-parameter surface invariant formulas.

![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue) ![License](https://img.shields.io/badge/license-MIT-green) ![Status](https://img.shields.io/badge/status-Active-brightgreen)

---

## 📋 Table of Contents

- [Features](#features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Configuration](#configuration)
- [Architecture](#architecture)
- [For Users](#for-users)
- [For Developers](#for-developers)
- [Theory](#theory)
- [References](#references)

---

## ✨ Features

- **Two system types**: Cemented and air-spaced doublets
- **Batch processing**: Automatically generate and rank designs
- **Flexible I/O**: JSON config files + CLI overrides
- **Multiple catalogs**: Support SCHOTT, OHARA, CDGM AGF files
- **Rich output**: CSV results, Markdown reports, JSON configs
- **Zero dependencies**: Standard library only (Python 3.10+)
- **Reusable core**: Use `src.synthesis.run(cfg)` from scripts, GUIs, or tools

---

## 📦 Installation

```bash
git clone https://github.com/Leyangf/AutoAchromat.git
cd AutoAchromat
```

**Requirements:** Python 3.10+ only (no pip packages needed)

---

## 🚀 Quick Start

### Basic Usage

```bash
# Run with defaults (cemented doublet, top 30 designs)
python run.py --system cemented --N 30

# Air-spaced system with custom focal length
python run.py --system air_spaced --f 150 --D 30 --N 50

# Using a config file
python run.py --config configs/example_config.json

# Override config file with CLI arguments
python run.py --config configs/example_config.json --N 50 --display 20

# Multiple glass catalogs
python run.py --agf glass_database/SCHOTT.AGF --agf glass_database/OHARA.agf --system air_spaced --N 20
```

### Output

After running, find results in `output/`:
- **`results_table.csv`** — All optical parameters (28 columns)
- **`report.md`** — Human-readable Markdown report with config + stats + results table
- **`resolved_config.json`** — Final merged configuration used for the run

---

## 🔧 Configuration

### Configuration Priority

**Lowest → Highest:**
1. **Hardcoded defaults** in `src/synthesis.py`
2. **JSON config file** (if `--config` provided)
3. **CLI arguments** (highest priority)

### Configuration File (JSON)

Create `configs/my_design.json`:

```json
{
  "agf_paths": ["glass_database/SCHOTT.AGF", "glass_database/OHARA.agf"],
  "system_type": "cemented",
  "N_keep": 30,
  "f": 100.0,
  "D": 25.0,
  "P0": 0.0,
  "W0": 0.0,
  "C0": 0.0,
  "lam0": 0.58756,
  "lam1": 0.48613,
  "lam2": 0.65627,
  "min_delta_nu": 10.0,
  "max_PE": 1e10,
  "allow_repeat": false,
  "d_air": 5.0,
  "out_dir": "output",
  "display": 10
}
```

Then run:
```bash
python run.py --config configs/my_design.json
```

### CLI Arguments Reference

| Argument | Type | Description | Default |
|----------|------|-------------|---------|
| `--config PATH` | str | JSON config file | - |
| `--agf PATH` | str | AGF catalog (repeatable) | `glass_database/SCHOTT.AGF` |
| `--system` | enum | `cemented` or `air_spaced` | `cemented` |
| `--N` | int | Top N candidates to keep | 30 |
| `--f` | float | Focal length (mm) | 100.0 |
| `--D` | float | Aperture diameter (mm) | 25.0 |
| `--P0, --W0, --C0` | float | Target aberrations | 0.0 |
| `--lam0, --lam1, --lam2` | float | Wavelengths (μm) | 0.58756, 0.48613, 0.65627 |
| `--min_delta_nu` | float | Min Abbe number difference | 10.0 |
| `--max_PE` | float | Max PE threshold | 1e10 |
| `--allow_repeat` / `--no_allow_repeat` | bool | Same glass both elements? | False |
| `--d_air` | float | Air gap for air-spaced (mm) | 5.0 |
| `--out_dir` | str | Output directory | `output` |
| `--display` | int | Rows to show in console | N_keep |

**Example:** Override config file's N_keep and display:
```bash
python run.py --config configs/example_config.json --N 50 --display 20
```

---

## 🏗️ Architecture

```
run.py (56 lines) ──────────→ Thin entrypoint
    │
    ├─→ src/cli.py ────────→ CLI parsing + config resolution
    │     • parse_args()
    │     • build_config()
    │     • print_config_summary()
    │
    ├─→ src/io.py ─────────→ Data extraction + export
    │     • extract_rows()
    │     • print_table()
    │     • save_csv()
    │     • save_report_md()
    │
    └─→ src/synthesis.py ──→ Core synthesis pipeline (1000+ lines)
          • run(cfg) → (results, stats)
          • solve_C() - achromatization
          • solve_PW() - aberration solving
          • compute_pe() - candidate scoring
```

### Project Structure

```
AutoAchromat/
├── run.py                    # Main CLI entry point
├── README.md                 # This file
├── configs/
│   └── example_config.json   # Example configuration
├── glass_database/           # AGF glass catalogs
│   ├── SCHOTT.AGF
│   ├── OHARA.agf
│   └── CDGM.AGF
├── src/                      # Python package
│   ├── __init__.py
│   ├── cli.py               # CLI argument parsing & config handling
│   ├── glass_reader.py      # AGF file parsing
│   ├── io.py                # Data extraction & export (CSV, MD, JSON)
│   ├── score.py             # PE scoring and candidate ranking
│   └── synthesis.py         # Core synthesis algorithm (1000+ lines)
├── tests/                   # Unit tests
│   ├── __init__.py
│   └── test_io.py
└── reference_paper/         # Academic references
    ├── Ivanova et al. (2017)
    └── Nguyen (2022)
```

---

## 👥 For Users

### Running Design Sweeps

Generate multiple designs by varying a parameter:

```bash
# Focal length sweep: 80mm, 100mm, 120mm
for f in 80 100 120; do
  python run.py --f $f --N 20 --out_dir output/f_${f}mm
done
```

### Comparing System Types

```bash
# Compare cemented vs air-spaced
python run.py --system cemented --N 30 --out_dir output/cemented
python run.py --system air_spaced --N 30 --out_dir output/air_spaced
```

### Custom Glass Combinations

Add your own AGF file:
```bash
python run.py --agf glass_database/SCHOTT.AGF --agf my_glasses.agf --N 50
```

---

## 👨‍💻 For Developers

### Using the Core API

```python
from src.synthesis import Config, default_cfg, run
from src import io

# Create or modify config
cfg = default_cfg()
cfg.system_type = "air_spaced"
cfg.N_keep = 100
cfg.f = 150.0

# Run synthesis
results, stats = run(cfg)

# Process results
rows = io.extract_rows(results)
io.save_csv(rows, "my_results.csv")
io.save_report_md(rows, stats, cfg, 10, "my_report.md")
```

### Config Resolution in Code

```python
from src.cli import parse_args, build_config

# Parse CLI args
args = parse_args()

# Build final config (defaults → config file → CLI overrides)
cfg, display_count = build_config(args)
if cfg is None:
    print(f"Error: {display_count}")
    exit(1)
```

### Testing

```bash
# Run unit tests
python tests/test_io.py
```

### Module Overview

| Module | Responsibility | Key Functions |
|--------|-----------------|---------------|
| `synthesis.py` | Core algorithm | `run()`, `solve_C()`, `solve_PW()`, `compute_pe()` |
| `glass_reader.py` | AGF parsing | `read_agf()`, `load_catalog()` |
| `score.py` | PE scoring | `Candidate`, `TopNContainer` |
| `io.py` | Data I/O | `extract_rows()`, `save_csv()`, `save_report_md()` |
| `cli.py` | CLI & config | `parse_args()`, `build_config()` |

---

## 📐 Theory

### Cemented Doublet

Two lens elements bonded at a shared interface (R₂):

```
         R₁        R₂        R₃
    ←─────────────────────────────→
    |      n₁      |      n₂      |
    | Lens 1       | Lens 2       |
    └──────────────┴──────────────┘
         cemented interface
```

**Scoring:** PE = |P₂| · 3^|W| / R₂²

Where:
- **P₂** = Spherical aberration surface invariant (cemented interface)
- **W** = Total system coma
- **R₂** = Cemented surface radius

### Air-Spaced Doublet

Two separate lenses with air gap:

```
     R₁   R₂  air  R₃   R₄
    ←────────┼───────────────→
    |    n₁    |    n₂      |
    | Lens 1 | | Lens 2     |
    └────────┘ └────────────┘
    4 surfaces instead of 3
```

**Scoring:** PE = mean(|Pᵢ|) for i = 1..4

Where:
- **Pᵢ** = Spherical aberration surface invariant (surface i)

### Achromatization

For perfect achromatic correction (C=0):

```
φ₁/ν₁ + φ₂/ν₂ = 0

Solution:
φ₁ = φ · ν₁ / (ν₁ - ν₂)
φ₂ = -φ · ν₂ / (ν₁ - ν₂)
```

Where **ν** is the Abbe number (dispersion).

---

## 📚 References

- **Ivanova et al. (2017)** — *Computer tool for achromatic and aplanatic cemented doublet design and analysis* [PDF in `reference_paper/`]
- **Nguyen (2022)** — *Automation of Synthesis and Ranking of Cemented and Air-Spaced Doublets* [PDF in `reference_paper/`]

---

## 🐛 Troubleshooting

### No AGF files found
**Error:** `ERROR: No AGF catalogs found.`

**Solution:** Ensure AGF files exist in `glass_database/` or specify explicit paths:
```bash
python run.py --agf /full/path/to/SCHOTT.AGF
```

### JSON config not found
**Error:** `ERROR: Config file not found: configs/my_config.json`

**Solution:** Check file path and ensure file exists:
```bash
ls -la configs/
```

### No valid candidates generated
**Problem:** Output shows `Candidates saved: 0`

**Solution:**
- Increase `--max_PE` threshold
- Reduce `--min_delta_nu` constraint
- Check that glass catalogs cover required wavelength range
- Try `--allow_repeat` for more combinations

---

## 📝 License

MIT License — See LICENSE file for details.

---

## 🤝 Contributing

Contributions welcome! Please ensure:
1. No changes to physics/math (synthesize.py algorithm)
2. Standard library only (no new pip packages)
3. Tests pass: `python tests/test_io.py`
4. Code is documented with docstrings

