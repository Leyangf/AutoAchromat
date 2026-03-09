# AutoAchromat

Automated achromatic doublet designer — from glass catalog to ray-traced evaluation in one pipeline.

Reads AGF glass catalogs (SCHOTT, OHARA, CDGM, …), enumerates two-element achromatic lens candidates via thin-lens theory, converts them to thick-lens prescriptions, builds ray-tracing models with [optiland](https://github.com/HarrisonKramer/optiland), and evaluates optical performance (spot size, Seidel aberrations, chromatic errors).

Supports two design modes:

| Mode | Surfaces | Ranking |
|------|----------|---------|
| `cemented` | R₁, R₂, R₃ (shared cement interface) | Top-N by \|W − W₀\|, then PE |
| `spaced` | R₁, R₂, R₃, R₄ (air gap between elements) | All accepted, sorted by PE |

## Requirements

- Python ≥ 3.10 (tested with 3.11)
- `numpy`, `optiland` (≥ 0.5.8), `scipy`

## Installation

### Option A: pip (recommended)

```bash
git clone <repo-url>
cd AutoAchromat
pip install -e .
```

This installs `numpy` and `optiland` (which pulls in `scipy`, `matplotlib`, `numba`, etc.).

### Option B: conda

```bash
git clone <repo-url>
cd AutoAchromat
conda env create -f environment.yml
conda activate autoachromat
```

## Glass Catalogs

The repo ships with **OHARA** and **CDGM** catalogs in `data/catalogs/`. The Schott catalog (`SCHOTT.AGF`) is not included due to redistribution restrictions.

To add Schott glasses:
1. Obtain `SCHOTT.AGF` from your Zemax/OpticStudio installation (typically under `Glasscat/`) or from the [Schott website](https://www.schott.com)
2. Copy it to `data/catalogs/SCHOTT.AGF`
3. Select it in the GUI catalog browser or reference it in your config JSON

The pipeline works with any combination of catalogs — you can run with OHARA or CDGM alone.

## Usage

### GUI

```bash
autoachromat-gui
```

Tkinter desktop application providing:

- **Input panel** — edit all optical parameters (f', D, wavelengths, C₀/P₀/W₀, constraints) and load/save JSON configs
- **Catalog browser** — add/remove AGF catalog files
- **One-click run** — executes the full pipeline (synthesis → thickening → build → evaluate) in a background thread with progress bar
- **Sortable results table** — columns for glass pair, EFL, F/#, RMS/GEO spot, SA, LchC, TchC, PE, thermal α_h; click any header to sort
- **Detail panel** — select a row to see full prescription, Seidel coefficients, thermal metrics, and 2D lens drawing
- **Export** — save results as JSON or CSV

### CLI

```bash
autoachromat --config config_example.json                  # full pipeline
autoachromat --config config_example.json --thin-only      # thin-lens only (no ray tracing)
autoachromat --config config_example.json --out results.json --max-n 50
```

| Flag | Description |
|------|-------------|
| `--config` | Path to JSON config (required) |
| `--out` | Export results to JSON |
| `--max-n` | Max candidates to evaluate (default: config `N`) |
| `--thin-only` | Skip thickening / ray tracing |

## Config File

See [config_example.json](config_example.json). Key fields:

| Field | Description |
|-------|-------------|
| `catalogs` | List of AGF file paths (relative to config dir) |
| `lam0`, `lam1`, `lam2` | Design wavelengths [µm] |
| `D` | Clear aperture diameter [mm] |
| `fprime` | Target effective focal length [mm] |
| `C0`, `P0`, `W0` | Aberration target constants |
| `min_delta_nu` | Minimum Abbe number difference for glass pairs |
| `max_PE` | Maximum PE threshold |
| `N` | Top-N candidates to keep (cemented) |
| `system_type` | `"cemented"` or `"spaced"` |
| `air_gap` | Air gap between elements [mm] (spaced only) |

## Pipeline

```
AGF files ──▶ glass_reader ──▶ filter ──▶ thin-lens synthesis ──▶ thickening ──▶ optiland build ──▶ evaluate
                                          (cemented | spaced)     (Table 10-2/3    (Optic model)    (spot, Seidel,
                                                                   + power corr.)                    chromatic)
```

1. **Glass loading** — parse AGF catalogs, compute n(λ) via 12 dispersion formulas, compute Abbe numbers, filter invalid/out-of-range glasses
2. **Thin-lens synthesis** — enumerate all glass pairs with |Δν| ≥ threshold, solve for shape parameters Q (cemented) or (Q₁, Q₂) (spaced), compute radii and aberration figures of merit
3. **Thickening** — assign center/edge thicknesses per standard tables, iteratively correct radii to preserve thin-lens power
4. **Build** — translate thick prescription into optiland `Optic` (aperture, field, wavelengths, image solve)
5. **Evaluate** — extract EFL, F/#, RMS/geometric spot radius, Seidel coefficients (SA, CC, AC, PC, DC), and chromatic aberrations (LchC, TchC)

## Running Tests

```bash
pytest
```

## Project Layout

```text
src/autoachromat/
  models.py           # data structures: Inputs, Candidate, ElementRx, ThickPrescription
  glass_reader.py     # AGF catalog parser (NM/CD/LD/ED/OD/TD tags)
  optics.py           # refractive index, Abbe number, achromat power split, radius checks
  cemented.py         # cemented doublet thin-lens solver
  spaced.py           # air-spaced doublet thin-lens solver
  thickening.py       # thin → thick lens conversion (sag, thickness tables, radius correction)
  thermal.py          # thermal analysis: dn/dT, thermo-optical coefficients, housing CTE
  pipeline.py         # unified thicken → build → evaluate pipeline
  cli.py              # CLI entry point
  gui.py              # Tkinter GUI
  optiland_bridge/
    builder.py        # ThickPrescription → optiland Optic
    evaluator.py      # ray tracing + aberration extraction → OpticMetrics
data/catalogs/        # AGF glass catalog files (OHARA, CDGM; add SCHOTT manually)
tests/                # pytest test suite
config_example.json   # example configuration
```
