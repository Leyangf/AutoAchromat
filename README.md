# AutoAchromat

Automated achromatic doublet lens designer. Generates initial cemented or air-spaced doublet structures from glass catalogs, evaluates their optical performance by ray tracing alongside manufacturing-feasibility and thermal checks, and exports designs for further optimization in optical design software.

![Achromatic doublet](example/doublet.png)

**TL;DR**

```bash
git clone https://github.com/Leyangf/AutoAchromat.git
cd AutoAchromat
uv sync
uv run autoachromat-gui
```

## Features

- **Glass catalogs** — reads industry-standard Zemax AGF catalogs (SCHOTT, OHARA, CDGM, …)
- **Exhaustive glass-pair search** — cemented and air-spaced doublets, from visible to LWIR
- **Custom aberration targets** — spherical aberration, coma, and longitudinal chromatic aberration (default: fully corrected)
- **Manufacturable prescriptions** — thin-lens solutions are converted to thick lenses that respect minimum edge/center thicknesses
- **Full evaluation** — spot size, Seidel aberrations, chromatic focal shift, secondary spectrum, thermal stability
- **Numerical refinement** (experimental) — least-squares optimization of the best designs
- **Export** — Zemax `.zmx`, JSON, or CSV

## How It Works

A two-stage pipeline:

- **Stage A — synthesis & screening.** For every glass pair: filter by Abbe-number difference → solve thin-lens powers for the chromatic target → solve lens shapes for the spherical-aberration and coma targets → reject designs with too-steep radii or high PE → thicken to a manufacturable prescription (with focal-length correction) → ray-trace and run thermal analysis.
- **Stage B — optimization (optional).** Refines radii, center thicknesses, and air gap of selected designs with a least-squares optimizer, minimizing spot size while holding the focal length.

![Pipeline overview](example/pipeline.png)

## Installation

This project uses [uv](https://docs.astral.sh/uv/) for environment and dependency management.

```bash
git clone https://github.com/Leyangf/AutoAchromat.git
cd AutoAchromat
uv sync
```

`uv sync` creates a `.venv/`, installs all dependencies from `uv.lock` (`numpy`, `optiland`, `scipy`, `matplotlib`), and installs the project in editable mode. Requires Python >= 3.10; uv downloads a matching interpreter if none is found.

Run commands through `uv run` (e.g. `uv run autoachromat-gui`), or activate the environment first:

```bash
.venv\Scripts\activate        # Windows
source .venv/bin/activate     # macOS / Linux
```

## Quick Start (GUI)

```bash
uv run autoachromat-gui
```

(Or, after activating `.venv`: `autoachromat-gui` / `python -m autoachromat.gui`.)

![Doublet Designer GUI](example/spaced_gui.png)

### Step 1: Set parameters

In the Global panel, configure:

- **f', D** — focal length and aperture diameter in mm (F/# is displayed automatically)
- **Wavelengths** — pick a preset (Visible, NIR, SWIR, MWIR, LWIR) or enter values in nm
- **Field [°]** — half-field angle, used for off-axis evaluation in Stage A and as an optimization target in Stage B
- **Air gap** — fixed during Stage A, becomes a free variable in Stage B (air-spaced doublets only)
- **Catalogs** — load one or more AGF glass catalog files

### Step 2: Run Stage A

Click **▶ Run Stage A** to search all glass-pair combinations. Results appear in a sortable table ranked by aberration performance.

Optional settings:

- **te min / tc min** — minimum edge thickness (positive lens) and center thickness (negative lens). These are auto-filled from the aperture diameter using standard manufacturing tables; override them for special requirements, or click **Auto** to restore the defaults.
- **▶ Advanced** — aberration targets (SA = spherical aberration, Coma, LCA = longitudinal chromatic aberration; default 0 = fully corrected) and screening thresholds (Δν min, PE max). See [Parameters](#parameters).
- **Hide extreme shapes** (air-spaced only) — filters out meniscus designs whose convex surface faces the air gap, which are difficult to manufacture and mount.

### Step 3: Analyze results

Select any row to see:

- Prescription and first-order data
- Seidel aberration coefficients
- Glass properties, transmittance, and thermal analysis
- 2D optical layout drawing
- Chromatic focal shift curve
- Per-surface Seidel bar chart
- Dispersion comparison of the two glasses

### Step 4: Optimize (Stage B, experimental)

Select one or more rows in the results table and click **⚡ Optimize**. This runs a numerical optimizer (optiland LeastSquares) that refines radii, thicknesses, and air gap to minimize spot size while preserving focal length.

> **Note:** Stage B is slow (~10–30 s per candidate). Merit-function weights (W_efl, W_rms, W_field), iteration count, and thickness/air-gap bounds can be adjusted in the Stage B panel.

### Step 5: Export

- **Export .zmx** — export the selected design to Zemax for further optimization
- **Export JSON / CSV** — export the full results table

## Command-Line (Batch) Use

The project also installs an `autoachromat` command for headless runs. It reads a JSON config, runs Stage A (thin-lens synthesis → thickening → ray tracing), prints a ranked results table, and can export the results to JSON.

```bash
uv run autoachromat --config config_example.json
uv run autoachromat --config config_example.json --out results.json   # also write JSON
uv run autoachromat --config config_example.json --thin-only          # synthesis only, skip ray tracing
uv run autoachromat --config config_example.json --max-n 20           # cap candidates evaluated
```

See [`config_example.json`](config_example.json) for a complete config. Catalog paths are resolved relative to the config file.

### Parameters

| Key | Meaning | Unit | GUI default |
|---|---|---|---|
| `system_type` | `"cemented"` or `"spaced"` | — | — |
| `fprime`, `D` | Focal length, aperture diameter | mm | — |
| `lam1`, `lam0`, `lam2` | Short, reference, and long wavelength | **µm** in the config (the GUI uses nm) | F, d, C lines |
| `P0` | Spherical aberration (SA) target | — | 0 |
| `W0` | Coma target | — | 0 |
| `C0` | Longitudinal chromatic aberration (LCA) target | — | 0 |
| `min_delta_nu` | Minimum Abbe-number difference Δν between the two glasses | — | 10 |
| `max_PE` | Upper limit on PE, the per-surface aberration load (a proxy for higher-order aberrations); higher-PE designs are discarded | — | 100 |
| `N` | Number of cemented candidates to keep (ignored for spaced) | — | — |
| `air_gap` | Air gap between the elements (spaced only) | mm | 1.0 |
| `catalogs` | List of AGF catalog paths | — | — |

> **Note:** The CLI covers Stage A only. Stage B is available exclusively in the GUI.

## Glass Catalogs

The repo includes SCHOTT, OHARA, and CDGM catalogs in `data/catalogs/`. To add other vendors, copy the corresponding `.AGF` file (e.g. `HOYA.AGF`) from your Zemax installation into `data/catalogs/`.

## Limitations

- **Stage B** is experimental, slow, and available only in the GUI.
- **`.zmx` export** writes glasses as refractive index / Abbe number only, without the full Sellmeier dispersion data. Re-assign the catalog glasses in Zemax for accurate dispersion.

## Project Structure

```
src/autoachromat/
├── gui.py              Tkinter GUI (autoachromat-gui)
├── cli.py              Command-line entry point (autoachromat)
├── pipeline.py         Shared synthesize → thicken → evaluate pipeline
├── cemented.py         Cemented doublet synthesis
├── spaced.py           Air-spaced doublet synthesis
├── optics.py           Core optics: refractive index, achromat powers, radius checks
├── thickening.py       Thin → thick lens conversion with focal-length correction
├── seidel_refine.py    Thickness-aware Seidel refinement
├── thermal.py          Thermal analysis of glass pairs
├── glass_reader.py     AGF catalog loader
├── zemax_export.py     .zmx export
├── models.py           Data classes (Inputs, Candidate, ThickPrescription, …)
└── optiland_bridge/    Ray tracing and optimization via optiland
    ├── builder.py      Prescription → optiland Optic
    ├── evaluator.py    Optical metrics extraction
    └── optimizer.py    Stage B least-squares optimizer
```

## Running Tests

```bash
uv run pytest
```

## Acknowledgements

This work was inspired by:

> D. H. Nguyen and A. V. Bakholdin, "Calculation automation of 'cemented doublet
> with a singlet' and 'singlet with a cemented doublet' objectives," *Journal of
> Optical Technology* **90**(1), 6–13 (2023). https://doi.org/10.1364/JOT.90.000006
