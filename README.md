# AutoAchromat

Automated achromatic doublet lens designer. Generates cemented and air-spaced doublet initial structures from glass catalogs, evaluates optical performance via ray tracing, and exports designs for further optimization in Zemax or other optical design software.

## What It Does

1. Reads industry-standard AGF glass catalogs (SCHOTT, OHARA, CDGM, etc.)
2. Automatically searches all valid glass pair combinations
3. Produces thick-lens prescriptions with manufacturing-feasible thicknesses
4. Evaluates each design: spot size, Seidel aberrations, chromatic errors, secondary spectrum, thermal stability
5. Exports the best candidates as `.zmx` (Zemax), JSON, or CSV

## Installation

```bash
git clone <repo-url>
cd AutoAchromat
pip install -e .
```

Requires Python >= 3.10. Dependencies (`numpy`, `optiland`, `scipy`, `matplotlib`) are installed automatically.

## Quick Start

```bash
python -m autoachromat.gui
```

### Step 1: Set Parameters

Configure focal length, aperture, wavelengths (visible or IR presets available), and glass catalogs in the Global panel.

### Step 2: Run Stage A

The system automatically sets minimum edge/center thicknesses (te min, tc min) based on the aperture diameter according to standard manufacturing tables. You can override these values for specific manufacturing requirements, or click **Auto** to reset to the recommended defaults.

Click **Run Stage A** to search all glass pair combinations. Results appear in a sortable table ranked by aberration performance.

For air-spaced doublets, check **Hide extreme shapes** to filter out meniscus designs with the convex surface facing the air gap, which are difficult to manufacture and mount.

### Step 3: Analyze Results

Select any row to see:
- Prescription & first-order data
- Seidel aberration coefficients
- Glass properties, transmittance, and thermal analysis
- 2D optical layout drawing
- Chromatic focal shift curve
- Per-surface Seidel bar chart
- Dispersion comparison of the two glasses

### Step 4: Optimize (Stage B, experimental)

Select one or more rows in the results table and click **Optimize**. This runs a numerical optimizer (optiland LeastSquares) that refines radii, thicknesses, and air gap to minimize spot size while preserving focal length.

> **Note:** Stage B is slow (~10-30s per candidate) and still experimental. Merit function weights (W_efl, W_rms, W_field) and iteration count can be adjusted in the Stage B panel.

### Step 5: Export

- **Export .zmx** — export selected design to Zemax for further optimization
- **Export JSON / CSV** — export full results table

## Examples

### Cemented Doublet

![Cemented doublet example](example/cemented_example.png)

### Air-Spaced Doublet

![Air-spaced doublet example](example/spaced_example.png)

## Glass Catalogs

The repo includes OHARA and CDGM catalogs in `data/catalogs/`. To add Schott glasses, copy `SCHOTT.AGF` from your Zemax installation to `data/catalogs/`.

## Running Tests

```bash
pytest
```
