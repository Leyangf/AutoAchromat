# AutoAchromat

AutoAchromat is a Python tool for enumerating achromatic two-element lens candidates from AGF glass catalogs.

It currently supports two design modes:

- `cemented`: cemented doublet (`R1, R2, R3`)
- `spaced`: separated doublet (`R1, R2, R3, R4`)

The solver reads glass data from AGF catalogs (for example SCHOTT, OHARA, CDGM), computes refractive index and Abbe number at configured wavelengths, applies constraints, and ranks accepted candidates.

## Requirements

- Python `>= 3.10`
- `numpy`

## Installation

```bash
pip install -e .
```

After installation, the CLI entry point is:

```bash
autoachromat
```

You can also run as a module:

```bash
python -m autoachromat.cli
```

## Quick Start

Run with catalogs declared in your JSON config:

```bash
autoachromat --config config_example.json
```

Write full accepted candidates to JSON:

```bash
autoachromat --config config_example.json --out results.json
```

## Config File

`--config` points to a JSON file. Example: [config_example.json](config_example.json).

Supported fields:

- `catalogs`: required list of AGF paths.
- `lam0`, `lam1`, `lam2`: wavelengths in micrometers.
- `D`: clear aperture (same unit system as `fprime`, typically mm).
- `fprime`: effective focal length target.
- `C0`, `P0`, `W0`: model target constants.
- `min_delta_nu`: minimum `|nu1 - nu2|` for a glass pair.
- `max_PE`: reject candidates with `PE > max_PE`.
- `N`: Top-N for cemented mode. (`spaced` currently does not use Top-N truncation.)
- `system_type`: `"cemented"` or `"spaced"`.
- `eps`: optional numeric tolerance, default `1e-12`.
- `root_imag_tol`: optional complex-root tolerance, default `1e-9`.

Path resolution behavior:

- If catalog paths are relative, they are resolved relative to the config file directory.
- If `catalogs` is missing or empty, CLI exits with an error.

## Output Behavior

CLI output includes:

- `Accepted: <count>`: number of candidates kept after filtering
- Top 20 formatted candidates printed to stdout
- Optional full JSON output via `--out`

Ranking differs by mode:

- `cemented`: keeps Top-N by smallest `|W - W0|`, then sorts by `(|W-W0|, PE)`
- `spaced`: keeps all accepted candidates, then sorts by `PE`

## Why Some Top Results Look Identical

This is expected in many runs:

- Different catalog entries can have effectively identical optical properties at (`lam0`, `lam1`, `lam2`).
- Those entries can produce the same `Q` (or `(Q1, Q2)`), radii, and `PE`.
- The CLI prints only the first 20 results, so repeated-equivalent solutions often appear grouped.

## Project Layout

```text
src/autoachromat/
  cli.py          # command-line entry
  glass_reader.py # AGF parser
  optics.py       # refractive index, Abbe, shared optics utilities
  cemented.py     # cemented solver
  spaced.py       # spaced solver
  models.py       # Inputs and Candidate data models
```

## Development Notes

- AGF parser supports common dispersion formula IDs used by catalogs.
- Glass filtering excludes invalid dispersion entries and out-of-range wavelength coverage.
- Candidate fields include optional glass relative costs from AGF `OD` records.
