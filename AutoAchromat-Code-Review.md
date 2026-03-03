---
title: AutoAchromat Comprehensive Code Review Report
date: 2026-03-02
tags: [optics, lens-design, code-review, achromat]
---

# AutoAchromat — Comprehensive Code Review Report

> Generated: 2026-03-02
> Based on line-by-line reading of actual source code, with all formulas and code cross-checked.

---

## Table of Contents

1. [[#1. Project Overview and Design Goals]]
2. [[#2. Repository Directory Structure]]
3. [[#3. Mathematical Background and Physical Foundations]]
4. [[#4. Module Details]]
5. [[#5. Data Flow and Call Chain]]
6. [[#6. Test Coverage]]
7. [[#7. Module Dependency Relationships]]
8. [[#8. Design Decisions and Review Points]]

---

## 1. Project Overview and Design Goals

**AutoAchromat** is an automated achromatic doublet lens design system, developed for master's thesis research.

### Functional Goals

1. Read industrial glass catalogs (AGF format, supporting SCHOTT / OHARA / CDGM)
2. Enumerate all glass pairs and identify achromatic doublet solutions via thin-lens Seidel aberration theory
3. Convert thin-lens solutions into manufacturable thick-lens prescriptions
4. Evaluate actual optical performance through the optiland ray tracing engine
5. Perform first-order thermal analysis (passive athermalization) for each design
6. Present results via CLI / Tkinter GUI

### Design Philosophy

| Principle | Implementation |
|-----------|---------------|
| Separation of pure math and optics | `thickening.py` has no external dependencies, only the Python standard library |
| Defensive programming | All functions return `None` on failure; exceptions are not propagated upward |
| Separation of concerns | Data structures, optical synthesis, thickening, and ray tracing are each independent modules |
| Testability | Unit tests + contract tests + integration tests, totaling 1,434 lines of test code |

---

## 2. Repository Directory Structure

```
c:\MasterThesis\AutoAchromat\
├── README.md
├── pyproject.toml              # Python project configuration
├── pytest.ini
├── config_example.json         # Example configuration file
├── results.json                # 30.6 KB historical output (committed to git)
├── data/
│   └── catalogs/
│       ├── SCHOTT.AGF          # 480 KB
│       ├── OHARA.agf           # 691 KB
│       └── CDGM.AGF            # 621 KB
├── reference paper/            # Reference literature
├── scripts/
│   └── smoke_test_bridge.py    # End-to-end integration test script
├── src/
│   └── autoachromat/
│       ├── __init__.py         # Public API exports (34 lines)
│       ├── models.py           # Dataclasses (161 lines)
│       ├── optics.py           # Thin-lens and aberration mathematics (244 lines)
│       ├── glass_reader.py     # AGF file parsing (304 lines)
│       ├── cemented.py         # Cemented doublet synthesis (185 lines)
│       ├── spaced.py           # Air-spaced doublet synthesis (435 lines)
│       ├── thickening.py       # Thin→thick conversion (658 lines)
│       ├── thermal.py          # Thermal analysis (256 lines)
│       ├── pipeline.py         # Pipeline coordination (279 lines)
│       ├── cli.py              # Command-line interface (221 lines)
│       ├── gui.py              # Tkinter desktop GUI (1,019 lines)
│       └── optiland_bridge/
│           ├── __init__.py     # Module exports (18 lines)
│           ├── builder.py      # Thick prescription → optiland Optic (~280 lines)
│           └── evaluator.py    # Optical metric extraction (~350 lines)
└── tests/
    ├── test_contracts.py       # API contract tests (237 lines)
    ├── test_refactoring.py     # Refactoring integrity tests (382 lines)
    ├── test_thermal.py         # Thermal analysis tests (393 lines)
    └── test_thickness.py       # Thick-lens calculation tests (422 lines)
```

**Code Size Statistics**

| Category | Files | Lines of Code |
|----------|-------|--------------|
| Main source | 13 | ~3,800 |
| optiland_bridge | 3 | ~650 |
| Tests | 4 | 1,434 |
| Scripts | 1 | ~150 |
| **Total** | **21** | **~6,034** |

---

## 3. Mathematical Background and Physical Foundations

### 3.1 Three-Color Design Wavelengths (Fraunhofer Spectral Lines)

$$\lambda_0 = 0.58756\ \mu m \quad \text{(d-line, yellow, primary design wavelength)}$$

$$\lambda_1 = 0.48613\ \mu m \quad \text{(F-line, blue, short wave)}$$

$$\lambda_2 = 0.65627\ \mu m \quad \text{(C-line, red, long wave)}$$

### 3.2 Refractive Index Dispersion Formulas (12 AGF Formats)

| formula_id | Name | Formula |
|:---:|---|---|
| 1 | Schott | $n^2 = A_0 + A_1\lambda^2 + A_2\lambda^{-2} + A_3\lambda^{-4} + A_4\lambda^{-6} + A_5\lambda^{-8}$ |
| 2 | Sellmeier 1 | $n^2 = 1 + \frac{K_1\lambda^2}{\lambda^2-L_1} + \frac{K_2\lambda^2}{\lambda^2-L_2} + \frac{K_3\lambda^2}{\lambda^2-L_3}$ |
| 3 | Herzberger | $L=\frac{1}{\lambda^2-0.028}$, $n = A_0+A_1L+A_2L^2+A_3\lambda^2+A_4\lambda^4[+A_5\lambda^6]$ |
| 4 | Sellmeier 2 | $n^2 = 1+A+\frac{B\lambda^2}{\lambda^2-C}+\frac{D\lambda^2}{\lambda^2-E}$ |
| 5 | Conrady | $n = A_0 + A_1/\lambda + A_2/\lambda^{3.5}$ |
| 6 | Sellmeier 3 | Sellmeier 1 extended, up to 4 terms |
| 7, 8 | Handbook 1&2 | $n^2 = 1+A_0+\frac{A_1}{\lambda^2-A_2}+\frac{A_3}{\lambda^2-A_4}[+A_5\lambda^2]$ |
| 9~12 | Extended/Sellmeier 5 | ⚠️ Code falls back to Sellmeier 1 format |

> [!warning] Review Point
> `formula_id` 9~12 all fall into the `else` branch and are computed using the Sellmeier 1 format. If the coefficient format for these IDs in the actual catalog differs, **silently incorrect refractive indices will be returned without any error**.

### 3.3 Abbe Number

$$\nu = \frac{n(\lambda_0) - 1}{n(\lambda_1) - n(\lambda_2)}$$

Returns $\infty$ when the denominator $< 10^{-10}$ (code is safe; physically meaningful, confirmation needed).

### 3.4 Achromatic Condition and Normalized Power

**General form used in this code** (`optics.py:180–191`):

$$\varphi_1 = \frac{\nu_1(1 - \nu_2 C_0)}{\nu_1 - \nu_2}, \qquad \varphi_2 = 1 - \varphi_1$$

- When $C_0 = 0$, reduces to the standard achromatic condition $\varphi_1/\nu_1 + \varphi_2/\nu_2 = 0$
- Always guarantees $\varphi_1 + \varphi_2 = 1$ (normalized total system power)

> [!note] Mathematical Review
> When $C_0 \neq 0$, $\varphi_1 + \varphi_2 = 1$ still holds (since $\varphi_2 = 1 - \varphi_1$), but the achromatic condition becomes the generalized form $\varphi_1/\nu_1 + \varphi_2/\nu_2 = C_0/(\nu_1\nu_2)$. The physical meaning should be confirmed against the reference literature.

### 3.5 Thin-Lens Seidel Third-Order Aberration Theory

#### Cemented Doublet (`cemented.py`)

**Spherical aberration minimization equation** (quadratic equation in shape parameter $Q$, `cemented.py:23–41`):

$$AQ^2 + BQ + (C - P_0) = 0$$

$$A = \left(1 + \frac{2}{n_1}\right)\varphi_1 + \left(1 + \frac{2}{n_2}\right)\varphi_2$$

$$B = \frac{3}{n_1-1}\varphi_1^2 - \frac{3}{n_2-1}\varphi_2^2 - 2\varphi_2$$

$$C = \frac{n_1}{(n_1-1)^2}\varphi_1^3 + \frac{n_2}{(n_2-1)^2}\varphi_2^3 + \left(\frac{1}{n_2-1}+1\right)\varphi_2^2$$

The two roots correspond to two bending forms (positive and negative Q values).

**Coma coefficient W** (`cemented.py:44–50`):

$$W = KQ + L, \qquad K = \frac{A+1}{2}, \qquad L = \frac{B - \varphi_2}{3}$$

**Radii of curvature for three surfaces** (`cemented.py:53–69`):

$$R_1 = \frac{f'}{\dfrac{n_1}{n_1-1}\varphi_1 + Q}, \quad R_2 = \frac{f'}{\varphi_1 + Q}, \quad R_3 = \frac{f'}{\dfrac{n_2}{n_2-1}\varphi_1 + Q - \dfrac{1}{n_2-1}}$$

**PE penalty function** (`cemented.py:72–85`):

$$u_2 = Q\!\left(1 - \frac{1}{n_1}\right) + \varphi_1, \quad u_3 = Q\!\left(1 - \frac{1}{n_2}\right) + \varphi_1$$

$$P_2 = \left[\frac{u_3-u_2}{1/n_2 - 1/n_1}\right]^2 \cdot \left(\frac{u_3}{n_2} - \frac{u_2}{n_1}\right)$$

$$PE = \frac{|P_2| \cdot 3^{|W|}}{R_2^2}$$

> [!warning] Review Point
> The physical origin and literature basis for the $3^{|W|}$ factor need to be confirmed. When $|W| > 3$, PE grows very rapidly ($3^3 = 27$, $3^4 = 81$), which has a dominant effect on candidate filtering.

**Candidate filtering strategy**: Maintain a Top-N heap, with primary sort key $|W - W_0|$ and secondary sort key $PE$.

#### Air-Spaced Doublet (`spaced.py`)

The two elements each have independent bending parameters $Q_1$ and $Q_2$.

$$\begin{aligned} A_1Q_1^2+B_1Q_1+A_2Q_2^2+B_2Q_2+C&=P_0\\ K_1Q_1+K_2Q_2+L&=W_0 \end{aligned}$$

**Seidel coefficients** (`spaced.py:54–70`):

$$A_1 = \varphi_1^3\!\left(1 + \frac{2}{n_1}\right), \quad B_1 = \varphi_1^3 \cdot \frac{3}{n_1-1}$$

$$A_2 = \varphi_2^3\!\left(1 + \frac{2}{n_2}\right), \quad B_2 = \varphi_2^3 \cdot \frac{3}{n_2-1} - 4\varphi_1\varphi_2^2\!\left(1 + \frac{1}{n_2}\right)$$

$$C=\varphi_1^3\cdot\frac{n_1}{(n_1-1)^2} +\varphi_2^3\cdot\frac{n_2}{(n_2-1)^2} -\varphi_1\varphi_2^2\left(4\cdot\frac{1}{n_2-1}+1\right) +\varphi_1^2\varphi_2\left(3+\frac{2}{n_2}\right)$$

$$K_1 = \varphi_1^2\!\left(1 + \frac{1}{n_1}\right), \quad K_2 = \varphi_2^2\!\left(1 + \frac{1}{n_2}\right)$$

$$L = \frac{\varphi_1^2}{n_1-1} + \frac{\varphi_2^2}{n_2-1} - \varphi_1\varphi_2\!\left(2 + \frac{1}{n_2}\right)$$

**Zero-coma constraint line** (linear equation):

$$K_1 Q_1 + K_2 Q_2 + (L - W_0) = 0 \implies Q_2 = -\frac{K_1}{K_2}Q_1 + \frac{W_0-L}{K_2}$$

**Finding the spherical aberration minimum on the constraint line** (substitute to eliminate $Q_2$, reducing to a quadratic in $Q_1$):

$$q_2 Q_1^2 + q_1 Q_1 + q_0 = 0, \quad \text{solved using numpy.roots()}$$

**Radii of curvature for each surface** (`spaced.py:252–274`):

$$f_1 = f'/\varphi_1, \quad f_2 = f'/\varphi_2$$

$$R_1 = \frac{f_1}{n_1/(n_1-1)+Q_1},\quad R_2 = \frac{f_1}{1+Q_1},\quad R_3 = \frac{f_2}{n_2/(n_2-1)+Q_2},\quad R_4 = \frac{f_2}{1+Q_2}$$

**Adjacent surface geometric interference check** (`spaced.py:374–375`):

$$\text{air\_gap} + \text{sag}(R_3, a) - \text{sag}(R_2, a) \geq 0, \quad \text{sag}(R, a) = R - \text{sign}(R)\sqrt{R^2-a^2}$$

**Root-finding search strategy** (algebraic + sweep):

| Path | Method | Result Characteristics |
|------|--------|----------------------|
| Algebraic path | `numpy.roots()` exact solution | $S_I = 0$ (spherical aberration strictly minimized), but may fall in a geometrically forbidden region |
| Sweep path | $Q_1 \in [-3,3]$ with 61 equally spaced points | Geometrically feasible; returns optimal point + two boundary points (up to 3) |

Two paths are deduplicated ($|\Delta Q_1| < 0.05$), then merged and sorted by $(|S_I\text{ residual}|,\ \sum|P_\text{surfaces}|)$.

> [!warning] Review Point
> The air-spaced doublet has **no Top-N filtering**. When using SCHOTT×SCHOTT (~400 glass pairs), potentially tens of thousands of candidates may be generated; `max_PE` is the only filtering layer.

### 3.6 Thin→Thick Lens Conversion (ABCD Matrix Method)

#### Spherical Surface Sag

$$\text{sag}(R, a) = R - \text{sign}(R)\sqrt{R^2 - a^2}, \quad a = D/2$$

Reports an error when $|R| < a + 10^{-6}$ mm (geometrically impossible).

#### Thick Lens Power (from matrix product)

$$\Phi_\text{thick} = (n-1)(c_1 - c_2) - \frac{(n-1)^2}{n} \cdot t \cdot c_1 c_2$$

#### Power-Preserving Curvature Uniform Scaling (Strategy A, `thickening.py:232–309`)

Let $c_1' = sc_1$, $c_2' = sc_2$, requiring $\Phi_\text{thick}(s) = \Phi_\text{thin}$:

$$B_\text{prod} \cdot s^2 + A \cdot s - A = 0, \quad A = c_1 - c_2,\quad B_\text{prod} = \frac{(n-1)tc_1c_2}{n}$$

$$\Delta = A^2 + 4AB_\text{prod}, \qquad s = \frac{-A \pm \sqrt{\Delta}}{2B_\text{prod}}$$

Take the positive root closest to 1 ($R' = R/s$, shape-preserving correction, shape factor $Q$ unchanged).

> [!warning] Review Point
> Returns `None` when $\Delta < 0$ (cannot be corrected; code is correct). However, when $B_\text{prod} < 0$ (negative lens), the sign consistency of the two roots needs to be verified.

#### ABCD Matrix EFL (`thickening.py:336–414`)

$$M_\text{refr}(R, n_1 \to n_2) = \begin{pmatrix}1 & 0 \\ -(n_2-n_1)/R & 1\end{pmatrix}, \qquad M_\text{tran}(t, n) = \begin{pmatrix}1 & t/n \\ 0 & 1\end{pmatrix}$$

**Cemented doublet (5 surfaces)**:

$$M = M_\text{refr}(R_3, n_2{\to}1) \cdot M_\text{tran}(t_2, n_2) \cdot M_\text{refr}(R_2, n_1{\to}n_2) \cdot M_\text{tran}(t_1, n_1) \cdot M_\text{refr}(R_1, 1{\to}n_1)$$

$$f' = -1/M_{10}$$

**EFL iterative correction** (up to 2 rounds, `_CORRECTION_ITER = 2`):

$$k = f'_\text{target} / f'_\text{actual}, \qquad R_i \leftarrow R_i \cdot k \quad \text{(uniform scaling, Q unchanged)}$$

Convergence condition: $|k - 1| < 10^{-9}$.

> [!note]
> In the air-spaced doublet, the air gap **does not participate in scaling** (fixed user constraint). For designs with a large air_gap, the air_gap/EFL ratio changes after scaling, introducing an error of order $O(d/f')$ (noted in code comments).

### 3.7 Thermal Analysis (Schott First-Order Thermal Dispersion, `thermal.py`)

**Refractive index temperature coefficient**:

$$\frac{dn}{dT} = \frac{n^2-1}{2n}\left[D_0 + \frac{E_0}{\lambda^2 - \lambda_{tk}^2}\right] \quad \text{[1/K]}$$

Parameters from the TD line in the AGF file ($D_0\ \text{[1/K]}$, $E_0\ \text{[µm²/K]}$, $\lambda_{tk}\ \text{[µm]}$).

**Thermo-optical coefficient**:

$$V = \frac{dn/dT}{n-1} - \alpha \quad \text{[1/K]}, \qquad \alpha = \text{CTE}_{(-40,+20)} \times 10^{-6}\ \text{[1/K]}$$

$V < 0$ (for most optical glasses) → focal length increases with temperature.

**System thermal power derivative** (normalized power, thick-lens approximation):

$$\frac{d\Phi}{dT}_\text{norm} = V_1\varphi_1 + V_2\varphi_2 \quad \text{[1/K]}$$

Exact for cemented doublets; error approximately $O(d/f') \approx 3\text{–}6\%$ for air-spaced doublets.

**Passive athermalization condition**:

$$\frac{df'}{dT} = \frac{dl}{dT} \implies \alpha_{h,\text{req}} = -\frac{d\Phi}{dT}_\text{norm} \quad \text{[1/K]}$$

**Thermal defocus**:

$$\delta f' = f' \cdot (\alpha_{h,\text{req}} - \alpha_{h,\text{actual}}) \cdot \Delta T \quad \text{[mm]}$$

> [!note] Approximation Note
> Only the $D_0$ and $E_0$ terms are used; higher-order corrections $D_1 \Delta T$ and $D_2 \Delta T^2$ are ignored. Valid for $\Delta T \leq 20\text{ K}$; the code documentation states this is a "first-order" approximation.

---

## 4. Module Details

### 4.1 models.py — Data Structures

#### `Inputs` (frozen dataclass)

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| lam0, lam1, lam2 | float | — | Three-color wavelengths [µm] |
| D | float | 50.0 | Entrance beam diameter [mm] |
| fprime | float | 200.0 | Target effective focal length [mm] |
| C0, P0, W0 | float | 0.0 | Seidel target values |
| min\_delta\_nu | float | 10.0 | Minimum Abbe number difference |
| max\_PE | float | 100.0 | PE penalty function upper limit |
| N | int | 20 | Top-N candidate count |
| system\_type | Literal | "cemented" | "cemented" or "spaced" |
| air\_gap | float | 1.0 | Air gap [mm] |
| eps | float | 1e-12 | Numerical zero threshold |
| root\_imag\_tol | float | 1e-9 | Imaginary root filter threshold |
| T\_ref | float | 20.0 | Reference temperature [°C] |
| T\_delta | float | 20.0 | Thermal analysis temperature difference [K] |
| alpha\_housing | Optional[float] | None | Housing thermal expansion CTE [1/K] |

Factory method `Inputs.with_defaults(**overrides)` for convenient testing and rapid development.

> [!warning] Review Point
> The `wavelengths` order in `ThickPrescription` is `(lam1, lam0, lam2)` (short/primary/long), not the intuitive `(lam0, lam1, lam2)`. This is correctly matched in the builder; please be careful during maintenance.

#### `Candidate` (mutable dataclass)

Core design data carrier:

- **Glass information**: glass1/2, catalog1/2, n1/2, ν1/2, φ1/2
- **Shape parameters**: Q (cemented) or Q1/Q2 (spaced)
- **Radii of curvature**: radii (3 or 4 values)
- **Aberration quantities**: P2, PE (cemented); P\_surfaces, PE (spaced)
- **Dispersion data**: formula\_id1/2 + cd1/2 (passed to builder)
- **Thermal analysis**: thermal: ThermalMetrics | None

#### `ElementRx` (single element thick-lens prescription)

```
R_front, R_back     Front and back radii of curvature [mm]
t_center, t_edge    Center/edge thickness [mm]
nd, vd              Refractive index and Abbe number
formula_id, cd      AGF dispersion data
```

#### `ThickPrescription` (complete thick-lens prescription)

```
system_type         "cemented" | "spaced"
elements            [ElementRx, ElementRx]  (always 2)
air_gap             float | None  (air-spaced doublet only)
back_focus_guess    0.9 × f'  (image_solve initial value)
D                   Entrance beam diameter [mm]
wavelengths         (lam1, lam0, lam2) [µm]
actual_efl          ABCD matrix actual EFL [mm]
efl_deviation       (actual − target) / target
```

---

### 4.2 glass_reader.py — AGF File Parsing

Supports the following record line types in the Zemax standard AGF format:

| Line Type | Key Fields | Description |
|-----------|-----------|-------------|
| NM | name, formula\_id, nd\_ref, vd\_ref, exclude\_sub | Glass identifier |
| CD | up to 8+ floating-point coefficients | Dispersion formula coefficients |
| LD | ld\_min\_um, ld\_max\_um | Valid wavelength range [µm] |
| ED | cte\_m40\_20, cte\_20\_300, density, dPgF | Thermal properties [ppm/K] |
| TD | D₀, D₁, D₂, E₀, E₁, λ\_tk, T\_ref | Schott thermal dispersion coefficients |
| OD | relative\_cost, CR, FR, SR, AR, PR | Cost and chemical stability |

Public functions:
- `read_agf(path) -> (comment, list[Glass])`: Single file parsing
- `load_catalog(paths) -> list[Glass]`: Multi-catalog merged loading

---

### 4.3 optics.py — Optical Mathematics

```python
refractive_index(glass, wavelength_um) -> float
    # 12 AGF formulas; formula_id 9~12 use Sellmeier 1 fallback

compute_abbe_number(glass, lam0, lam1, lam2) -> float
    # ν = (n₀−1)/(n₁−n₂); returns inf when denominator < 1e-10

filter_glasses(glasses, cfg) -> list[Glass]
    # Filters: exclude_sub / wavelength range / no formula / n(lam0) computation failure

achromat_power(nu1, nu2, C0) -> (phi1, phi2)
    # φ₁ = ν₁(1−ν₂C₀)/(ν₁−ν₂), φ₂ = 1−φ₁

check_min_radius(radii, D, margin=0.05) -> bool
    # |R| ≥ (D/2)(1+0.05) for all surfaces

prepare_glass_data(glasses, lam0, lam1, lam2) -> list[(Glass, n0, nu)]
    # Calls filter + compute, pre-computes (n₀, ν), used by cemented/spaced
```

---

### 4.4 cemented.py — Cemented Doublet Synthesis

```
prepare_glass_data()
    ↓
for (g1,n1,ν1) × (g2,n2,ν2):
    ── Filter same glass / Δν < min_delta_nu
    ── achromat_power()  →  (φ₁, φ₂)
    ── _solve_Q_roots()  →  [Q₁, Q₂]  ← quadratic equation, up to 2 roots
    ── for Q in roots:
        ── _coma_W()     →  W
        ── _radii()      →  (R₁, R₂, R₃)
        ── check_min_radius()
        ── _P2_and_PE()  →  (P₂, PE)
        ── PE > max_PE   →  skip
        ── compute_thermal_metrics()
        ── Maintain Top-N heap (sorted by |W−W₀|)
    ↓
sort by (|W−W₀|, PE)
return list[Candidate]
```

---

### 4.5 spaced.py — Air-Spaced Doublet Synthesis

```
for (g1,g2) glass pair:
    _coeffs()        →  Seidel coefficient dict {A1,B1,A2,B2,K1,K2,L}
    _C_const()       →  constant term C_const

    [Algebraic path]  _solve_Q_pairs()      →  pairs_alg   (≤ 2 pairs)
    [Sweep path]      _sweep_Q1_coma_line() →  pairs_sweep (≤ 3 pairs)

    Deduplication (|ΔQ₁| < 0.05) → all_pairs

    for (Q₁, Q₂) in all_pairs:
        _radii()              →  (R₁, R₂, R₃, R₄)
        check_min_radius()
        Adjacent surface geometric interference check
        _seidel_SA_residual() →  PE = |S_I residual|
        PE > max_PE           →  skip
        _Ps_and_PE()          →  P_surfaces
        ↓
        Candidate(...)

sort by (PE, Σ|P_surfaces|)
return list[Candidate]   ← No Top-N limit
```

---

### 4.6 thickening.py — Thin→Thick Conversion

**Manufacturing Standard Reference Tables (Chinese Optical Manufacturing Standards)**

Table 10-2 (Single-element machining allowance Δ):

| Aperture D ≤ (mm) | Δ (mm) |
|:---:|:---:|
| 6 | — (no requirement) |
| 10 | 1.0 |
| 18 | 1.5 |
| 30 | 2.0 |
| 50 | 2.5 |
| 80 | 3.0 |
| 120 | 3.5 |
| ∞ | 4.5 |

Table 10-3 (Minimum thickness):

| Aperture D ≤ (mm) | t_edge_min (positive lens) | t_center_min (negative lens) |
|:---:|:---:|:---:|
| 6 | 0.4 | 0.6 |
| 10 | 0.6 | 0.8 |
| 18 | 0.8 | 1.0 |
| 30 | 1.2 | 1.5 |
| 50 | 1.8 | 2.2 |
| 80 | 2.4 | 3.5 |
| 120 | 3.0 | 5.0 |
| ∞ | 4.0 | 8.0 |

**Thickness allocation logic**:

```python
t_from_edge = te_min + sag_front - sag_back   # Satisfies minimum edge thickness

if phi >= 0:   # Positive lens
    t_center = max(t_from_edge, 0.5)
else:          # Negative lens
    t_center = max(t_from_edge, tc_min_neg, 0.5)

t_edge = t_center - sag_front + sag_back
```

**EFL iterative correction process** (3 loop iterations = 1 initial + 2 corrections):

```
for _iter in range(3):
    Compute thickness element_thickness()
    Compute actual_efl  ← ABCD matrix
    if _iter < 2:
        k = fprime / actual_efl
        Check k validity (0 < k ≤ 5, |k−1| > 1e-9)
        R₁,R₂,R₃[,R₄] *= k   ← Uniform scaling, Q unchanged
```

---

### 4.7 thermal.py — Thermal Analysis

Call chain (never throws exceptions; the entire body is wrapped in `try/except Exception`):

```
compute_thermal_metrics(g1, g2, n1, n2, φ1, φ2, λ)
    ── dn_dT(g1, λ, n1)                → dndt1 | None
    ── dn_dT(g2, λ, n2)                → dndt2 | None
    ── thermo_optical_coeff(g1, λ, n1) → V1 | None
    ── thermo_optical_coeff(g2, λ, n2) → V2 | None
    ── system_thermal_power_derivative(V1,φ1,V2,φ2) → dphi | None
    ── required_housing_cte(dphi)      → alpha_req | None
    ── return ThermalMetrics(...)
```

`ThermalMetrics.thermal_data_available = True` if and only if both glasses have complete TD and ED data.

---

### 4.8 pipeline.py — Unified Pipeline

Three public entry functions:

| Function | Input | Output | Description |
|----------|-------|--------|-------------|
| `process_candidate(cand, inputs)` | Single Candidate | PipelineResult | Complete single-candidate pipeline |
| `run_pipeline(candidates, inputs, max_n, on_progress)` | Candidate list | list[PipelineResult] | Batch processing with optional progress callback |
| `run_design(inputs, catalog_paths, on_progress)` | Config + catalog paths | DesignResult | High-level entry point, end-to-end |

`PipelineResult.to_dict()` serialization is guaranteed to include 25+ fixed keys, enforced by contract tests in `test_contracts.py`.

---

### 4.9 optiland_bridge/builder.py — Optical Model Construction

**Custom material `AGFMaterial`**:

```python
class AGFMaterial(BaseMaterial):
    def _calculate_n(self, wavelength):
        return refractive_index(fake_glass, wavelength)  # Any wavelength
    def _calculate_k(self, wavelength):
        return 0.0   # Absorption coefficient ignored
```

Breaks through optiland's default 0.38–0.75 µm visible-range restriction.

**Surface layout**:

| Surface | Cemented Doublet (5 surfaces) | Air-Spaced Doublet (6 surfaces) |
|:---:|---|---|
| 0 | Object plane (∞) | Object plane (∞) |
| 1 | R₁ (entrance surface), material n₁, thickness t₁ | R₁ (entrance surface), material n₁, thickness t₁ |
| 2 | R₂ (cemented surface), material n₂, thickness t₂ | R₂, air, thickness air\_gap |
| 3 | R₃ (air), thickness BFD\_guess | R₃, material n₂, thickness t₂ |
| 4 | Image plane | R₄ (air), thickness BFD\_guess |
| 5 | — | Image plane |

System configuration: EPD = D, field angles 0° + 1°, `image_solve()` moves the image plane to the paraxial focus.

**Numerical stability**: Radii of curvature are rounded to 6 decimal places (nanometer scale) to avoid NaN values inside optiland.

---

### 4.10 optiland_bridge/evaluator.py — Ray Tracing Evaluation

**`OpticMetrics` field overview**:

| Category | Fields |
|----------|--------|
| Paraxial quantities | efl, fno, bfd |
| Spot diagrams | rms\_spot\_radius, geo\_spot\_radius (per wavelength individually + mean) |
| Seidel third-order | seidel\_sa (spherical), seidel\_cc (coma), seidel\_ac (astigmatism), seidel\_pc (field curvature), seidel\_dc (distortion) |
| Lateral Seidel | TSC, TCC, TAC, TPC, TDC |
| Chromatic aberrations | lch\_c (longitudinal chromatic), tch\_c (transverse chromatic) |
| Thick-lens prescription | PE, Q/Q1/Q2, cost1/2 (from Candidate) |
| Status | success, error\_msg, build\_time\_ms, eval\_time\_ms |

**Chromatic aberration computation patch** (`_patched_dispersion` context manager):

optiland hardcodes F/C lines (0.4861/0.6563 µm) for chromatic aberration calculation. This context manager temporarily replaces them with the design wavelengths (lam1, lam2), which is critical for non-visible-range designs.

Spot diagram settings: 6 rings, hexapolar distribution, covering all three design wavelengths.

---

### 4.11 cli.py / gui.py — User Interface

**CLI** entry point:

```bash
autoachromat --config config_example.json [--thin-only] [--out results.json] [--max-n 50]
```

**GUI** (Tkinter, 1,019 lines, 5-panel layout):

```
───────────────────────────────────────────────────────────────────
│ ① Input Parameter Panel                                           │
│    Optical parameters (f', D, λ) / Design targets (C₀,P₀,W₀)    │
│    / Glass catalog management                                      │
───────────────────────────────────────────────────────────────────
│ ② Control Bar                                                      │
│    Run Pipeline / Load Config / Export JSON/CSV / Progress bar     │
───────────────────────────────────────────────────────────────────
│ ③ Results Table (16-column Treeview)                               │
│    #·Glass1·Glass2·Type·EFL·FNO·t₁·t₂·gap·RMS·GEO·SA·LchC·PE·α_h │
───────────────────────────────────────────┬───────────────────────
│ ④-Left: Surface list  │ ④-Mid: Aberrations │ ④-Right: Thin→thick  │
│                       │                    │ comparison + thermal  │
───────────────────────────────────────────┴───────────────────────
│ ⑤ 2D Ray Diagram (Matplotlib embedded, 5 rays, figsize=(12,3))    │
───────────────────────────────────────────────────────────────────
```

Pipeline runs in a background thread; `self.after()` updates the UI in a thread-safe manner.

---

## 5. Data Flow and Call Chain

```
User configuration (JSON / GUI)
        ↓
  run_design(inputs, catalog_paths)
        ↓
  load_catalog()  →  list[Glass]
        ↓
  run_cemented() or run_spaced()
        ↓  ── prepare_glass_data()
        ↓  ── Enumerate glass pairs × Seidel solving
        ↓  ── compute_thermal_metrics()  (per pair)
        ↓
  list[Candidate]  (thick prescription + thermal analysis)
        ↓
  run_pipeline(candidates[:N], inputs)
        ↓  per candidate:
        ↓  ── thicken()
        ↓  │    ── element_thickness() × 2
        ↓  │    ── _system_efl_cemented/spaced()
        ↓  │    ── Iterative EFL correction (up to 2 rounds)
        ↓  │    → ThickPrescription | None
        ↓  ── build_optic_from_prescription()
        ↓  │    ── Construct AGFMaterial × 2
        ↓  │    ── Assemble surface layout
        ↓  │    ── image_solve()
        ↓  │    → optiland.Optic | None
        ↓  ── evaluate(optic, candidate, inputs)
        ↓       ── Paraxial quantities
        ↓       ── SpotDiagram (6-ring hexapolar)
        ↓       ── Seidel coefficients
        ↓       ── _patched_dispersion → chromatic aberrations
        ↓       → OpticMetrics
        ↓
  list[PipelineResult]
        ↓
  GUI display / JSON / CSV export
```

---

## 6. Test Coverage

| Test File | Main Test Classes | Key Verification Content |
|-----------|------------------|--------------------------|
| `test_contracts.py` | `TestToDictKeyStability` | `to_dict()` always contains all 25+ fixed keys |
| `test_contracts.py` | `TestFailureSemantics` | Pipeline never throws exceptions |
| `test_contracts.py` | `TestInputsWithDefaults` | Frozen dataclass is immutable; `with_defaults` is correct |
| `test_refactoring.py` | `TestPrepareGlassData` | Filters exclude\_sub / wavelength range handling |
| `test_refactoring.py` | `TestRunPipeline` | Mock verifies batch pipeline behavior |
| `test_thermal.py` | `TestDnDT` | N-BK7: −1.39×10⁻⁶ /K, N-SF5: −1.44×10⁻⁶ /K |
| `test_thermal.py` | `TestRequiredHousingCte` | $\alpha_h = -d\Phi/dT$, BK7+SF5 approximately 1–10 ppm/K |
| `test_thermal.py` | `TestThermalDefocus` | Linear verification of thermal defocus formula |
| `test_thickness.py` | `TestSag` | Positive/negative curvature; error when $|R| < a+\text{margin}$ |
| `test_thickness.py` | `TestCorrectRadiiForThickness` | After correction, $\Phi_\text{thick} \approx \Phi_\text{thin}$ |
| `test_thickness.py` | `TestReconcileCementedRadius` | Average taken when cemented surface radius mismatches |

Tests use real SCHOTT AGF data (N-BK7, N-SF5); results depend on the catalog file being present.

**Run commands**:

```bash
cd c:\MasterThesis\AutoAchromat
pytest tests/ -v
python scripts/smoke_test_bridge.py --config config_example.json --max-n 5
```

---

## 7. Module Dependency Relationships

```
models.py          ← No external dependencies
glass_reader.py    ← No external dependencies
optics.py          ← glass_reader
thermal.py         ← glass_reader
cemented.py        ← glass_reader, models, optics, thermal
spaced.py          ← glass_reader, models, optics, thermal, numpy
thickening.py      ← models (only)
pipeline.py        ← glass_reader, models, cemented, spaced,
                      thickening, optiland_bridge.*
cli.py             ← pipeline, models, glass_reader
gui.py             ← pipeline, models, glass_reader, optiland_bridge.*
optiland_bridge/
  builder.py       ← models, optics, glass_reader, thickening, optiland
  evaluator.py     ← models, builder, optiland
```

**External Dependencies**:

| Dependency | Scope | Purpose |
|------------|-------|---------|
| `numpy` | `spaced.py` only | `numpy.roots()` for solving quadratic equations |
| `optiland` | `optiland_bridge/` only | Ray tracing engine |
| Python standard library | All modules | math, dataclasses, logging, typing, json, tkinter |

---

## 8. Design Decisions and Review Points

### 8.1 Known Design Trade-offs

| Item | Decision | Rationale |
|------|----------|-----------|
| Normalized power | $\varphi_1 + \varphi_2 = 1$ | Simplifies formulas; system EFL is controlled by the fprime parameter |
| Seidel third-order theory | Seidel approximation throughout | Fast screening; accurate evaluation performed by optiland ray tracing |
| First-order thermal analysis | Only $D_0$ and $E_0$ terms | Valid for $\Delta T \leq 20\text{ K}$; stated in documentation |
| 2-round curvature correction | Does not pursue machine-precision convergence | Error typically < 0.1% after 2 rounds; meets optical design precision requirements |
| Air gap not scaled | Fixed user constraint | Introduces sub-order error (stated); does not practically affect thermal approximation error |

### 8.2 Mathematical Points Requiring Manual Verification

> [!warning] Priority Review — PE Penalty Function Exponential Term
> `cemented.py:84`: `PE = |P₂| · 3^|W| / R₂²`
>
> - The physical origin and literature basis of the `3^|W|` exponential factor need to be clarified
> - When $|W| > 3$, PE is maximized ($3^3 = 27$, $3^4 = 81$), exerting a dominant effect on candidate filtering
> - Dimensional consistency: the dimensions of $P_2$, and the normalization dimensions of $R_2^2$

> [!warning] Priority Review — Air-Spaced Doublet Refraction Angle Formula
> `spaced.py:288`: `u4 = phi1 + n2 + Q2 * (n2 - 1.0)`
>
> This expression mixes dimensionless quantities (phi1) with refractive indices (n2). The dimensions and physical meaning need to be verified against the derivation in the reference literature.

> [!warning] Priority Review — Cemented Doublet Radius of Curvature Formula Origin
> `cemented.py:57–68` comments are labeled "Your spec" (a custom-derived formula, not from standard literature). The expressions for the three radii $R_1, R_2, R_3$ as functions of $Q$, $\varphi_1$, $n_1$, $n_2$ need to be verified against the derivation document.

> [!warning] Priority Review — formula_id 9~12 Fallback
> `optics.py:111–119`: All unknown formula\_ids are computed using the Sellmeier 1 format (silent handling, no warning). If glasses with these IDs exist in the catalog, incorrect refractive indices will be returned without any error.

> [!note] Potential Improvement — Generalized Achromatic Parameter C₀
> `optics.py:183`: When $C_0 \neq 0$, the generalized achromatic condition gives $\varphi_1/\nu_1 + \varphi_2/\nu_2 = C_0/(\nu_1\nu_2)$. The physical definition needs to be confirmed as consistent with the reference literature.

### 8.3 Code Quality Observations

**Strengths**:

- None-propagation defensive programming is consistently applied throughout
- ABCD matrix implementation has no numpy dependency (pure Python tuples), good readability
- The `AGFMaterial` design for breaking through optiland's wavelength restriction is elegant
- Contract tests guarantee the stability of API serialization
- `thermal.py` and `thickening.py` have detailed module-level docstrings

**Potential Issues**:

| Issue | Location | Severity |
|-------|----------|:---:|
| Air-spaced doublet has no Top-N filtering | `spaced.py` | Medium |
| GUI is 1,019 lines in a single file with logic and UI mixed | `gui.py` | Low |
| `results.json` is tracked in git | Root directory | Low |
| `test.ipynb` is in tests/ but is not a pytest test | `tests/` | Low |
| formula_id 9~12 silent fallback | `optics.py:111` | High |
