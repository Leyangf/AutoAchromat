# Architectural Degeneracy in Automated Achromatic Doublet Synthesis and Its Implications for Opto-Mechanical-Thermal Co-Design

## Abstract

This report examines the structural limitations of classical thin-lens Seidel synthesis as a generator for opto-mechanical-thermal (OMT) co-design exploration. An automated pipeline was implemented that enumerates achromatic doublet candidates from optical glass catalogs, converts them to thick-lens prescriptions, and evaluates their performance via ray tracing. During systematic evaluation, a significant architectural degeneracy was observed: candidates from diverse glass pairs converge to a narrow family of geometrically similar designs. Mathematical analysis traces this degeneracy to three sources — the zero residual degrees of freedom in the cemented doublet constraint system, the correlation structure of the glass catalog in the Abbe diagram, and the premature enforcement of the spherical aberration constraint in the spaced doublet synthesis. These findings have direct implications for the viability of thin-lens synthesis as a starting-point generator for multi-domain design exploration.

---

## 1. Introduction

### 1.1 Motivation

The design of optical systems increasingly requires simultaneous consideration of optical performance, mechanical packaging, and thermal stability. In aerospace, defense, and precision instrumentation, the lens system must maintain performance across temperature excursions, survive mechanical loads, and fit within strict volume and mass budgets. These constraints interact: glass selection affects both chromatic correction and thermal sensitivity; element thickness affects both aberration balance and mechanical stiffness; housing material selection couples thermal expansion to optical focus shift.

Traditional optical design workflows address these domains sequentially. An optical designer selects glasses and optimizes surface curvatures for image quality; a mechanical engineer then designs the housing, mounts, and retention mechanisms; a thermal analyst evaluates focus shift under temperature changes. If the thermal or mechanical analysis reveals infeasibility, the process iterates — often requiring a restart from glass selection.

An alternative approach is to generate a diverse set of optically viable architectures at the earliest design stage, each characterized not only by optical merit but also by mechanical and thermal properties. This would enable concurrent evaluation across all three domains, allowing system-level trade studies before committing to a single optical configuration.

### 1.2 Scope

This work investigates whether classical thin-lens achromatic doublet synthesis — automated to enumerate all viable glass pairs from commercial catalogs — can serve as this early-stage architecture generator. The investigation is based on a complete implementation of the synthesis-to-evaluation pipeline, covering cemented and air-spaced doublets, with thermal analysis integrated at the synthesis stage.

### 1.3 Reference Baseline

The implemented pipeline follows the algorithmic framework established by two prior works from the ITMO Applied Optics group:

- **Nguyen & Bakholdin (2022)**: Automated synthesis and ranking of cemented and air-spaced doublets using a preliminary evaluation criterion based on third-order aberration theory.
- **Ivanova, Romanova, Zhukova & Kalinkina (2019)**: Extension of the cemented doublet synthesis to include passive athermalization via thermo-optical coefficient computation.

The current implementation reproduces and extends both works by combining optical synthesis with thermal analysis for both doublet types, and by adding a thick-lens verification stage via numerical ray tracing.

---

## 2. Methodology

### 2.1 Pipeline Architecture

The automated pipeline proceeds through five stages:

1. **Glass catalog loading**: AGF-format catalogs (SCHOTT, OHARA, CDGM) are parsed. For each glass, the refractive index at the design wavelength and the Abbe number are computed from the glass dispersion formula. Glasses outside the valid wavelength range or with insufficient dispersion data are excluded.

2. **Thin-lens synthesis**: All ordered glass pairs satisfying a minimum Abbe number difference are enumerated. For each pair, the achromatic power split is computed from the generalized achromatic condition. The shape parameter is then solved from the third-order spherical aberration equation (a quadratic yielding 0–2 real roots). Coma is evaluated but not constrained in the cemented case; both coma and spherical aberration are simultaneously constrained in the spaced case.

3. **Candidate screening**: Candidates are filtered by minimum radius of curvature (to exclude near-hemispherical surfaces), maximum power error metric, and geometric overlap constraints. A ranking metric combines the spherical aberration residual, coma deviation, and surface curvature.

4. **Thick-lens conversion**: Surviving thin-lens candidates are assigned physical center and edge thicknesses based on manufacturing standard tables, with iterative radius correction to preserve the thin-lens focal length.

5. **Ray-tracing evaluation**: Each thick-lens prescription is converted to a numerical optical model. Paraxial properties (EFL, F/#, back focal distance), spot diagram metrics (RMS and geometric spot radii), Seidel aberration coefficients, and chromatic aberration values are extracted.

### 2.2 Thermal Analysis

For each glass pair with available thermal dispersion data, the pipeline computes:

- The thermo-optical coefficients $V_1$ and $V_2$ for each element, incorporating the temperature-dependent refractive index change $dn/dT$ and the coefficient of thermal expansion (CTE) of the glass.
- The normalized thermal power derivative: $\frac{d\Phi}{dT} = V_1 \varphi_1 + V_2 \varphi_2$
- The required housing CTE for passive athermalization: $\alpha_h = -\frac{d\Phi/dT}{\Phi}$

These metrics are computed at the thin-lens stage and reported alongside the optical performance.

---

## 3. Observations

### 3.1 Convergence of Candidate Geometries

Systematic runs across the full SCHOTT catalog (approximately 200 usable glasses after filtering) revealed an unexpected pattern: despite enumerating thousands of glass pairs, the resulting designs cluster into a narrow family of geometrically similar configurations. Specifically:

- **Power split convergence**: The element powers $\varphi_1$ and $\varphi_2$ span a bounded range (approximately 1.7 to 3.5 and −0.7 to −2.5, respectively) for the overwhelming majority of crown-flint pairs.
- **Radius convergence**: The front surface radius $R_1$, cemented/inner radius $R_2$, and back radius $R_3$ cluster within factors of 2–3 across hundreds of distinct glass pairs, for a fixed focal length.
- **Post-thickening collapse**: Many candidates that are distinguishable at the thin-lens stage become geometrically indistinguishable after thick-lens conversion, as the iterative radius correction further homogenizes the designs.

A deduplication mechanism based on optical fingerprinting (grouping by rounded $n_1, n_2, \nu_1, \nu_2$) was introduced to address this redundancy. This confirmed that optically distinct glass pairs frequently produce identical or near-identical geometries.

### 3.2 Geometric Infeasibility After Thickening

A secondary observation is that some thin-lens candidates become geometrically infeasible when physical thicknesses are assigned:

- Surface sag at the aperture edge exceeds the available center thickness.
- Edge thickness falls below the manufacturing minimum.
- The iterative radius correction diverges or produces unphysical values.

Conversely, some candidates pruned early by the power error or minimum radius criteria might become feasible designs if the synthesis accounted for realistic element thicknesses from the outset.

### 3.3 Thermal Metric as a Passive Output

The required housing CTE $\alpha_h$ varies substantially across glass pairs (from negative values to values exceeding 30 ppm/K). However, this metric is currently computed and reported without being used as a selection criterion. The pipeline outputs designs requiring physically unrealizable housing materials (e.g., negative CTE) alongside thermally feasible designs, without distinguishing them.

---

## 4. Mathematical Analysis of the Degeneracy

### 4.1 Degrees of Freedom in the Cemented Doublet

The cemented doublet has three surface curvatures $(R_1, R_2, R_3)$ as geometric degrees of freedom. The synthesis imposes three constraints:

1. **Total power**: $\varphi_1 + \varphi_2 = \Phi$ (normalization to unit focal length)
2. **Achromatic condition**: $\varphi_1 / \nu_1 + \varphi_2 / \nu_2 = C_0$ (chromatic correction)
3. **Spherical aberration**: $A Q^2 + B Q + C = P_0$ (third-order SA target)

Constraints (1) and (2) form a linear system that uniquely determines the power split $(\varphi_1, \varphi_2)$ as a function of the Abbe numbers $(\nu_1, \nu_2)$ alone. Constraint (3) is a quadratic in the shape parameter $Q$, whose coefficients $A$, $B$, $C$ are explicit functions of $(n_1, n_2, \varphi_1, \varphi_2)$.

The result is that the design — all three surface curvatures — is **uniquely determined** (up to two roots of the quadratic) by the four-parameter tuple $(n_1, n_2, \nu_1, \nu_2)$. There are **zero residual continuous degrees of freedom**. Each glass pair produces at most two discrete designs. No bending, packaging, or thermal optimization is possible within the cemented synthesis framework without changing the glass pair.

### 4.2 Glass Catalog Clustering in the Abbe Diagram

The Abbe diagram for standard optical glasses exhibits a well-known correlation between refractive index and Abbe number. Crown glasses cluster around $(\nu \approx 55\text{–}70,\; n \approx 1.50\text{–}1.60)$; flint glasses cluster around $(\nu \approx 25\text{–}40,\; n \approx 1.60\text{–}1.85)$. Along the main glass families, knowing $\nu$ approximately determines $n$ to within $\Delta n \approx 0.05\text{–}0.10$.

This means the nominally 4-dimensional glass pair space $(n_1, n_2, \nu_1, \nu_2)$ is effectively reduced to approximately 2 dimensions (parameterized by $\nu_1$ and $\nu_2$), because the refractive indices are correlated with the Abbe numbers along the glass line.

Consequently:
- The power split $\varphi_1 = \nu_1 / (\nu_1 - \nu_2)$ is a smooth function of two variables.
- The shape parameter $Q$ varies continuously but slowly across the catalog, because its coefficients depend on $(n_1, n_2)$ which co-vary with $(\nu_1, \nu_2)$.
- The resulting surface curvatures trace a 2-dimensional surface in the 3-dimensional $(R_1, R_2, R_3)$ space.

Glasses that deviate from the main glass line — lanthanum crowns, fluorocrowns, anomalous partial dispersion glasses — do produce genuinely distinct designs, but these constitute a small minority of the catalog.

### 4.3 The Spaced Doublet: a Discarded Continuous Freedom

For the air-spaced doublet, the synthesis solves two simultaneous constraints:

1. **Coma**: $K_1 Q_1 + K_2 Q_2 = W_0 - L$ — a linear constraint defining a **line** in $(Q_1, Q_2)$ space.
2. **Spherical aberration**: a quadratic condition — defining a **conic section** in $(Q_1, Q_2)$ space.

The intersection of a line with a conic section yields at most 2 points. These are the only designs retained by the synthesis.

However, every point on the coma-zero line satisfies the coma constraint with a varying spherical aberration residual. Different positions along this line correspond to different element bendings and therefore different physical layouts:

- Different sag distributions → different center thickness requirements
- Different curvature ratios → different total track lengths
- Different bending balance → different sensitivity to manufacturing tolerances

This 1-dimensional manifold of valid (coma-corrected) designs is entirely discarded by the synthesis, which retains only the 2 points where spherical aberration also vanishes. If a subsequent thick-lens optimization stage is available to correct spherical aberration, the SA=0 constraint at the synthesis stage is redundant — and its enforcement collapses a continuous design space into two discrete points per glass pair.

### 4.4 The Discriminant as a Hard Filter

For a substantial fraction of crown-flint combinations, the discriminant $B^2 - 4AC$ of the shape parameter quadratic is negative, meaning no real solution exists. This is not a limitation of the algorithm but a genuine physical constraint: those glass pairs cannot form an achromatic doublet with simultaneously corrected spherical aberration in the thin-lens limit. This further concentrates the surviving designs into a subset of the catalog.

---

## 5. Implications for Opto-Mechanical-Thermal Co-Design

### 5.1 Architectural Diversity

The observed degeneracy has direct consequences for the goal of generating diverse starting points for OMT exploration:

| Diversity dimension | Current state | Required for OMT |
|---|---|---|
| Glass pair selection | O(n²) enumeration, but clustered | Sufficient if thermal filtering is added |
| Element bending (shape factor) | 1–2 discrete values per pair | Continuous sampling needed |
| Doublet orientation | Single orientation only | Crown-first and flint-first variants |
| Total track length | Not parameterized | Explicit packaging constraint |
| Field angle | Fixed at 1° | Configurable per application |
| Architecture type | Cemented and spaced only | Reversed, meniscus-first, triplet |

The pipeline currently generates **glass-pair diversity** but not **architectural diversity**. For OMT co-design, the relevant diversity dimensions are packaging geometry (track length, housing diameter), thermal behavior (housing CTE compatibility), and mechanical properties (element mass, gravity sensitivity) — none of which are varied by the current synthesis.

### 5.2 The Housing CTE Gap

The thermal analysis computes the required housing CTE $\alpha_h$ for each design, but this value is never used as a constraint or selection criterion. For practical OMT co-design, the housing material is often predetermined (aluminum: $\alpha_h \approx 23.6$ ppm/K; titanium: $\alpha_h \approx 8.6$ ppm/K; Invar: $\alpha_h \approx 1.3$ ppm/K). Designs requiring housing CTE values far from the available materials are thermally infeasible and should be filtered or penalized.

Furthermore, the thermal metric is computed at the thin-lens stage and not recomputed after thick-lens conversion. The thick-lens thermal focus shift includes a mechanical spacer coupling term that is absent from the thin-lens model:

$$\Delta f' / \Delta T \approx -f^2 \left( V_1 \varphi_1 + V_2 \varphi_2 - \alpha_h \left(1 + d_{\text{gap}} / f' \right) \right)$$

where $d_{\text{gap}}$ is the physical air gap. This correction is first-order in the spacer length and is non-negligible for spaced doublets.

### 5.3 Mechanical Properties Not Propagated

The glass catalogs contain density data (from the AGF ED line), but the pipeline does not propagate this to the output. Element mass, moment of inertia, and total system mass cannot be computed from the exported results without re-accessing the glass catalog. For OMT co-design, these are essential selection criteria.

---

## 6. Comparison with Reference Works

### 6.1 Nguyen & Bakholdin (2022)

The 2022 paper presents the synthesis and ranking algorithm that this pipeline reproduces. The ranking metric (preliminary evaluation criterion, termed PE in this work) combines the Petzval sum residual, coma deviation, and surface curvature into a scalar figure of merit. The paper validates the algorithm by demonstrating automated glass selection and ranking for both cemented and air-spaced doublets.

The paper does not characterize the architectural degeneracy described in this report. The synthesis is presented as producing "different variants" without quantifying how many of these variants are geometrically distinct. The deduplication analysis performed in this work suggests that the effective diversity is substantially lower than the nominal candidate count.

### 6.2 Ivanova et al. (2019)

The 2019 paper extends the cemented doublet synthesis to include passive athermalization. Its core contribution is the integration of the thermo-optical coefficient $V$ into the automated glass selection loop, enabling simultaneous optical and thermal screening. The paper uses a nomogram method for visual glass selection on a $V_1\varphi_1 + V_2\varphi_2 = 0$ diagram.

The current implementation computes the same thermal metrics but does not enforce the athermalization condition as a synthesis constraint. The 2019 paper's approach — using thermal feasibility as a filter during synthesis rather than as a post-hoc display metric — is closer to the OMT co-design objective, but is limited to cemented doublets and does not consider housing material constraints.

### 6.3 Contributions Beyond Both References

The current work extends both references in three directions:

1. **Thermal analysis for spaced doublets**: Neither reference addresses thermal behavior of air-spaced doublets. The current pipeline computes $V_1$, $V_2$, and $\alpha_h$ for both types.
2. **Thick-lens verification**: Neither reference includes a thick-lens evaluation stage. The ray-tracing pipeline provides quantitative validation of the thin-lens predictions, revealing the geometric convergence and infeasibility issues described above.
3. **Degeneracy characterization**: The mathematical analysis of why thin-lens synthesis produces a narrow family of solutions (zero residual DOF, glass clustering, discarded coma-line freedom) is a new finding not present in either reference.

---

## 7. Proposed Direction

The analysis suggests a specific algorithmic modification that would directly address the degeneracy:

**Relax the spherical aberration constraint at synthesis and sample the coma-zero line as a continuous parameter.**

For the spaced doublet, instead of solving $\text{SA} = 0 \cap \text{Coma} = 0$ (yielding 2 discrete points per glass pair), the synthesis would:

1. Parameterize the coma-zero line $K_1 Q_1 + K_2 Q_2 = W_0 - L$ by a free variable (e.g., $Q_1$).
2. Sample $Q_1$ at regular intervals along the feasible range.
3. Accept the resulting SA residual as a starting condition for thick-lens optimization.
4. Filter by packaging constraints (total track length, housing CTE compatibility) rather than by thin-lens SA.

Each sampled point on the coma line corresponds to a distinct physical layout — different element bendings, different sag distributions, different total track lengths — while maintaining zero coma. The thick-lens optimizer would then correct the SA residual while respecting edge thickness, radius bounds, and focal length constraints.

This transforms the synthesis from a **constraint satisfaction problem** (producing 0–2 discrete solutions) into a **design space sampling problem** (producing a continuous family of starting points parameterized by packaging geometry). The resulting candidate set would be diverse along the dimensions that matter for OMT co-design: track length, element thickness distribution, and thermal sensitivity.

For the cemented doublet, where no continuous freedom exists within the standard doublet form, diversity would come from:

- Adding the reversed orientation (flint-first) as a separate architecture type.
- Using the housing CTE target as a glass-pair selection criterion (following the 2019 reference more closely).
- Including anomalous partial dispersion glasses that lie off the main glass line.

---

## 8. Conclusion

The classical thin-lens Seidel synthesis for achromatic doublets, when automated across a full glass catalog, produces a substantially narrower design space than the nominal O(n²) glass pair enumeration suggests. The architectural degeneracy arises from three compounding factors: the constraint system leaves zero continuous degrees of freedom per glass pair, the glass catalog clusters in the Abbe diagram reducing the effective parameter dimensionality, and the spaced doublet synthesis discards a 1-dimensional manifold of coma-corrected designs by prematurely enforcing the spherical aberration constraint.

These findings do not invalidate the thin-lens synthesis approach. Rather, they identify a specific algorithmic modification — relaxing the SA constraint and sampling the coma-zero line — that would transform the synthesis from a discrete enumeration into a continuous design space generator suitable for opto-mechanical-thermal co-design exploration. The implemented pipeline, with its integrated thermal analysis and thick-lens verification, provides the infrastructure needed to evaluate this modified approach.

---

## References

1. Nguyen DH, Bakholdin AB. Automation of synthesis and ranking of cemented and air-spaced doublets. *Computer Optics*. 2022;46(1):83–89. doi:10.18287/2412-6179-CO-923.

2. Ivanova TV, Romanova GE, Zhukova TI, Kalinkina OS. Automation of athermal cemented doublet synthesis. *Scientific and Technical Journal of Information Technologies, Mechanics and Optics*. 2019;19(4):594–601.

3. Tamagawa Y, Tajime T. Dual-element optical system for passive athermalization. *Applied Optics*. 1996;35(10):1751–1756.

4. Schott AG. TIE-19: Temperature coefficient of the refractive index. *Technical Information*. 2016.
