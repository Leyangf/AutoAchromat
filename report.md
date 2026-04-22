
## Chapter 2 — Theoretical Foundations

The two tools developed in this thesis draw on overlapping but not identical theoretical foundations. Glass dispersion models (§2.1), ray tracing procedures (§2.2), and chromatic aberration quantities (§2.4) are shared by both tools. The Seidel aberration framework (§2.3) is specific to the doublet synthesis of Chapter 3, and the imaging simulation architecture (§2.5) is specific to the wave-optics image prediction of Chapter 4. This chapter develops all five foundations in a unified treatment, so that Chapters 3 and 4 can reference them without redundancy.

### 2.1 Optical Glass Dispersion

The refractive index $n(\lambda)$ of an optical glass varies with wavelength according to empirical dispersion formulas. A common form is the Sellmeier equation:

$$n^2(\lambda) = 1 + \sum_{i} \frac{K_i \lambda^2}{\lambda^2 - L_i}$$

Other common dispersion representations (Schott polynomial, Herzberger, Conrady, etc.) are functionally equivalent for the present purposes; the specific catalog formats and formula handling are documented in §3.2.1 and Appendix A.

The Abbe number quantifies the dispersion strength:

$$\nu = \frac{n(\lambda_0) - 1}{n(\lambda_1) - n(\lambda_2)}$$

where $\lambda_0$ is the primary design wavelength and $\lambda_1$, $\lambda_2$ are the short and long wavelengths, respectively. High $\nu$ indicates low dispersion (crown glasses, typically $\nu > 50$); low $\nu$ indicates high dispersion (flint glasses, typically $\nu < 35$).

![[Abbe_diagram.png]]
*Figure: Abbe diagram of 544 Sellmeier-described glasses from the CDGM, Ohara, and Schott catalogs. Crown (high $\nu$, left) and flint (low $\nu$, right) regions are clearly separated. The full catalogs include additional dispersion representations accepted by the pipeline (§3.2.1).*

For a thin doublet with normalized element powers $\varphi_1$ and $\varphi_2$ (satisfying $\varphi_1 + \varphi_2 = 1$), the achromatic condition requires the chromatic power contributions to cancel:

$$\frac{\varphi_1}{\nu_1} + \frac{\varphi_2}{\nu_2} = 0$$

This uniquely determines the power split:

$$\varphi_1 = \frac{\nu_1}{\nu_1 - \nu_2}, \qquad \varphi_2 = \frac{-\nu_2}{\nu_1 - \nu_2}$$

A large difference $|\nu_1 - \nu_2|$ yields milder element powers, reducing higher-order aberrations; in practice a minimum Abbe number difference is imposed to discard pairs whose element powers become unmanageably large (the specific threshold used in the pipeline is documented in §3.2.1). This two-wavelength achromatic condition generalizes to higher-order color correction by imposing additional constraints on partial dispersions: apochromatic (3 wavelengths) and superachromatic (4+ wavelengths) designs require matching partial dispersions across the glass pair, as illustrated by Mikš and Novák [8].

### 2.2 Ray Tracing Fundamentals

Both AutoAchromat (§3) and ChromeFringe (§4) rely on ray tracing to connect lens prescriptions to aberration quantities. This section develops the paraxial and real-ray tracing procedures that underpin both tools. The notation throughout this chapter follows [39].

#### 2.2.1 Paraxial Ray Tracing

Paraxial theory assumes that the angle between rays and the optical axis is extremely small ($\sin\theta \approx \theta$), linearizing Snell's law. For an optical system consisting of $K$ refracting surfaces, the recurrence relations for paraxial marginal ray tracing are:

$$n'_k u'_k = n_k u_k - h_k \phi_k$$

$$h_{k+1} = h_k + u'_k \cdot d_k$$

where:

| Symbol                        | Meaning                                   |
| ----------------------------- | ----------------------------------------- |
| $h_k$                         | Ray height at surface $k$                 |
| $u_k$, $u'_k$                 | Ray slope before/after refraction         |
| $n_k$, $n'_k$                 | Refractive index before/after surface $k$ |
| $\phi_k = (n'_k - n_k) / R_k$ | Optical power of surface $k$              |
| $R_k$                         | Radius of curvature of surface $k$        |
| $d_k$                         | Spacing from surface $k$ to surface $k+1$ |

Tracing a unit-height, zero-slope marginal ray ($h_0 = 1$, $u_0 = 0$) through all surfaces yields the back focal length:

$$\text{BFL}(\lambda) = -\frac{h_K}{u'_K} \quad (\text{mm})$$

The BFL is wavelength-dependent because each refractive index $n_k(\lambda)$ varies with wavelength through the dispersion relations of §2.1. This dependence is the basis of chromatic aberration (§2.4).

#### 2.2.2 ABCD (Ray-Transfer) Matrix Method

With the reduced state vector $\mathbf{r} = (h,\; nu)^\top$, the paraxial ray trace of §2.2.1 takes the matrix form $\mathbf{r}' = M\,\mathbf{r}$, where the system matrix $M$ is the ordered product of refraction and transfer matrices.

**Refraction matrix** at surface $k$:

$$M_{\text{refr},k} = \begin{pmatrix}1 & 0 \\ -\phi_k & 1\end{pmatrix}, \qquad \phi_k = \frac{n'_k - n_k}{R_k}$$

**Transfer matrix** from surface $k$ to surface $k+1$:

$$M_{\text{tran},k} = \begin{pmatrix}1 & d_k/n'_k \\ 0 & 1\end{pmatrix}$$

For a system of $K$ surfaces, the composite matrix $M = M_{\text{refr},K}\,M_{\text{tran},K-1}\cdots M_{\text{tran},1}\,M_{\text{refr},1} = \begin{pmatrix} A & B \\ C & D \end{pmatrix}$ maps the input state to the output state in a single multiplication. For a parallel input ray $(h_0, 0)^\top$, the output is $(Ah_0,\; Ch_0)^\top$. The EFL (measured from the rear principal plane) and BFL (measured from the last surface) in air are then:

$$f' = -\frac{1}{C}, \qquad \text{BFL} = -\frac{A}{C}$$

This method is used in §3.3 for the thin-to-thick lens EFL correction.

#### 2.2.3 Real Ray Tracing

When the paraxial approximation is insufficient — for example, when evaluating spherical aberration or tracing rays at large aperture — Snell's law must be applied without the small-angle linearisation. Real-ray tracing applies the exact refraction law at each surface and propagates rays through finite thicknesses, capturing all orders of aberration. It provides the exact transverse aberration and back-focal intercept at each pupil height, which are essential for the RoRi chromatic model (§4.3) and the ray-fan ESF method (§4.6).

### 2.3 Seidel Aberration Framework

The doublet synthesis equations in §3.2 are built from the third-order aberration contributions of individual refracting surfaces. The underlying theory is well established; this section recapitulates the derivation of the normalized per-surface coefficients $P$ and $W$ that underpin the synthesis, following the Slyusarev notation as systematized by Romanova et al. [2].

#### 2.3.1 Seidel Surface Contributions

The paraxial ray tracing of §2.2.1 defines the marginal ray quantities $h_k$, $u_k$ at each surface. The Seidel framework additionally requires the chief ray, which obeys the same recurrence relations with quantities $\bar{u}$, $\bar{h}$, $\bar{i}$. The angle of incidence at a surface with curvature $c = 1/R$ is $i = u + hc$ for the marginal ray and $\bar{i} = \bar{u} + \bar{h}c$ for the chief ray. The Lagrange invariant $H = n(u\bar{h} - \bar{u}h)$ is preserved across all surfaces.

For on-axis analysis, the relevant Seidel coefficients are $S_I$ (spherical aberration) and $S_{II}$ (coma). The third-order contributions of a single refracting surface are:

$$S_I = -(ni)^2 \, h \, \Delta\!\left(\frac{u}{n}\right), \qquad S_{II} = -(ni)(n\bar{i}) \, h \, \Delta\!\left(\frac{u}{n}\right)$$

where $ni = n'i'$ is the refraction invariant (Snell's law), $\bar{i}$ is the chief ray incidence angle, and $\Delta(u/n) \equiv u'/n' - u/n$ is the change in reduced ray slope across the surface. The ratio $S_{II}/S_I = \bar{i}/i$ links coma to spherical aberration through the chief-to-marginal incidence angle ratio.

The above expressions contain the incidence angle $i$, which depends on surface curvature. To obtain forms that use only quantities directly available from paraxial ray tracing — ray slope changes $\Delta u$ and refractive index changes $\Delta(1/n)$ — the refraction invariant $ni$ is rewritten as follows. From the refraction equation:

$$\Delta u \equiv u' - u = \frac{nu - (n'-n)hc}{n'} - u = -\frac{(n'-n)\,i}{n'}$$

Solving for $i$ and multiplying by $n$ gives $ni = -nn'\Delta u/(n'-n)$. Since $\Delta(1/n) \equiv 1/n' - 1/n = -(n'-n)/(nn')$, this simplifies to:

$$ni = \frac{\Delta u}{\Delta(1/n)}$$

By the same argument for the chief ray: $n\bar{i} = \Delta\bar{u}\,/\,\Delta(1/n)$. Substituting into the expressions above:

$$S_I = -\left(\frac{\Delta u}{\Delta(1/n)}\right)^{\!2} h \, \Delta\!\left(\frac{u}{n}\right), \qquad S_{II} = -\frac{\Delta u}{\Delta(1/n)} \cdot \frac{\Delta\bar{u}}{\Delta(1/n)} \cdot h \, \Delta\!\left(\frac{u}{n}\right)$$

Both expressions share a common structure: a product of difference quotients in marginal ray slopes and refractive indices, multiplied by the ray height $h$ at that surface. The factor $\Delta\bar{u}/\Delta(1/n)$ in $S_{II}$ is the chief ray refraction invariant, which depends on stop position and field angle rather than surface curvatures.

#### 2.3.2 Normalized Coefficients $P$ and $W$

The Seidel expressions of §2.3.1 factor into a ray-height term and a core that depends only on ray slopes and refractive indices. The normalized per-surface coefficients isolate this core:

$$P_k \equiv \left(\frac{\Delta u}{\Delta(1/n)}\right)_{\!k}^{\!2} \Delta\!\left(\frac{u}{n}\right)_{\!k}, \qquad W_k \equiv -\left(\frac{\Delta u}{\Delta(1/n)}\right)_{\!k}\,\Delta\!\left(\frac{u}{n}\right)_{\!k}$$

so that the per-surface Seidel contributions become:

$$S_{I,k} = -P_k \cdot h_k, \qquad S_{II,k} = W_k \cdot \frac{\Delta\bar{u}}{\Delta(1/n)}\bigg|_k \cdot h_k$$

The chief ray factor $\Delta\bar{u}/\Delta(1/n)$ in $S_{II}$ depends on stop position and field angle, not on surface curvatures or glass choice. In the doublet synthesis (§3.2), the aperture stop is fixed at the front surface, so this factor is a known constant that separates from the design variables — hence $W_k$ alone suffices as the per-surface coma coefficient.

In the **thin-lens approximation**, all surfaces of an element share the same ray height $h$, so:

$$S_{I,\text{sys}} = -h \sum_{k=1}^{K} P_k = -h \cdot P_{\text{sys}}$$

The common factor $h$ divides out, and the design equations reduce to:

$$P_{\text{sys}} = \sum_{k=1}^{K} P_k, \qquad W_{\text{sys}} = \sum_{k=1}^{K} W_k$$

Each $P_k$ and $W_k$ depends only on the ray slopes and refractive indices at surface $k$.


## Chapter 3 — Automated Achromatic Doublet Design

This chapter presents AutoAchromat, an automated pipeline that takes a doublet specification (focal length, aperture, wavelength band) and produces a ranked set of manufacturable thick-lens prescriptions. The pipeline proceeds in two stages: an analytical synthesis (Stage A) that enumerates glass pairs from industrial catalogs and solves the Seidel equations derived in Chapter 2, followed by a constrained numerical optimization (Stage B) that refines the thick-lens model. Thermal stability analysis is integrated throughout. The chapter concludes with a characterization of the inherent design space degeneracy of thin-lens synthesis — a finding that motivates the coma-line sampling strategy discussed as future work in §5.4.

### 3.1 Two-Stage Pipeline Overview

The pipeline implements a two-stage design process:

**Stage A — Analytical Synthesis and Evaluation:**
1. Load glass catalogs (any AGF-format catalog from the Zemax glass database can be imported).
2. Filter glasses by wavelength-range coverage, dispersion-formula validity, and catalog exclusion flags; entries whose computed $n(\lambda_0)$ falls outside the physically plausible range $1.3 < n < 2.6$ are rejected as corrupted.
3. Enumerate all glass pairs satisfying $|\nu_1 - \nu_2| \geq \Delta\nu_{\min}$.
4. For each pair, compute the achromatic power split and solve the Seidel equations.
5. Filter solutions by geometric feasibility (minimum radius, non-overlap for spaced doublets) and by the preliminary evaluation (PE) threshold. Cemented candidates are then Top-$N$ ranked using $|W - W_0|$ as the primary metric and PE as the secondary metric. Air-spaced candidates are not Top-$N$ culled at this stage; all surviving algebraic solutions are sorted by PE, and the subsequent thickening/evaluation stages may process only the first $N$.
6. Compute the first-order thermal metrics ($V_1$, $V_2$, $d\Phi/dT$, $\alpha_{h,\text{required}}$) for each retained candidate, so that thermal information accompanies the candidate through all subsequent stages.
7. Convert thin-lens solutions to thick-lens prescriptions with iterative EFL correction.
8. Apply a thickness-aware Seidel refinement (§3.3.5) so that the aberration balance reflects the actual surface ray heights, not the thin-lens $h \equiv 1$ assumption.
9. Build ray-trace models and evaluate via spot diagram analysis and Seidel coefficient extraction.

**Stage B — Thick-Lens Optimisation:**
1. Deduplicate user-selected Stage A candidates by the rounded optical fingerprint defined in §3.4.1; one representative per fingerprint group is optimized.
2. For each unique candidate, formulate a constrained nonlinear least-squares problem with user-specified thickness bounds per element.
3. Optimize radii, center thicknesses, and air gap (for spaced doublets) using Trust-Region Reflective least squares. Revert to the Stage A starting values if the on-axis polychromatic RMS spot radius worsens by more than 10 %.

### 3.2 Stage A — Thin-Lens Seidel Synthesis

The analytical core of Stage A follows the third-order synthesis methodology originated by Slyusarev and systematized through Trubko's lookup tables. Ivanova et al. [1] formalized the calculation and analysis of cemented components within this framework, and Romanova et al. [2] automated the glass search and radii computation, replacing the original static tables with software. Nguyen and Bakholdin [3] extended the method to air-spaced doublets and introduced the preliminary evaluation (PE) ranking criterion. However, neither the source code nor the detailed algorithmic implementation of these works has been published, so the present pipeline re-derives and re-implements the synthesis from the published equations, adding thick-lens conversion (§3.3), constrained optimization (§3.4), and thermal analysis (§3.5). The cemented and air-spaced synthesis tracks are implemented as independent solvers, because the algebraic variables, feasibility checks, and ranking criteria differ between the two doublet types.

Alternative analytical formulations of the same thin-lens Seidel equations exist in different notation traditions. Mikš and Novák [8] parameterise the shape via an auxiliary variable $\rho_i$ and derive closed-form solutions for superachromatic doublets and triplets; Banerjee and Hazra [7] encode the Coddington shape factor $X$ and glass indices as chromosomes in a genetic algorithm. As noted in §2.2, the shape parameterisations ($Q$, $X$, $\rho$) are related by linear transformations and yield mathematically identical quadratic equations for the lens shape when $S_I = S_{II} = 0$ [7, 8]. The present work retains the Slyusarev $Q$-parameterisation for consistency with the ITMO lineage [1–3] and because the exhaustive enumeration strategy eliminates the need for stochastic search.

#### 3.2.1 Glass Pair Enumeration

The pipeline begins by loading glass catalogs in the industry-standard AGF (ANSI Glass Format) format, supporting multiple manufacturers (e.g. SCHOTT, OHARA, CDGM). The AGF catalogs include several dispersion representations — among them the Sellmeier, Schott polynomial, Herzberger, Conrady, and Handbook of Optics forms; the implementation explicitly evaluates the formula classes encountered in the present catalog set (fallback handling is documented in Appendix A). For each glass entry, the reader extracts the dispersion formula type and coefficients, transmission wavelength limits, thermo-optic (TD) coefficients and coefficients of thermal expansion (CTE) for thermal analysis, relative cost, and a status flag indicating whether the glass is preferred, standard, or obsolete.

Before pairing, individual glasses are filtered by the following criteria: (1) glasses flagged as excluded by the manufacturer are removed, (2) the dispersion formula must exist with valid coefficients and cover the user-specified wavelength range, and (3) the refractive index at the design wavelength must pass a sanity check ($1.3 < n < 2.6$) to guard against corrupted catalog entries.

All ordered pairs of surviving glasses are then tested against the minimum Abbe number difference constraint $|\nu_1 - \nu_2| \geq \Delta\nu_{\min}$. AutoAchromat uses $\Delta\nu_{\min} = 10$ as the default, exposed as a user-configurable screening parameter rather than a theoretical hard limit. This threshold ensures that the achromatic power split produces element powers that are not excessively large. The implementation uses the general form with a chromatic correction target offset $C_0$:

$$\varphi_1 = \frac{\nu_1(1 - \nu_2 C_0)}{\nu_1 - \nu_2}, \qquad \varphi_2 = 1 - \varphi_1$$

The simplified expressions in §2.1 correspond to the $C_0 = 0$ case. A large $|\nu_1 - \nu_2|$ yields milder element powers, reducing higher-order aberrations and improving manufacturing feasibility. Pairs passing this filter proceed to the Seidel synthesis stages (§3.2.2 for cemented doublets, §3.2.3 for air-spaced doublets).

The enumeration strategy — testing all ordered pairs from the full catalog — contrasts with targeted glass selection methods in the literature. Rayces and Rosete-Aguilar [13, 14] select glass pairs for reduced secondary spectrum by ranking the secondary spectrum magnitude across all pairs and applying tolerance conditions on spherochromatism and fifth-order spherical aberration. Sun et al. [6] propose an "illustration method" based on plotting glasses in a combined $P_{AB}$–$P_{dC}$ partial dispersion diagram, selecting outlier glasses (e.g., N-FK58) that maximize the triangle area in the $V_d$–$P_{dC}$ plane (large $G$) while minimizing the triangle area in the $P_{AB}$–$P_{dC}$ plane (small $E$). Dai et al. [15] treat glass selection as a discontinuous optimization problem, embedding discrete glass switching into the merit function of a continuous optimizer. The present pipeline avoids such pre-filtering: by enumerating exhaustively and deferring selection to the ranking stage, it ensures that no viable glass pair is missed due to assumptions about which partial dispersion region is optimal.

#### 3.2.2 Cemented Doublet Synthesis

The cemented doublet synthesis follows the Slyusarev methodology as automated by Romanova et al. [2]. The central variable is the signed Abbe bending invariant $Q$, defined as:

$$Q \equiv \frac{f'}{R_2} - \varphi_1$$

where $R_2$ is the radius of the cemented interface and $f'$ is the system focal length. The quantity $Q$ is the signed Abbe bending invariant used in the Slyusarev synthesis notation. Let

$$\rho_2 \equiv \frac{f'}{R_2}.$$

In the normalized thin-lens coordinates of the synthesis ($h = 1$, $f' = 1$), $\rho_2$ is the normalized curvature of the cemented interface. Substituting the surface-2 slopes used below gives

$$n_1(\rho_2 - u_2) \;=\; n_2(\rho_2 - u_3) \;=\; Q.$$

Thus $Q$ is not an arbitrary bending parameter: it is the Slyusarev cemented-interface bending invariant written in the normalized variables of the synthesis. Relative to the refraction-invariant convention of §2.3.1,

$$\left.\frac{\Delta u}{\Delta(1/n)}\right|_{\text{surface 2}} \;=\; \frac{u_3 - u_2}{1/n_2 - 1/n_1} \;=\; -Q,$$

so the sign difference is only a consequence of the slope/incidence convention. In the cemented doublet this signed invariant also serves as the shape variable, since once $Q$ is chosen all three radii follow directly. It is related to the standard Coddington shape factor $X_1$ of the first element by:

$$Q = \frac{(X_1 - (2n_1 - 1))\,\varphi_1}{2(n_1 - 1)}$$

so that varying $Q$ is equivalent to varying the bending of both elements simultaneously.

The total third-order spherical aberration is quadratic in $Q$:

$$AQ^2 + BQ + C = P$$

with coefficients:

$$A = \left(1 + \frac{2}{n_1}\right)\varphi_1 + \left(1 + \frac{2}{n_2}\right)\varphi_2$$

$$B = \frac{3}{n_1 - 1}\varphi_1^2 - \frac{3}{n_2 - 1}\varphi_2^2 - 2\varphi_2$$

$$C = \frac{n_1}{(n_1 - 1)^2}\varphi_1^3 + \frac{n_2}{(n_2 - 1)^2}\varphi_2^3 + \left(\frac{1}{n_2 - 1} + 1\right)\varphi_2^2$$

For a general target $P=P_0$, this quadratic yields up to two real roots in $Q$ — the two bending solutions that achieve the specified third-order spherical aberration. The special case $P_0 = 0$ corresponds to zero spherical aberration.

The third-order coma is linear in $Q$:

$$W = KQ + L, \qquad K = \frac{A + 1}{2}, \qquad L = \frac{B - \varphi_2}{3}$$

The user-specified coma target $W_0$ enters only through $|W - W_0|$ in the ranking step; setting $W_0 = 0$ corresponds to the aplanatic (coma-free) condition. For each glass pair, both real bending solutions of the quadratic in $Q$ are retained whenever they pass the radius and PE filters. The cemented Top-$N$ selection is then applied globally across all valid (glass pair, $Q$-root) combinations, using $|W - W_0|$ as the primary ranking metric and PE as the secondary metric.

Given $Q$ and the power split, the three surface radii are determined:

$$R_1 = \frac{f'}{\frac{n_1}{n_1 - 1}\varphi_1 + Q}, \quad R_2 = \frac{f'}{\varphi_1 + Q}, \quad R_3 = \frac{f'}{\frac{n_2}{n_2 - 1}\varphi_1 + Q - \frac{1}{n_2 - 1}}$$

**Per-surface Seidel contribution at the cemented interface.** The cemented doublet has three surfaces with the following refractive index transitions:

| Surface | Before ($n$) | After ($n'$) | Slopes |
|:-------:|:---:|:---:|:---:|
| 1 | $1$ | $n_1$ | $u_1 \to u_2$ |
| 2 (cemented) | $n_1$ | $n_2$ | $u_2 \to u_3$ |
| 3 | $n_2$ | $1$ | $u_3 \to u_4$ |

The marginal ray slopes before and after refraction at the cemented surface (surface 2) are:

$$u_2 = Q\!\left(1 - \frac{1}{n_1}\right) + \varphi_1, \qquad u_3 = Q\!\left(1 - \frac{1}{n_2}\right) + \varphi_1$$

Applying the §2.3.2 formula $P_k = (\Delta u / \Delta(1/n))^2\,\Delta(u/n)$ at the cemented surface:

$$P_2 = \left(\frac{u_3 - u_2}{1/n_2 - 1/n_1}\right)^{\!2}\!\left(\frac{u_3}{n_2} - \frac{u_2}{n_1}\right)$$

This per-surface contribution is used in the preliminary evaluation criterion (§3.2.4).

#### 3.2.3 Air-Spaced Doublet Synthesis

The same invariant construction also applies to a single thin element. Physically, a singlet may be regarded as the limiting case in which the second cemented component is replaced by air, but the cemented-doublet formulas above are not evaluated by setting $n_2 = 1$ because those expressions contain terms singular in $n_2 - 1$. Instead, the singlet is written directly as an air–glass–air element with one back-surface bending invariant. An air-spaced doublet is then modeled as two such singlets separated by an air gap, each with its own shape parameter:

$$Q_1 \equiv \frac{f_1}{R_2} - 1, \qquad Q_2 \equiv \frac{f_2}{R_4} - 1,$$

where $f_i = f' / \varphi_i$ is the element focal length. 

Although the two elements are physically separated, their Seidel contributions are not fully independent: the ray state at element 2 depends on the refraction by element 1, introducing cross-terms in $\varphi_1$ and $\varphi_2$. Furthermore, the shape parameter $Q_i$ is defined relative to the element focal length $f_i$ (rather than the system focal length $f'$ as in the cemented case), which shifts the powers of $\varphi_i$ in the coefficients.

The total spherical aberration is quadratic in each shape parameter:

$$P = A_1 Q_1^2 + B_1 Q_1 + A_2 Q_2^2 + B_2 Q_2 + C$$

with coefficients:

$$A_1 = \varphi_1^3\!\left(1 + \frac{2}{n_1}\right), \qquad B_1 = \frac{3\varphi_1^3}{n_1 - 1}$$

$$A_2 = \varphi_2^3\!\left(1 + \frac{2}{n_2}\right), \qquad B_2 = \frac{3\varphi_2^3}{n_2 - 1} - 4\varphi_1\varphi_2^2\!\left(1 + \frac{1}{n_2}\right)$$

$$C = \frac{n_1\varphi_1^3}{(n_1 - 1)^2} + \frac{n_2\varphi_2^3}{(n_2 - 1)^2} - \varphi_1\varphi_2^2\!\left(\frac{4}{n_2 - 1} + 1\right) + \varphi_1^2\varphi_2\!\left(3 + \frac{2}{n_2}\right)$$

The total coma is bilinear in $(Q_1, Q_2)$:

$$W = K_1 Q_1 + K_2 Q_2 + L$$

with:

$$K_1 = \varphi_1^2\!\left(1 + \frac{1}{n_1}\right), \qquad K_2 = \varphi_2^2\!\left(1 + \frac{1}{n_2}\right)$$

$$L = \frac{\varphi_1^2}{n_1 - 1} + \frac{\varphi_2^2}{n_2 - 1} - \varphi_1\varphi_2\!\left(2 + \frac{1}{n_2}\right)$$

Setting $W = W_0$ defines the **coma-zero line** in the $(Q_1, Q_2)$ plane:

$$Q_2 = -\frac{K_1}{K_2} Q_1 + \frac{W_0 - L}{K_2} \equiv aQ_1 + b$$

This reduces the two-dimensional problem to one dimension. Substituting $Q_2 = aQ_1 + b$ into $P = P_0$ yields a quadratic in $Q_1$:

$$(A_1 + A_2 a^2)\,Q_1^2 + (B_1 + 2A_2 ab + B_2 a)\,Q_1 + (A_2 b^2 + B_2 b + C - P_0) = 0$$

Up to two real roots are obtained, each corresponding to a distinct doublet geometry.

The four surface radii are:

$$R_1 = \frac{f_1}{\frac{n_1}{n_1 - 1} + Q_1}, \quad R_2 = \frac{f_1}{1 + Q_1}, \quad R_3 = \frac{f_2}{\frac{n_2}{n_2 - 1} + Q_2}, \quad R_4 = \frac{f_2}{1 + Q_2}$$

A geometric feasibility check ensures that the inner surfaces do not overlap at the aperture edge:

$$d_{\text{air}} + \text{sag}(R_3, a) - \text{sag}(R_2, a) \geq 0$$

where $a = D/2$ is the semi-aperture and $\text{sag}(R, a) = R - \text{sign}(R)\sqrt{R^2 - a^2}$. If both algebraic roots fall into the overlap-forbidden region for the specified air gap, the current implementation rejects the candidate rather than performing a numerical sweep along the coma-zero line; the latter is discussed as future work in §5.4.

**Per-surface Seidel contributions.** The marginal ray is traced through the four surfaces (object at infinity, normalized total power $\Phi = 1$), giving successive slopes:

$$u_1 = 0, \quad u_2 = \frac{(n_1 - 1)\,\rho_1\,\varphi_1}{n_1}, \quad u_3 = \varphi_1, \quad u_4 = \frac{\varphi_1 + (n_2 - 1)\,\rho_3\,\varphi_2}{n_2}, \quad u_5 = 1$$

where $\rho_1 = \frac{n_1}{n_1 - 1} + Q_1$ and $\rho_3 = \frac{n_2}{n_2 - 1} + Q_2$ are the normalized curvatures of the first surfaces of each element. The refractive index transitions at each surface are:

|   Surface    | Before ($n$) | After ($n'$) |    Slopes     |
| :----------: | :----------: | :----------: | :-----------: |
|      1       |     $1$      |    $n_1$     | $u_1 \to u_2$ |
| 2(air space) |    $n_1$     |     $1$      | $u_2 \to u_3$ |
|      3       |     $1$      |    $n_2$     | $u_3 \to u_4$ |
|      4       |    $n_2$     |     $1$      | $u_4 \to u_5$ |

Applying the §2.3.2 formula $P_k = (\Delta u / \Delta(1/n))^2\,\Delta(u/n)$ at each surface:

$$P_1 = \left(\frac{u_2 - u_1}{1/n_1 - 1}\right)^{\!2}\!\left(\frac{u_2}{n_1} - u_1\right), \qquad P_2 = \left(\frac{u_3 - u_2}{1 - 1/n_1}\right)^{\!2}\!\left(u_3 - \frac{u_2}{n_1}\right)$$

$$P_3 = \left(\frac{u_4 - u_3}{1/n_2 - 1}\right)^{\!2}\!\left(\frac{u_4}{n_2} - u_3\right), \qquad P_4 = \left(\frac{u_5 - u_4}{1 - 1/n_2}\right)^{\!2}\!\left(u_5 - \frac{u_4}{n_2}\right)$$

These per-surface contributions are used in the preliminary evaluation criterion (§3.2.4).

#### 3.2.4 Preliminary Evaluation Criterion (PE)

The PE criterion, introduced by Nguyen and Bakholdin [3], provides a scalar figure of merit for ranking candidates beyond the primary aberration targets. While the algebraic synthesis targets the system-level totals ($P_{\text{sys}}$ and $W_{\text{sys}}$), PE examines how the aberration is distributed across individual surfaces.

For the **cemented doublet**, PE combines the per-surface Seidel contribution $P_2$ at the cemented interface (derived in §3.2.2) with a coma penalty:

$$PE = \frac{|P_2| \cdot 3^{|W - W_0|}}{R_2^2}$$

The exponential factor $3^{|W - W_0|}$ amplifies the penalty for designs deviating from the coma target.

For the **air-spaced doublet**, PE is the mean absolute per-surface Seidel contribution across all four surfaces (derived in §3.2.3):

$$PE = \frac{1}{4}\sum_{i=1}^{4} |P_i|$$

Since the algebraic solutions already satisfy $P_{\text{sys}} = P_0$ and $W = W_0$ by construction, PE quantifies how unevenly the aberration is distributed across surfaces — a proxy for higher-order aberration susceptibility and manufacturing sensitivity.

### 3.3 Thin-to-Thick Lens Conversion

Thin-lens theory neglects element thickness. For a 50 mm aperture system at f/4, center thicknesses of 5–10 mm represent 2.5–5% of the focal length. The thick-lens effect systematically shortens the EFL because more strongly curved surfaces are brought closer together by the finite element thickness. A correction procedure is therefore required to convert thin-lens solutions into manufacturable thick-lens prescriptions while preserving the target focal length.

Two principal strategies for thin-to-thick conversion appear in the literature. Wang et al. [5] assign thicknesses from manufacturing standards (the Chinese national standard GB/T 1205–1975 in their case), then adjust one radius per element — chosen to produce the largest possible radius — to hold each element's focal power constant. Sun et al. [6] propose an $\alpha$-preserved strategy: defining $\alpha_i = h_i C_i$ (the product of ray height and curvature at each surface), they show that maintaining the numerical values of all $\alpha_i$ during thickening exactly preserves every element's power contribution at all wavelengths, including chromatic aberration, without iteration. The present pipeline adopts a third approach — uniform radius scaling with iterative EFL correction via the ABCD matrix. Uniform scaling leaves the curvature ratios (and thus the Coddington shape factor $X = (c_1 + c_2)/(c_1 - c_2)$) unchanged, preserving the relative bending of the element surfaces established during synthesis, at the cost of requiring a convergence loop, as detailed below.

The conversion proceeds in three stages: (1) compute the outside diameter, (2) assign physical thicknesses from manufacturing tables, and (3) iteratively correct the radii so the thick-lens EFL matches the target.

#### 3.3.1 Outside Diameter and Thickness Limits

The physical outside diameter of each element exceeds the clear aperture $D$ by an increment $\Delta(D)$ that depends on the mounting method; for retaining-ring mounting, the outside diameter is $\phi = D + \Delta(D)$. Element thicknesses must in turn satisfy manufacturing constraints that depend on lens sign (positive or negative power) and clear aperture: the minimum edge thickness $t_{e,\min}(D)$ is binding for positive elements, and the minimum center thickness $t_{c,\min}(D)$ is binding for negative elements. The implementation uses two retaining-ring mounting lookup tables, labelled as Table 10-2 and Table 10-3 in the source notes; their numerical values are reproduced in Appendix B.

#### 3.3.2 Thickness Allocation

The relationship between center and edge thickness is governed by the sag geometry:

$$t_{\text{edge}} = t_{\text{center}} - \text{sag}(R_{\text{front}}, a) + \text{sag}(R_{\text{back}}, a)$$

where $a = D/2$ is the semi-aperture and $\text{sag}(R, a) = R - \text{sign}(R)\sqrt{R^2 - a^2}$ is the signed spherical sagitta. A surface with $|R| < a$ is geometrically impossible and the candidate is rejected.

The thickness allocation depends on the lens sign, determined by the thin-lens paraxial power $\Phi = (n-1)(1/R_{\text{front}} - 1/R_{\text{back}})$:

- **Positive lens** ($\Phi \geq 0$): the edge is thinner than the center. The minimum edge thickness constraint dominates:

$$t_{\text{center}} = \max\left(t_{e,\min} + \text{sag}(R_{\text{front}}, a) - \text{sag}(R_{\text{back}}, a),\; 0.5\text{ mm}\right)$$

- **Negative lens** ($\Phi < 0$): the center is thinner than the edge. Both the edge and center thickness constraints must be satisfied:

$$t_{\text{center}} = \max\left(t_{e,\min} + \text{sag}(R_{\text{front}}, a) - \text{sag}(R_{\text{back}}, a),\; t_{c,\min,\text{neg}},\; 0.5\text{ mm}\right)$$

The absolute minimum of 0.5 mm is always enforced regardless of lens sign or aperture. Candidates producing non-finite or non-positive thickness values are rejected as unmanufacturable.

#### 3.3.3 EFL Correction via ABCD Matrix

Assigning finite thicknesses to elements that were designed under the thin-lens approximation shifts the system EFL. The paraxial EFL of the thick-lens system is computed using the ABCD (ray-transfer) matrix method.

Using the ABCD refraction and transfer matrices defined in §2.2.2, the composite system matrix is formed by multiplying these matrices in surface order. For a **cemented doublet** (3 surfaces, 2 glass thicknesses):

$$M = M_{\text{refr}}(R_3, n_2{\to}1) \cdot M_{\text{tran}}(t_2, n_2) \cdot M_{\text{refr}}(R_2, n_1{\to}n_2) \cdot M_{\text{tran}}(t_1, n_1) \cdot M_{\text{refr}}(R_1, 1{\to}n_1)$$

For an **air-spaced doublet** (4 surfaces, 2 glass thicknesses, 1 air gap $d$):

$$M = M_{\text{refr}}(R_4, n_2{\to}1) \cdot M_{\text{tran}}(t_2, n_2) \cdot M_{\text{refr}}(R_3, 1{\to}n_2) \cdot M_{\text{tran}}(d, 1) \cdot M_{\text{refr}}(R_2, n_1{\to}1) \cdot M_{\text{tran}}(t_1, n_1) \cdot M_{\text{refr}}(R_1, 1{\to}n_1)$$

The system EFL is extracted from the lower-left element: $f' = -1/C$ where $C = M_{10}$.

#### 3.3.4 Iterative Radius Scaling

To restore the target EFL, all radii are scaled by a uniform factor $k = f'_{\text{target}} / f'_{\text{actual}}$ after each thickness assignment. This procedure is iterated because scaling the radii changes the sag values, which in turn changes the thicknesses, which again affects the EFL. The iteration loop proceeds as:

1. Assign thicknesses to both elements using the current radii.
2. Compute the thick-lens EFL via the ABCD matrix.
3. Compute the correction factor $k = f'_{\text{target}} / f'_{\text{actual}}$.
4. If $|k - 1| < 10^{-9}$ (convergence), stop.
5. Scale all radii: $R_i \leftarrow k \cdot R_i$.
6. Repeat from step 1.

A maximum of 4 correction rounds is applied, corresponding to 5 total evaluations including the final one.

**Key property**: Uniform scaling of all radii leaves curvature ratios unchanged, preserving the Coddington shape factor $X = (c_1 + c_2)/(c_1 - c_2)$, so the relative bending established during thin-lens synthesis is maintained. In contrast, the $\alpha$-preserved strategy of Sun et al. [6] adjusts each surface curvature independently to maintain $\alpha_i = h_i C_i$, preserving per-surface power contribution at all wavelengths (including chromatic aberration) without iteration, but altering $X$ and therefore the bending configuration established during synthesis. The present method trades exact per-surface power preservation for exact shape preservation; the thickness-aware Seidel refinement in §3.3.5 then restores the spherical-aberration target, and any remaining higher-order and finite-aperture residuals are left for the ray-trace evaluation and the optional Stage B optimizer. For air-spaced doublets, the air gap is intentionally **not scaled** — it is a user-specified design constraint, and the achromatic condition is only weakly sensitive to the gap when $d \ll f'$.

**Pathological guards**: The correction is abandoned if $k \leq 0$, $k > 5$, or the EFL is non-finite. These conditions indicate configurations where the thick-lens model diverges qualitatively from the thin-lens design (e.g., a surface becomes nearly hemispherical).

#### 3.3.5 Thickness-Aware Seidel Refinement

After thickness assignment and EFL correction, the thin-lens assumption that all surfaces of an element share the same marginal ray height no longer holds: the ray height $h_k$ at each surface depends on the actual element thickness. Because the synthesis in §3.2 was derived under $h \equiv 1$, the Seidel spherical aberration of the thick prescription no longer matches the target $P_0$. AutoAchromat therefore performs a thickness-aware Seidel refinement before ray-trace evaluation.

In each iteration, the current finite-thickness state is approximated by tracing a paraxial marginal ray through the radii implied by the current bending parameter together with the current element thicknesses. The per-surface Seidel sum

$$P_{\text{thick}} = \sum_{k} h_k \left(\frac{\Delta u}{\Delta(1/n)}\right)_{\!k}^{\!2}\Delta\!\left(\frac{u}{n}\right)_{\!k}$$

is evaluated on this trace, and a single damped Newton step is applied to the bending variable to drive the residual $P_{\text{thick}} - P_0$ toward zero. For cemented doublets, the bending variable is $Q$; for air-spaced doublets, $Q_1$ is refined along the coma-zero line and $Q_2$ follows from $Q_2 = a Q_1 + b$, so the coma target $W_0$ is preserved exactly. After each Newton update, the prescription is re-thickened (new thicknesses, new EFL-corrected radii) so that the iterate remains a valid manufacturable geometry.

The Newton step on $Q$ (or $Q_1$) uses a finite-difference estimate of $dP_{\text{thick}}/dQ$ and is clamped to $|\Delta Q| \leq 0.5\,\max(|Q|, 1)$ to prevent overshoot. The refinement stops when $|P_{\text{thick}} - P_0| < 10^{-6}$, when 10 iterations are exhausted, or when any intermediate step fails (e.g. a degenerate radius denominator, a non-finite trace, or an infeasible thickening); in the latter case the last valid prescription is retained.

### 3.4 Stage B — Thick-Lens Optimization

The Stage B optimization minimizes a weighted sum of squared residuals:

$$\min_{\mathbf{x}} \sum_{k} w_k^2 \left(f_k(\mathbf{x}) - t_k\right)^2$$

**Variables** $\mathbf{x}$:
- Surface radii of curvature (3 for cemented, 4 for spaced), bounded sign-preservingly to preserve lens topology. For each starting radius $R_0$, the maximum allowed magnitude is $R_{\max} = \max(10 f',\; 5|R_0|)$, and the minimum allowed magnitude is $D/2 + 1\text{ mm}$.
- Element center thicknesses (2), with lower bounds incorporating sag correction so that the edge-thickness manufacturing limit of Appendix B is still feasible.
- Air gap (spaced doublets only).

**Operands** $f_k(\mathbf{x})$ and weights:

| Operand | Target $t_k$ | Weight $w_k$ |
|---------|:---:|:---:|
| Effective focal length | $f'$ | 10 |
| On-axis polychromatic RMS spot radius | 0 | 4 |
| Off-axis polychromatic RMS spot radius | 0 | 1 |
| Edge thickness per element (soft inequality) | $\geq t_{e,\min}$ | 5 |

The high weight on EFL preserves the focal length specification. The on-axis spot receives higher weight than the off-axis spot, reflecting the typical priority in doublet design. Edge thickness is enforced through a soft inequality operand evaluated at the specified semi-aperture; the operand's interaction with optiland's surface-aperture state is an implementation detail and is described in Appendix C.

**Solver**: The Trust-Region Reflective (TRF) algorithm is used, with a default maximum of 200 iterations and convergence tolerance $\text{tol} = 10^{-6}$.

**Revert mechanism**: If the on-axis polychromatic RMS spot radius after optimisation exceeds 1.10 times its Stage A starting value, the optimization result is discarded and the Stage A variable values are restored.

#### 3.4.1 Deduplication

Many glass pairs from different manufacturers share nearly identical optical properties (e.g., SCHOTT N-BK7 $\approx$ OHARA S-BSL7 $\approx$ CDGM H-K9L). Candidates promoted to Stage B are grouped by the optical fingerprint

$$\text{fp}(c) = \left(\text{round}(n_1, 3),\; \text{round}(n_2, 3),\; \text{round}(\nu_1, 1),\; \text{round}(\nu_2, 1)\right),$$

i.e. $n$ rounded to three decimals and $\nu$ to one decimal. Only the highest-ranked Stage A representative in each fingerprint group is optimized, which avoids repeated optimization of optically equivalent glass substitutions.

### 3.5 Thermal Analysis

Temperature affects lens performance through two coupled mechanisms that shift the focal length in opposite or reinforcing directions:

1. **Thermo-optic effect** ($dn/dT$): the refractive index of each glass changes with temperature, altering the optical power of every surface.
2. **Mechanical expansion** (CTE $\alpha$): the glass elements expand, changing surface radii and element thicknesses; simultaneously, the housing expands, shifting the image plane position.

A design that is insensitive to uniform temperature changes — **passive athermalization** — requires these effects to cancel at the system level [9, 10]. The present pipeline computes thermal metrics for each retained candidate but deliberately separates them from the optical ranking: it reports $V_1$, $V_2$, the normalised system power drift $d\Phi/dT$, and the required housing CTE $\alpha_{h,\text{required}}$ alongside optical merit without filtering. These exported quantities are sufficient to compute the residual thermal defocus in post-processing once a housing CTE is chosen, using the expression derived in §3.5.2. This preserves designer flexibility over the performance–thermal–cost trade-off, whose priorities are project-specific.

#### 3.5.1 Thermal Model

The Schott first-order thermal dispersion model computes the rate of change of refractive index with temperature for a given glass at wavelength $\lambda$:

$$\frac{dn}{dT} = \frac{n^2 - 1}{2n} \left[D_0 + \frac{E_0}{\lambda^2 - \lambda_{tk}^2}\right]$$

where $D_0$, $E_0$, and $\lambda_{tk}$ are material-specific parameters from the glass catalog. This is the starting point: knowing how fast each glass's index drifts with temperature.

For a single thin lens element, two effects change its focal length simultaneously: the index change ($dn/dT$) alters the optical power, and the thermal expansion ($\alpha$) scales the surface radii. These combine into a single scalar — the **thermo-optical coefficient**:

$$V = \frac{dn/dT}{n - 1} - \alpha$$

The first term $\frac{dn/dT}{n-1}$ is the fractional power change due to the index shift (since power $\varphi \propto (n-1)/R$). The second term $\alpha$ is the fractional power change due to radii scaling ($R \to R(1+\alpha\Delta T)$, so $\varphi$ decreases). Thus $V > 0$ means that the element's power increases with temperature, while $V < 0$ means it decreases.

For a doublet with normalized powers $\varphi_1$, $\varphi_2$ ($\varphi_1 + \varphi_2 = 1$), the normalized system power drift is the weighted sum of each element's thermo-optical coefficient:

$$\frac{d\Phi}{dT}\bigg|_{\text{norm}} = V_1\varphi_1 + V_2\varphi_2$$

This is an element-power thermal model. Each coefficient $V_i$ is computed for one glass element using the air–glass–air thin-lens expression, and the system drift is then formed from the normalized element powers $\varphi_i$. For a cemented doublet, the paraxial power decomposition

$$\Phi = (n_1 - 1)(c_1 - c_2) + (n_2 - 1)(c_2 - c_3) = \varphi_1 + \varphi_2$$

is exact in the thin-lens limit $d = 0$, so the refractive contribution of the glass–glass interface is absorbed into the shared $c_2 = 1/R_2$ terms of $\varphi_1$ and $\varphi_2$. The model does not, however, assign an independent thermo-mechanical deformation law to the physical cemented interface: differential expansion of the two glasses, adhesive stress, and interface deformation are mechanical reliability effects outside the scope of this first-order passive-athermal analysis.

A complementary surface-resolved decomposition that makes the glass–glass interface contribution explicit, including the role of the shared cemented curvature $R_2$, is derived in Appendix D. It is presented as a theoretical diagnostic framework rather than a replacement for the element-power model used by the current thermal metric, and is not currently implemented in the pipeline.

Since $f' \approx 1/\Phi$, a positive system power drift shortens the focal length as temperature rises, while a negative drift lengthens it. In the sign convention of the implementation, many glass pairs give a negative system drift, which leads to a positive required housing CTE (derived in §3.5.2).

#### 3.5.2 Athermalization Assessment

The lens is mounted in a housing that also expands with temperature (CTE $\alpha_h$), shifting the image plane away from the lens. Passive athermalization is achieved when the housing expansion exactly compensates the focal length drift, requiring:

$$\alpha_{h,\text{required}} = -\frac{d\Phi}{dT}\bigg|_{\text{norm}}$$

If $\alpha_{h,\text{required}}$ matches a practical housing material (e.g., aluminium at $\sim$23.6 ppm/K), the design is passively athermalized without any active compensation. If no match is found, the designer must either accept a thermal defocus budget or choose a different glass pair.

When the actual housing CTE $\alpha_h$ does not match $\alpha_{h,\text{required}}$, the residual focal shift over a temperature change $\Delta T$ is:

$$\delta f' = f' \cdot (\alpha_{h,\text{required}} - \alpha_h) \cdot \Delta T$$

This is the quantity that ultimately matters for system performance: it determines whether the thermal defocus falls within the depth of focus budget (typically $\pm 2\lambda F_\#^2$ for diffraction-limited systems). For air-spaced doublets, the finite air gap introduces an additional degree of freedom that expands the athermal solution space [10].
