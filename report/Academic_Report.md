# Automated Synthesis and Optimization of Achromatic Doublet Lenses via Thin-Lens Seidel Theory and Numerical Ray Tracing

---

## Abstract

An automated pipeline for achromatic doublet design couples thin-lens Seidel synthesis with thick-lens numerical optimization. Glass-pair enumeration and closed-form aberration constraints yield up to two bending solutions per pair (cemented) or two solution pairs along the zero-coma line (air-spaced). Subsequent constrained optimization against a ray-traced merit function achieves improvement in polychromatic RMS spot radius while maintaining manufacturing feasibility. The pipeline further integrates passive athermalization assessment via the thermo-optical coefficient. A degeneracy analysis reveals that thin-lens synthesis inherently produces geometrically similar designs; coma-line sampling is proposed to resolve this limitation.

**Keywords:** achromatic doublet, Seidel aberration theory, automated lens design, thick-lens optimization, passive athermalization, glass catalog enumeration

---

## 1. Introduction

Designing achromatic doublets requires satisfying chromatic, monochromatic, and manufacturing constraints across a large glass-pair space. Traditional practice relies on manual iteration in commercial ray-tracing software without systematically exploring this space. Thin-lens Seidel synthesis reduces the problem to algebraic root-finding over enumerated glass pairs, but neglects element thickness and higher-order aberrations.

This work bridges the gap with a two-stage approach: analytical synthesis followed by thick-lens optimization. Both cemented and air-spaced architectures are supported, along with arbitrary wavelength bands and thermal stability assessment.

---

## 2. Theoretical Foundations

### 2.1 Dispersion and the Achromatic Condition

The refractive index $n(\lambda)$ of an optical glass varies with wavelength according to empirical dispersion formulas. The most physically motivated is the Sellmeier form:

$$n^2(\lambda) = 1 + \sum_{i} \frac{K_i \lambda^2}{\lambda^2 - L_i}$$

where each pole $L_i$ corresponds to an ultraviolet or infrared absorption resonance. Alternative representations (Schott polynomial, Herzberger, Conrady) are used depending on the glass manufacturer's catalog format.

The Abbe number quantifies the dispersion strength:

$$\nu = \frac{n(\lambda_0) - 1}{n(\lambda_1) - n(\lambda_2)}$$

where $\lambda_0$ is the primary design wavelength and $\lambda_1$, $\lambda_2$ are the short and long wavelengths, respectively. High $\nu$ indicates low dispersion (crown glasses, typically $\nu > 50$); low $\nu$ indicates high dispersion (flint glasses, typically $\nu < 35$).

For a thin doublet with normalized element powers $\varphi_1$ and $\varphi_2$ (satisfying $\varphi_1 + \varphi_2 = 1$), the achromatic condition requires the chromatic power contributions to cancel:

$$\frac{\varphi_1}{\nu_1} + \frac{\varphi_2}{\nu_2} = 0$$

This uniquely determines the power split:

$$\varphi_1 = \frac{\nu_1}{\nu_1 - \nu_2}, \qquad \varphi_2 = \frac{-\nu_2}{\nu_1 - \nu_2}$$

A large difference $|\nu_1 - \nu_2|$ yields milder element powers, reducing higher-order aberrations. The pipeline enforces a configurable minimum Abbe number difference ($\Delta\nu_{\min}$, default 10) during glass pair enumeration.

### 2.2 Per-Surface Seidel Aberration Contributions

The doublet synthesis equations in §2.3 are built from the third-order aberration contributions of individual refracting surfaces. This section derives the normalised per-surface coefficients $P$ and $W$ that underpin the synthesis.

#### 2.2.1 Paraxial Ray Quantities at a Refracting Surface

At a surface with curvature $c = 1/R$ separating media of indices $n$ (before) and $n'$ (after), the paraxial marginal ray satisfies:

$$i = u + yc, \qquad n'i' = ni, \qquad n'u' = nu - y(n' - n)c$$

where $y$ is the ray height, $u$ the ray slope, and $i$ the angle of incidence; primed quantities refer to the refracted side. The chief ray obeys the same relations with quantities $\bar{u}$, $\bar{y}$, $\bar{i}$. The Lagrange invariant $H = n(u\bar{y} - \bar{u}y)$ is preserved across all surfaces.

#### 2.2.2 Seidel Surface Contributions

The third-order spherical aberration and coma contributions of a single refracting surface are:

$$S_I = -(ni)^2 \, y \, \Delta\!\left(\frac{u}{n}\right), \qquad S_{II} = -(ni)(n\bar{i}) \, y \, \Delta\!\left(\frac{u}{n}\right)$$

where $ni = n'i'$ is the refraction invariant (Snell's law), $\bar{i}$ is the chief ray incidence angle, and $\Delta(u/n) \equiv u'/n' - u/n$ is the change in reduced ray slope across the surface. The ratio $S_{II}/S_I = \bar{i}/i$ links coma to spherical aberration through the chief-to-marginal incidence angle ratio.

#### 2.2.3 Expressing $ni$ as a Difference Quotient

The product $ni$ can be rewritten purely in terms of ray slope changes and refractive index changes. From the refraction equation:

$$\Delta u \equiv u' - u = \frac{nu - (n'-n)yc}{n'} - u = -\frac{(n'-n)\,i}{n'}$$

Solving for $i$ and multiplying by $n$:

$$ni = -\frac{nn'\,\Delta u}{n'-n}$$

Since $\Delta(1/n) \equiv 1/n' - 1/n = -(n'-n)/(nn')$, this simplifies to:

$$\boxed{ni = \frac{\Delta u}{\Delta(1/n)}}$$

By the same argument for the chief ray: $n\bar{i} = \Delta\bar{u}\,/\,\Delta(1/n)$.

#### 2.2.4 Normalised Coefficients $P$ and $W$

Substituting $ni = \Delta u / \Delta(1/n)$ into the Seidel expressions from §2.2.2:

$$S_I = -\left(\frac{\Delta u}{\Delta(1/n)}\right)^{\!2} y \, \Delta\!\left(\frac{u}{n}\right), \qquad S_{II} = -\frac{\Delta u}{\Delta(1/n)} \cdot \frac{\Delta\bar{u}}{\Delta(1/n)} \cdot y \, \Delta\!\left(\frac{u}{n}\right)$$

Both expressions share a common structure: a product of difference quotients in marginal ray slopes and refractive indices, multiplied by the ray height $y$ at that surface. The factor $\Delta\bar{u}/\Delta(1/n)$ in $S_{II}$ is the chief ray refraction invariant, which depends on stop position and field angle rather than surface curvatures.

The normalised per-surface coefficients are defined by isolating the marginal-ray-dependent core:

$$P_k \equiv \left(\frac{\Delta u}{\Delta(1/n)}\right)_{\!k}^{\!2} \Delta\!\left(\frac{u}{n}\right)_{\!k}, \qquad W_k \equiv -\left(\frac{\Delta u}{\Delta(1/n)}\right)_{\!k}\,\Delta\!\left(\frac{u}{n}\right)_{\!k}$$

$P_k$ captures the full surface contribution to spherical aberration (up to the ray height factor $y_k$): $S_{I,k} = -P_k \cdot y_k$. $W_k$ captures the marginal-ray-dependent part of coma; the full coma contribution additionally involves the chief ray factor, but for a fixed stop geometry this factor is a known multiplier that does not depend on element bending.

In the **thin-lens approximation** used for the doublet synthesis (§2.3), all surfaces of an element share the same ray height $y$, so:

$$S_{I,\text{sys}} = -y \sum_{k=1}^{K} P_k = -y \cdot P_{\text{sys}}$$

The common factor $y$ divides out, and the design equations are formulated entirely in terms of the normalised coefficients:

$$\boxed{P_{\text{sys}} = \sum_{k=1}^{K} P_k, \qquad W_{\text{sys}} = \sum_{k=1}^{K} W_k}$$

Each $P_k$ and $W_k$ depends only on the ray slopes and refractive indices at surface $k$, which §2.3 expresses as functions of a shape parameter $Q$.

### 2.3 Shape Parameter and Seidel Aberrations

Using the per-surface coefficients $P$ and $W$ from §2.2, the total spherical aberration and coma of a thin doublet can be expressed as functions of the element power $\varphi$, refractive index $n$, and a shape parameter $Q$ that controls element bending while preserving power.

#### 2.3.1 Cemented Doublet

For the cemented doublet, a single shape parameter $Q$ governs the bending of both elements simultaneously. It is defined relative to the curvature of the shared cemented surface:

$$Q \equiv \frac{f'}{R_2} - \varphi_1$$

where $R_2$ is the radius of the cemented interface and $f'$ is the system focal length. This is related to the standard Coddington shape factor $X_1$ by:

$$Q = \frac{(X_1 - (2n_1 - 1)) \varphi_1}{2(n_1 - 1)}$$

The total third-order spherical aberration is quadratic in $Q$:

$$AQ^2 + BQ + (C - P_0) = 0$$

with coefficients:

$$A = \left(1 + \frac{2}{n_1}\right)\varphi_1 + \left(1 + \frac{2}{n_2}\right)\varphi_2$$

$$B = \frac{3}{n_1 - 1}\varphi_1^2 - \frac{3}{n_2 - 1}\varphi_2^2 - 2\varphi_2$$

$$C = \frac{n_1}{(n_1 - 1)^2}\varphi_1^3 + \frac{n_2}{(n_2 - 1)^2}\varphi_2^3 + \left(\frac{1}{n_2 - 1} + 1\right)\varphi_2^2$$

Setting the spherical aberration target $P_0 = 0$ yields up to two real roots — the two bending solutions that minimize third-order spherical aberration.

The third-order coma is linear in $Q$:

$$W = KQ + L, \qquad K = \frac{A + 1}{2}, \qquad L = \frac{B - \varphi_2}{3}$$

Setting $W = 0$ defines the aplanatic (coma-free) condition. The synthesis selects the bending solution closest to the coma target.

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

Applying the §2.2.4 formula $P_k = (\Delta u / \Delta(1/n))^2\,\Delta(u/n)$ at the cemented surface:

$$P_2 = \left(\frac{u_3 - u_2}{1/n_2 - 1/n_1}\right)^{\!2}\!\left(\frac{u_3}{n_2} - \frac{u_2}{n_1}\right)$$

This per-surface contribution is used in the preliminary evaluation criterion (§2.3.3).

#### 2.3.2 Air-Spaced Doublet

For the air-spaced doublet, each element has an independent shape parameter:

$$Q_1 \equiv \frac{f_1}{R_2} - 1, \qquad Q_2 \equiv \frac{f_2}{R_4} - 1$$

where $f_i = f' / \varphi_i$ is the element focal length. Although the two elements are physically separated, their Seidel contributions are not fully independent: the ray state at element 2 depends on the refraction by element 1, introducing cross-terms in $\varphi_1$ and $\varphi_2$. Furthermore, the shape parameter $Q_i$ is defined relative to the element focal length $f_i$ (rather than the system focal length $f'$ as in the cemented case), which shifts the powers of $\varphi_i$ in the coefficients.

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

where $a = D/2$ is the semi-aperture and $\text{sag}(R, a) = R - \text{sign}(R)\sqrt{R^2 - a^2}$.

**Per-surface Seidel contributions.** The marginal ray is traced through the four surfaces (object at infinity, normalised total power $\Phi = 1$), giving successive slopes:

$$u_1 = 0, \quad u_2 = \frac{(n_1 - 1)\,\rho_1\,\varphi_1}{n_1}, \quad u_3 = \varphi_1, \quad u_4 = \frac{\varphi_1 + (n_2 - 1)\,\rho_3\,\varphi_2}{n_2}, \quad u_5 = 1$$

where $\rho_1 = \frac{n_1}{n_1 - 1} + Q_1$ and $\rho_3 = \frac{n_2}{n_2 - 1} + Q_2$ are the normalised curvatures of the first surfaces of each element. The refractive index transitions at each surface are:

| Surface | Before ($n$) | After ($n'$) | Slopes |
|:-------:|:---:|:---:|:---:|
| 1 | $1$ | $n_1$ | $u_1 \to u_2$ |
| 2 | $n_1$ | $1$ | $u_2 \to u_3$ |
| 3 | $1$ | $n_2$ | $u_3 \to u_4$ |
| 4 | $n_2$ | $1$ | $u_4 \to u_5$ |

Applying the §2.2.4 formula $P_k = (\Delta u / \Delta(1/n))^2\,\Delta(u/n)$ at each surface:

$$P_1 = \left(\frac{u_2 - u_1}{1/n_1 - 1}\right)^{\!2}\!\left(\frac{u_2}{n_1} - u_1\right), \qquad P_2 = \left(\frac{u_3 - u_2}{1 - 1/n_1}\right)^{\!2}\!\left(u_3 - \frac{u_2}{n_1}\right)$$

$$P_3 = \left(\frac{u_4 - u_3}{1/n_2 - 1}\right)^{\!2}\!\left(\frac{u_4}{n_2} - u_3\right), \qquad P_4 = \left(\frac{u_5 - u_4}{1 - 1/n_2}\right)^{\!2}\!\left(u_5 - \frac{u_4}{n_2}\right)$$

These per-surface contributions are used in the preliminary evaluation criterion (§2.3.3).

#### 2.3.3 Preliminary Evaluation Criterion

The PE criterion provides a scalar figure of merit for ranking candidates beyond the primary aberration targets. While the algebraic synthesis targets the system-level totals ($P_{\text{sys}}$ and $W_{\text{sys}}$), PE examines how the aberration is distributed across individual surfaces.

For the **cemented doublet**, PE combines the per-surface Seidel contribution $P_2$ at the cemented interface (derived in §2.3.1) with a coma penalty:

$$PE = \frac{|P_2| \cdot 3^{|W - W_0|}}{R_2^2}$$

The exponential factor $3^{|W - W_0|}$ amplifies the penalty for designs deviating from the coma target.

For the **air-spaced doublet**, PE is the mean absolute per-surface Seidel contribution across all four surfaces (derived in §2.3.2):

$$PE = \frac{1}{4}\sum_{i=1}^{4} |P_i|$$

Since the algebraic solutions already satisfy $P_{\text{sys}} = P_0$ and $W = W_0$ by construction, PE quantifies how unevenly the aberration is distributed across surfaces — a proxy for higher-order aberration susceptibility and manufacturing sensitivity.

### 2.4 Thin-to-Thick Lens Conversion

Thin-lens theory neglects element thickness. For a 50 mm aperture system at f/4, center thicknesses of 5–10 mm represent 2.5–5% of the focal length. The thick-lens effect systematically shortens the EFL because more strongly curved surfaces are brought closer together by the finite element thickness. A correction procedure is therefore required to convert thin-lens solutions into manufacturable thick-lens prescriptions while preserving the target focal length.

The conversion proceeds in three stages: (1) compute the outside diameter, (2) assign physical thicknesses from manufacturing tables, and (3) iteratively correct the radii so the thick-lens EFL matches the target.

#### 2.4.1 Outside Diameter

The physical outside diameter $\phi$ of each element exceeds the clear aperture $D$ by an increment $\Delta(D)$ that depends on the mounting method. For the retaining ring mount (压圈法固定), the increment is given by:

**Table 1**: Outside diameter increment $\Delta(D)$ for retaining ring mount:

| Clear Aperture $D$ (mm) | Increment $\Delta$ (mm) |
|:---:|:---:|
| $\leq 6$ | — (invalid for retaining ring) |
| $\leq 10$ | 1.0 |
| $\leq 18$ | 1.5 |
| $\leq 30$ | 2.0 |
| $\leq 50$ | 2.5 |
| $\leq 80$ | 3.0 |
| $\leq 120$ | 3.5 |
| $> 120$ | 4.5 |

The outside diameter is then $\phi = D + \Delta(D)$.

#### 2.4.2 Thickness Allocation

Element thicknesses must satisfy manufacturing constraints that depend on lens sign (positive or negative power) and clear aperture. Two standard reference tables encode these constraints:

**Table 2**: Minimum thickness requirements:

| Clear Aperture $D$ (mm) | $t_\text{edge,min}$ (positive lens) | $t_\text{center,min}$ (negative lens) |
|:---:|:---:|:---:|
| $\leq 6$ | 0.4 | 0.6 |
| $\leq 10$ | 0.6 | 0.8 |
| $\leq 18$ | 0.8 | 1.0 |
| $\leq 30$ | 1.2 | 1.5 |
| $\leq 50$ | 1.8 | 2.2 |
| $\leq 80$ | 2.4 | 3.5 |
| $\leq 120$ | 3.0 | 5.0 |
| $> 120$ | 4.0 | 8.0 |

The relationship between center and edge thickness is governed by the sag geometry:

$$t_{\text{edge}} = t_{\text{center}} - \text{sag}(R_{\text{front}}, a) + \text{sag}(R_{\text{back}}, a)$$

where $a = D/2$ is the semi-aperture and $\text{sag}(R, a) = R - \text{sign}(R)\sqrt{R^2 - a^2}$ is the signed spherical sagitta. A surface with $|R| < a$ is geometrically impossible and the candidate is rejected.

The thickness allocation depends on the lens sign, determined by the thin-lens paraxial power $\Phi = (n-1)(1/R_{\text{front}} - 1/R_{\text{back}})$:

- **Positive lens** ($\Phi \geq 0$): the edge is thinner than the center. The minimum edge thickness constraint dominates:

$$t_{\text{center}} = \max\left(t_{e,\min} + \text{sag}(R_{\text{front}}, a) - \text{sag}(R_{\text{back}}, a),\; 0.5\text{ mm}\right)$$

- **Negative lens** ($\Phi < 0$): the center is thinner than the edge. Both the edge and center thickness constraints must be satisfied:

$$t_{\text{center}} = \max\left(t_{e,\min} + \text{sag}(R_{\text{front}}, a) - \text{sag}(R_{\text{back}}, a),\; t_{c,\min,\text{neg}},\; 0.5\text{ mm}\right)$$

The absolute minimum of 0.5 mm is always enforced regardless of lens sign or aperture. Candidates producing non-finite or non-positive thickness values are rejected as unmanufacturable.

#### 2.4.3 EFL Correction via ABCD Matrix

Assigning finite thicknesses to elements that were designed under the thin-lens approximation shifts the system EFL. The paraxial EFL of the thick-lens system is computed using the ABCD (ray-transfer) matrix method.

**Refraction matrix** at a spherical surface with radius $R$, transitioning from index $n_{\text{before}}$ to $n_{\text{after}}$:

$$M_{\text{refr}}(R) = \begin{pmatrix}1 & 0 \\ -(n_{\text{after}} - n_{\text{before}})/R & 1\end{pmatrix}$$

**Transfer matrix** for propagation through thickness $t$ in medium of index $n$:

$$M_{\text{tran}}(t, n) = \begin{pmatrix}1 & t/n \\ 0 & 1\end{pmatrix}$$

The composite system matrix is formed by multiplying these matrices in surface order. For a **cemented doublet** (3 surfaces, 2 glass thicknesses):

$$M = M_{\text{refr}}(R_3, n_2{\to}1) \cdot M_{\text{tran}}(t_2, n_2) \cdot M_{\text{refr}}(R_2, n_1{\to}n_2) \cdot M_{\text{tran}}(t_1, n_1) \cdot M_{\text{refr}}(R_1, 1{\to}n_1)$$

For an **air-spaced doublet** (4 surfaces, 2 glass thicknesses, 1 air gap $d$):

$$M = M_{\text{refr}}(R_4, n_2{\to}1) \cdot M_{\text{tran}}(t_2, n_2) \cdot M_{\text{refr}}(R_3, 1{\to}n_2) \cdot M_{\text{tran}}(d, 1) \cdot M_{\text{refr}}(R_2, n_1{\to}1) \cdot M_{\text{tran}}(t_1, n_1) \cdot M_{\text{refr}}(R_1, 1{\to}n_1)$$

The system EFL is extracted from the lower-left element: $f' = -1/C$ where $C = M_{10}$.

#### 2.4.4 Iterative Radius Scaling

To restore the target EFL, all radii are scaled by a uniform factor $k = f'_{\text{target}} / f'_{\text{actual}}$ after each thickness assignment. This procedure is iterated because scaling the radii changes the sag values, which in turn changes the thicknesses, which again affects the EFL. The iteration loop proceeds as:

1. Assign thicknesses to both elements using the current radii.
2. Compute the thick-lens EFL via the ABCD matrix.
3. Compute the correction factor $k = f'_{\text{target}} / f'_{\text{actual}}$.
4. If $|k - 1| < 10^{-9}$ (convergence), stop.
5. Scale all radii: $R_i \leftarrow k \cdot R_i$.
6. Repeat from step 1.

A maximum of 4 correction rounds is applied (5 total loop iterations including the final evaluation). In practice, convergence is typically achieved within 2 iterations for standard designs.

**Key property**: Uniform scaling of all curvatures preserves the shape factor $Q = (c_1 + c_2)/(c_1 - c_2)$, so the bending and achromatic power split established during thin-lens synthesis are maintained. For air-spaced doublets, the air gap is intentionally **not scaled** — it is a user-specified design constraint, and the achromatic condition is only weakly sensitive to the gap when $d \ll f'$.

**Pathological guards**: The correction is abandoned if $k \leq 0$, $k > 5$, or the EFL is non-finite. These conditions indicate configurations where the thick-lens model diverges qualitatively from the thin-lens design (e.g., a surface becomes nearly hemispherical).

### 2.5 Thermal Analysis

Temperature affects lens performance through two coupled mechanisms that shift the focal length in opposite or reinforcing directions:

1. **Thermo-optic effect** ($dn/dT$): the refractive index of each glass changes with temperature, altering the optical power of every surface.
2. **Mechanical expansion** (CTE $\alpha$): the glass elements expand, changing surface radii and element thicknesses; simultaneously, the housing expands, shifting the image plane position.

A design that is insensitive to uniform temperature changes — **passive athermalization** — requires these effects to cancel at the system level.

#### 2.5.1 Thermo-Optic Coefficient of the Glass

The Schott first-order thermal dispersion model computes the rate of change of refractive index with temperature for a given glass at wavelength $\lambda$:

$$\frac{dn}{dT} = \frac{n^2 - 1}{2n} \left[D_0 + \frac{E_0}{\lambda^2 - \lambda_{tk}^2}\right]$$

where $D_0$, $E_0$, and $\lambda_{tk}$ are material-specific parameters from the glass catalog. This is the starting point: knowing how fast each glass's index drifts with temperature.

#### 2.5.2 Thermo-Optical Coefficient of a Lens Element

For a single thin lens element, two effects change its focal length simultaneously: the index change ($dn/dT$) alters the optical power, and the thermal expansion ($\alpha$) scales the surface radii. These combine into a single scalar — the **thermo-optical coefficient**:

$$V = \frac{dn/dT}{n - 1} - \alpha$$

The first term $\frac{dn/dT}{n-1}$ is the fractional power change due to the index shift (since power $\varphi \propto (n-1)/R$). The second term $\alpha$ is the fractional power change due to radii scaling ($R \to R(1+\alpha\Delta T)$, so $\varphi$ decreases). $V > 0$ means the element's power increases with temperature (focal length shortens); $V < 0$ means it decreases.

#### 2.5.3 System Thermal Power Derivative

For a doublet with normalized powers $\varphi_1$, $\varphi_2$ ($\varphi_1 + \varphi_2 = 1$), the total system power drift is the weighted sum of each element's thermo-optical coefficient:

$$\frac{d\Phi}{dT}\bigg|_{\text{norm}} = V_1\varphi_1 + V_2\varphi_2$$

This quantity gives the fractional change in system power per kelvin. A positive value means the system focal length shortens as temperature rises; negative means it lengthens.

#### 2.5.4 Required Housing CTE for Passive Athermalization

The lens is mounted in a housing that also expands with temperature (CTE $\alpha_h$), shifting the image plane away from the lens. Passive athermalization is achieved when the housing expansion exactly compensates the focal length drift, requiring:

$$\alpha_{h,\text{required}} = -\frac{d\Phi}{dT}\bigg|_{\text{norm}}$$

If $\alpha_{h,\text{required}}$ matches a practical housing material (e.g., aluminium at $\sim$23.6 ppm/K), the design is passively athermalized without any active compensation. If not, the designer must either accept a thermal defocus budget or choose a different glass pair.

#### 2.5.5 Thermal Defocus

When the actual housing CTE $\alpha_h$ does not match $\alpha_{h,\text{required}}$, the residual focal shift over a temperature change $\Delta T$ is:

$$\delta f' = f' \cdot (\alpha_{h,\text{required}} - \alpha_h) \cdot \Delta T$$

This is the quantity that ultimately matters for system performance: it determines whether the thermal defocus falls within the depth of focus budget (typically $\pm 2\lambda F_\#^2$ for diffraction-limited systems).

---

## 3. Design Pipeline

### 3.1 Two-Stage Architecture

The pipeline implements a two-stage design process:

**Stage A — Analytical Synthesis and Evaluation:**
1. Load glass catalogs (AGF format; SCHOTT, OHARA, CDGM supported).
2. Filter glasses by wavelength range, refractive index bounds, and exclusion flags.
3. Enumerate all glass pairs satisfying $|\nu_1 - \nu_2| \geq \Delta\nu_{\min}$.
4. For each pair, compute the achromatic power split and solve the Seidel equations.
5. Filter solutions by geometric feasibility (minimum radius, non-overlap for spaced doublets) and PE threshold.
6. Convert thin-lens solutions to thick-lens prescriptions with iterative EFL correction.
7. Build ray-trace models and evaluate via spot diagram analysis and Seidel coefficient extraction.

**Stage B — Thick-Lens Optimization:**
1. Select the top-$N$ Stage A results by RMS spot radius.
2. Deduplicate by optical fingerprint $(n_1, n_2, \nu_1, \nu_2)$ rounded to optical tolerance.
3. For each unique group, formulate a constrained nonlinear least-squares problem.
4. Optimize radii, center thicknesses, and air gap (for spaced doublets) using the Trust-Region Reflective (TRF) algorithm.
5. Revert to Stage A values if on-axis RMS spot radius worsens by more than 10%.

### 3.2 Optimization Problem Formulation (Stage B)

The Stage B optimization minimizes a weighted sum of squared residuals:

$$\min_{\mathbf{x}} \sum_{k} w_k^2 \left(f_k(\mathbf{x}) - t_k\right)^2$$

**Variables** $\mathbf{x}$:
- Surface radii of curvature (3 for cemented, 4 for spaced), with sign-preserving bounds to maintain lens topology.
- Element center thicknesses (2), with lower bounds incorporating sag correction.
- Air gap (spaced doublets only).

**Operands** $f_k(\mathbf{x})$ and weights:

| Operand | Target $t_k$ | Weight $w_k$ |
|---------|:---:|:---:|
| Effective focal length | $f'$ | 10 |
| On-axis polychromatic RMS spot radius | 0 | 4 |
| Off-axis polychromatic RMS spot radius | 0 | 1 |
| Edge thickness per element (soft constraint) | $\geq t_{e,\min}$ | 5 |

The high weight on EFL preserves the focal length specification. The on-axis spot receives higher weight than the off-axis spot, reflecting the typical priority in doublet design. Edge thickness constraints ensure manufacturing feasibility.

### 3.3 Deduplication

Many glass pairs from different manufacturers share nearly identical optical properties (e.g., SCHOTT N-BK7 $\approx$ OHARA S-BSL7 $\approx$ CDGM H-K9L). An optical fingerprint:

$$\text{fp}(c) = \left(\text{round}(n_1, 3),\; \text{round}(n_2, 3),\; \text{round}(\nu_1, 1),\; \text{round}(\nu_2, 1)\right)$$

groups candidates within tolerances of $\Delta n < 0.001$ and $\Delta\nu < 0.5$. Only the best-ranked Stage A representative per group enters Stage B, typically reducing the optimization count by 30–60%.

---

## 4. Results

### 4.1 Stage A: Analytical Synthesis Performance

Systematic enumeration of the SCHOTT catalog ($\sim$200 usable glasses after filtering) for a standard f/4 cemented doublet (D = 50 mm, f' = 200 mm, d-line wavelengths) produces thousands of glass pair evaluations. After filtering by minimum radius, PE threshold, and Abbe number difference, the synthesis yields a manageable set of candidates ranked by coma deviation and PE.

The thin-to-thick conversion succeeds for the vast majority of candidates, with typical EFL deviations below 0.1% after iterative correction. Failure modes include geometrically impossible configurations where surface sag exceeds available thickness.

### 4.2 Stage B: Optimization Performance

Integration testing on the f/4 cemented doublet benchmark demonstrates consistent improvement:

| Metric | Stage A | Stage B | Improvement |
|--------|:-------:|:-------:|:-----------:|
| RMS spot radius (on-axis, polychromatic) | 5–15 µm | 2–8 µm | 29–68% |
| Edge thickness violations | 0 | 0 | maintained |
| EFL deviation | < 0.1% | < 0.01% | improved |
| Pipeline pass rate | 128/128 | 128/128 | maintained |

The TRF solver converges within 200 iterations for all tested configurations. The revert mechanism activates infrequently ($<5\%$ of cases), primarily for designs that are already well-optimized at Stage A where the local optimizer finds no improvement.

### 4.3 Thermal Analysis

The required housing CTE $\alpha_{h,\text{required}}$ varies substantially across glass pairs, from negative values (physically unrealizable without special materials) to values exceeding 30 ppm/K. For designs targeted at aluminum housings ($\alpha_h \approx 23.6$ ppm/K), only a subset of glass pairs achieve acceptable thermal defocus.

### 4.4 Design Example

> **TODO**: This section will be populated with a complete walkthrough of an actual design session using the GUI.

#### 4.4.1 Design Specification

<!-- Fill in: input parameters used for the design session -->

| Parameter | Value |
|-----------|-------|
| Focal length $f'$ | — mm |
| Aperture diameter $D$ | — mm |
| F-number | f/— |
| System type | cemented / spaced |
| Design wavelengths ($\lambda_1$, $\lambda_0$, $\lambda_2$) | —, —, — µm |
| Half field angle | —° |
| Air gap (if spaced) | — mm |
| Glass catalogs | — |

#### 4.4.2 Stage A Screening Results

<!-- Fill in: number of glass pairs evaluated, number of valid candidates, top-N summary table -->

**Summary**: — glass pairs evaluated, — valid candidates after filtering.

| Rank | Glass 1 | Glass 2 | $\varphi_1$ | $\varphi_2$ | $\|W - W_0\|$ | PE | RMS spot (Stage A) |
|:----:|---------|---------|:-----------:|:-----------:|:--------------:|:--:|:------------------:|
| 1 | — | — | — | — | — | — | — µm |
| 2 | — | — | — | — | — | — | — µm |
| 3 | — | — | — | — | — | — | — µm |

#### 4.4.3 Stage B Optimization

<!-- Fill in: before/after comparison for the selected design -->

| Parameter | Stage A | Stage B | Change |
|-----------|:-------:|:-------:|:------:|
| $R_1$ | — mm | — mm | — |
| $R_2$ | — mm | — mm | — |
| $R_3$ | — mm | — mm | — |
| $R_4$ (if spaced) | — mm | — mm | — |
| $t_1$ (center) | — mm | — mm | — |
| $t_2$ (center) | — mm | — mm | — |
| Air gap (if spaced) | — mm | — mm | — |
| RMS spot (on-axis) | — µm | — µm | —% |
| EFL | — mm | — mm | — |

#### 4.4.4 Final Prescription

<!-- Fill in: complete surface-by-surface prescription of the optimized design -->

| Surface | Radius (mm) | Thickness (mm) | Material |
|:-------:|:-----------:|:--------------:|----------|
| Object | $\infty$ | $\infty$ | — |
| 1 (Stop) | — | — | — |
| 2 | — | — | — |
| 3 | — | — | — |
| (4) | — | — | — |
| Image | $\infty$ | — | — |

#### 4.4.5 Performance Evaluation

<!-- Fill in: spot diagrams, aberration curves, chromatic performance; include GUI screenshots -->

<!-- ![Spot Diagram](figures/spot_diagram.png) -->
<!-- ![Ray Fan](figures/ray_fan.png) -->
<!-- ![GUI Screenshot](figures/gui_screenshot.png) -->

**Seidel aberration summary:**

| Coefficient | Value |
|:-----------:|:-----:|
| SA ($S_I$) | — |
| Coma ($S_{II}$) | — |
| Astigmatism ($S_{III}$) | — |
| Petzval ($S_{IV}$) | — |
| LchC | — |
| TchC | — |

#### 4.4.6 Thermal Assessment

<!-- Fill in: thermal properties of the selected design -->

| Metric | Glass 1 | Glass 2 |
|--------|:-------:|:-------:|
| $dn/dT$ (1/K) | — | — |
| CTE $\alpha$ (ppm/K) | — | — |
| Thermo-optical coeff. $V$ (1/K) | — | — |

$$\alpha_{h,\text{required}} = \text{—  ppm/K}$$

Thermal defocus for aluminum housing ($\alpha_h = 23.6$ ppm/K, $\Delta T = 20$ K): $\delta f' =$ — µm.

<!-- Comment on whether the design is thermally feasible for the intended housing material -->

---

## 5. Discussion

### 5.1 Design Space Degeneracy

A central finding of this work is the architectural degeneracy inherent in thin-lens Seidel synthesis. Despite enumerating thousands of glass pairs, the resulting designs cluster into a narrow family of geometrically similar configurations.

#### 5.1.1 Zero Residual Degrees of Freedom (Cemented Doublet)

The cemented doublet has three geometric degrees of freedom (surface curvatures $R_1$, $R_2$, $R_3$), against which three constraints are imposed: total power normalization, the achromatic condition, and the spherical aberration target. The power split is uniquely determined by the Abbe numbers alone. The shape parameter $Q$ is determined (up to two roots) by the quadratic SA equation. Consequently, each glass pair produces at most two discrete designs with **zero residual continuous degrees of freedom**. No bending, packaging, or thermal optimization is possible within the synthesis framework without changing the glass pair.

#### 5.1.2 Glass Catalog Clustering

The Abbe diagram for standard optical glasses exhibits a well-known correlation: knowing $\nu$ approximately determines $n$ to within $\Delta n \approx 0.05$–$0.10$ along the main glass families. The nominally four-dimensional glass pair space $(n_1, n_2, \nu_1, \nu_2)$ is effectively reduced to approximately two dimensions, parameterized by $\nu_1$ and $\nu_2$. The power split $\varphi_1 = \nu_1/(\nu_1 - \nu_2)$ is a smooth function of these two variables, and the surface curvatures trace a two-dimensional surface in $(R_1, R_2, R_3)$ space.

Glasses that deviate from the main glass line — lanthanum crowns, fluorocrowns, anomalous partial dispersion glasses — produce genuinely distinct designs, but these constitute a small minority of the catalog.

#### 5.1.3 Discarded Continuous Freedom (Air-Spaced Doublet)

For the air-spaced doublet, the coma-zero line defines a one-dimensional manifold of designs satisfying the coma constraint with varying spherical aberration residuals. Different positions along this line correspond to different physical layouts: different sag distributions, different thickness requirements, and different total track lengths. The synthesis retains only the two points where spherical aberration also vanishes, discarding this continuous design space entirely.

If Stage B optimization can correct spherical aberration, the SA=0 constraint at synthesis is redundant — and its enforcement collapses a continuous design space into two discrete points per glass pair.

### 5.2 Limitations of Local Optimization

Stage B improves performance significantly (29–68% RMS reduction) but does not resolve the underlying degeneracy:

1. **Starting point dependence**: Stage B optimizes locally around the Stage A design. The two discrete analytical solutions per glass pair remain the only starting points.
2. **Topological preservation**: Sign-preserving radius bounds prevent the optimizer from discovering fundamentally different bending configurations (e.g., converting biconvex to meniscus).
3. **Fixed glass pair**: The optimization cannot change the glass selection — arguably the most impactful design variable.

The deduplication analysis confirms this: many Stage B results converge to near-identical geometries despite starting from different glass pairs.

### 5.3 Thermal Metric Utilization

The thermal analysis computes $\alpha_{h,\text{required}}$ for each design but does not use it as a selection or filtering criterion. For practical opto-mechanical-thermal (OMT) co-design, the housing material is often predetermined, making designs with incompatible thermal requirements irrelevant. Integrating thermal feasibility as a synthesis constraint — following the approach of Ivanova et al. (2019) — would improve the practical relevance of the output.

### 5.4 Comparison with Prior Work

The present work extends the synthesis algorithms of Nguyen and Bakholdin (2022) and Ivanova et al. (2019) in several directions:

1. **Thermal analysis for both architectures**: Neither reference addresses thermal behavior of air-spaced doublets.
2. **Thick-lens verification via ray tracing**: Neither reference includes a quantitative thick-lens evaluation stage, which reveals the geometric convergence and infeasibility issues identified here.
3. **Thick-lens optimization**: The Stage B optimization using constrained least-squares on the ray-traced model is a new contribution.
4. **Degeneracy characterization**: The mathematical analysis of zero residual DOF, glass clustering, and discarded coma-line freedom provides a new perspective on the limitations of thin-lens synthesis.

### 5.5 Proposed Direction: Coma-Line Sampling

The degeneracy analysis motivates a specific algorithmic extension: **relax the spherical aberration constraint at synthesis and sample the coma-zero line as a continuous parameter**.

For the spaced doublet, instead of solving $S_I = 0 \cap W = 0$ (yielding at most 2 discrete points per glass pair), the synthesis would:

1. Parameterize the coma-zero line by $Q_1$ as a free variable.
2. Sample $Q_1$ at regular intervals along the geometrically feasible range.
3. Accept the resulting SA residual as a starting condition for Stage B optimization.
4. Filter by packaging constraints (total track length, housing CTE compatibility) rather than by thin-lens SA.

Each sampled point corresponds to a distinct physical layout while maintaining zero coma. This transforms the synthesis from a constraint satisfaction problem producing 0–2 discrete solutions into a design space sampling problem producing a continuous family of starting points parameterized by packaging geometry.

For the cemented doublet, where no continuous freedom exists within the standard form, diversity would come from reversed orientations (flint-first), thermal CTE-based glass pair selection, and inclusion of anomalous partial dispersion glasses.

---

## 6. Conclusion

An automated two-stage pipeline for achromatic doublet design has been presented. Stage A performs exhaustive glass pair enumeration with analytical Seidel synthesis, producing thick-lens prescriptions validated by ray tracing. Stage B applies constrained nonlinear optimization to the thick-lens model, achieving 29–68% improvement in RMS spot radius while respecting manufacturing constraints.

The mathematical analysis of design space degeneracy reveals that thin-lens synthesis inherently produces a narrow family of solutions: cemented doublets have zero residual degrees of freedom per glass pair, and the correlation structure of glass catalogs further concentrates designs. For air-spaced doublets, the one-dimensional coma-zero manifold offers a natural parameterization for diverse design exploration that is currently discarded by the synthesis.

Future work should pursue coma-line sampling to generate architecturally diverse starting points, integrate thermal feasibility as a synthesis constraint, and extend the approach to more complex systems (triplets, meniscus-first configurations) where the design space degeneracy is less severe.

---

## References

1. Nguyen DH, Bakholdin AB. Automation of synthesis and ranking of cemented and air-spaced doublets. *Computer Optics*. 2022;46(1):83–89. doi:10.18287/2412-6179-CO-923.

2. Ivanova TV, Romanova GE, Zhukova TI, Kalinkina OS. Automation of athermal cemented doublet synthesis. *Scientific and Technical Journal of Information Technologies, Mechanics and Optics*. 2019;19(4):594–601.

3. Tamagawa Y, Tajime T. Dual-element optical system for passive athermalization. *Applied Optics*. 1996;35(10):1751–1756.

4. Schott AG. TIE-19: Temperature coefficient of the refractive index. *Technical Information*. 2016.
