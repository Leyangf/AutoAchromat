# Pseudocode Blocks

本文件按照在 `thesisEN.tex` 中出现的顺序汇总了所有伪代码（algorithm 环境）块，并改写为 Markdown 友好的伪代码格式。

约定：
- `←` 表示赋值
- `⊥` 表示拒绝/无效返回
- `∅` 表示空集，`∪` 集合并集，`∈` 属于
- `▷` 表示行尾注释
- 控制结构（`for`、`if/else`、`return` 等）以加粗英文关键字呈现
- 数学符号尽量使用 Unicode；下标采用 `x_y` 形式以便阅读

---

## 1. Algorithm: Complete Stage A synthesis and thickening


```text
Algorithm 1: Complete Stage A synthesis and thickening
─────────────────────────────────────────────────────────────────────────────
Input : glass list G; specification
        (f', D, λ_0, λ_1, λ_2, Δν_min, C_0, P_0, W_0, max_PE,
         N, system_type, ...)
Output: ranked Stage A candidate list S_A with refined thick prescriptions
─────────────────────────────────────────────────────────────────────────────
 1:  G' ← FilterGlasses(G, [λ_0, λ_1, λ_2])                      ▷ Sec. glass-enum
 2:  for each g ∈ G' do
 3:      compute n(λ_0) and ν(λ_0, λ_1, λ_2)
 4:  end for
 5:  S_cand ← ∅
 6:  for each ordered pair (g_1, g_2) ∈ G' × G' with g_1 ≠ g_2 do
 7:      if |ν_1 − ν_2| < Δν_min then continue
 8:      (φ_1, φ_2) ← AchromatPower(ν_1, ν_2, C_0)
 9:      if system_type = cemented then
10:          Λ ← CementedSolve(n_1, n_2, φ_1, φ_2, P_0, W_0,
                               f', D, max_PE)                    ▷ Alg. 2
11:      else
12:          Λ ← SpacedSolve(n_1, n_2, φ_1, φ_2, P_0, W_0,
                             f', D, d_air, max_PE)               ▷ Alg. 3
13:      end if
14:      if Λ = ∅ then continue                                  ▷ skip thermal when no candidate
15:      θ ← ThermalMetrics(g_1, g_2, n_1, n_2, φ_1, φ_2, λ_0)   ▷ depends only on pair
16:      for each cand ∈ Λ do
17:          cand.glasses ← (g_1, g_2);  cand.thermal ← θ
18:          if system_type = cemented then
19:              UpdateTopN(S_cand, cand, key = (|W − W_0|, PE),
                            capacity = N)                        ▷ online: replace worst
20:          else
21:              S_cand ← S_cand ∪ {cand}                        ▷ defer top-N to line 28
22:          end if
23:      end for
24:  end for
25:  if system_type = cemented then
26:      sort S_cand by (|W − W_0|, PE)                          ▷ already capacity-N
27:  else
28:      sort S_cand ascending by PE;  truncate to first N
29:  end if
30:  S_A ← ∅
31:  for each cand ∈ S_cand do
32:      rx_0 ← Thicken(cand, inputs)                            ▷ Alg. 4
33:      if rx_0 = ⊥ then continue
34:      rx ← RefineSeidel(rx_0, cand, inputs)                   ▷ Alg. 5
35:      op ← BuildOptic(rx)                                     ▷ Optiland model
36:      if op = ⊥ then continue
37:      cand.metrics ← Evaluate(op, cand, inputs, rx)           ▷ ray trace
38:      cand.rx ← rx
39:      S_A ← S_A ∪ {cand}
40:  end for
41:  return S_A
```

---

## 2. Algorithm: Cemented doublet synthesis


```text
Algorithm 2: Cemented doublet synthesis
─────────────────────────────────────────────────────────────────────────────
Input : n_1, n_2, φ_1, φ_2, P_0, W_0, f', D, max_PE
Output: list of feasible cemented candidates
─────────────────────────────────────────────────────────────────────────────
 1:  compute A, B, C from (n_1, n_2, φ_1, φ_2)                   ▷ Eq. cemented-quadratic
 2:  roots ← real roots of  A·Q² + B·Q + (C − P_0) = 0
 3:  out ← ∅
 4:  for each Q ∈ roots do
 5:      W ← K·Q + L                                             ▷ Eq. cemented-coma
 6:      (R_1, R_2, R_3) ← RadiiCemented(f', n_1, n_2, φ_1, Q)   ▷ Eq. cemented-radii
 7:      if min(|R_1|, |R_2|, |R_3|) < 1.05 · D/2 then continue
 8:      compute u_2, u_3, P_2 at the cemented interface         ▷ App. per-surface-seidel
 9:      PE ← ComputePE(P_2, W, W_0, R_2)                        ▷ Sec. pe-criterion
10:      if PE > max_PE then continue
11:      out ← out ∪ { (Q, W, R_1, R_2, R_3, P_2, PE) }
12:  end for
13:  return out
```

---

## 3. Algorithm: Air-spaced doublet synthesis

```text
Algorithm 3: Air-spaced doublet synthesis
─────────────────────────────────────────────────────────────────────────────
Input : n_1, n_2, φ_1, φ_2, P_0, W_0, f', D, d_air, max_PE
Output: list of feasible air-spaced candidates
─────────────────────────────────────────────────────────────────────────────
 1:  compute A_1, B_1, A_2, B_2, C, K_1, K_2, L                  ▷ Sec. spaced-synthesis coefficients
 2:  if |K_2| ≈ 0 then return ∅
 3:  compute coma-target line (a, b) such that  Q_2 = a·Q_1 + b  ▷ Sec. spaced-synthesis
 4:  compute quadratic coefficients q_0, q_1, q_2 along this line ▷ Sec. spaced-synthesis
 5:  roots ← real roots of  q_2·Q_1² + q_1·Q_1 + q_0 = 0
 6:  out ← ∅;   a_sag ← D/2
 7:  for each Q_1 ∈ roots do
 8:      Q_2 ← a·Q_1 + b
 9:      (R_1, R_2, R_3, R_4) ← RadiiSpaced(f', n_1, n_2,
                                            φ_1, φ_2, Q_1, Q_2)  ▷ App. per-surface-seidel
10:      if min_i |R_i| < 1.05 · D/2 then continue
11:      if d_air + sag(R_3, a_sag) − sag(R_2, a_sag) < 0
            then continue                                        ▷ inner-surface overlap
12:      (P_1, ..., P_4) ← PerSurfacePk(n_1, n_2, φ_1, φ_2,
                                        Q_1, Q_2)                ▷ App. per-surface-seidel
13:      PE ← ComputePE(P_1, P_2, P_3, P_4)
14:      if PE > max_PE then continue
15:      out ← out ∪ { (Q_1, Q_2, R_1, ..., R_4,
                        P_1, ..., P_4, PE) }
16:  end for
17:  return out
```

---

## 4. Algorithm: Inner thin-to-thick EFL correction


```text
Algorithm 4: Inner thin-to-thick EFL correction
─────────────────────────────────────────────────────────────────────────────
Input : candidate (radii, glasses, n_1, n_2, signs);
        f'_target, D, d_air, system_type
Output: ThickPrescription, or ⊥ (rejection)
─────────────────────────────────────────────────────────────────────────────
 1:  R ← candidate.radii
 2:  for k ← 0 to 4 do                                           ▷ up to 5 evaluations
 3:      for each element i ∈ {1, 2} do
 4:          s_f ← sag(R_front,i, D/2)
 5:          s_b ← sag(R_back,i,  D/2)
 6:          if element i is positive (φ_i ≥ 0) then
 7:              t_i ← max( t_e,min(D) + s_f − s_b,
                            0.5 mm )
 8:          else
 9:              t_i ← max( t_e,min(D) + s_f − s_b,
                            t_c,min,neg(D),
                            0.5 mm )
10:          end if
11:          if t_i is not finite or t_i ≤ 0 then return ⊥
12:      end for
13:      if system_type = spaced  AND
            d_air + sag(R_3, D/2) − sag(R_2, D/2) < 0 then
14:          return ⊥                                            ▷ inner-surface overlap (rechecked after scaling)
15:      end if
16:      f'_actual ← EflABCD(R, t_1, t_2, d_air, n_1, n_2)       ▷ ABCD matrix
17:      if f'_actual is not finite or |f'_actual| < 10⁻⁹
            then break                                           ▷ keep last valid
18:      κ ← f'_target / f'_actual
19:      if κ ≤ 0 or κ > 5 then break                            ▷ pathological, keep last valid
20:      if |κ − 1| < 10⁻⁹ then break                            ▷ converged
21:      R ← κ · R                                               ▷ uniform scaling
22:  end for
23:  return ThickPrescription(R, t_1, t_2, d_air, f'_actual)
```

---

## 5. Algorithm: Outer thickness-aware Seidel refinement


```text
Algorithm 5: Outer thickness-aware Seidel refinement
─────────────────────────────────────────────────────────────────────────────
Input : ThickPrescription rx_0, synthesis candidate c, inputs
Output: refined ThickPrescription rx
─────────────────────────────────────────────────────────────────────────────
 1:  Q ← c.Q   (or Q_1 for air-spaced)
 2:  rx ← rx_0
 3:  for k ← 1 to 10 do
 4:      P_thick ← ThickSeidelP(rx)                              ▷ paraxial marginal trace + Σ h_k·A_k²·Δ(u/n)_k
 5:      r ← P_thick − P_0                                       ▷ residual
 6:      if |r| < 10⁻⁶ then return rx
 7:      estimate dP/dQ ← (ThickSeidelP(rx | Q+δ) − P_thick) / δ,
            δ = 10⁻⁷                                             ▷ t held fixed; thickening's t-dependence on Q treated as 2nd-order
 8:      if |dP/dQ| < 10⁻¹⁵ then return rx                       ▷ singular Jacobian
 9:      ΔQ ← − r / (dP/dQ)
10:      clamp:  |ΔQ| ≤ 0.5 · max(|Q|, 1)                        ▷ damping
11:      Q_new ← Q + ΔQ
12:      rebuild radii from Q_new                                ▷ Q_2 ← a·Q_1 + b for air-spaced
13:      rx_new ← Thicken(candidate(Q_new), inputs)              ▷ Alg. 4
14:      if rx_new = ⊥  or trace fails  then return rx           ▷ keep last valid
15:      Q ← Q_new;   rx ← rx_new
16:  end for
17:  return rx                                                   ▷ iteration budget exhausted
```

---

## 6. Algorithm: Stage B thick-lens optimization


```text
Algorithm 6: Stage B thick-lens optimization
─────────────────────────────────────────────────────────────────────────────
Input : Stage A optic op_A, prescription rx_A, inputs
Output: optimized optic op_B (or op_A on revert)
─────────────────────────────────────────────────────────────────────────────
 1:  problem ← OptimizationProblem()
 2:  for each surface i with finite R_i do
 3:      R_max ← max( 10·f',  5·|R_i| )
 4:      add variable  R_i ∈ sign(R_i) · [ D/2 + 1 mm,  R_max ]  ▷ sign-preserving
 5:  end for
 6:  for each element j ∈ {1, 2} do
 7:      mfg_min ← t_e,min + sag(R_f,A) − sag(R_b,A)             ▷ at Stage A radii
 8:      (t_min, t_max) ← clamp( user_bounds_j,  mfg_min )
 9:      add variable  t_j ∈ [ t_min, t_max ]
10:  end for
11:  if rx_A.system_type = spaced then
12:      add variable  d_air ∈ [ gap_min, gap_max ]
13:  end if
                                                                 ▷ operands (target, weight)
14:  add operand  f'           = f'_target,    w = 10
15:  add operand  RMS_onaxis   = 0,            w = 4             ▷ polychromatic
16:  add operand  RMS_offaxis  = 0,            w = 1             ▷ polychromatic
17:  for each element j do
18:      add operand  edge_thickness_j ≥ t_e,min,  w = 5         ▷ soft inequality
19:  end for
20:  RMS_0 ← Evaluate(RMS_onaxis)
21:  x_0   ← Snapshot(problem.variables)
22:  run TRF least-squares,  maxiter = 200,  tol = 10⁻⁶
23:  RMS*  ← Evaluate(RMS_onaxis)
24:  if RMS* > 1.10 · RMS_0 then restore variables ← x_0          ▷ Sec. stage-b revert rule
25:  recompute image plane  (op.ImageSolve())
26:  return op
```
