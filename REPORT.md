# NACA 2412 vs NACA 0012: Thin-Airfoil-Theory Analysis

All numbers below are **model estimates** produced by `main.py`. None are experimental measurements.

## 1. Problem

How does camber change the lift and drag-polar behavior of an airfoil? We compare a cambered airfoil (NACA 2412) with a symmetric one (NACA 0012) of equal thickness (12 %).

## 2. Geometry (`geometry.py`)

NACA 4-digit code MPTT: m = M/100 (max camber), p = P/10 (its chordwise position), t = TT/100 (max thickness), chord c = 1.

- Thickness (empirical NACA fit): y_t = 5t(0.2969√x − 0.1260x − 0.3516x² + 0.2843x³ − 0.1036x⁴)
- Camber line, two parabolic arcs joined at x = p:
  - x < p: y_c = (m/p²)(2px − x²)
  - x ≥ p: y_c = (m/(1−p)²)((1−2p) + 2px − x²)
- Surfaces are offset perpendicular to the camber line, with θ = arctan(dy_c/dx): x_u = x − y_t sinθ, y_u = y_c + y_t cosθ (lower surface with opposite signs).

Verification (`check_geometry`): NACA 2412 gives camber 0.0200 at x = 0.400 and thickness 0.1200 at x = 0.300, with a closed trailing edge. These match the definition of the code.

## 3. Lift (`aerodynamics.py`)

**Physical principle.** In thin airfoil theory the airfoil is replaced by a vortex sheet on the camber line. The flow-tangency condition fixes the vortex strength, and the Kutta condition fixes the circulation.

**Assumptions.** 2D, inviscid, incompressible, small angles, thin airfoil.

**Result.** With x = (1 − cosθ)/2, the tangency condition is a Fourier problem. Its first two coefficients give

- C_l = 2π(α − α_L0)
- α_L0 = −(1/π) ∫₀^π (dy_c/dx)(cosθ − 1) dθ

**Interpretation.** The lift slope 2π per radian does not depend on camber. Camber only shifts α_L0. For a symmetric airfoil dy_c/dx = 0, so α_L0 = 0.

**Finite wing (Prandtl lifting line).** a = a₀ / (1 + a₀/(π e AR)) with a₀ = 2π. For AR = 8 and e = 0.85 this gives a = 4.855 per radian. Here AR and e are assumed inputs.

**Numerical result.** NACA 2412: α_L0 = −2.077°, so C_L(α = 0°) = 0.176. NACA 0012: α_L0 = 0, so C_L(0°) = 0. Quoted values for the 2412 are about −2°; check this against a cited source before using it in a paper.

## 4. Drag

**Assumed model.** C_D = C_D0 + C_L²/(π e AR).

- The second term is induced drag. It exists only for a finite wing, not for a 2D airfoil.
- **C_D0 = 0.010 is an assumed placeholder for both airfoils.** It is not computed by this model and not measured.

## 5. Efficiency

Setting d(C_L/C_D)/dC_L = 0 gives C_L* = √(C_D0 π e AR) and (L/D)max = 1/(2√(C_D0 k)), with k = 1/(π e AR).

Numerical check: both airfoils give C_L* = 0.462 and (L/D)max = 23.11. The simulation and the closed form agree.

| | α_L0 | C_L(0°) | (L/D)max | α at (L/D)max |
|---|---|---|---|---|
| NACA 2412 | −2.08° | 0.176 | 23.11 | 3.4° |
| NACA 0012 | 0° | 0 | 23.11 | 5.5° |

## 6. Findings

1. **Supported by the model:** a cambered airfoil produces positive lift at α = 0°. Its lift curve is shifted by |α_L0| ≈ 2.1° relative to the symmetric one.
2. **Supported by the model:** the best-glide point occurs at the same C_L for both airfoils, but at a smaller angle of attack for the cambered one (3.4° versus 5.5°).
3. **Not supported by this model:** that the NACA 2412 has a *higher* peak L/D. In a parabolic polar with equal C_D0, (L/D)max is independent of camber. Any real difference would come from airfoil-specific C_D0, which depends on Reynolds number and surface condition.
4. **Not modeled:** stall. The linear lift curve has no maximum, so the drag rise and efficiency loss at high angle of attack seen in real wings do not appear here.

## 7. Limitations

- Linear regime only (α about −4° to 12°).
- Thickness does not affect C_L in thin airfoil theory.
- No viscosity, Reynolds number, Mach number or ground effect.
- C_D0, AR and e are assumptions.

## 8. Next steps (requires data not yet available)

- Obtain viscous polars for both airfoils (XFOIL or a wind-tunnel dataset, with cited Reynolds number) and replace the assumed C_D0 with airfoil-specific values.
- Compare the thin-airfoil C_L(α) against that data and quantify the deviation.
- Add a stall model or restrict conclusions to the pre-stall range.

## 9. Viscous Benchmark Comparison (XFOIL) - TEMPLATE, TO BE COMPLETED AFTER RUNNING `validate_viscous.py`

**Status:** not yet run. No XFOIL data have been added to `data/`. Every number below is a placeholder.

**What is compared.** The v1 thin-airfoil model against XFOIL polars from Airfoil Tools at Re = [RE, from the downloaded file] and Ncrit = [NCRIT, from the page]. XFOIL is itself a model (viscous panel method), so this is a model-to-model comparison, not an experimental validation.

**2D to 3D mapping.** Each XFOIL point (alpha_2D, c_l, c_d) is mapped to a finite wing of aspect ratio AR = [AR] with span efficiency e = [e] using alpha_3D = alpha_2D + c_l/(pi e AR), C_L = c_l, C_D = c_d + c_l^2/(pi e AR). This assumes elliptic spanwise loading and neglects 3D viscous effects, sweep and taper.

**Results.** Paste the contents of `section8_table.md` here.

| Metric | v1 inviscid theory | XFOIL viscous | Error / divergence |
|---|---|---|---|
| [from section8_table.md] | | | |

![Lift Curve Validation](figures/cl_alpha_validation.png)

![L/D and Drag Polar Comparison](figures/ld_polar_comparison.png)

**Interpretation (write after reading the table).** Thin-airfoil theory assumes attached, inviscid flow, so agreement is expected only at small and moderate angles of attack. At higher angles the adverse pressure gradient on the suction side thickens the boundary layer and can cause separation, which limits c_l and raises c_d; the inviscid model contains neither effect. Whether the NACA 2412 has a higher (L/D)max than the NACA 0012 must be taken from the computed table, not assumed. State the size of the difference relative to the angle step of the data and note that XFOIL results near stall may not be fully converged.

## 10. Numerical Lifting-Line, Experimental Comparison and C++ Benchmark

### 10.1 Numerical lifting-line solver (`lifting_line.py`) - computed results

**Principle.** Prandtl lifting-line theory represents the wing by a bound vortex whose circulation varies along the span. **Assumptions:** incompressible, inviscid, small angles, unswept untwisted wing, section lift slope a0 = 2 pi, symmetric loading.

**Equation.** With Gamma = 2 b V sum A_n sin(n theta), y = -(b/2) cos(theta), and the local effective angle alpha_eff = alpha - alpha_L0 - alpha_i, where alpha_i = sum n A_n sin(n theta)/sin(theta), the section relation Gamma = (1/2) V c a0 alpha_eff gives the monoplane equation at N stations theta_i = i pi/(2N):

sum over odd n of A_n sin(n theta_i) [ 4b/(a0 c(theta_i)) + n/sin(theta_i) ] = alpha - alpha_L0.

The solution gives C_L = pi AR A_1, delta = sum_{n>=3} n (A_n/A_1)^2, e = 1/(1 + delta), C_Di = C_L^2/(pi e AR).

**Verification (computed).** An elliptic planform returns e = 1.000000000000 and C_L equal to the analytic value a alpha with a = a0/(1 + a0/(pi AR)). For N = 10, 20, 40, 80, 160 the rectangular-wing e at AR = 8 is 0.936737, 0.936674, 0.936671, 0.936670, 0.936670, so N = 40 is converged to about 1e-6.

**Span efficiency, AR = 8, N = 40, alpha_L0 = 0 (model results):**

| Planform | e | C_L at 5 deg | C_Di at 5 deg |
|---|---|---|---|
| Elliptic | 1.0000 | 0.4386 | 0.00766 |
| Taper 1.0 (rectangular) | 0.9367 | 0.4222 | 0.00757 |
| Taper 0.8 | 0.9575 | 0.4269 | 0.00757 |
| Taper 0.6 | 0.9761 | 0.4313 | 0.00758 |
| Taper 0.5 | 0.9831 | 0.4332 | 0.00760 |
| Taper 0.4 | 0.9872 | 0.4346 | 0.00761 |
| Taper 0.2 | 0.9755 | 0.4339 | 0.00768 |

![Spanwise loading](figures/spanwise_loading.png)

**Consequence for v1/v2.** The value e = 0.85 used there was an assumption. For this planform family the solver gives e between 0.937 and 1.000 at AR = 8. Note also that v1 used a single e in both the lift-slope formula and the drag formula. For non-elliptic wings the lift-slope correction and the induced-drag factor are different quantities, so v1's slope for such wings is approximate; the solver's C_L is the direct result. These conclusions apply only to the stated model (no twist, no sweep, a0 = 2 pi).

### 10.2 Wind-tunnel comparison (`validate_experimental.py`) - NOT YET RUN

**Status:** no wind-tunnel data have been added; all values below are placeholders.

**Required conditions.** Wind-tunnel data (source: [CITATION, e.g. NACA Report 824, figure/table numbers]) and XFOIL polars must share the same Reynolds number [RE], and ideally comparable Ncrit/turbulence level and surface condition. Airfoil Tools polars stop at Re = 1e6, so a match with higher-Re tunnel data requires running XFOIL yourself at that Re. If the Reynolds numbers differ, report it and do not attribute the RMSE to model error alone. State the digitization uncertainty of any values read from plots: [UNCERTAINTY].

Paste `section10_experimental_table.md` here.

| Airfoil | Region | N pts | RMSE(XFOIL vs tunnel) | RMSE(theory vs tunnel) | c_l,max XFOIL | c_l,max tunnel |
|---|---|---|---|---|---|---|
| [from script] | | | | | | |

![Experimental validation](figures/experimental_validation.png)

### 10.3 Python vs C++ benchmark - RE-RUN ON YOUR MACHINE

`python benchmark.py` compiles `cpp/lifting_line.cpp`, checks that C_L, e and C_Di agree with `lifting_line.py`, and times 10,000 solves (100 angles x 100 taper ratios, N = 40).

| Quantity | Value |
|---|---|
| Max difference, C++ vs Python (rect., AR = 8, 5 deg) | [from output] |
| Python runtime (numpy.linalg.solve) | [PYTHON_S] s |
| C++ -O3 runtime (hand-written Gaussian elimination) | [CPP_S] s |
| Speedup | [SPEEDUP] x |
| Machine / compiler | [from cpp_benchmark_results.csv] |

**Interpretation.** Timings are machine- and implementation-dependent. NumPy's solver calls optimized LAPACK, whereas the C++ code is a plain elimination, so the speedup measures these two implementations and not the languages in general. Because e depends only on the planform, the sweep could be done with 100 solves; the benchmark repeats all 10,000 deliberately.
