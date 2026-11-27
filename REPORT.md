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

![NACA 2412 vs NACA 0012 Comparison](naca_comparison.png)

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

## 8. Next steps

- Validate against published wind-tunnel experimental measurements (e.g., Abbott & von Doenhoff) to separate XFOIL panel-method transition modeling from physical boundary-layer behavior.
- Port the Prandtl lifting-line and polar sweep routines to C++ (`solver.cpp`) and benchmark execution speed against the Python implementation.

## 9. Viscous Benchmark Comparison (XFOIL)

**What is compared.** The v1 thin-airfoil model against XFOIL polars from Airfoil Tools (`naca2412-il` and `naca0012-il`) at $\text{Re} = 1{,}000{,}000$, $M = 0$, and $N_{\text{crit}} = 9$. Because XFOIL is itself a viscous panel/boundary-layer model rather than wind-tunnel data, this serves as a model-to-model benchmark.

**2D to 3D mapping.** Each XFOIL point $(\alpha_{2D}, c_l, c_d)$ is mapped to a finite wing of aspect ratio $AR = 8.0$ with span efficiency $e = 0.85$ using $\alpha_{3D} = \alpha_{2D} + c_l/(\pi e AR)$, $C_L = c_l$, and $C_D = c_d + c_l^2/(\pi e AR)$. This assumes elliptic spanwise loading and neglects 3D crossflow, sweep, and taper effects.

**Results** ($\text{Re} = 1\text{e}6$, $AR = 8.0$, $e = 0.85$; attached-flow window $-2^\circ$ to $8^\circ$; v1 assumed $C_{D0} = 0.010$).

**NACA 2412**

| Metric | v1 inviscid theory | XFOIL viscous | Error / divergence |
|---|---|---|---|
| 2D lift slope [1/rad] | 6.283 | 6.224 | +0.9 % |
| α_L0 [deg] | −2.08 | −2.30 | +0.22° |
| RMSE of c_l, linear window | n/a | n/a | 0.0256 (4.4 % of mean \|c_l\|) |
| c_l at α_stall = 16.5° | 2.037 | 1.582 (max) | 0.455 |
| (L/D)max, 3D | 23.11 | 31.24 | −8.13 |
| α at (L/D)max [deg] | 3.38 | 1.93 | +1.45° |

**NACA 0012**

| Metric | v1 inviscid theory | XFOIL viscous | Error / divergence |
|---|---|---|---|
| 2D lift slope [1/rad] | 6.283 | 6.566 | −4.3 % |
| α_L0 [deg] | 0.00 | 0.04 | −0.04° |
| RMSE of c_l, linear window | n/a | n/a | 0.0237 (6.2 % of mean \|c_l\|) |
| c_l at α_stall = 15.5° | 1.700 | 1.389 (max) | 0.311 |
| (L/D)max, 3D | 23.11 | 28.61 | −5.50 |
| α at (L/D)max [deg] | 5.45 | 3.54 | +1.91° |

![Lift Curve Validation](figures/cl_alpha_validation.png)

![L/D and Drag Polar Comparison](figures/ld_polar_comparison.png)

**Interpretation.** Thin-airfoil theory accurately captures attached-flow lift behavior between $-2^\circ$ and $8^\circ$, keeping $c_l$ RMSE to $0.0256$ (4.4%) for NACA 2412 and $0.0237$ (6.2%) for NACA 0012, with 2D lift slopes within 0.9% and 4.3% of $2\pi\text{ rad}^{-1}$. At high angles of attack, adverse pressure gradients thicken and separate the suction-side boundary layer—capping $c_{l,\max}$ at $1.582$ ($\alpha_{\text{stall}} = 16.5^\circ$) for NACA 2412 and $1.389$ ($\alpha_{\text{stall}} = 15.5^\circ$) for NACA 0012, where the unseparated inviscid model overpredicts lift by $0.455$ and $0.311$, respectively. Finally, replacing v1's identical placeholder $C_{D0} = 0.010$ with viscous drag polars resolves the v1 efficiency tie: NACA 2412 achieves a 3D $(L/D)_{\max}$ of $31.24$ at $\alpha_{3D} = 1.93^\circ$ compared to $28.61$ at $\alpha_{3D} = 3.54^\circ$ for NACA 0012 (a $+2.64$ advantage, well above data-step resolution), while confirming that actual profile drag at $\text{Re} = 10^6$ is substantially lower than the conservative $0.010$ assumption. Note that XFOIL predictions near stall remain sensitive to boundary-layer transition modeling ($N_{\text{crit}} = 9$) and convergence limits.