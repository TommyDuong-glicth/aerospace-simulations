"""main.py - compare NACA 2412 and NACA 0012; writes naca_polars.csv and naca_comparison.png."""
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from geometry import generate_naca4, check_geometry
from aerodynamics import aerodynamic_polars, best_glide, lift_curve_slope

# ---- INPUTS (assumptions, not measurements) -------------------------------------------
AIRFOILS = ["2412", "0012"]
AR, E = 8.0, 0.85
CD0 = {"2412": 0.010, "0012": 0.010}   # SAME on purpose: no airfoil-specific drag data yet
ALPHA = np.linspace(-4, 12, 161)       # stay in the linear (pre-stall) range
# ---------------------------------------------------------------------------------------

res = {c: aerodynamic_polars(c, ALPHA, AR, E, CD0[c]) for c in AIRFOILS}

print("Geometry checks:")
for c in AIRFOILS:
    g = check_geometry(c)
    print(f"  NACA {c}: camber {g['max_camber']:.4f} @ x={g['x_max_camber']:.3f}, "
          f"thickness {g['max_thickness']:.4f} @ x={g['x_max_thickness']:.3f}, TE gap {g['te_gap']:.1e}")

print(f"Finite-wing slope: {lift_curve_slope(AR, E):.3f} /rad (2D: {2*np.pi:.3f})")
print("Results (model estimates):")
for c in AIRFOILS:
    r = res[c]; i = r["LD"].argmax()
    CLs, LDs = best_glide(CD0[c], AR, E)
    print(f"  NACA {c}: alpha_L0={r['alpha_L0_deg']:.3f} deg | CL(0)={np.interp(0, ALPHA, r['CL']):.3f} | "
          f"(L/D)max={r['LD'][i]:.2f} at {ALPHA[i]:.2f} deg | closed form: {LDs:.2f} at CL={CLs:.3f}")

with open("naca_polars.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["airfoil", "alpha_deg", "CL", "CD", "L_over_D"])
    for c in AIRFOILS:
        for a, cl, cd, ld in zip(ALPHA, res[c]["CL"], res[c]["CD"], res[c]["LD"]):
            w.writerow([f"NACA {c}", f"{a:.2f}", f"{cl:.5f}", f"{cd:.5f}", f"{ld:.4f}"])

fig, ax = plt.subplots(2, 2, figsize=(11, 8))
for c in AIRFOILS:
    xu, yu, xl, yl = generate_naca4(c)
    ax[0, 0].plot(np.r_[xu, xl[::-1]], np.r_[yu, yl[::-1]], label=f"NACA {c}")
    r = res[c]; i = r["LD"].argmax()
    ax[0, 1].plot(ALPHA, r["CL"], label=f"NACA {c} (α_L0 = {r['alpha_L0_deg']:.2f}°)")
    ax[1, 0].plot(r["CD"], r["CL"], label=f"NACA {c}")
    ax[1, 1].plot(ALPHA, r["LD"], label=f"NACA {c}: max {r['LD'][i]:.1f} at {ALPHA[i]:.1f}°")
ax[0, 0].set_aspect("equal"); ax[0, 0].set(title="Geometry (chord = 1)", xlabel="x/c", ylabel="y/c")
ax[0, 1].set(title="Lift curve", xlabel="α [deg]", ylabel="$C_L$")
ax[1, 0].set(title="Drag polar", xlabel="$C_D$", ylabel="$C_L$")
ax[1, 1].set(title="Efficiency", xlabel="α [deg]", ylabel="$C_L/C_D$")
for a in ax.flat:
    a.grid(alpha=.3); a.legend(fontsize=8)
fig.suptitle(f"Thin-airfoil model: AR={AR}, e={E}, $C_{{D0}}$=0.010 (assumed, equal), no stall")
fig.tight_layout(); fig.savefig("naca_comparison.png", dpi=150)
print("Wrote naca_polars.csv and naca_comparison.png")
