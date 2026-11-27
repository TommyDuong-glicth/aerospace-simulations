"""
validate_viscous.py - compare the v1 inviscid model with XFOIL polars (from Airfoil Tools).

WHAT THIS IS: a comparison between two MODELS (thin-airfoil theory vs XFOIL, a viscous panel
code with an e^N transition model). XFOIL output is simulation, not experimental measurement.

2D -> 3D MAPPING (elliptic-loading lifting line, consistent with aerodynamics.py)
  A section on a finite wing flies at an effective angle alpha_2D = alpha - alpha_i,
  alpha_i = C_L/(pi*e*AR). So each XFOIL point (alpha_2D, cl, cd) maps to
      alpha_3D = alpha_2D + cl/(pi*e*AR),  C_L = cl,  C_D = cd + cl^2/(pi*e*AR).
  (Evaluating cd at alpha_3D instead would use the wrong section angle.)

Usage:  python validate_viscous.py --re 1e6
Files:  data/xfoil_naca0012_Re1e6.csv, data/xfoil_naca2412_Re1e6.csv  (columns Alpha,Cl,Cd[,...])
"""
import argparse, os, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from aerodynamics import zero_lift_angle, lift_curve_slope, aerodynamic_polars

AIRFOILS = ["2412", "0012"]
COLORS = {"2412": "tab:blue", "0012": "tab:red"}
WINDOW = (-2.0, 8.0)          # attached-flow window [deg], from the plan


def load_polar(path, re_label):
    """Read an Airfoil Tools CSV (metadata lines allowed) -> sorted arrays alpha, cl, cd."""
    pre, header, rows = [], None, []
    with open(path, encoding="utf-8-sig") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if header is None:
                if s.lower().replace(" ", "").startswith("alpha,"):
                    header = [h.strip().lower() for h in s.split(",")]
                else:
                    pre.append(s)
                continue
            try:
                rows.append([float(v) for v in s.split(",")[:len(header)]])
            except ValueError:
                continue
    if header is None or not {"alpha", "cl", "cd"} <= set(header):
        raise ValueError(f"{path}: need a header with Alpha, Cl, Cd")
    d = np.array(rows)
    a, cl, cd = (d[:, header.index(k)] for k in ("alpha", "cl", "cd"))
    o = np.argsort(a)
    meta = " ".join(pre)
    m = re.search(r"Reynolds number\s*=\s*([0-9.eE+]+)", meta)
    if m and abs(float(m.group(1)) - float(re_label)) > 1e-6 * float(re_label):
        print(f"  WARNING: file says Re = {m.group(1)} but --re = {re_label}. Fix the label or the file.")
    return a[o], cl[o], cd[o], meta[:160]


def rmse(x, y):
    return float(np.sqrt(np.mean((x - y) ** 2)))


def analyse(code, a, cl, cd, AR, e, CD0_v1):
    aL0 = zero_lift_angle(code)
    th = 2 * np.pi * (np.radians(a) - aL0)                       # thin-airfoil 2D cl
    w = (a >= WINDOW[0]) & (a <= WINDOW[1])
    out = {"a": a, "cl": cl, "cd": cd, "cl_theory": th, "w": w}
    out["rmse"] = rmse(th[w], cl[w])
    out["nrmse_pct"] = 100 * out["rmse"] / np.mean(np.abs(cl[w]))
    slope, icpt = np.polyfit(np.radians(a[w]), cl[w], 1)
    out["slope_fit"], out["aL0_fit_deg"] = slope, np.degrees(-icpt / slope)
    out["aL0_th_deg"] = np.degrees(aL0)
    i = int(np.argmax(cl))
    out.update(i_stall=i, a_stall=a[i], clmax=cl[i], cl_th_at_stall=th[i],
               gap=abs(th[i] - cl[i]), stall_at_edge=(i in (0, len(a) - 1)))
    k = 1 / (np.pi * e * AR)
    a3 = a + np.degrees(cl * k)
    CD3 = cd + cl**2 * k
    LD = cl / CD3
    pos = cl > 0
    j = int(np.argmax(np.where(pos, LD, -np.inf)))
    out.update(a3=a3, CD3=CD3, LD=LD, LDmax=LD[j], a_opt=a3[j], j=j)
    v1 = aerodynamic_polars(code, np.linspace(-4, 18, 2201), AR, e, CD0_v1)
    out["v1"] = v1
    jj = int(np.argmax(v1["LD"]))
    out["v1_LDmax"], out["v1_aopt"] = v1["LD"][jj], np.linspace(-4, 18, 2201)[jj]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--re", default="1e6", help="Reynolds number label used in file names, e.g. 1e6")
    ap.add_argument("--ar", type=float, default=8.0)
    ap.add_argument("--e", type=float, default=0.85)
    ap.add_argument("--cd0-v1", type=float, default=0.010, help="assumed C_D0 used by v1 (placeholder)")
    ap.add_argument("--data-dir", default="data")
    ap.add_argument("--out-dir", default=".")
    A = ap.parse_args()
    figdir = os.path.join(A.out_dir, "figures"); os.makedirs(figdir, exist_ok=True)

    R = {}
    for c in AIRFOILS:
        path = os.path.join(A.data_dir, f"xfoil_naca{c}_Re{A.re}.csv")
        if not os.path.exists(path):
            raise SystemExit(f"Missing {path}. Download the polar from Airfoil Tools (see data/README.md).")
        a, cl, cd, meta = load_polar(path, A.re)
        print(f"NACA {c}: {len(a)} points, alpha {a.min():.2f}..{a.max():.2f} deg | {meta}")
        R[c] = analyse(c, a, cl, cd, A.ar, A.e, A.cd0_v1)

    # ---------------- terminal report ----------------
    print(f"\nWindow {WINDOW} deg, AR={A.ar}, e={A.e}, Re label={A.re}, v1 assumed CD0={A.cd0_v1}")
    for c, r in R.items():
        print(f"\nNACA {c}")
        print(f"  1) RMSE(cl, linear window) = {r['rmse']:.4f}  ({r['nrmse_pct']:.1f}% of mean |cl|)")
        print(f"     fitted slope = {r['slope_fit']:.3f}/rad (theory {2*np.pi:.3f}); "
              f"alpha_L0 fit = {r['aL0_fit_deg']:.2f} deg (theory {r['aL0_th_deg']:.2f})")
        edge = "  [NO STALL INSIDE DATA RANGE: max is at the edge]" if r["stall_at_edge"] else ""
        print(f"  2) alpha_stall = {r['a_stall']:.2f} deg, cl_max = {r['clmax']:.3f}, "
              f"theory there = {r['cl_th_at_stall']:.3f}, gap = {r['gap']:.3f}{edge}")
        print(f"  3) (L/D)max XFOIL-derived 3D = {r['LDmax']:.2f} at alpha_3D = {r['a_opt']:.2f} deg | "
              f"v1 (assumed CD0) = {r['v1_LDmax']:.2f} at {r['v1_aopt']:.2f} deg")
    d = R["2412"]["LDmax"] - R["0012"]["LDmax"]
    verdict = ("no meaningful difference (|diff| < 0.1)" if abs(d) < 0.1
               else f"higher in this dataset: NACA {'2412' if d > 0 else '0012'}")
    print(f"\n4) (L/D)max difference 2412 - 0012 = {d:+.2f} -> {verdict}")
    print("   (data step size limits resolution; interpret small differences with care)")

    # ---------------- figure 1: lift curves ----------------
    fig, ax = plt.subplots(figsize=(8, 5.5))
    ax.axvspan(*WINDOW, color="gray", alpha=0.12, label=f"attached-flow window {WINDOW[0]:.0f}..{WINDOW[1]:.0f}°")
    for n, (c, r) in enumerate(R.items()):
        ax.plot(r["a"], r["cl"], "-", color=COLORS[c], label=f"NACA {c} XFOIL (Re={A.re})")
        ax.plot(r["a"], r["cl_theory"], "--", color=COLORS[c], label=f"NACA {c} thin-airfoil theory")
        if not r["stall_at_edge"]:
            offset = (10, 18) if c == "2412" else (12, -28)
            ax.annotate(
                f"max $c_l$={r['clmax']:.2f}\nat {r['a_stall']:.1f}°",
                xy=(r["a_stall"], r["clmax"]),
                xytext=offset,
                textcoords="offset points",
                arrowprops=dict(arrowstyle="->", color=COLORS[c]),
                color=COLORS[c],
                fontsize=8,
            )
    ax.set(xlabel="α [deg]", ylabel="$c_l$ (2D)", title="Lift curve: inviscid theory vs XFOIL (both are models)")
    ax.grid(alpha=.3); ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(figdir, "cl_alpha_validation.png"), dpi=300); plt.close(fig)

    # ---------------- figure 2: polar + L/D ----------------
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    for c, r in R.items():
        col = COLORS[c]
        ax[0].plot(r["CD3"], r["cl"], "-", color=col, label=f"NACA {c} XFOIL-derived")
        ax[0].plot(r["v1"]["CD"], r["v1"]["CL"], "--", color=col, label=f"NACA {c} v1 parabola")
        ax[1].plot(r["a3"], r["LD"], "-", color=col, label=f"NACA {c} XFOIL-derived")
        ax[1].plot(np.linspace(-4, 18, 2201), r["v1"]["LD"], "--", color=col, label=f"NACA {c} v1")
        ax[1].plot(r["a_opt"], r["LDmax"], "o", color=col)
        ax[1].annotate(f"{r['LDmax']:.1f} @ {r['a_opt']:.1f}°", xy=(r["a_opt"], r["LDmax"]),
                       xytext=(5, 6), textcoords="offset points", color=col, fontsize=8)
    ax[0].set(xlabel="$C_D$ (3D, AR=%g)" % A.ar, ylabel="$C_L$", title="Drag polar")
    ax[1].set(xlabel="α₃D [deg]", ylabel="$C_L/C_D$", title="Efficiency", xlim=(-4, 18))
    for a_ in ax:
        a_.grid(alpha=.3); a_.legend(fontsize=8)
    fig.suptitle(f"Finite wing AR={A.ar}, e={A.e}; XFOIL sections at Re={A.re}; v1 uses assumed $C_{{D0}}$={A.cd0_v1}", fontsize=10)
    fig.tight_layout(); fig.savefig(os.path.join(figdir, "ld_polar_comparison.png"), dpi=300); plt.close(fig)

    # ---------------- Section 8 table (markdown) ----------------
    L = [f"Re = {A.re}, AR = {A.ar}, e = {A.e}; window {WINDOW[0]:.0f}° to {WINDOW[1]:.0f}°; v1 C_D0 = {A.cd0_v1} (assumed).\n"]
    for c, r in R.items():
        L += [f"**NACA {c}**\n", "| Metric | v1 inviscid theory | XFOIL viscous | Error / divergence |", "|---|---|---|---|",
              f"| 2D lift slope [1/rad] | {2*np.pi:.3f} | {r['slope_fit']:.3f} | {100*(2*np.pi-r['slope_fit'])/r['slope_fit']:+.1f} % |",
              f"| α_L0 [deg] | {r['aL0_th_deg']:.2f} | {r['aL0_fit_deg']:.2f} | {r['aL0_th_deg']-r['aL0_fit_deg']:+.2f}° |",
              f"| RMSE of c_l, linear window | n/a | n/a | {r['rmse']:.4f} ({r['nrmse_pct']:.1f} % of mean abs c_l) |",
              f"| c_l at α_stall = {r['a_stall']:.1f}° | {r['cl_th_at_stall']:.3f} | {r['clmax']:.3f} (max) | {r['gap']:.3f} |",
              f"| (L/D)max, 3D | {r['v1_LDmax']:.2f} | {r['LDmax']:.2f} | {r['v1_LDmax']-r['LDmax']:+.2f} |",
              f"| α at (L/D)max [deg] | {r['v1_aopt']:.2f} | {r['a_opt']:.2f} | {r['v1_aopt']-r['a_opt']:+.2f}° |", ""]
    with open(os.path.join(A.out_dir, "section8_table.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(L))
    print("\nWrote figures/cl_alpha_validation.png, figures/ld_polar_comparison.png, section8_table.md")


if __name__ == "__main__":
    main()
