"""
validate_experimental.py - thin-airfoil theory (v1) vs XFOIL (v2) vs wind-tunnel points (v3).

Inputs (you supply; nothing here is typed from memory):
  data/xfoil_naca<code>_Re<label>.csv      XFOIL polar (see data/README.md)
  data/windtunnel_naca<code>.csv           columns alpha,cl,cd (cd optional); '#' comment lines allowed
IMPORTANT: wind-tunnel and XFOIL data must be at the SAME Reynolds number (and comparable
Ncrit / turbulence level / surface roughness) or the RMSE mixes Reynolds effects with model error.
Values read off published plots carry digitization uncertainty; record it in the CSV header.
"""
import argparse, csv, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from aerodynamics import zero_lift_angle
from validate_viscous import load_polar, WINDOW

COLORS = {"2412": "tab:blue", "0012": "tab:red"}


def load_wt(path):
    rows = []
    with open(path, encoding="utf-8-sig") as f:
        lines = [l for l in f if l.strip() and not l.lstrip().startswith("#")]
    for r in csv.DictReader(lines):
        r = {k.strip().lower(): v for k, v in r.items()}
        rows.append((float(r["alpha"]), float(r["cl"]), float(r["cd"]) if r.get("cd") not in (None, "") else np.nan))
    d = np.array(sorted(rows))
    return d[:, 0], d[:, 1], d[:, 2]


def rmse(x, y):
    return float(np.sqrt(np.mean((x - y) ** 2))) if len(x) else float("nan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--re", default="1e6", help="Re label of the XFOIL files")
    ap.add_argument("--re-tunnel", default="UNSPECIFIED", help="Re of the wind-tunnel data (for the plot label)")
    ap.add_argument("--near-stall-deg", type=float, default=4.0,
                    help="'near cl_max' = tunnel points with alpha >= alpha(cl_max,exp) - this value")
    ap.add_argument("--data-dir", default="data"); ap.add_argument("--out-dir", default=".")
    A = ap.parse_args()
    os.makedirs(os.path.join(A.out_dir, "figures"), exist_ok=True)

    if A.re_tunnel != A.re:
        print(f"WARNING: tunnel Re ({A.re_tunnel}) != XFOIL Re ({A.re}). Differences will include Reynolds-number effects.")
    fig, ax = plt.subplots(2, 2, figsize=(11, 8.5))
    L = [f"XFOIL Re = {A.re}; wind-tunnel Re = {A.re_tunnel}; linear window {WINDOW[0]:.0f}..{WINDOW[1]:.0f} deg; "
         f"near-stall = alpha >= alpha(cl_max,exp) - {A.near_stall_deg:g} deg.\n",
         "| Airfoil | Region | N pts | RMSE(XFOIL vs tunnel) | RMSE(theory vs tunnel) | cl_max XFOIL | cl_max tunnel |",
         "|---|---|---|---|---|---|---|"]
    for k, c in enumerate(("2412", "0012")):
        xa, xcl, xcd, _ = load_polar(os.path.join(A.data_dir, f"xfoil_naca{c}_Re{A.re}.csv"), A.re)
        ta, tcl, tcd = load_wt(os.path.join(A.data_dir, f"windtunnel_naca{c}.csv"))
        th = lambda a: 2 * np.pi * (np.radians(a) - zero_lift_angle(c))
        inr = (ta >= xa.min()) & (ta <= xa.max())
        if not inr.all():
            print(f"NACA {c}: {(~inr).sum()} tunnel points lie outside the XFOIL alpha range and are skipped.")
        i_exp = int(np.argmax(tcl))
        regions = {"linear": inr & (ta >= WINDOW[0]) & (ta <= WINDOW[1]),
                   "near cl_max": inr & (ta >= ta[i_exp] - A.near_stall_deg)}
        xi = np.interp(ta, xa, xcl)
        for name, m in regions.items():
            L.append(f"| NACA {c} | {name} | {int(m.sum())} | {rmse(xi[m], tcl[m]):.4f} | {rmse(th(ta[m]), tcl[m]):.4f} "
                     f"| {xcl.max():.3f} | {tcl.max():.3f} |")
            print(f"NACA {c} [{name}, N={int(m.sum())}]: RMSE XFOIL-tunnel = {rmse(xi[m], tcl[m]):.4f}, "
                  f"theory-tunnel = {rmse(th(ta[m]), tcl[m]):.4f}; cl_max XFOIL {xcl.max():.3f} vs tunnel {tcl.max():.3f}")
        a_ = ax[0, k]; col = COLORS[c]
        grid = np.linspace(xa.min(), xa.max(), 400)
        a_.plot(grid, th(grid), "--", color="k", label="thin-airfoil theory (v1)")
        a_.plot(xa, xcl, "-", color=col, label=f"XFOIL (v2), Re={A.re}")
        a_.plot(ta, tcl, "o", color=col, mfc="none", label=f"wind tunnel (v3), Re={A.re_tunnel}")
        a_.set(title=f"NACA {c}", xlabel="α [deg]", ylabel="$c_l$"); a_.grid(alpha=.3); a_.legend(fontsize=8)
        b_ = ax[1, k]
        b_.plot(xcd, xcl, "-", color=col, label="XFOIL")
        if np.isfinite(tcd).any():
            b_.plot(tcd, tcl, "o", color=col, mfc="none", label="wind tunnel")
        b_.set(title=f"NACA {c} drag polar", xlabel="$c_d$", ylabel="$c_l$"); b_.grid(alpha=.3); b_.legend(fontsize=8)
    fig.suptitle("Theory vs XFOIL vs wind-tunnel data (2D section)", fontsize=11)
    fig.tight_layout(); fig.savefig(os.path.join(A.out_dir, "figures", "experimental_validation.png"), dpi=300)
    open(os.path.join(A.out_dir, "section10_experimental_table.md"), "w", encoding="utf-8").write("\n".join(L) + "\n")
    print("Wrote figures/experimental_validation.png and section10_experimental_table.md")


if __name__ == "__main__":
    main()
