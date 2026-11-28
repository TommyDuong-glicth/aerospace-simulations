"""
lifting_line.py - numerical Prandtl lifting-line (Glauert Fourier method).

[THEORY]  Prandtl lifting-line, incompressible, small angles, 2D section slope a0 = 2*pi,
          straight (unswept) wing, no twist, symmetric loading -> odd Fourier terms only.
Monoplane equation at station theta_i = i*pi/(2N), i = 1..N (b = 1, so c means c/b):
    sum_{n=1,3,..,2N-1} A_n sin(n th) [ 4/(a0 c(th)/b) + n/sin(th) ] = alpha - alpha_L0
Outputs:  C_L = pi*AR*A1,  delta = sum_{n>=3} n (A_n/A_1)^2,  e = 1/(1+delta),  C_Di = C_L^2/(pi e AR)
Planform: c(y) = c_root*(1 - (1-lambda)*|2y/b|), c_root/b = 2/(AR(1+lambda)); elliptic: c/b = 4 sin(th)/(pi AR).
e depends only on planform (not on alpha) for an untwisted wing with constant alpha_L0.
"""
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

A0 = 2 * np.pi


def chord_over_b(theta, AR, taper):
    if taper == "elliptic":
        return 4 / (np.pi * AR) * np.sin(theta)
    return 2 / (AR * (1 + taper)) * (1 - (1 - taper) * np.cos(theta))


def solve_lifting_line(AR, taper, alpha_deg, alpha_L0_deg=0.0, N=40, a0=A0):
    i = np.arange(1, N + 1)
    th = i * np.pi / (2 * N)
    n = 2 * i - 1
    c = chord_over_b(th, AR, taper)
    M = np.sin(np.outer(th, n)) * (4 / (a0 * c)[:, None] + n[None, :] / np.sin(th)[:, None])
    A = np.linalg.solve(M, np.full(N, np.radians(alpha_deg - alpha_L0_deg)))
    CL = np.pi * AR * A[0]
    if abs(A[0]) < 1e-14:                       # alpha = alpha_L0: e undefined from this ratio
        return dict(CL=CL, e=np.nan, CDi=np.nan, delta=np.nan, A=A, n=n)
    delta = np.sum(n[1:] * (A[1:] / A[0]) ** 2)
    e = 1 / (1 + delta)
    return dict(CL=CL, e=e, CDi=CL**2 / (np.pi * e * AR), delta=delta, A=A, n=n)


def gamma_distribution(A, n, npts=401):
    """Gamma(y)/Gamma_root over the full span. Returns (2y/b, normalized Gamma)."""
    th = np.linspace(1e-6, np.pi - 1e-6, npts)
    g = np.sin(np.outer(th, n)) @ A
    return -np.cos(th), g / (np.sin(np.outer([np.pi / 2], n)) @ A)[0]


if __name__ == "__main__":
    AR = 8.0
    # ---- check 1: elliptic planform must give e = 1 ----
    el = solve_lifting_line(AR, "elliptic", 5.0)
    print(f"Check: elliptic AR={AR}: e = {el['e']:.12f} (expect 1), "
          f"CL = {el['CL']:.6f}, analytic a*alpha = {A0/(1+A0/(np.pi*AR))*np.radians(5):.6f}")
    # ---- check 2: convergence in N ----
    print("Convergence of e with N:")
    for N in (10, 20, 40, 80, 160):
        print(f"  N={N:4d}  rect e={solve_lifting_line(AR,1.0,5,N=N)['e']:.8f}  "
              f"taper0.4 e={solve_lifting_line(AR,0.4,5,N=N)['e']:.8f}")
    # ---- e table ----
    rows = [("elliptic", "elliptic")] + [(f"taper {l}", l) for l in (1.0, 0.8, 0.6, 0.5, 0.4, 0.2)]
    with open("lifting_line_efficiency.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["planform", "AR", "N", "e", "CL_at_5deg_alphaL0_0", "CDi_at_5deg"])
        print(f"\nSpan efficiency, AR={AR}, N=40:")
        for name, t in rows:
            r = solve_lifting_line(AR, t, 5.0)
            w.writerow([name, AR, 40, f"{r['e']:.6f}", f"{r['CL']:.6f}", f"{r['CDi']:.6f}"])
            print(f"  {name:12s} e = {r['e']:.4f}   CL = {r['CL']:.4f}   CDi = {r['CDi']:.5f}")
    print("\nCompare with the v1/v2 assumed e = 0.85 (see REPORT for the caveat on non-elliptic wings).")
    # ---- figure ----
    fig, ax = plt.subplots(figsize=(7.5, 5))
    for name, t, st in (("elliptic", "elliptic", "-"), ("taper λ=0.4", 0.4, "--"), ("rectangular λ=1.0", 1.0, ":")):
        r = solve_lifting_line(AR, t, 5.0)
        y, g = gamma_distribution(r["A"], r["n"])
        ax.plot(y, g, st, lw=2, label=f"{name}  (e = {r['e']:.4f})")
    ax.set(xlabel="2y/b", ylabel="Γ(y) / Γ(root)", title=f"Spanwise circulation, AR = {AR:g} (lifting-line model)")
    ax.grid(alpha=.3); ax.legend()
    fig.tight_layout(); fig.savefig("figures/spanwise_loading.png", dpi=300)
    print("Wrote figures/spanwise_loading.png and lifting_line_efficiency.csv")
