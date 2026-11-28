"""
benchmark.py - verify C++ == Python and time a 10,000-case sweep (100 alphas x 100 taper ratios).
Writes cpp_benchmark_results.csv. Timings depend on THE MACHINE; re-run on your own computer.
Note: numpy.linalg.solve calls optimized LAPACK, the C++ code is a hand-written Gaussian
elimination, and this Python loop pays interpreter overhead per case (matrix build + solve).
Physics note: e depends only on the planform, so a smart sweep needs only 100 solves, not 10,000.
This benchmark deliberately repeats the full solve per case to compare solver implementations.
"""
import csv, os, platform, subprocess, sys, tempfile, time
import numpy as np
from lifting_line import solve_lifting_line

N, REPS_PY = 40, 3
ALPHAS = np.linspace(0.1, 10.0, 100)
TAPERS = np.linspace(0.1, 1.0, 100)

BIN = os.path.join(tempfile.mkdtemp(), "lifting_line" + (".exe" if os.name == "nt" else ""))  # no binary kept in repo
subprocess.run(["g++", "-O3", "-std=c++17", "cpp/lifting_line.cpp", "-o", BIN], check=True)

# ---- accuracy check: rectangular wing, AR=8, alpha=5 deg, alpha_L0 = 0 ----
py = solve_lifting_line(8.0, 1.0, 5.0, 0.0, N)
out = subprocess.run([BIN, "verify"], capture_output=True, text=True, check=True).stdout.split()
cpp = {out[i]: float(out[i + 1]) for i in range(0, 6, 2)}
print("Verification (rectangular, AR=8, alpha=5 deg):")
worst = 0.0
for k, kp in (("CL", "CL"), ("e", "e"), ("CDi", "CDi")):
    d = abs(py[kp] - cpp[k]); worst = max(worst, d)
    print(f"  {k:3s} python {py[kp]:.12f}  c++ {cpp[k]:.12f}  |diff| = {d:.2e}")
print(f"  match to 6 decimals: {worst < 5e-7}")

# ---- timing ----
best_py, chk_py = 1e300, 0.0
for _ in range(REPS_PY):
    s = 0.0; t0 = time.perf_counter()
    for a in ALPHAS:
        for t in TAPERS:
            r = solve_lifting_line(8.0, t, a, 0.0, N)
            s += r["CL"] + r["e"]
    best_py = min(best_py, time.perf_counter() - t0); chk_py = s
b = subprocess.run([BIN, "bench", "5"], capture_output=True, text=True, check=True).stdout.split()
best_cpp, chk_cpp = float(b[1]), float(b[2])
agree = abs(chk_py - chk_cpp) / abs(chk_py) < 1e-9
print(f"\nPython best of {REPS_PY}: {best_py:.4f} s | C++ -O3 best of 5: {best_cpp:.4f} s | "
      f"speedup x{best_py / best_cpp:.1f} | checksums agree (rel 1e-9): {agree}")
cxx = subprocess.run(["g++", "--version"], capture_output=True, text=True).stdout.splitlines()[0]
with open("cpp_benchmark_results.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["cases", "N", "python_s", "cpp_O3_s", "speedup", "checksum_python", "checksum_cpp",
                "checksums_agree", "max_abs_diff_verify", "platform", "python", "numpy", "compiler"])
    w.writerow([10000, N, f"{best_py:.5f}", f"{best_cpp:.5f}", f"{best_py/best_cpp:.2f}", f"{chk_py:.12f}",
                f"{chk_cpp:.12f}", agree, f"{worst:.2e}", platform.platform(), sys.version.split()[0],
                np.__version__, cxx])
print("Wrote cpp_benchmark_results.csv")
