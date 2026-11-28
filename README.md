# aerospace-simulations

Python simulations in aerodynamics, flight dynamics, and orbital mechanics.

## Projects

| Folder | Status | Description |
| :--- | :---: | :--- |
| `naca-airfoil-analysis/` | Done (v3) | NACA 2412 vs NACA 0012: geometry, thin-airfoil theory (v1); XFOIL viscous polar comparison at $`Re = 10^6`$ (v2); numerical Prandtl lifting-line solver in Python and C++ with a runtime benchmark, and comparison with published wind-tunnel data (v3) |

## Run

```bash
pip install numpy matplotlib
cd naca-airfoil-analysis

python main.py                          # v1: thin-airfoil theory, parabolic drag polar
python validate_viscous.py --re 1e6     # v2: needs data/xfoil_*.csv
python lifting_line.py                  # v3: span efficiency table and spanwise loading
python benchmark.py                     # v3: Python vs C++ (requires g++)
python validate_experimental.py --re 1e6 --re-tunnel <Re>   # v3: needs data/windtunnel_*.csv
```

## Repository layout

```
naca-airfoil-analysis/
├── data/             XFOIL and wind-tunnel CSV files (see data/README.md)
├── cpp/              C++ lifting-line solver
├── figures/          generated plots
├── geometry.py       NACA 4-digit geometry
├── aerodynamics.py   thin-airfoil lift, drag polar
├── lifting_line.py   numerical lifting-line solver
└── REPORT.md         methods, assumptions, results, limitations
```

## Data and scope

- Thin-airfoil and lifting-line results are model estimates under the assumptions stated in `REPORT.md`.
- XFOIL polars are simulation results (Airfoil Tools, Ncrit as stated in the report), not measurements.
- Wind-tunnel data are from the sources cited in `REPORT.md`; check that the Reynolds number matches the XFOIL data before comparing.
- Benchmark timings depend on the machine and compiler; see `cpp_benchmark_results.csv`.
