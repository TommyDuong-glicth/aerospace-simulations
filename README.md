# aerospace-simulations

Python simulations in aerodynamics, flight dynamics, and orbital mechanics.

## Projects

| Folder | Status | Description |
|---|---|---|
| `naca-airfoil-analysis/` | Done (v2) | NACA 2412 vs NACA 0012: geometry, thin-airfoil theory, finite-wing lifting line, and XFOIL viscous polar benchmarking ($Re = 10^6$) |

## Run

```bash
pip install numpy matplotlib
cd naca-airfoil-analysis
python main.py
python validate_viscous.py --re 1e6
