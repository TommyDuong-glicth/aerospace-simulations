"""
aerodynamics.py - lift, drag and L/D from thin airfoil theory + a parabolic drag polar.

MODEL LABELS
  [THEORY]      thin airfoil theory (2D, inviscid, incompressible, small angles):
                  C_l = 2*pi*(alpha - alpha_L0)
                  alpha_L0 = -(1/pi) * Int_0^pi (dy_c/dx)(cos(th) - 1) dth,  x = (1 - cos th)/2
  [THEORY]      Prandtl lifting-line slope correction for finite aspect ratio AR, efficiency e
  [ASSUMPTION]  drag polar C_D = C_D0 + C_L^2/(pi*e*AR); C_D0 is an INPUT (not computed)
  [NOT MODELED] viscosity, stall, Reynolds/Mach effects, thickness effect on C_L
"""
import numpy as np
from geometry import parse_naca4, camber_line

_trapz = getattr(np, "trapezoid", None) or np.trapz   # NumPy 1.x and 2.x


def zero_lift_angle(code: str, n: int = 4001) -> float:
    """alpha_L0 in radians (thin airfoil theory)."""
    m, p, _ = parse_naca4(code)
    if m == 0:
        return 0.0
    th = np.linspace(0, np.pi, n)
    x = 0.5 * (1 - np.cos(th))
    _, dyc = camber_line(x, m, p)
    return -_trapz(dyc * (np.cos(th) - 1), th) / np.pi


def lift_curve_slope(AR: float, e: float) -> float:
    """Finite-wing slope dC_L/dalpha [1/rad]; a0 = 2*pi is the 2D value."""
    a0 = 2*np.pi
    return a0 / (1 + a0 / (np.pi * e * AR))


def aerodynamic_polars(code, alpha_deg, AR=8.0, e=0.85, CD0=0.010):
    """Return dict with CL, CD, LD arrays and alpha_L0 [deg]."""
    alpha = np.radians(np.asarray(alpha_deg, dtype=float))
    aL0 = zero_lift_angle(code)
    CL = lift_curve_slope(AR, e) * (alpha - aL0)
    CD = CD0 + CL**2 / (np.pi * e * AR)
    return {"CL": CL, "CD": CD, "LD": CL / CD, "alpha_L0_deg": np.degrees(aL0)}


def best_glide(CD0: float, AR: float, e: float):
    """Closed form for the parabolic polar: CL* = sqrt(CD0*pi*e*AR), (L/D)max = 1/(2*sqrt(CD0*k))."""
    k = 1 / (np.pi * e * AR)
    return np.sqrt(CD0 / k), 1 / (2*np.sqrt(CD0*k))
