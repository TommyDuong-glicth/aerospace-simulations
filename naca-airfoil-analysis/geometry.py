"""
geometry.py - NACA 4-digit airfoil geometry (chord = 1).

Code "MPTT": m = max camber (MP/100 of chord), p = its location (tenths of chord),
TT = max thickness (percent of chord).  NACA 2412 -> m=0.02, p=0.4, t=0.12.
"""
import numpy as np


def parse_naca4(code: str):
    if len(code) != 4 or not code.isdigit():
        raise ValueError(f"Not a NACA 4-digit code: {code!r}")
    return int(code[0]) / 100, int(code[1]) / 10, int(code[2:]) / 100


def thickness(x, t):
    """Half-thickness y_t(x); the -0.1036 coefficient gives a closed trailing edge."""
    return 5 * t * (0.2969*np.sqrt(x) - 0.1260*x - 0.3516*x**2
                    + 0.2843*x**3 - 0.1036*x**4)


def camber_line(x, m, p):
    """Mean camber line y_c(x) and slope dy_c/dx (two parabolic arcs joined at x = p)."""
    x = np.asarray(x, dtype=float)
    yc, dyc = np.zeros_like(x), np.zeros_like(x)
    if m > 0:
        f = x < p
        yc[f] = m / p**2 * (2*p*x[f] - x[f]**2)
        dyc[f] = 2*m / p**2 * (p - x[f])
        yc[~f] = m / (1-p)**2 * ((1 - 2*p) + 2*p*x[~f] - x[~f]**2)
        dyc[~f] = 2*m / (1-p)**2 * (p - x[~f])
    return yc, dyc


def generate_naca4(code: str, n: int = 200):
    """Return (xu, yu, xl, yl). Cosine spacing clusters points at leading/trailing edge."""
    m, p, t = parse_naca4(code)
    beta = np.linspace(0, np.pi, n)
    x = 0.5 * (1 - np.cos(beta))
    yt = thickness(x, t)
    yc, dyc = camber_line(x, m, p)
    th = np.arctan(dyc)
    xu, yu = x - yt*np.sin(th), yc + yt*np.cos(th)
    xl, yl = x + yt*np.sin(th), yc - yt*np.cos(th)
    return xu, yu, xl, yl


def check_geometry(code: str):
    """Sanity checks against the definition of the code. Returns a dict of measured values."""
    m, p, t = parse_naca4(code)
    x = np.linspace(0, 1, 20001)
    yc, _ = camber_line(x, m, p)
    yt = thickness(x, t)
    return {"max_camber": yc.max(), "x_max_camber": x[yc.argmax()],
            "max_thickness": 2*yt.max(), "x_max_thickness": x[yt.argmax()],
            "te_gap": 2*thickness(np.array([1.0]), t)[0]}
