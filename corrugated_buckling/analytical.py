"""Closed-form buckling loads. All loads are line loads in N per metre of loaded edge."""
from __future__ import annotations

import numpy as np
from .homogenise import Sheet, Stiffness
from .profiles import Profile

# Euler effective-length factor K (L_e = K L) for the loaded ends
END_FACTORS = {"pinned": 1.0, "fixed": 0.5, "pinned-fixed": 0.699156, "fixed-free": 2.0}


def wide_strut(D: float, length: float, ends: str = "pinned") -> float:
    """Euler load of a wide plate strut bending cylindrically: pi^2 D / (K L)^2."""
    return np.pi ** 2 * D / (END_FACTORS[ends] * length) ** 2


def orthotropic_plate(st: Stiffness, length: float, width: float, direction: str = "x", max_m: int = 12):
    """Simply supported orthotropic plate on all four edges, uniaxial compression.

    direction 'x' loads across the corrugations, 'y' along them.
    Returns (N_cr, m), m being the number of half-waves in the load direction.
    """
    Dl, Dt = (st.D11, st.D22) if direction == "x" else (st.D22, st.D11)
    H = st.D12 + 2 * st.D66
    best = (np.inf, 1)
    for m in range(1, max_m + 1):
        N = np.pi ** 2 * (Dl * (m / length) ** 2 + 2 * H / width ** 2
                          + Dt * length ** 2 / (m ** 2 * width ** 4))
        if N < best[0]:
            best = (N, m)
    return best


def squash_load(sheet: Sheet, p: Profile, direction: str = "y") -> float:
    """Line load at which the section yields in pure compression.

    Along the corrugations the whole developed section carries the load.
    Across them the sheet is bent by the offset of the thrust line, so first
    yield is reached when membrane plus bending stress equals the yield stress.
    """
    if direction == "y":
        return sheet.yield_stress * sheet.t * p.stretch
    lever = p.depth / 2                       # largest offset from the thrust line
    return sheet.yield_stress / (1 / sheet.t + 6 * lever / sheet.t ** 2)


def design_resistance(N_cr: float, N_y: float, alpha: float = 0.49) -> float:
    """Indicative resistance with imperfections, EN 1993-1-1 column curve (default curve c)."""
    lam = np.sqrt(N_y / N_cr)
    phi = 0.5 * (1 + alpha * (lam - 0.2) + lam ** 2)
    chi = min(1.0, 1 / (phi + np.sqrt(max(phi ** 2 - lam ** 2, 0.0))))
    return chi * N_y


def mass_per_area(sheet: Sheet, p: Profile) -> float:
    """Panel mass per unit plan area [kg/m^2]."""
    return sheet.density * sheet.t * p.stretch
