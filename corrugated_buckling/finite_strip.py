"""Semi-analytical finite strip buckling of the panel compressed ALONG the corrugations.

The cross-section is divided into flat strips. Along the corrugations the
displacements vary as a single sine half-wave of length `half_wave`, which is
exact for simply supported loaded ends. Each strip carries plane-stress
membrane action and Kirchhoff plate bending, so local buckling of the flats
and webs, and global buckling of the whole panel, come out of the same model.
"""
from __future__ import annotations

import numpy as np
from scipy.linalg import eigh
from .homogenise import Sheet
from .profiles import Profile

_GP = np.array([-0.8611363116, -0.3399810436, 0.3399810436, 0.8611363116])
_GW = np.array([0.3478548451, 0.6521451549, 0.6521451549, 0.3478548451])


def _strip(sheet: Sheet, b: float, k: float, half_wave: float):
    """Local stiffness and geometric matrices, DOF order [u1 v1 w1 th1 u2 v2 w2 th2]."""
    E1 = sheet.E / (1 - sheet.nu ** 2); t = sheet.t
    Dm = t * np.array([[E1, sheet.nu * E1, 0], [sheet.nu * E1, E1, 0], [0, 0, sheet.G]])
    D = sheet.D11
    Db = np.array([[D, sheet.nu * D, 0], [sheet.nu * D, D, 0], [0, 0, D * (1 - sheet.nu) / 2]])
    Ke = np.zeros((8, 8)); Kg = np.zeros((8, 8))
    iu, iv, iw = [0, 4], [1, 5], [2, 3, 6, 7]
    for gp, gw in zip(_GP, _GW):
        s = (gp + 1) / 2; wt = gw * b / 2 * half_wave / 2
        N = np.array([1 - s, s]); dN = np.array([-1, 1]) / b
        H = np.array([1 - 3 * s ** 2 + 2 * s ** 3, b * (s - 2 * s ** 2 + s ** 3),
                      3 * s ** 2 - 2 * s ** 3, b * (s ** 3 - s ** 2)])
        dH = np.array([-6 * s + 6 * s ** 2, b * (1 - 4 * s + 3 * s ** 2),
                       6 * s - 6 * s ** 2, b * (3 * s ** 2 - 2 * s)]) / b
        ddH = np.array([-6 + 12 * s, b * (-4 + 6 * s), 6 - 12 * s, b * (6 * s - 2)]) / b ** 2
        Bm = np.zeros((3, 8))
        Bm[0, iu] = dN; Bm[1, iv] = -k * N; Bm[2, iu] = k * N; Bm[2, iv] = dN
        Bb = np.zeros((3, 8))
        Bb[0, iw] = -ddH; Bb[1, iw] = k ** 2 * H; Bb[2, iw] = 2 * k * dH
        Ke += wt * (Bm.T @ Dm @ Bm + Bb.T @ Db @ Bb)
        G = np.zeros((3, 8)); G[0, iu] = N; G[1, iv] = N; G[2, iw] = H
        Kg += wt * t * k ** 2 * (G.T @ G)          # unit compressive stress
    return Ke, Kg


def buckle(sheet: Sheet, profile: Profile, width: float, half_wave: float,
           edges: str = "supported"):
    """Critical line load [N/m] for one buckle half-wavelength, with its mode.

    edges  'supported' holds the two unloaded edges against out-of-plane
           movement, 'free' leaves them free.
    Returns (N_cr, nodes, mode) with mode = (U, V, W, theta) per node.
    """
    periods = max(1, int(round(width / (2 * profile.c))))
    xz = profile.polyline(periods)
    xz = xz * [width / xz[-1, 0], 1.0]
    n = len(xz); ndof = 4 * n
    K = np.zeros((ndof, ndof)); Kg = np.zeros((ndof, ndof))
    k = np.pi / half_wave
    for e in range(n - 1):
        d = xz[e + 1] - xz[e]; b = np.hypot(*d); cs, sn = d / b
        ke, kg = _strip(sheet, b, k, half_wave)
        R = np.eye(8)
        for o in (0, 4):
            R[o, o], R[o, o + 2], R[o + 2, o], R[o + 2, o + 2] = cs, sn, -sn, cs
        idx = np.r_[4 * e:4 * e + 8]
        K[np.ix_(idx, idx)] += R.T @ ke @ R; Kg[np.ix_(idx, idx)] += R.T @ kg @ R
    fixed = [2, 4 * (n - 1) + 2] if edges == "supported" else []
    free = np.setdiff1d(np.arange(ndof), fixed)
    mu, vec = eigh(Kg[np.ix_(free, free)], K[np.ix_(free, free)])
    j = np.argmax(mu)
    full = np.zeros(ndof); full[free] = vec[:, j]
    mode = full.reshape(n, 4)
    mode = mode / np.abs(mode[:, [0, 2]]).max()
    stress = 1 / mu[j]
    return stress * sheet.t * profile.stretch, xz, mode


def signature_curve(sheet, profile, width, half_waves, edges="supported"):
    """Critical line load against buckle half-wavelength."""
    return np.array([buckle(sheet, profile, width, hw, edges)[0] for hw in half_waves])


def panel_load(sheet, profile, width, length, ends="pinned", edges="supported", max_m=40):
    """Lowest critical load of a panel of given length, and its half-wave count.

    Pinned ends admit half-wavelengths length/m. Fixed ends are approximated
    by the effective length length/2, the usual column result.
    """
    eff = length * {"pinned": 1.0, "fixed": 0.5, "pinned-fixed": 0.699156, "fixed-free": 2.0}[ends]
    best = (np.inf, 1)
    for m in range(1, max_m + 1):
        N = buckle(sheet, profile, width, eff / m, edges)[0]
        if N < best[0]:
            best = (N, m)
        elif m > best[1] + 6:
            break
    return best
