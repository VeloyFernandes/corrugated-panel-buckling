"""Eigenvalue buckling of the panel compressed ACROSS the corrugations.

A wide panel loaded across its corrugations bends cylindrically, so a strip of
unit width behaves as a folded plane-strain beam. The profile is meshed with
Euler-Bernoulli frame elements (axial stiffness A11, bending stiffness D11 of
the sheet) and the linear buckling problem K phi = lambda Kg phi is solved.
The geometric stiffness uses the axial force each element carries under a unit
end thrust, N = cos(alpha), alpha being the slope of the element.
"""
from __future__ import annotations

import numpy as np
from scipy.linalg import eigh
from .homogenise import Sheet
from .profiles import Profile


def _element(EA, EI, L, cs, sn, N):
    k = np.zeros((6, 6)); g = np.zeros((6, 6))
    a = EA / L; b = EI / L ** 3
    k[np.ix_([0, 3], [0, 3])] = a * np.array([[1, -1], [-1, 1]])
    k[np.ix_([1, 2, 4, 5], [1, 2, 4, 5])] = b * np.array([
        [12, 6 * L, -12, 6 * L],
        [6 * L, 4 * L ** 2, -6 * L, 2 * L ** 2],
        [-12, -6 * L, 12, -6 * L],
        [6 * L, 2 * L ** 2, -6 * L, 4 * L ** 2]])
    g[np.ix_([1, 2, 4, 5], [1, 2, 4, 5])] = N / (30 * L) * np.array([
        [36, 3 * L, -36, 3 * L],
        [3 * L, 4 * L ** 2, -3 * L, -L ** 2],
        [-36, -3 * L, 36, -3 * L],
        [3 * L, -L ** 2, -3 * L, 4 * L ** 2]])
    T = np.zeros((6, 6))
    for o in (0, 3):
        T[o:o + 3, o:o + 3] = [[cs, sn, 0], [-sn, cs, 0], [0, 0, 1]]
    return T.T @ k @ T, T.T @ g @ T


def buckle(sheet: Sheet, profile: Profile, length: float, ends: str = "pinned", n_modes: int = 1):
    """Critical line load [N/m] and mode shape for compression across the corrugations.

    Returns (N_cr, nodes, mode) for n_modes == 1, where `mode` holds the (u, w)
    displacement of every node, or (loads, nodes, modes) otherwise.
    """
    periods = max(1, int(round(length / (2 * profile.c))))
    xz = profile.polyline(periods)
    xz = xz * [length / xz[-1, 0], 1.0]            # fit a whole number of periods
    n = len(xz); ndof = 3 * n
    K = np.zeros((ndof, ndof)); Kg = np.zeros((ndof, ndof))
    for e in range(n - 1):
        d = xz[e + 1] - xz[e]; L = np.hypot(*d); cs, sn = d / L
        ke, ge = _element(sheet.A11, sheet.D11, L, cs, sn, cs)
        idx = np.r_[3 * e:3 * e + 6]
        K[np.ix_(idx, idx)] += ke; Kg[np.ix_(idx, idx)] += ge
    last = 3 * (n - 1)
    fixed = {"pinned": [0, 1, last + 1],
             "fixed": [0, 1, 2, last + 1, last + 2],
             "pinned-fixed": [0, 1, last + 1, last + 2],
             "fixed-free": [0, 1, 2]}[ends]
    free = np.setdiff1d(np.arange(ndof), fixed)
    # K phi = lam Kg phi  <=>  Kg phi = (1/lam) K phi, K positive definite
    mu, vec = eigh(Kg[np.ix_(free, free)], K[np.ix_(free, free)])
    order = np.argsort(-mu)[:n_modes]
    loads = 1 / mu[order]
    modes = []
    for j in order:
        full = np.zeros(ndof); full[free] = vec[:, j]
        m = full.reshape(n, 3)[:, :2]
        modes.append(m / np.abs(m).max())
    if n_modes == 1:
        return loads[0], xz, modes[0]
    return loads, xz, modes
