"""Cross-section geometry of flat, trapezoidal and round corrugated panels.

Coordinates: x runs across the corrugations, y along them, z is out of plane.
One period has width 2c and developed (unfolded) length 2l. Every profile
starts and ends on the mid-plane z = 0, so the line of thrust of an in-plane
load passes through the centroid of the section.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class Profile:
    """One period of a corrugation, described by straight segments."""

    kind: str
    c: float            # half period [m]
    l: float            # half developed length [m]
    I1: float           # integral of (dx/ds)^2 over one period [m]
    I2: float           # integral of z^2 over one period [m^3]
    depth: float        # total profile depth, peak to peak [m]
    points: np.ndarray  # (n, 2) polyline of one period, x from 0 to 2c

    @property
    def stretch(self) -> float:
        """Developed length per unit panel width, l / c (1 for a flat plate)."""
        return self.l / self.c

    def polyline(self, periods: int) -> np.ndarray:
        """Polyline of `periods` consecutive corrugations."""
        pts = [self.points]
        for k in range(1, periods):
            shifted = self.points[1:].copy()
            shifted[:, 0] += 2 * self.c * k
            pts.append(shifted)
        return np.vstack(pts)


def _subdivide(corners: np.ndarray, n: int) -> np.ndarray:
    out = [corners[0]]
    for p, q in zip(corners[:-1], corners[1:]):
        for i in range(1, n + 1):
            out.append(p + (q - p) * i / n)
    return np.array(out)


def flat(width: float = 0.1016, n: int = 8) -> Profile:
    """A flat plate, treated as a 'corrugation' of zero depth."""
    c = width / 2
    pts = _subdivide(np.array([[0.0, 0.0], [width, 0.0]]), n)
    return Profile("flat", c, c, 2 * c, 0.0, 0.0, pts)


def trapezoidal(f: float, c: float, theta: float, n: int = 4) -> Profile:
    """Trapezoidal corrugation.

    f      half depth (the profile runs from z = -f to z = +f) [m]
    c      half period [m]
    theta  web angle measured from the panel plane [rad]
    """
    run = 2 * f / np.tan(theta)          # horizontal projection of one web
    flat_w = c - run                     # width of each flat
    if flat_w < -1e-12:
        raise ValueError("webs overlap: need c >= 2 f / tan(theta)")
    flat_w = max(flat_w, 0.0)
    l = 2 * f / np.sin(theta) + flat_w
    I1 = 2 * flat_w + 4 * f * np.cos(theta) ** 2 / np.sin(theta)
    I2 = 4 * f ** 3 / (3 * np.sin(theta)) + 2 * f ** 2 * flat_w
    corners = [[0.0, 0.0], [run / 2, f]]
    if flat_w > 0:
        corners.append([run / 2 + flat_w, f])
    corners.append([run / 2 + flat_w + run, -f])
    if flat_w > 0:
        corners.append([run / 2 + 2 * flat_w + run, -f])
    corners.append([2 * c, 0.0])
    return Profile("trapezoidal", c, l, I1, I2, 2 * f, _subdivide(np.array(corners), n))


def round_(R: float, L: float, n_arc: int = 10, n_straight: int = 3) -> Profile:
    """Round corrugation: semicircles of radius R joined by vertical straights.

    R  radius of the semicircular crests [m]
    L  half length of each vertical straight [m]
    The half period is c = 2R and the profile depth is 2(R + L).
    """
    c = 2 * R
    l = np.pi * R + 2 * L
    I1 = np.pi * R
    I2 = 4 * L ** 3 / 3 + 2 * np.pi * L ** 2 * R + 8 * L * R ** 2 + np.pi * R ** 3
    pts = [[0.0, 0.0]]
    if L > 0:
        for i in range(1, n_straight + 1):
            pts.append([0.0, L * i / n_straight])
    for i in range(1, n_arc + 1):                       # crest
        a = np.pi * i / n_arc
        pts.append([R - R * np.cos(a), L + R * np.sin(a)])
    if L > 0:
        for i in range(1, 2 * n_straight + 1):
            pts.append([2 * R, L - L * i / n_straight])
    for i in range(1, n_arc + 1):                       # trough
        a = np.pi * i / n_arc
        pts.append([3 * R - R * np.cos(a), -L - R * np.sin(a)])
    if L > 0:
        for i in range(1, n_straight + 1):
            pts.append([4 * R, -L + L * i / n_straight])
    return Profile("round", c, l, I1, I2, 2 * (R + L), np.array(pts))
