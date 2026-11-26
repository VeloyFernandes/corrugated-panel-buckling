"""Equivalent orthotropic plate properties of a corrugated panel.

Closed-form homogenisation of Xia, Friswell and Saavedra Flores (2012),
"Equivalent models of corrugated panels", Int. J. Solids Struct. 49, 1453-1462.
Direction 1 (x) is across the corrugations, direction 2 (y) along them.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
from .profiles import Profile


@dataclass(frozen=True)
class Sheet:
    """Isotropic sheet the panel is formed from."""

    E: float = 2.1e11        # Young's modulus [Pa]
    nu: float = 0.3          # Poisson's ratio [-]
    t: float = 0.00635       # thickness [m]
    yield_stress: float = 250e6   # [Pa]
    density: float = 7850.0  # [kg/m^3]

    @property
    def G(self): return self.E / (2 * (1 + self.nu))
    @property
    def A11(self): return self.E * self.t / (1 - self.nu ** 2)
    @property
    def A12(self): return self.nu * self.A11
    @property
    def A66(self): return self.G * self.t
    @property
    def D11(self): return self.E * self.t ** 3 / (12 * (1 - self.nu ** 2))
    @property
    def D12(self): return self.nu * self.D11
    @property
    def D66(self): return self.G * self.t ** 3 / 12


@dataclass(frozen=True)
class Stiffness:
    A11: float; A12: float; A22: float; A66: float   # [N/m]
    D11: float; D12: float; D22: float; D66: float   # [N m]

    def as_dict(self): return asdict(self)


def equivalent_plate(sheet: Sheet, p: Profile) -> Stiffness:
    """Membrane (A) and bending (D) stiffness of the equivalent flat plate."""
    A11, A12, A22, A66 = sheet.A11, sheet.A12, sheet.A11, sheet.A66
    D11, D12, D22, D66 = sheet.D11, sheet.D12, sheet.D11, sheet.D66
    c, l = p.c, p.l
    eA11 = 2 * c / (p.I1 / A11 + p.I2 / D11)
    eA12 = A12 / A11 * eA11
    eA22 = A12 / A11 * eA12 + (l / c) * (A11 * A22 - A12 ** 2) / A11
    eA66 = (c / l) * A66
    eD11 = (c / l) * D11
    eD12 = D12 / D11 * eD11
    eD22 = (p.I2 * A22 + p.I1 * D22) / (2 * c)
    eD66 = (l / c) * D66
    return Stiffness(eA11, eA12, eA22, eA66, eD11, eD12, eD22, eD66)
