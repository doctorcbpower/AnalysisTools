#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
analysistools.orbits.cosmology
------------------------------
Flat LambdaCDM background and unit constants for orbit tracking.

Internal units throughout ``analysistools.orbits``: kpc, km/s, Gyr, Msun.
Radiation is neglected (fine for a > ~0.01), flat curvature assumed.
"""
from __future__ import annotations

import numpy as np

# G in kpc (km/s)^2 / Msun
G_KPC_KMS2_MSUN = 4.30091727e-6
# 1 km/s expressed in kpc/Gyr
KMS_TO_KPC_PER_GYR = 1.0227121650537077
# 100 km/s/Mpc expressed in 1/Gyr
H100_PER_GYR = 0.10227121650537077


class FlatLCDM:
    """Flat LambdaCDM: age, H(a) and the exact inverse a(t)."""

    def __init__(self, h: float, omega_m: float, omega_lambda: float = None):
        self.h = float(h)
        self.omega_m = float(omega_m)
        self.omega_lambda = (1.0 - self.omega_m if omega_lambda is None
                             else float(omega_lambda))
        if abs(self.omega_m + self.omega_lambda - 1.0) > 1e-3:
            raise ValueError(
                "FlatLCDM needs Omega_m + Omega_Lambda = 1 (got "
                f"{self.omega_m + self.omega_lambda:.4f}).")
        self.H0 = H100_PER_GYR * self.h            # 1/Gyr

    def E(self, a):
        a = np.asarray(a, dtype=float)
        return np.sqrt(self.omega_m * a ** -3 + self.omega_lambda)

    def hubble(self, a):
        """H(a) in 1/Gyr."""
        return self.H0 * self.E(a)

    def age(self, a):
        """Cosmic time at scale factor ``a``, Gyr."""
        a = np.asarray(a, dtype=float)
        return (2.0 / (3.0 * self.H0 * np.sqrt(self.omega_lambda))
                * np.arcsinh(np.sqrt(self.omega_lambda / self.omega_m)
                             * a ** 1.5))

    def a_of_t(self, t):
        """Inverse of :meth:`age`."""
        t = np.asarray(t, dtype=float)
        x = 1.5 * np.sqrt(self.omega_lambda) * self.H0 * t
        return (self.omega_m / self.omega_lambda) ** (1.0 / 3.0) \
            * np.sinh(x) ** (2.0 / 3.0)

    @classmethod
    def from_tree(cls, tree) -> "FlatLCDM":
        """Build from a ``MergerTreeTools`` instance carrying SubFind-HBT
        ``Parameters`` (``Omega0``/``OmegaLambda``/``HubbleParam``)."""
        d = tree.data
        return cls(h=float(d.HubbleParam), omega_m=float(d.Omega0),
                   omega_lambda=float(d.OmegaLambda))

    def __repr__(self) -> str:
        return (f"FlatLCDM(h={self.h:g}, Om={self.omega_m:g}, "
                f"OL={self.omega_lambda:g})")
