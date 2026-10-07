#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
analysistools.orbits.potential
------------------------------
Spherical host potentials for orbit analysis: closed-form NFW, a tabulated
spherical potential built from an enclosed-mass profile (typically straight
from snapshot particles), and a time-dependent stack of the latter.

Everything is spherical by construction -- a host's triaxiality, its disc
and the (non-inertial) motion of its centre are *not* represented. That is
an assumption of this module, not a property of the data; it is adequate for
a first-order E, L, pericentre/apocentre and orphan-orbit estimate and is
the thing to revisit if those results turn out to matter at the
tens-of-percent level.

Units: kpc, km/s, Msun, Gyr (see ``cosmology.py``). ``Phi`` -> 0 at infinity.
Every method takes ``r`` (physical kpc, array-like) and an optional cosmic
time ``t`` (Gyr); static potentials ignore ``t``.
"""
from __future__ import annotations

from typing import Optional, Sequence

import numpy as np
from scipy.optimize import brentq

from .cosmology import G_KPC_KMS2_MSUN as G


class Potential:
    """Interface: subclasses provide ``menc`` and ``phi``."""

    def menc(self, r, t=None):
        raise NotImplementedError

    def phi(self, r, t=None):
        raise NotImplementedError

    # -- derived -----------------------------------------------------
    def accel(self, r, t=None):
        """Inward gravitational acceleration magnitude, (km/s)^2/kpc."""
        r = np.asarray(r, dtype=float)
        return G * self.menc(r, t) / np.maximum(r, 1e-12) ** 2

    def vcirc(self, r, t=None):
        r = np.asarray(r, dtype=float)
        return np.sqrt(G * self.menc(r, t) / np.maximum(r, 1e-12))

    def density(self, r, t=None, eps: float = 1e-3):
        """rho(r) from dM/dr (central difference in ln r), Msun/kpc^3."""
        r = np.asarray(r, dtype=float)
        lo, hi = r * (1 - eps), r * (1 + eps)
        dm = self.menc(hi, t) - self.menc(lo, t)
        return dm / (4.0 * np.pi * r ** 2 * (hi - lo))

    def specific_energy(self, r, v, t=None):
        """0.5 v^2 + Phi(r), (km/s)^2."""
        return 0.5 * np.asarray(v, dtype=float) ** 2 + self.phi(r, t)

    # -- E, L -> turning points ---------------------------------------
    def turning_points(self, E, L, r0, t=None, r_search_max: float = 1e4):
        """Pericentre and apocentre of the orbit through radius ``r0`` with
        specific energy ``E`` and angular momentum ``L`` (kpc km/s), in this
        (frozen at ``t``) spherical potential.

        Returns ``(r_peri, r_apo)``; ``r_apo = inf`` if the orbit is
        unbound / reaches ``r_search_max``, ``r_peri = 0`` for a purely
        radial orbit. NaN if the inputs are not self-consistent
        (f(r0) < 0 beyond roundoff).
        """
        E, L, r0 = float(E), abs(float(L)), float(r0)

        def f(r):
            return 2.0 * (E - float(self.phi(r, t))) - L * L / (r * r)

        f0 = f(r0)
        if not np.isfinite(f0) or f0 < -1e-6 * max(abs(2 * E), 1.0):
            return (np.nan, np.nan)

        # apocentre: walk outward until f < 0
        r_hi = r0
        while True:
            r_next = r_hi * 1.2
            if r_next > r_search_max:
                r_apo = np.inf
                break
            if f(r_next) < 0:
                r_apo = brentq(f, r_hi, r_next, xtol=1e-8 * r_hi)
                break
            r_hi = r_next

        # pericentre: walk inward until f < 0
        if L == 0.0:
            return (0.0, r_apo)
        r_lo = r0
        r_peri = 0.0
        while True:
            r_next = r_lo / 1.2
            if r_next < 1e-6:
                r_peri = 0.0
                break
            if f(r_next) < 0:
                r_peri = brentq(f, r_next, r_lo, xtol=1e-8 * r_lo)
                break
            r_lo = r_next
        return (r_peri, r_apo)

    def radial_period(self, r_peri, r_apo, E, L, t=None, n: int = 64):
        """Radial period 2 * integral dr / v_r between the turning points,
        Gyr (Gauss-Legendre in theta after r = r_m - d cos(theta), which
        removes the inverse-sqrt endpoint singularities)."""
        if not (np.isfinite(r_peri) and np.isfinite(r_apo)) or r_apo <= r_peri:
            return np.nan
        x, w = np.polynomial.legendre.leggauss(n)
        theta = 0.5 * np.pi * (x + 1.0)
        wt = 0.5 * np.pi * w
        rm, d = 0.5 * (r_apo + r_peri), 0.5 * (r_apo - r_peri)
        r = rm - d * np.cos(theta)
        f = 2.0 * (E - self.phi(r, t)) - (L / r) ** 2
        f = np.maximum(f, 1e-30)
        integral = np.sum(wt * d * np.sin(theta) / np.sqrt(f))   # kpc/(km/s)
        return 2.0 * integral / 1.0227121650537077               # -> Gyr


# ---------------------------------------------------------------------
# NFW
# ---------------------------------------------------------------------

def _nfw_m(x):
    x = np.asarray(x, dtype=float)
    return np.log1p(x) - x / (1.0 + x)


class NFWPotential(Potential):
    """NFW halo normalised so that M(<r200) = m200."""

    def __init__(self, m200: float, r200: float, rs: float):
        self.m200, self.r200, self.rs = float(m200), float(r200), float(rs)
        self.c = self.r200 / self.rs
        self._norm = self.m200 / float(_nfw_m(self.c))      # = 4 pi rho0 rs^3

    @classmethod
    def from_vmax_rmax(cls, m200: float, r200: float, vmax: float = None,
                       rmax: float = None, vmax_tolerance: float = 0.3):
        """rs from the Vmax radius (``Rmax = 2.163 rs`` for an NFW),
        normalisation from (M200, R200). ``vmax`` is accepted for a
        consistency check only (a mismatch of more than ~30 per cent
        between the NFW's implied Vmax and ``vmax`` (``vmax_tolerance``) raises ValueError --
        the halo is then not NFW-like and the caller should use a
        tabulated potential instead)."""
        rs = float(rmax) / 2.1626
        obj = cls(m200, r200, rs)
        if vmax is not None:
            v = float(obj.vcirc(rmax))
            if abs(v - vmax) > vmax_tolerance * vmax:
                raise ValueError(
                    f"NFW(M200={m200:.3g}, R200={r200:.3g}, rs={rs:.3g}) "
                    f"implies Vmax={v:.1f} at Rmax, tree says {vmax:.1f}; "
                    "halo is not NFW-like.")
        return obj

    def menc(self, r, t=None):
        return self._norm * _nfw_m(np.asarray(r, dtype=float) / self.rs)

    def phi(self, r, t=None):
        r = np.maximum(np.asarray(r, dtype=float), 1e-12)
        return -G * self._norm * np.log1p(r / self.rs) / r

    def __repr__(self):
        return (f"NFWPotential(m200={self.m200:.3e}, r200={self.r200:.3g}, "
                f"rs={self.rs:.3g}, c={self.c:.2f})")


# ---------------------------------------------------------------------
# Tabulated spherical potential
# ---------------------------------------------------------------------

class TabulatedSpherical(Potential):
    """Spherical potential from M(<r) tabulated on a radial grid.

    Between grid points M is interpolated linearly in (ln r, ln M); inside
    the first point M ~ r^3 (constant density), beyond the last M is
    constant (point mass), and Phi follows by integrating G M / r^2 inwards
    from -G M_tot / r_max.
    """

    def __init__(self, r_grid: Sequence[float], m_enc: Sequence[float]):
        r = np.asarray(r_grid, dtype=float)
        m = np.asarray(m_enc, dtype=float)
        if r.ndim != 1 or r.shape != m.shape or r.size < 3:
            raise ValueError("r_grid/m_enc must be matching 1-D arrays "
                             "with at least 3 points.")
        if np.any(np.diff(r) <= 0) or np.any(r <= 0):
            raise ValueError("r_grid must be positive and strictly increasing.")
        if np.any(m <= 0):
            raise ValueError("m_enc must be positive everywhere on the grid "
                             "(raise rmin so the first shell is populated).")
        m = np.maximum.accumulate(m)            # enclosed mass is monotonic
        self.r, self.m = r, m
        self._lr, self._lm = np.log(r), np.log(m)
        # Phi on the grid: integrate G M / r^2 dr = G M / r dln r inwards
        integrand = G * m / r
        seg = 0.5 * (integrand[1:] + integrand[:-1]) * np.diff(self._lr)
        tail = np.concatenate([np.cumsum(seg[::-1])[::-1], [0.0]])
        # Phi increases outward (dPhi/dr = +g), so Phi(r) = Phi(rmax) - int_r^rmax g
        self._phi_grid = -G * m[-1] / r[-1] - tail

    # -- constructors ---------------------------------------------------
    @classmethod
    def from_particles(cls, pos, mass, centre, boxsize: float = None,
                       rmin: float = 1.0, rmax: float = None, n_grid: int = 80):
        """Spherical M(<r) about ``centre`` from particles.

        ``pos`` and ``centre`` in physical kpc (caller converts), ``mass``
        scalar or per-particle, Msun. ``boxsize`` (physical kpc) enables
        periodic wrapping. M(<r) is evaluated exactly at ``n_grid``
        log-spaced radii by sorting particle radii -- no histogram
        binning -- so the cost is O(N log N) in the particles inside
        ``rmax``. Choose ``rmin`` of a few softening lengths: M(<r) below
        that is not trustworthy.
        """
        pos = np.asarray(pos, dtype=float)
        d = pos - np.asarray(centre, dtype=float)
        if boxsize:
            d = (d + 0.5 * boxsize) % boxsize - 0.5 * boxsize
        rp = np.linalg.norm(d, axis=1)
        mass = np.broadcast_to(np.asarray(mass, dtype=float), rp.shape)
        if rmax is None:
            rmax = float(rp.max())
        keep = rp <= rmax
        rp, mp = rp[keep], mass[keep]
        order = np.argsort(rp)
        rp, cm = rp[order], np.cumsum(mp[order])
        grid = np.geomspace(rmin, rmax, n_grid)
        idx = np.searchsorted(rp, grid, side="right")
        m_grid = np.where(idx > 0, cm[np.maximum(idx - 1, 0)], 0.0)
        return cls(grid, m_grid)

    @classmethod
    def from_potential(cls, other: Potential, r_grid, t=None):
        """Tabulate any other potential on ``r_grid``."""
        r = np.asarray(r_grid, dtype=float)
        return cls(r, other.menc(r, t))

    # -- evaluation -----------------------------------------------------
    def menc(self, r, t=None):
        r = np.asarray(r, dtype=float)
        lr = np.log(np.maximum(r, 1e-12))
        inner = self.m[0] * (r / self.r[0]) ** 3
        mid = np.exp(np.interp(lr, self._lr, self._lm))
        return np.where(r < self.r[0], inner,
                        np.where(r > self.r[-1], self.m[-1], mid))

    def phi(self, r, t=None):
        r = np.asarray(r, dtype=float)
        lr = np.log(np.maximum(r, 1e-12))
        mid = np.interp(lr, self._lr, self._phi_grid)
        inner = self._phi_grid[0] - (G * self.m[0] / (2 * self.r[0] ** 3)) \
            * (self.r[0] ** 2 - r ** 2)
        outer = -G * self.m[-1] / np.maximum(r, 1e-12)
        return np.where(r < self.r[0], inner,
                        np.where(r > self.r[-1], outer, mid))


# ---------------------------------------------------------------------
# Time-dependent stack
# ---------------------------------------------------------------------

class TimeDependentPotential(Potential):
    """A sequence of spherical potentials at cosmic times ``times`` (Gyr),
    resampled on one common radial grid; M(<r, t) is linearly interpolated
    between neighbouring times (and held fixed outside the range). Phi is
    interpolated the same way (not re-derived from the interpolated M, which
    differs only at second order in the time step).
    """

    def __init__(self, times: Sequence[float],
                 potentials: Sequence[Potential],
                 r_grid: Optional[Sequence[float]] = None):
        t = np.asarray(times, dtype=float)
        if t.ndim != 1 or t.size != len(potentials) or t.size < 1:
            raise ValueError("times and potentials must match in length.")
        order = np.argsort(t)
        self.t = t[order]
        pots = [potentials[i] for i in order]
        if r_grid is None:
            r_grid = np.geomspace(1.0, 1.0e3, 120)
        self.r = np.asarray(r_grid, dtype=float)
        self._lr = np.log(self.r)
        self._m = np.array([p.menc(self.r) for p in pots])
        self._phi = np.array([p.phi(self.r) for p in pots])

    def _weights(self, t):
        if t is None:
            raise ValueError("TimeDependentPotential needs a time `t` (Gyr).")
        t = float(t)
        if self.t.size == 1 or t <= self.t[0]:
            return 0, 0, 0.0
        if t >= self.t[-1]:
            return self.t.size - 1, self.t.size - 1, 0.0
        k = int(np.searchsorted(self.t, t, side="right") - 1)
        return k, k + 1, (t - self.t[k]) / (self.t[k + 1] - self.t[k])

    def _eval(self, table, r, t, log=False):
        r = np.asarray(r, dtype=float)
        lr = np.log(np.maximum(r, 1e-12))
        k0, k1, w = self._weights(t)
        a = np.interp(lr, self._lr, table[k0])
        if k1 == k0 or w == 0.0:
            return a
        b = np.interp(lr, self._lr, table[k1])
        return (1 - w) * a + w * b

    def menc(self, r, t=None):
        r = np.asarray(r, dtype=float)
        v = self._eval(self._m, r, t)       # linear in t, log-r in space
        v = np.where(r < self.r[0], v * (r / self.r[0]) ** 3, v)
        return v

    def phi(self, r, t=None):
        r = np.asarray(r, dtype=float)
        v = self._eval(self._phi, r, t)
        k0, k1, w = self._weights(t)
        mtot = (1 - w) * self._m[k0, -1] + w * self._m[k1, -1]
        return np.where(r > self.r[-1], -G * mtot / np.maximum(r, 1e-12), v)
