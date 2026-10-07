#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
analysistools.orbits.benson
---------------------------
Infall orbital parameters in the conventions of Benson (2005, MNRAS 358, 551):

* radii in units of the host's virial radius r_vir and velocities in units of
  V_vir, with r_vir the radius enclosing the mean density of a just-collapsed
  spherical top-hat (here via the Bryan & Norman 1998 fit for flat LCDM),
  *not* R200c;
* the Hubble flow is added to the relative peculiar velocity (Benson 2005,
  footnote 3), which matters: at r_vir it shifts V_r by sqrt(2/Delta_vir)
  (~0.14 V_vir at z=0);
* V_r > 0 is *inward*; V_theta is the tangential speed; the infall angle is
  phi = arctan(V_theta / V_r);
* two-body quantities for point masses with f2 = 1 + M_sat/M_host:
  specific energy (eqs 1/5), eccentricity (eq 7), pericentre/apocentre
  (eqs 10/11);
* the fitting function of eqs (2)-(4) with the parameters of Table 2
  (z = 0, 0.5, 1.0).

Equation/table numbers refer to the paper.
"""
from __future__ import annotations

import numpy as np

from .cosmology import FlatLCDM, G_KPC_KMS2_MSUN as G


def omega_m_of_z(cosmology: FlatLCDM, z):
    a = 1.0 / (1.0 + np.asarray(z, dtype=float))
    return cosmology.omega_m * a ** -3 / cosmology.E(a) ** 2


def bryan_norman_delta(cosmology: FlatLCDM, z):
    """Bryan & Norman (1998) virial overdensity relative to the *critical*
    density for flat LCDM: 18 pi^2 + 82 x - 39 x^2, x = Omega_m(z) - 1.
    (~102 at z=0 for Omega_m=0.3; 18 pi^2 for Einstein-de Sitter.)"""
    x = omega_m_of_z(cosmology, z) - 1.0
    return 18.0 * np.pi ** 2 + 82.0 * x - 39.0 * x ** 2


def v_vir_from_r_vir(cosmology: FlatLCDM, a, r_vir, delta):
    """V_vir = sqrt(G M_vir / r_vir) = H r_vir sqrt(Delta/2)  (km/s) for
    M_vir = Delta (4 pi / 3) rho_crit r_vir^3; r_vir in physical kpc."""
    H = cosmology.hubble(a) / 1.0227121650537077          # 1/Gyr -> km/s/kpc
    return H * np.asarray(r_vir) * np.sqrt(np.asarray(delta) / 2.0)


def nfw_virial_radius(m200, r200, rs, delta, cosmology: FlatLCDM, a):
    """Radius (physical kpc) where an NFW halo of (M200, R200, rs) encloses
    ``delta`` times the critical density at scale factor ``a``. Vectorised
    over host snapshots. This needs a density profile, which the tree does
    not carry: ``rs`` (hence the concentration) comes from the host's Rmax
    and is therefore uncertain -- r_vir/R200c varies only slowly with
    concentration, which the caller should quantify."""
    from scipy.optimize import brentq
    m200, r200, rs = (np.atleast_1d(np.asarray(v, dtype=float))
                      for v in (m200, r200, rs))
    delta = np.broadcast_to(np.asarray(delta, dtype=float), m200.shape)
    a = np.broadcast_to(np.asarray(a, dtype=float), m200.shape)
    out = np.empty_like(m200)
    mnfw = lambda x: np.log1p(x) - x / (1.0 + x)
    for i in range(m200.size):
        rho_c = 3.0 * (100.0 * cosmology.h * float(cosmology.E(a[i])) / 1e3) ** 2 \
            / (8.0 * np.pi * G)
        norm = m200[i] / mnfw(r200[i] / rs[i])
        f = lambda r: norm * mnfw(r / rs[i]) - delta[i] * (4 * np.pi / 3) * rho_c * r ** 3
        out[i] = brentq(f, 1e-3 * r200[i], 10.0 * r200[i])
    return out


# ---------------------------------------------------------------------
# two-body quantities in Benson's units (r_vir = V_vir = 1)
# ---------------------------------------------------------------------

def f2_of(m_sat, m_host):
    """f2 = 1 + M2/M1 (eq. after 1)."""
    return 1.0 + np.asarray(m_sat, dtype=float) / np.asarray(m_host, dtype=float)


def specific_energy(vr, vt, f2=1.0):
    """Eq. (5): E = -1 + [Vr^2 + ((2-f2)^2/f2^2) Vt^2] / (2 f2)."""
    vr, vt, f2 = (np.asarray(v, dtype=float) for v in (vr, vt, f2))
    return -1.0 + (vr ** 2 + ((2.0 - f2) ** 2 / f2 ** 2) * vt ** 2) / (2.0 * f2)


def _root(vr, vt, f2):
    return np.sqrt((1.0 - f2 / vt ** 2) ** 2 + (vr / vt) ** 2)


def eccentricity(vr, vt, f2=1.0):
    """Eq. (7)."""
    vr, vt, f2 = (np.asarray(v, dtype=float) for v in (vr, vt, f2))
    return (vt ** 2 / f2) * _root(vr, vt, f2)


def pericentre(vr, vt, f2=1.0):
    """Eq. (10), in units of r_vir."""
    vr, vt, f2 = (np.asarray(v, dtype=float) for v in (vr, vt, f2))
    return 1.0 / (f2 / vt ** 2 + _root(vr, vt, f2))


def apocentre(vr, vt, f2=1.0):
    """Eq. (11), in units of r_vir (negative/inf: unbound)."""
    vr, vt, f2 = (np.asarray(v, dtype=float) for v in (vr, vt, f2))
    with np.errstate(divide="ignore"):
        return 1.0 / (f2 / vt ** 2 - _root(vr, vt, f2))


# ---------------------------------------------------------------------
# fitting function, eqs (2)-(4), Table 2
# ---------------------------------------------------------------------

#: Table 2: a1..a9 at z = 0.0, 0.5, 1.0 (combined VLS + VIRGO LCDM)
BENSON05_TABLE2 = {
    0.0: (3.90, 2.49, 10.2, 0.684, 0.354, 1.08, 0.510, 0.206, 0.315),
    0.5: (4.46, 2.98, 11.0, 1.11, 0.494, 1.16, 0.261, -0.279, 0.331),
    1.0: (6.38, 2.30, 18.8, 0.506, -0.0934, 1.05, 0.267, -0.154, 0.157),
}


def benson05_fit(vr, vt, z=0.0):
    """Eq. (2): f(Vr, Vt) = a1 Vt exp{-a2 (Vt-a9)^2 - b1(Vt) [Vr - b2(Vt)]^2},
    with b1 = a3 exp[-a4 (Vt-a5)^2], b2 = a6 exp[-a7 (Vt-a8)^2]. ``z`` must be
    one of the tabulated redshifts (0, 0.5, 1.0). Not normalised by the paper
    to exactly 1; see :func:`benson05_grid`."""
    if z not in BENSON05_TABLE2:
        raise ValueError(f"Table 2 has z in {sorted(BENSON05_TABLE2)}, got {z}.")
    a1, a2, a3, a4, a5, a6, a7, a8, a9 = BENSON05_TABLE2[z]
    vr, vt = np.asarray(vr, dtype=float), np.asarray(vt, dtype=float)
    b1 = a3 * np.exp(-a4 * (vt - a5) ** 2)
    b2 = a6 * np.exp(-a7 * (vt - a8) ** 2)
    return a1 * vt * np.exp(-a2 * (vt - a9) ** 2 - b1 * (vr - b2) ** 2)


def benson05_grid(z=0.0, vmax=2.5, n=501):
    """The fit on a regular grid, normalised to unit integral.

    Returns ``(vr, vt, pdf, norm)`` with ``pdf[i, j]`` at (vr[i], vt[j]) and
    ``norm`` the integral of the *unnormalised* eq. (2) over the grid (close
    to 1 if the paper's fit is itself a normalised density)."""
    v = np.linspace(0.0, vmax, n)
    VR, VT = np.meshgrid(v, v, indexing="ij")
    f = benson05_fit(VR, VT, z)
    norm = np.trapezoid(np.trapezoid(f, v, axis=1), v)
    return v, v, f / norm, float(norm)


def benson05_marginals(z=0.0, vmax=2.5, n=501):
    """(v, df/dVr, df/dVt, phi_deg, df/dphi): marginal densities implied by
    the fit, and the infall-angle density (phi = arctan(Vt/Vr), degrees)."""
    v, _, pdf, _ = benson05_grid(z, vmax, n)
    dvr = np.trapezoid(pdf, v, axis=1)            # integrate out Vt
    dvt = np.trapezoid(pdf, v, axis=0)            # integrate out Vr
    VR, VT = np.meshgrid(v, v, indexing="ij")
    phi = np.degrees(np.arctan2(VT, VR)).ravel()
    bins = np.linspace(0.0, 90.0, 91)
    w = (pdf * (v[1] - v[0]) ** 2).ravel()
    h, _ = np.histogram(phi, bins=bins, weights=w)
    centres = 0.5 * (bins[1:] + bins[:-1])
    return v, dvr, dvt, centres, h / np.diff(bins)
