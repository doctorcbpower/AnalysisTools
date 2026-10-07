#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
analysistools.orbits.orphans
----------------------------
Orphan-galaxy orbits: continue a satellite whose subhalo the tree has lost
by integrating it in the host's (time-dependent, spherical) potential from
its last resolved state.

Assumptions of this module -- not facts about the data:
  * orbit in a *spherical* potential about the host centre; the host's own
    motion (a non-inertial frame) and the cosmological expansion terms are
    neglected, which is reasonable for objects bound inside ~R200 over a few
    Gyr and degrades for long integrations from early times;
  * no mass loss: the orphan keeps its last-resolved mass (matters only if
    dynamical friction is on);
  * dynamical friction is *off* by default. When on, it is Chandrasekhar's
    formula with an isotropic-isothermal velocity distribution
    (sigma = Vc / sqrt(2)) and a Coulomb logarithm ln(1 + M_host(<r)/M_sat).
    That is a modelling choice with a material effect on orphan radii
    (it makes them sink), so it is never silently enabled.
"""
from __future__ import annotations

from typing import Optional, Sequence

import numpy as np
from scipy.special import erf

from .cosmology import G_KPC_KMS2_MSUN as G, KMS_TO_KPC_PER_GYR
from .potential import Potential


def _acceleration(pot: Potential, x, v, t, m_sat, dynamical_friction):
    """Acceleration (km/s)^2/kpc, shape (N,3)."""
    r = np.linalg.norm(x, axis=1)
    g = -pot.accel(r, t)[:, None] * x / np.maximum(r, 1e-12)[:, None]
    if dynamical_friction:
        vm = np.maximum(np.linalg.norm(v, axis=1), 1e-6)
        sigma = pot.vcirc(r, t) / np.sqrt(2.0)
        X = vm / (np.sqrt(2.0) * sigma)
        rho = pot.density(r, t)
        lnl = np.log1p(pot.menc(r, t) / m_sat)
        bracket = erf(X) - 2.0 * X / np.sqrt(np.pi) * np.exp(-X * X)
        mag = 4.0 * np.pi * G * G * m_sat * rho * lnl * bracket / vm ** 2
        g = g - mag[:, None] * v / vm[:, None]
    return g


def integrate_orphans(potential: Potential, x0, v0, t_out: Sequence[float],
                      m_sat=None, dynamical_friction: bool = False,
                      dt_max: float = 0.005):
    """Integrate ``N`` orphans together in a spherical potential.

    Parameters
    ----------
    x0, v0 : (N,3) host-centric *physical* position (kpc) and velocity
        (km/s) at ``t_out[0]``.
    t_out : cosmic times (Gyr), ascending; the state is returned at each.
        Step sizes are subdivided so that no step exceeds ``dt_max`` (Gyr)
        and every output time is hit exactly.
    m_sat : (N,) satellite masses in Msun (only used by dynamical friction).

    Returns ``(x, v)``, each ``(N, len(t_out), 3)``. Kick-drift-kick leapfrog
    (second order, symplectic when friction is off).
    """
    x = np.array(x0, dtype=float, copy=True)
    v = np.array(v0, dtype=float, copy=True)
    t_out = np.asarray(t_out, dtype=float)
    n = x.shape[0]
    if dynamical_friction and m_sat is None:
        raise ValueError("dynamical_friction=True needs m_sat.")
    m_sat = None if m_sat is None else np.asarray(m_sat, dtype=float)
    xs = np.empty((n, t_out.size, 3))
    vs = np.empty_like(xs)
    xs[:, 0], vs[:, 0] = x, v
    # (km/s)^2/kpc * Gyr -> km/s: 1 km/s = 1.0227 kpc/Gyr, so
    # (km/s)*Gyr/kpc = 1.0227
    inv = KMS_TO_KPC_PER_GYR
    t = t_out[0]
    acc = _acceleration(potential, x, v, t, m_sat, dynamical_friction)
    for k in range(1, t_out.size):
        span = t_out[k] - t_out[k - 1]
        nstep = max(1, int(np.ceil(span / dt_max)))
        dt = span / nstep
        for _ in range(nstep):
            v += 0.5 * dt * acc * inv
            x += dt * v * KMS_TO_KPC_PER_GYR
            t += dt
            acc = _acceleration(potential, x, v, t, m_sat, dynamical_friction)
            v += 0.5 * dt * acc * inv
        t = t_out[k]
        xs[:, k], vs[:, k] = x, v
    return xs, vs
