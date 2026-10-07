#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
analysistools.orbits.kinematics
-------------------------------
Host-centric orbit kinematics from tree tracks: unit handling, Hermite
interpolation between snapshots, and event detection (crossings of a
reference radius, pericentres, apocentres).

Why interpolate. At ~100-200 Myr snapshot spacing a satellite moving at a few
hundred km/s travels tens of kpc between snapshots, so the minimum of the
*sampled* distance systematically over-estimates pericentre. The tree gives
positions *and* velocities at every snapshot, so the separation is
interpolated with a cubic Hermite spline per interval.

Coordinate choice. Interpolation is done on the *comoving* host-centric
separation ``dx_c(t)``, whose time derivative is simply
``v_pec_rel / a`` -- no Hubble-flow term, no need for H(a) -- and the
physical separation is ``a(t) |dx_c|`` with ``a(t)`` exact from the
cosmology. Periodic wrapping is applied to ``dx_c`` before interpolation.

Assumptions (all inherited by every number derived here):
  * the host centre is the position of the host's central subhalo;
  * the tree's per-snapshot subhalo position/velocity are the quantities to
    interpolate (no sub-snapshot information is used);
  * a satellite's orbit is that of its subhalo, so stripped-to-nothing
    subhaloes are only followed as long as the tree follows them (the tracker
    extends them as orphans separately).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np

from .cosmology import FlatLCDM, KMS_TO_KPC_PER_GYR


# ---------------------------------------------------------------------
# Units
# ---------------------------------------------------------------------

@dataclass
class TrackUnits:
    """How a ``HaloTrack``'s numbers are expressed.

    Mirrors the ``comoving``/``little_h`` switches of ``MergerTreeTools``.
    ``velocity`` says how the *stored* velocity relates to the physical
    peculiar velocity: ``"physical"`` (stored = v_pec), ``"sqrt_a"``
    (stored = v_pec / sqrt(a), the GADGET-internal-style comoving velocity)
    or ``"inv_sqrt_a"``. ``"auto"`` guesses from ``comoving`` alone and is
    only right when the reader had to convert a *physical* native tree;
    ``OrbitTracker.from_epoch`` instead derives it from what the tree
    reader actually did (``native_is_comoving`` x ``comoving``). Verified
    directly for the Dorcha SubFind-HBT tree: positions comoving Mpc/h,
    velocities physical peculiar km/s (finite-difference dx/dt over 40
    snapshots matches v/a to 0.1 per cent, and v/sqrt(a) is off by 13 per
    cent).

    ``length_kpc`` is the size of one length unit in kpc (1000 for Mpc).
    ``mass_msun`` likewise for mass (1e10 for GADGET-style).
    """
    comoving: bool = True
    little_h: bool = False
    h: float = 0.6774
    length_kpc: float = 1000.0
    mass_msun: float = 1.0e10
    velocity: str = "auto"            # "auto" | "physical" | "sqrt_a" | "inv_sqrt_a"

    def _vel_mode(self) -> str:
        if self.velocity == "auto":
            return "sqrt_a" if self.comoving else "physical"
        if self.velocity not in ("physical", "sqrt_a", "inv_sqrt_a"):
            raise ValueError("velocity must be 'auto', 'physical', 'sqrt_a' "
                             "or 'inv_sqrt_a'.")
        return self.velocity

    def _hfac(self) -> float:
        return 1.0 / self.h if self.little_h else 1.0

    def comoving_kpc(self, length, a):
        """Convert length-like values to comoving kpc (h-free)."""
        x = np.asarray(length, dtype=float) * self.length_kpc * self._hfac()
        if not self.comoving:
            a = np.asarray(a, dtype=float)
            x = x / (a[..., None] if x.ndim > a.ndim else a)
        return x

    def physical_kpc(self, length, a):
        a = np.asarray(a, dtype=float)
        x = np.asarray(length, dtype=float) * self.length_kpc * self._hfac()
        if self.comoving:
            x = x * (a[..., None] if x.ndim > a.ndim else a)
        return x

    def peculiar_kms(self, vel, a):
        v = np.asarray(vel, dtype=float)
        if self._vel_mode() == "sqrt_a":
            a = np.asarray(a, dtype=float)
            v = v * np.sqrt(a)[..., None]
        elif self._vel_mode() == "inv_sqrt_a":
            a = np.asarray(a, dtype=float)
            v = v / np.sqrt(a)[..., None]
        return v

    def mass(self, m):
        return np.asarray(m, dtype=float) * self.mass_msun * self._hfac()


# ---------------------------------------------------------------------
# Hermite interpolation
# ---------------------------------------------------------------------

def hermite(p0, p1, m0, m1, h, s):
    """Cubic Hermite interpolant. p*, m*: (..., 3); h: (...,) interval
    length; s: (..., ns) fractions in [0, 1]. Returns (..., ns, 3)."""
    s = np.asarray(s, dtype=float)[..., None]
    h = np.asarray(h, dtype=float)[..., None, None]
    s2, s3 = s * s, s * s * s
    h00 = 2 * s3 - 3 * s2 + 1
    h10 = s3 - 2 * s2 + s
    h01 = -2 * s3 + 3 * s2
    h11 = s3 - s2
    return (h00 * p0[..., None, :] + h10 * h * m0[..., None, :]
            + h01 * p1[..., None, :] + h11 * h * m1[..., None, :])


# ---------------------------------------------------------------------
# Host-frame series
# ---------------------------------------------------------------------

@dataclass
class HostFrameSeries:
    """A satellite's orbit in the host frame on the *snapshot* grid.

    r (physical kpc), vr / vt (km/s, from the peculiar relative velocity;
    vr excludes the Hubble flow), L (kpc km/s, specific, physical), dx_c and
    v_pec_rel (the inputs the interpolation uses), and the reference radius
    ``r_ref`` (physical kpc) at each snapshot. ``t`` is cosmic time (Gyr).
    """
    snapnum: np.ndarray
    a: np.ndarray
    t: np.ndarray
    redshift: np.ndarray
    dx_c: np.ndarray            # (n,3) comoving kpc, host-centric
    v_rel: np.ndarray           # (n,3) km/s peculiar, relative to host
    r: np.ndarray
    vr: np.ndarray
    vt: np.ndarray
    L: np.ndarray               # (n,3)
    r_ref: np.ndarray           # physical kpc
    mass: np.ndarray            # Msun
    vmax: Optional[np.ndarray]  # km/s, if the track carries it
    is_subhalo: np.ndarray
    host_id_sat: np.ndarray     # FOF-central id the satellite sits in
    host_id_ref: np.ndarray     # id of the reference host's own central
    n_gaps: int = 0             # snapshots skipped between consecutive rows
    npart: Optional[np.ndarray] = None   # bound particle count, if the track has it

    def __len__(self):
        return len(self.snapnum)


R200C_FROM_M200 = "R200c_from_M200"


def r200c_from_m200(m200_msun, a, cosmology: FlatLCDM, overdensity=200.0):
    """Physical radius (kpc) enclosing ``overdensity`` x the *critical*
    density at scale factor ``a``, for mass ``m200_msun`` (Msun)."""
    from .cosmology import G_KPC_KMS2_MSUN as G
    h_kms_kpc = 100.0 * cosmology.h * cosmology.E(a) / 1000.0   # km/s/kpc
    rho_c = 3.0 * h_kms_kpc ** 2 / (8.0 * np.pi * G)            # Msun/kpc^3
    return (3.0 * np.asarray(m200_msun, dtype=float)
            / (4.0 * np.pi * overdensity * rho_c)) ** (1.0 / 3.0)


def host_reference(host_track, units: TrackUnits,
                   radius_key: str = "GroupR200", radius_factor: float = 1.0,
                   fixed_radius_kpc: Optional[float] = None,
                   cosmology: Optional[FlatLCDM] = None) -> np.ndarray:
    """Reference ("inside the host") radius at every snapshot of
    ``host_track``, physical kpc. See :func:`host_frame_series`.

    ``radius_key="R200c_from_M200"`` derives R200c from the track's
    ``extra["GroupM200"]`` and the cosmology (needs ``cosmology``) -- for
    trees that carry M200 but no R200 (the Dorcha SubFind-HBT trees).
    """
    n = len(host_track.SnapNum)
    if fixed_radius_kpc is not None:
        return np.full(n, float(fixed_radius_kpc))
    a = np.asarray(host_track.Time, dtype=float)
    if radius_key == R200C_FROM_M200:
        if cosmology is None:
            raise ValueError("radius_key='R200c_from_M200' needs cosmology.")
        m = host_track.extra.get("GroupM200")
        if m is None:
            raise KeyError("host track has no extra['GroupM200'].")
        return r200c_from_m200(units.mass(m), a, cosmology) * radius_factor
    rad = host_track.extra.get(radius_key)
    if rad is None:
        raise KeyError(
            f"host track has no extra['{radius_key}'] "
            f"(have {sorted(host_track.extra)}); pass radius_key= for a "
            f"field it carries (or 'R200c_from_M200'), or fixed_radius_kpc=.")
    return units.physical_kpc(np.asarray(rad), a) * radius_factor


def host_frame_series(track, host_track, cosmology: FlatLCDM,
                      units: TrackUnits, boxsize: Optional[float] = None,
                      radius_key: str = "GroupR200",
                      radius_factor: float = 1.0,
                      fixed_radius_kpc: Optional[float] = None
                      ) -> HostFrameSeries:
    """Express ``track`` in the frame of ``host_track`` on their common
    snapshots.

    The reference radius defining "inside the host" is, in order of
    precedence: ``fixed_radius_kpc`` (physical, constant); else
    ``radius_factor * host_track.extra[radius_key]`` (a length in the
    tree's own units, e.g. ``GroupR200``). Pass ``radius_factor=2`` for
    2 R200c, or any other ``radius_key`` the host track carries (e.g. a
    spherical-overdensity radius you added to ``extra``).
    """
    snaps, i_s, i_h = np.intersect1d(track.SnapNum, host_track.SnapNum,
                                     return_indices=True)
    if snaps.size == 0:
        raise ValueError("track and host_track share no snapshots.")
    a = np.asarray(host_track.Time, dtype=float)[i_h]          # scale factor
    t = cosmology.age(a)
    r_ref_all = host_reference(host_track, units, radius_key, radius_factor,
                               fixed_radius_kpc, cosmology)

    # comoving separation, periodic-wrapped in comoving kpc
    xs = units.comoving_kpc(track.Pos[i_s], a)
    xh = units.comoving_kpc(host_track.Pos[i_h], a)
    dx = xs - xh
    if boxsize:
        # boxsize is the *comoving* box side in the tree's length units
        box_c = float(boxsize) * units.length_kpc * units._hfac()
        dx = (dx + 0.5 * box_c) % box_c - 0.5 * box_c
    v_rel = (units.peculiar_kms(track.Vel[i_s], a)
             - units.peculiar_kms(host_track.Vel[i_h], a))

    r_phys_vec = dx * a[:, None]
    r = np.linalg.norm(r_phys_vec, axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        rhat = np.where(r[:, None] > 0, r_phys_vec / r[:, None], 0.0)
    vr = np.sum(v_rel * rhat, axis=1)
    vt = np.linalg.norm(v_rel - vr[:, None] * rhat, axis=1)
    L = np.cross(r_phys_vec, v_rel)

    r_ref = r_ref_all[i_h]

    vmax = track.extra.get("Vmax")
    gaps = int(np.sum(np.diff(snaps) > 1))
    host_ids = np.asarray(host_track.HostID)[i_h]
    return HostFrameSeries(
        snapnum=snaps, a=a, t=t, redshift=1.0 / a - 1.0, dx_c=dx, v_rel=v_rel,
        r=r, vr=vr, vt=vt, L=L, r_ref=r_ref,
        mass=units.mass(track.Mass[i_s]),
        vmax=None if vmax is None else np.asarray(vmax)[i_s],
        is_subhalo=np.asarray(track.IsSubhalo)[i_s],
        host_id_sat=np.asarray(track.HostID)[i_s],
        host_id_ref=host_ids, n_gaps=gaps,
        npart=(np.asarray(track.extra["Len"])[i_s]
               if "Len" in track.extra else None))


@dataclass
class FineOrbit:
    """Sub-snapshot-resolution orbit from Hermite interpolation."""
    t: np.ndarray
    a: np.ndarray
    r: np.ndarray
    r_ref: np.ndarray
    snap_index: np.ndarray      # interval index each fine point belongs to


def interpolate_orbit(series: HostFrameSeries, cosmology: FlatLCDM,
                      n_sub: int = 16) -> FineOrbit:
    """Interpolate ``series`` to ``n_sub`` points per snapshot interval."""
    n = len(series)
    if n < 2:
        return FineOrbit(series.t.copy(), series.a.copy(), series.r.copy(),
                         series.r_ref.copy(), np.zeros(n, dtype=int))
    t0, t1 = series.t[:-1], series.t[1:]
    h = t1 - t0
    a0, a1 = series.a[:-1], series.a[1:]
    # comoving slope = peculiar velocity / a, in kpc_comoving / Gyr
    m = series.v_rel * KMS_TO_KPC_PER_GYR / series.a[:, None]
    s = np.linspace(0.0, 1.0, n_sub, endpoint=False)
    pts = hermite(series.dx_c[:-1], series.dx_c[1:], m[:-1], m[1:], h,
                  np.broadcast_to(s, (n - 1, n_sub)))          # (n-1, ns, 3)
    t_f = t0[:, None] + s[None, :] * h[:, None]
    a_f = cosmology.a_of_t(t_f)
    r_f = np.linalg.norm(pts, axis=2) * a_f
    # reference radius: linear in t between snapshots
    rr = series.r_ref
    ref_f = rr[:-1, None] + s[None, :] * (rr[1:] - rr[:-1])[:, None]
    idx = np.repeat(np.arange(n - 1), n_sub)
    # close with the final snapshot itself
    t_all = np.concatenate([t_f.ravel(), series.t[-1:]])
    a_all = np.concatenate([a_f.ravel(), series.a[-1:]])
    r_all = np.concatenate([r_f.ravel(), series.r[-1:]])
    ref_all = np.concatenate([ref_f.ravel(), series.r_ref[-1:]])
    idx = np.concatenate([idx, [n - 2]])
    return FineOrbit(t_all, a_all, r_all, ref_all, idx)


# ---------------------------------------------------------------------
# Events
# ---------------------------------------------------------------------

CROSS_IN, CROSS_OUT, PERICENTRE, APOCENTRE = 0, 1, 2, 3
EVENT_NAMES = {CROSS_IN: "cross_in", CROSS_OUT: "cross_out",
               PERICENTRE: "pericentre", APOCENTRE: "apocentre"}


@dataclass
class Event:
    kind: int
    t: float            # cosmic time, Gyr
    a: float
    r: float            # physical kpc
    r_over_ref: float


def find_crossings(fine: FineOrbit, hysteresis: float = 0.05
                   ) -> Tuple[List[Event], bool]:
    """Crossings of ``r = r_ref`` with a hysteresis band.

    A crossing is only counted when the orbit subsequently moves a fraction
    ``hysteresis`` beyond the reference radius on the far side, so a
    satellite skimming the boundary (or noise in the host's R200c) cannot
    produce a spurious infall/backsplash pair.

    Returns ``(events, born_inside)``; ``born_inside`` means the track begins
    inside the reference radius, so the first infall is only an upper
    limit (the object had already fallen in before the tree starts
    following it).
    """
    f = fine.r / fine.r_ref
    inside = bool(f[0] < 1.0)
    born_inside = inside
    events: List[Event] = []
    last_origin = 0                      # last index on the origin side of f=1
    for i in range(1, f.size):
        if inside:
            if f[i] < 1.0:
                last_origin = i
            elif f[i] > 1.0 + hysteresis:
                events.append(_crossing(fine, f, last_origin, CROSS_OUT))
                inside, last_origin = False, i
        else:
            if f[i] > 1.0:
                last_origin = i
            elif f[i] < 1.0 - hysteresis:
                events.append(_crossing(fine, f, last_origin, CROSS_IN))
                inside, last_origin = True, i
    return events, born_inside


def _crossing(fine, f, i0, kind):
    """Linear interpolation of f = 1 between fine index i0 and i0 + 1."""
    i1 = i0 + 1
    denom = f[i1] - f[i0]
    w = 0.0 if denom == 0 else (1.0 - f[i0]) / denom
    w = min(max(w, 0.0), 1.0)
    t = fine.t[i0] + w * (fine.t[i1] - fine.t[i0])
    a = fine.a[i0] + w * (fine.a[i1] - fine.a[i0])
    r = fine.r[i0] + w * (fine.r[i1] - fine.r[i0])
    return Event(kind, float(t), float(a), float(r), 1.0)


def zigzag_extrema(y: np.ndarray, prominence: float) -> List[Tuple[int, int]]:
    """Alternating confirmed extrema of ``y`` (> 0), as
    ``(index, +1 for maximum / -1 for minimum)``.

    An extremum is confirmed only when ``y`` has subsequently moved a
    relative fraction ``prominence`` away from it. The first and last
    samples are never reported (an endpoint is not a turning point).
    """
    out: List[Tuple[int, int]] = []
    n = y.size
    if n < 3:
        return out
    direction = 0
    hi = lo = 0
    for i in range(1, n):
        if y[i] > y[hi]:
            hi = i
        if y[i] < y[lo]:
            lo = i
        if direction <= 0 and y[i] >= y[lo] * (1.0 + prominence):
            if lo > 0:
                out.append((lo, -1))
            direction, hi = 1, i
        elif direction >= 0 and y[i] <= y[hi] * (1.0 - prominence):
            if hi > 0:
                out.append((hi, +1))
            direction, lo = -1, i
    return out


def _refine_extremum(fine: FineOrbit, i: int) -> Tuple[float, float, float]:
    """Parabolic refinement of a fine-grid extremum: (t, a, r)."""
    if i <= 0 or i >= fine.t.size - 1:
        return float(fine.t[i]), float(fine.a[i]), float(fine.r[i])
    t = fine.t[i - 1:i + 2]
    r = fine.r[i - 1:i + 2]
    c = np.polyfit(t - t[1], r, 2)
    if c[0] == 0:
        return float(t[1]), float(fine.a[i]), float(r[1])
    dt = -c[1] / (2 * c[0])
    dt = float(np.clip(dt, t[0] - t[1], t[2] - t[1]))
    r_ext = float(np.polyval(c, dt))
    t_ext = float(t[1] + dt)
    a_ext = float(np.interp(t_ext, fine.t[i - 1:i + 2], fine.a[i - 1:i + 2]))
    return t_ext, a_ext, r_ext


def find_turning_points(fine: FineOrbit, prominence: float = 0.05
                        ) -> List[Event]:
    """Pericentres and apocentres of the interpolated ``r(t)``."""
    events = []
    for i, sign in zigzag_extrema(fine.r, prominence):
        t, a, r = _refine_extremum(fine, i)
        ref = float(np.interp(t, fine.t, fine.r_ref))
        events.append(Event(APOCENTRE if sign > 0 else PERICENTRE,
                            t, a, r, r / ref))
    return events
