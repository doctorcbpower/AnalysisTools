#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
analysistools.orbits.tracker
----------------------------
``OrbitTracker``: turn merger-tree tracks into host-frame orbits, events and
per-satellite orbital scalars; optionally continue orphans in the host
potential; write/read a standalone HDF5 sidecar (``OrbitCatalogue``) that the
master-catalogue pipeline consumes.

Definitions (each is a modelling choice, recorded in the sidecar's
attributes):

``infall``           first inward crossing of the reference radius
                     ``radius_factor * host[radius_key]`` (default 1.0 x
                     GroupR200 == R200c), found on the Hermite-interpolated
                     orbit, with a hysteresis band (``hysteresis``) so
                     boundary-skimming cannot count as two events.
``pericentre/apo``   confirmed local minima/maxima of the interpolated
                     physical separation after first infall (an extremum
                     counts only once r has moved a fraction ``prominence``
                     away from it); the first/last samples never count.
``eccentricity``     (r_apo - r_peri) / (r_apo + r_peri) from the two most
                     recent consecutive confirmed turning points -- NaN until
                     the orbit has completed half a radial period.
``*_EL``             turning points / radial period from the infall (E, L) in
                     the *frozen-at-infall* spherical host potential; only if
                     a potential was supplied. This is a different estimator
                     from the measured events and is reported separately.
``is_backsplash``    has crossed inward at least once and is now outside.
``preprocessed``     was a subhalo of some *other* FOF group before first
                     infall into the reference host.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field, asdict
from typing import Dict, List, Optional, Sequence

import h5py
import numpy as np

from .cosmology import FlatLCDM
from .kinematics import (R200C_FROM_M200, CROSS_IN, CROSS_OUT, PERICENTRE, APOCENTRE,
                         EVENT_NAMES, Event, FineOrbit, HostFrameSeries,
                         TrackUnits, find_crossings, find_turning_points,
                         hermite, host_frame_series, host_reference,
                         interpolate_orbit)
from .orphans import integrate_orphans
from .potential import NFWPotential, Potential, TimeDependentPotential

logger = logging.getLogger(__name__)

SIDECAR_VERSION = "0.1"

# name -> (dtype, fill, units, description)
SCALARS: Dict[str, tuple] = {
    "born_inside": ("bool", False, "dimensionless",
                    "Track begins inside the reference radius; infall quantities are then upper limits in time."),
    "t_infall": ("float64", np.nan, "Gyr", "Cosmic time of first inward crossing of the reference radius."),
    "z_infall": ("float64", np.nan, "dimensionless", "Redshift of first inward crossing."),
    "snap_infall": ("int32", -1, "dimensionless", "First snapshot at/after the first inward crossing."),
    "r_ref_infall": ("float64", np.nan, "kpc", "Physical reference radius at first crossing."),
    "mass_infall": ("float64", np.nan, "Msun", "Subhalo mass at the last snapshot before first crossing."),
    "vmax_infall": ("float64", np.nan, "km/s", "Subhalo Vmax at the last snapshot before first crossing."),
    "vr_infall": ("float64", np.nan, "km/s", "Radial velocity relative to host at first crossing."),
    "vt_infall": ("float64", np.nan, "km/s", "Tangential velocity relative to host at first crossing."),
    "L_infall": ("float64", np.nan, "kpc km/s", "Specific angular momentum vector at first crossing."),
    "t_fof_infall": ("float64", np.nan, "Gyr", "Cosmic time of the first snapshot with IsSubhalo (FOF membership)."),
    "snap_fof_infall": ("int32", -1, "dimensionless", "First snapshot with IsSubhalo."),
    "preprocessed": ("bool", False, "dimensionless", "Subhalo of a different FOF group before first infall."),
    "n_infalls": ("int16", 0, "dimensionless", "Inward crossings of the reference radius."),
    "n_outfalls": ("int16", 0, "dimensionless", "Outward crossings of the reference radius."),
    "is_backsplash": ("bool", False, "dimensionless", "Crossed inward at least once and is now outside."),
    "inside_now": ("bool", False, "dimensionless", "Inside the reference radius (default R200c) at the final snapshot: the Dorcha satellite definition."),
    "n_pericentres": ("int16", 0, "dimensionless", "Confirmed pericentres after first infall."),
    "t_peri_first": ("float64", np.nan, "Gyr", "Time of first pericentre after infall."),
    "r_peri_first": ("float64", np.nan, "kpc", "First pericentre after infall."),
    "r_peri_min": ("float64", np.nan, "kpc", "Smallest pericentre after infall."),
    "t_peri_last": ("float64", np.nan, "Gyr", "Time of most recent pericentre."),
    "r_peri_last": ("float64", np.nan, "kpc", "Most recent pericentre."),
    "r_apo_first": ("float64", np.nan, "kpc", "First apocentre after infall."),
    "r_apo_last": ("float64", np.nan, "kpc", "Most recent apocentre."),
    "ecc_last": ("float64", np.nan, "dimensionless", "(apo-peri)/(apo+peri) from the two latest consecutive turning points."),
    "period_peri": ("float64", np.nan, "Gyr", "Mean interval between successive pericentres (needs >= 2)."),
    "E_infall": ("float64", np.nan, "(km/s)^2", "Specific energy at first crossing in the host potential (needs potential)."),
    "r_peri_EL": ("float64", np.nan, "kpc", "Pericentre from infall (E, L) in the frozen host potential."),
    "r_apo_EL": ("float64", np.nan, "kpc", "Apocentre from infall (E, L) in the frozen host potential."),
    "period_EL": ("float64", np.nan, "Gyr", "Radial period from infall (E, L) in the frozen host potential."),
    "r_final": ("float64", np.nan, "kpc", "Physical host-centric distance at the last snapshot (orphan-extended)."),
    "orphan_continued": ("bool", False, "dimensionless", "Position after the last resolved snapshot comes from orbit integration (only if inside orphan_max_r_over_ref x R_ref when lost)."),
    "mass_final": ("float64", np.nan, "Msun", "Bound subhalo mass at the last resolved snapshot."),
    "mass_peak": ("float64", np.nan, "Msun", "Peak bound subhalo mass over the resolved history."),
    "npart_final": ("int32", -1, "dimensionless", "Bound particle count at the last resolved snapshot (-1: not in the tree)."),
    "npart_peak": ("int32", -1, "dimensionless", "Peak bound particle count over the resolved history (-1: not in the tree)."),
    "in_count_sample": ("bool", False, "dimensionless", "Passes the resolution cuts on the counts basis (default: peak). Combine with inside_now for the satellite counts sample."),
    "in_orbit_sample": ("bool", False, "dimensionless", "Passes the resolution cuts on the orbits basis (default: final). Combine with inside_now for the satellite orbit sample."),
    "last_resolved_snap": ("int32", -1, "dimensionless", "Last snapshot at which the tree resolves the subhalo."),
    "is_orphan_final": ("bool", False, "dimensionless", "Tree track ends before the host's last snapshot (see orphan_continued for whether its position was continued)."),
    "n_snapshot_gaps": ("int16", 0, "dimensionless", "Snapshots the track skips (interpolated across)."),
}


@dataclass
class OrbitConfig:
    radius_key: str = "GroupR200"
    radius_factor: float = 1.0
    fixed_radius_kpc: Optional[float] = None
    hysteresis: float = 0.05
    prominence: float = 0.05
    n_sub: int = 16
    # Resolution cuts defining the two samples (all optional, AND-ed; the
    # mass cut is max(min_mass, min_mass_fraction * M200c_host(z=0))).
    # "counts" sample (abundances, luminosity functions, occupation) uses
    # ``counts_basis`` -- peak by default, the usual completeness criterion;
    # "orbits" sample (pericentres, eccentricities, positions) uses
    # ``orbits_basis`` -- final by default, since an orbit measurement needs
    # the object to be resolved *now*, not once.
    min_mass: Optional[float] = None              # Msun
    min_mass_fraction: Optional[float] = None     # of host M200c at z=0
    min_particles: Optional[int] = None
    counts_basis: str = "peak"
    orbits_basis: str = "final"
    orphans: bool = True
    # Continue a lost subhalo only if it was within this multiple of the
    # reference radius at its last resolved snapshot. A bound spherical
    # potential with no Hubble flow is the wrong model for fringe objects
    # (FOF membership extends well beyond R200c): integrating them sent
    # 68 per cent of Dorcha's lost fringe subhaloes beyond R200c, out to
    # ~1 Mpc. Set to inf to continue everything.
    orphan_max_r_over_ref: float = 1.0
    dynamical_friction: bool = False
    orphan_dt_max: float = 0.005       # Gyr


@dataclass
class OrbitRecord:
    """One satellite's host-frame orbit. ``series`` arrays are on the
    snapshot grid (resolved rows then orphan-integrated rows,
    ``extrapolated`` marks the latter)."""
    halo_id: int
    series: Dict[str, np.ndarray]
    events: List[Event]
    scalars: Dict[str, object]


def tree_nfw_potential(host_track, cosmology: FlatLCDM, units: TrackUnits,
                       radius_key: str = "GroupR200",
                       check_vmax: bool = True,
                       vmax_tolerance: float = 0.5,
                       normalise: str = "m200") -> TimeDependentPotential:
    """Time-dependent NFW host potential straight from the host's tree
    track: M200 and R200 from ``extra['GroupM200']``/``extra[radius_key]``,
    scale radius from ``extra['VmaxRad']`` (Rmax = 2.163 rs).

    Snapshots where that is not NFW-like (the NFW's Vmax at Rmax differs
    from the tree's by more than ``vmax_tolerance``, or the halo is not yet
    resolved) are skipped with a warning rather than patched with an
    assumed concentration; before the first accepted snapshot the
    potential is held fixed at that first epoch.

    ``normalise="m200"`` (default) fixes the NFW mass by (M200, R200) and
    the scale radius by Rmax; ``"vmax"`` instead fixes it by the tree's
    (Vmax, Rmax) -- the bound inner subhalo -- and lets M200 follow.
    Neither is "right": in a truncation test on the real Dorcha tree
    (satellites the tree follows to z=0, cut at z~0.3-0.8, continued as
    orphans, compared with the tree's own z=0 radius) the m200 version puts
    ~7 points too many orphans inside R200c (median -0.16 dex in r), the
    vmax version ~16 points too few (median -0.03 dex); both scatter by
    ~0.4 dex. See docs/orbits.md.

    Known systematic (checked on the Dorcha host, z < 2.7): (M200, Rmax)
    implies a Vmax 12-30 per cent *above* the tree's bound-subhalo Vmax
    (M200 includes unbound/diffuse mass), i.e. ~25-70 per cent more inner
    mass than the bound subhalo carries. The default tolerance (0.5) accepts
    that; a potential measured from snapshot particles
    (``TabulatedSpherical.from_particles``) avoids having to choose.
    """
    need = ["GroupM200", "VmaxRad"] + ([] if radius_key == R200C_FROM_M200
                                        else [radius_key])
    for key in need:
        if key not in host_track.extra:
            raise KeyError(f"host track lacks extra['{key}'] needed for an "
                           f"NFW potential (have {sorted(host_track.extra)}).")
    a = np.asarray(host_track.Time, dtype=float)
    # M200 and R200 must be mutually consistent for the NFW normalisation,
    # so R200 is *not* scaled by any radius_factor here
    r200 = host_reference(host_track, units, radius_key, 1.0, None, cosmology)
    rmax = units.physical_kpc(np.asarray(host_track.extra["VmaxRad"]), a)
    m200 = units.mass(host_track.extra["GroupM200"])
    vmax = host_track.extra.get("Vmax")
    if normalise not in ("m200", "vmax"):
        raise ValueError("normalise must be 'm200' or 'vmax'.")
    times, pots, skipped, ratios = [], [], 0, []
    for i in range(a.size):
        if not (m200[i] > 0 and r200[i] > 0 and rmax[i] > 0):
            skipped += 1
            continue
        try:
            p = NFWPotential.from_vmax_rmax(
                m200[i], r200[i],
                vmax=(float(vmax[i]) if (check_vmax and vmax is not None)
                      else None),
                rmax=rmax[i], vmax_tolerance=vmax_tolerance)
        except ValueError:
            skipped += 1
            continue
        if vmax is not None and vmax[i] > 0:
            ratios.append(float(p.vcirc(rmax[i])) / float(vmax[i]))
            if normalise == "vmax":
                # rescale to the bound inner subhalo: M(<Rmax) = Vmax^2 Rmax / G
                from .cosmology import G_KPC_KMS2_MSUN as _G
                m_rmax = float(vmax[i]) ** 2 * rmax[i] / _G
                scale = m_rmax / float(p.menc(rmax[i]))
                p = NFWPotential(p.m200 * scale, p.r200, p.rs)
        times.append(float(cosmology.age(a[i])))
        pots.append(p)
    if not pots:
        raise ValueError("no host snapshot yields an NFW-like potential.")
    if skipped:
        logger.warning("tree_nfw_potential: %d/%d host snapshots skipped "
                       "(unresolved or Vmax mismatch > %.0f%%); potential "
                       "held fixed before the first accepted epoch (z=%.2f).",
                       skipped, a.size, 100 * vmax_tolerance,
                       1.0 / cosmology.a_of_t(times[0]) - 1.0)
    if ratios:
        logger.info("tree_nfw_potential: median NFW-implied/tree Vmax = %.2f "
                    "over %d accepted snapshots.", np.median(ratios),
                    len(ratios))
    return TimeDependentPotential(times, pots,
                                  r_grid=np.geomspace(0.5, 3.0e3, 160))


class OrbitTracker:
    """Compute orbits for many satellites relative to one reference host.

    Parameters
    ----------
    cosmology : FlatLCDM
    units : TrackUnits
        How the *tracks* are expressed. Prefer tracks loaded with the
        tree reader's defaults and say so here, rather than converting
        twice.
    config : OrbitConfig
    potential : Potential, optional
        Host potential (a ``TimeDependentPotential`` for cosmological use).
        Required for orphan continuation and the ``*_EL``/``E_infall``
        quantities. ``tree_nfw_potential(host_track, ...)`` builds one from
        the tree alone; ``TabulatedSpherical.from_particles`` from snapshot
        particles.
    boxsize : float, optional
        Comoving box side in the tracks' length units (periodic wrapping).
    """

    def __init__(self, cosmology: FlatLCDM, units: TrackUnits,
                 config: Optional[OrbitConfig] = None,
                 potential: Optional[Potential] = None,
                 boxsize: Optional[float] = None):
        self.cosmology = cosmology
        self.units = units
        self.config = config or OrbitConfig()
        self.potential = potential
        self.boxsize = boxsize

    @classmethod
    def from_epoch(cls, epoch, host_track, length_kpc: float, mass_msun: float,
                   potential="tree_nfw", velocity: str = "auto",
                   cosmology=None, potential_options=None,
                   **config) -> "OrbitTracker":
        """Build a tracker for ``epoch``'s merger tree.

        ``length_kpc``/``mass_msun`` -- the size of one length/mass unit of
        the tree's tracks -- are required and never guessed; check them
        against your tree file. Comoving/little-h state is read from the
        tree reader (``epoch.tree.backend``). ``cosmology`` may be a
        ``FlatLCDM``, ``{"h":.., "omega_m":..}``, or None to read it from a
        SubFind-HBT tree's ``Parameters``. ``potential`` is ``"tree_nfw"``
        (NFW from the host's tree track), a ``Potential`` instance, or None
        (no orphan continuation / E, L quantities); ``potential_options``
        is passed to ``tree_nfw_potential`` (e.g. ``{"normalise": "vmax"}``).
        Remaining kwargs are
        ``OrbitConfig`` fields.
        """
        backend = epoch.tree.backend
        h = backend._get_hubble_param()
        if isinstance(cosmology, dict):
            cosmo = FlatLCDM(cosmology["h"], cosmology["omega_m"],
                             cosmology.get("omega_lambda"))
        elif cosmology is None:
            cosmo = FlatLCDM.from_tree(backend)
        else:
            cosmo = cosmology
        units = TrackUnits(comoving=bool(backend.comoving),
                           little_h=bool(backend.little_h),
                           h=float(h) if h else cosmo.h,
                           length_kpc=float(length_kpc),
                           mass_msun=float(mass_msun), velocity=velocity)
        if velocity == "auto":
            from ..merger_tree_tools import TREE_NATIVE_IS_COMOVING
            native_c = backend.native_is_comoving
            if native_c is None:
                native_c = TREE_NATIVE_IS_COMOVING.get(
                    backend.treefileformat, True)
            # the reader multiplies velocities by 1/sqrt(a) only when it
            # converts a *physical* native tree to comoving, and by sqrt(a)
            # only in the opposite case; a native-comoving tree viewed
            # comoving is left untouched (verified physical-peculiar for the
            # Dorcha SubFind-HBT tree)
            velocity = {(True, True): "physical", (False, True): "sqrt_a",
                        (False, False): "physical",
                        (True, False): "inv_sqrt_a"}[
                (bool(native_c), bool(backend.comoving))]
            units.velocity = velocity
        if "radius_key" not in config and "fixed_radius_kpc" not in config \
                and "GroupR200" not in host_track.extra:
            config["radius_key"] = R200C_FROM_M200
            logger.info("tree has no GroupR200: using R200c derived from "
                        "GroupM200 and the cosmology.")
        cfg = OrbitConfig(**config)
        if isinstance(potential, str):
            if potential != "tree_nfw":
                raise ValueError("potential must be 'tree_nfw', a Potential "
                                 "instance or None.")
            potential = tree_nfw_potential(host_track, cosmo, units,
                                           radius_key=cfg.radius_key,
                                           **(potential_options or {}))
        return cls(cosmo, units, cfg, potential, boxsize=epoch.boxsize)

    # -----------------------------------------------------------------
    def run(self, tracks: Sequence, host_track,
            ids: Optional[Sequence[int]] = None) -> "OrbitCatalogue":
        """Orbit records for ``tracks`` (``HaloTrack`` or ``None`` for an
        unresolved satellite -- those get a fill-valued record)."""
        cfg = self.config
        ids = list(range(len(tracks))) if ids is None else list(ids)
        host_last_snap = int(np.max(host_track.SnapNum))
        self._set_resolution(host_track)
        h_t = self.cosmology.age(np.asarray(host_track.Time, dtype=float))
        h_ref = host_reference(host_track, self.units, cfg.radius_key,
                               cfg.radius_factor, cfg.fixed_radius_kpc,
                               self.cosmology)

        parts: List[Optional[dict]] = []
        for tr in tracks:
            if tr is None or len(tr) == 0:
                parts.append(None)
                continue
            try:
                ser = host_frame_series(
                    tr, host_track, self.cosmology, self.units, self.boxsize,
                    cfg.radius_key, cfg.radius_factor, cfg.fixed_radius_kpc)
            except ValueError:
                parts.append(None)
                continue
            parts.append({"series": ser,
                          "fine": interpolate_orbit(ser, self.cosmology,
                                                    cfg.n_sub),
                          "orphan": None})

        if cfg.orphans and self.potential is not None:
            self._continue_orphans(parts, host_track, host_last_snap,
                                   h_t, h_ref)
        elif cfg.orphans and self.potential is None:
            n_orph = sum(1 for p in parts if p is not None
                         and p["series"].snapnum[-1] < host_last_snap)
            if n_orph:
                logger.warning(
                    "OrbitTracker: %d satellites end before the host's last "
                    "snapshot and no potential was given, so their orbits "
                    "stop at the last resolved snapshot (no orphan "
                    "positions). Pass potential= to continue them.", n_orph)

        records = [self._finalise(hid, p, host_last_snap)
                   for hid, p in zip(ids, parts)]
        return OrbitCatalogue(records, self._metadata())

    # -----------------------------------------------------------------
    def _set_resolution(self, host_track):
        cfg = self.config
        for b in (cfg.counts_basis, cfg.orbits_basis):
            if b not in ("peak", "final"):
                raise ValueError("counts_basis/orbits_basis must be 'peak' "
                                 "or 'final'.")
        self._host_m200 = None
        m200 = host_track.extra.get("GroupM200")
        if m200 is not None:
            self._host_m200 = float(self.units.mass(m200[-1]))
        mass_cuts = []
        if cfg.min_mass is not None:
            mass_cuts.append(float(cfg.min_mass))
        if cfg.min_mass_fraction is not None:
            if self._host_m200 is None:
                raise ValueError("min_mass_fraction needs the host track's "
                                 "extra['GroupM200'].")
            mass_cuts.append(float(cfg.min_mass_fraction) * self._host_m200)
        self._min_mass_eff = max(mass_cuts) if mass_cuts else None

    def _passes(self, sc, basis):
        cfg = self.config
        ok = True
        if self._min_mass_eff is not None:
            ok &= sc[f"mass_{basis}"] >= self._min_mass_eff
        if cfg.min_particles is not None:
            ok &= sc[f"npart_{basis}"] >= int(cfg.min_particles)
        return bool(ok)

    # -----------------------------------------------------------------
    def _continue_orphans(self, parts, host_track, host_last_snap, h_t, h_ref):
        cfg = self.config
        host_snaps = np.asarray(host_track.SnapNum)
        groups: Dict[int, List[int]] = {}
        for i, p in enumerate(parts):
            if p is None:
                continue
            last = int(p["series"].snapnum[-1])
            ser = p["series"]
            if last < host_last_snap and (
                    ser.r[-1] <= cfg.orphan_max_r_over_ref * ser.r_ref[-1]):
                groups.setdefault(last, []).append(i)
        for last, members in groups.items():
            later = np.flatnonzero(host_snaps >= last)
            later = later[np.argsort(host_snaps[later])]
            ts = h_t[later]
            # dense grid: n_sub points per interval, plus the final snapshot
            s = np.linspace(0.0, 1.0, cfg.n_sub, endpoint=False)
            t_dense = np.concatenate(
                [ts[j] + s * (ts[j + 1] - ts[j]) for j in range(ts.size - 1)]
                + [ts[-1:]])
            x0 = np.array([parts[i]["series"].dx_c[-1]
                           * parts[i]["series"].a[-1] for i in members])
            v0 = np.array([parts[i]["series"].v_rel[-1] for i in members])
            m_sat = np.array([parts[i]["series"].mass[-1] for i in members])
            xs, vs = integrate_orphans(
                self.potential, x0, v0, t_dense, m_sat=m_sat,
                dynamical_friction=cfg.dynamical_friction,
                dt_max=cfg.orphan_dt_max)
            a_dense = self.cosmology.a_of_t(t_dense)
            ref_dense = np.interp(t_dense, h_t[later], h_ref[later])
            for k, i in enumerate(members):
                parts[i]["orphan"] = {
                    "t": t_dense, "a": a_dense, "x": xs[k], "v": vs[k],
                    "r_ref": ref_dense, "snaps": host_snaps[later],
                    "n_sub": cfg.n_sub}

    # -----------------------------------------------------------------
    def _finalise(self, halo_id, part, host_last_snap) -> OrbitRecord:
        scal = {k: v[1] for k, v in SCALARS.items()}
        if part is None:
            return OrbitRecord(int(halo_id), {}, [], scal)
        cfg = self.config
        ser: HostFrameSeries = part["series"]
        fine: FineOrbit = part["fine"]
        orph = part["orphan"]

        # ---- assemble snapshot-grid series (+ orphan rows) -------------
        S = {"snapnum": ser.snapnum.copy(), "t": ser.t.copy(),
             "a": ser.a.copy(), "r": ser.r.copy(), "vr": ser.vr.copy(),
             "vt": ser.vt.copy(), "L": ser.L.copy(),
             "r_ref": ser.r_ref.copy(), "mass": ser.mass.copy(),
             "extrapolated": np.zeros(len(ser), dtype=bool)}
        f_t, f_a, f_r, f_ref = fine.t, fine.a, fine.r, fine.r_ref
        if orph is not None:
            x, v = orph["x"], orph["v"]
            r_o = np.linalg.norm(x, axis=1)
            rhat = x / np.maximum(r_o, 1e-12)[:, None]
            vr_o = np.sum(v * rhat, axis=1)
            vt_o = np.linalg.norm(v - vr_o[:, None] * rhat, axis=1)
            L_o = np.cross(x, v)
            step = orph["n_sub"]
            rows = np.concatenate([np.arange(step, x.shape[0], step),
                                   [x.shape[0] - 1]]).astype(int)
            rows = np.unique(rows)
            n_old = len(ser)
            snaps_after = orph["snaps"][1:]            # drops the shared start
            if rows.size == snaps_after.size:
                S["snapnum"] = np.concatenate([S["snapnum"], snaps_after])
                S["t"] = np.concatenate([S["t"], orph["t"][rows]])
                S["a"] = np.concatenate([S["a"], orph["a"][rows]])
                S["r"] = np.concatenate([S["r"], r_o[rows]])
                S["vr"] = np.concatenate([S["vr"], vr_o[rows]])
                S["vt"] = np.concatenate([S["vt"], vt_o[rows]])
                S["L"] = np.concatenate([S["L"], L_o[rows]])
                S["r_ref"] = np.concatenate([S["r_ref"], orph["r_ref"][rows]])
                S["mass"] = np.concatenate(
                    [S["mass"], np.full(rows.size, ser.mass[-1])])
                S["extrapolated"] = np.concatenate(
                    [S["extrapolated"], np.ones(rows.size, dtype=bool)])
            # dense orbit for event detection: drop the shared first point
            f_t = np.concatenate([f_t, orph["t"][1:]])
            f_a = np.concatenate([f_a, orph["a"][1:]])
            f_r = np.concatenate([f_r, r_o[1:]])
            f_ref = np.concatenate([f_ref, orph["r_ref"][1:]])
        dense = FineOrbit(f_t, f_a, f_r, f_ref,
                          np.zeros(f_t.size, dtype=int))

        # ---- events ----------------------------------------------------
        crossings, born_inside = find_crossings(dense, cfg.hysteresis)
        turning = find_turning_points(dense, cfg.prominence)
        events = sorted(crossings + turning, key=lambda e: e.t)
        self._fill_scalars(scal, ser, S, events, born_inside,
                           orphan=orph is not None, host_last=host_last_snap)
        scal["in_count_sample"] = self._passes(scal, cfg.counts_basis)
        scal["in_orbit_sample"] = self._passes(scal, cfg.orbits_basis)
        return OrbitRecord(int(halo_id), S, events, scal)

    # -----------------------------------------------------------------
    def _state_at(self, ser: HostFrameSeries, t: float):
        """Host-frame (r_vec_phys, v_rel) at cosmic time t on the resolved
        series: Hermite position, linear velocity."""
        j = int(np.clip(np.searchsorted(ser.t, t) - 1, 0, len(ser) - 2)) \
            if len(ser) > 1 else 0
        if len(ser) == 1 or t >= ser.t[-1]:
            j = max(len(ser) - 2, 0)
        if len(ser) == 1:
            return ser.dx_c[0] * ser.a[0], ser.v_rel[0]
        t0, t1 = ser.t[j], ser.t[j + 1]
        s = float(np.clip((t - t0) / (t1 - t0), 0.0, 1.0))
        m = ser.v_rel * 1.0227121650537077 / ser.a[:, None]
        pos = hermite(ser.dx_c[j][None], ser.dx_c[j + 1][None], m[j][None],
                      m[j + 1][None], np.array([t1 - t0]),
                      np.array([[s]]))[0, 0]
        a = float(self.cosmology.a_of_t(t))
        v = (1 - s) * ser.v_rel[j] + s * ser.v_rel[j + 1]
        return pos * a, v

    def _fill_scalars(self, sc, ser, S, events, born_inside, orphan,
                      host_last):
        sc["born_inside"] = bool(born_inside)
        sc["n_snapshot_gaps"] = ser.n_gaps
        sc["mass_final"] = float(ser.mass[-1])
        sc["mass_peak"] = float(np.max(ser.mass))
        if ser.npart is not None:
            sc["npart_final"] = int(ser.npart[-1])
            sc["npart_peak"] = int(np.max(ser.npart))
        sc["last_resolved_snap"] = int(ser.snapnum[-1])
        sc["is_orphan_final"] = bool(ser.snapnum[-1] < host_last)
        sc["orphan_continued"] = bool(orphan)
        sc["r_final"] = float(S["r"][-1])
        ins = [e for e in events if e.kind == CROSS_IN]
        outs = [e for e in events if e.kind == CROSS_OUT]
        sc["n_infalls"] = len(ins) + (1 if born_inside else 0)
        sc["n_outfalls"] = len(outs)
        sc["inside_now"] = bool(S["r"][-1] < S["r_ref"][-1])
        sc["is_backsplash"] = bool(sc["n_infalls"] > 0 and not sc["inside_now"])

        # FOF-based infall, kept alongside as a second definition
        fof = np.flatnonzero(ser.is_subhalo)
        if fof.size:
            sc["snap_fof_infall"] = int(ser.snapnum[fof[0]])
            sc["t_fof_infall"] = float(ser.t[fof[0]])

        # first infall (R_ref crossing)
        t_inf = None
        if born_inside:
            t_inf = float(S["t"][0])
            j_in = 0
            sc["snap_infall"] = int(S["snapnum"][0])
            sc["t_infall"], sc["z_infall"] = t_inf, float(1 / S["a"][0] - 1)
            sc["r_ref_infall"] = float(S["r_ref"][0])
        elif ins:
            e0 = ins[0]
            t_inf = e0.t
            sc["t_infall"], sc["z_infall"] = e0.t, float(1 / e0.a - 1)
            j_in = int(np.searchsorted(S["t"], e0.t))
            j_in = min(j_in, len(S["t"]) - 1)
            sc["snap_infall"] = int(S["snapnum"][j_in])
            sc["r_ref_infall"] = float(np.interp(e0.t, S["t"], S["r_ref"]))
            jb = max(j_in - 1, 0)
            sc["mass_infall"] = float(S["mass"][jb])
            if ser.vmax is not None and jb < len(ser):
                sc["vmax_infall"] = float(ser.vmax[jb])
        if born_inside:
            sc["mass_infall"] = float(S["mass"][0])
            if ser.vmax is not None:
                sc["vmax_infall"] = float(ser.vmax[0])

        if t_inf is not None and t_inf <= ser.t[-1]:
            rvec, vrel = self._state_at(ser, t_inf)
            r = float(np.linalg.norm(rvec))
            if r > 0:
                vr = float(np.dot(vrel, rvec) / r)
                vt = float(np.linalg.norm(vrel - vr * rvec / r))
                L = np.cross(rvec, vrel)
                sc["vr_infall"], sc["vt_infall"] = vr, vt
                sc["L_infall"] = L
                if self.potential is not None:
                    E = float(self.potential.specific_energy(
                        r, np.linalg.norm(vrel), t_inf))
                    sc["E_infall"] = E
                    rp, ra = self.potential.turning_points(
                        E, np.linalg.norm(L), r, t_inf)
                    sc["r_peri_EL"], sc["r_apo_EL"] = rp, ra
                    sc["period_EL"] = self.potential.radial_period(
                        rp, ra, E, np.linalg.norm(L), t_inf)
            # pre-processing: in some other FOF group before infall
            before = ser.t < t_inf
            sc["preprocessed"] = bool(np.any(
                before & ser.is_subhalo
                & (ser.host_id_sat != ser.host_id_ref)))

        # turning points after infall
        t0 = -np.inf if t_inf is None else t_inf
        peri = [e for e in events if e.kind == PERICENTRE and e.t >= t0]
        apo = [e for e in events if e.kind == APOCENTRE and e.t >= t0]
        sc["n_pericentres"] = len(peri)
        if peri:
            sc["t_peri_first"], sc["r_peri_first"] = peri[0].t, peri[0].r
            sc["r_peri_min"] = min(e.r for e in peri)
            sc["t_peri_last"], sc["r_peri_last"] = peri[-1].t, peri[-1].r
        if apo:
            sc["r_apo_first"], sc["r_apo_last"] = apo[0].r, apo[-1].r
        if len(peri) >= 2:
            sc["period_peri"] = float(np.mean(
                np.diff([e.t for e in peri])))
        tp = [e for e in events if e.kind in (PERICENTRE, APOCENTRE)
              and e.t >= t0]
        if len(tp) >= 2:
            ra = max(tp[-1].r, tp[-2].r)
            rp = min(tp[-1].r, tp[-2].r)
            sc["ecc_last"] = float((ra - rp) / (ra + rp))

    # -----------------------------------------------------------------
    def _metadata(self) -> Dict[str, object]:
        md = {"sidecar_version": SIDECAR_VERSION,
              "cosmology": repr(self.cosmology),
              "cosmo_h": self.cosmology.h, "cosmo_Om": self.cosmology.omega_m,
              "cosmo_OL": self.cosmology.omega_lambda,
              "potential": repr(self.potential) if self.potential is not None
              else "none"}
        md["host_m200_final_msun"] = (-1.0 if getattr(self, "_host_m200", None)
                                      is None else self._host_m200)
        md["min_mass_effective_msun"] = (-1.0 if getattr(self, "_min_mass_eff",
                                                         None) is None
                                         else self._min_mass_eff)
        md.update({f"config_{k}": ("none" if v is None else v)
                   for k, v in asdict(self.config).items()})
        md.update({f"units_{k}": v for k, v in asdict(self.units).items()})
        return md


# ---------------------------------------------------------------------
# Sidecar I/O
# ---------------------------------------------------------------------

class OrbitCatalogue:
    """A set of ``OrbitRecord`` plus provenance, with HDF5 round-tripping.

    ``scalar(name)`` returns an array over satellites (in the order given to
    ``OrbitTracker.run``); ``record(i)`` returns the full ``OrbitRecord``.
    """

    def __init__(self, records: List[OrbitRecord], metadata: Dict[str, object]):
        self.records = records
        self.metadata = dict(metadata)

    def __len__(self):
        return len(self.records)

    def record(self, i: int) -> OrbitRecord:
        return self.records[i]

    @property
    def ids(self) -> np.ndarray:
        return np.array([r.halo_id for r in self.records], dtype=np.int64)

    def satellites(self, sample: str = "counts") -> np.ndarray:
        """Boolean mask of the configured satellite sample: inside the
        reference radius at the final snapshot AND passing the resolution
        cuts on that sample's basis (``OrbitConfig``: counts -> peak, orbits
        -> final by default). ``sample`` is ``"counts"`` or ``"orbits"``."""
        if sample not in ("counts", "orbits"):
            raise ValueError("sample must be 'counts' or 'orbits'.")
        flag = self.scalar("in_count_sample" if sample == "counts"
                           else "in_orbit_sample")
        return self.scalar("inside_now") & flag

    def select(self, min_mass: Optional[float] = None,
               min_particles: Optional[int] = None, basis: str = "peak",
               inside_only: bool = False,
               min_mass_fraction: Optional[float] = None) -> np.ndarray:
        """Boolean mask over satellites passing user-defined resolution
        thresholds (all given criteria must hold).

        ``min_mass`` in Msun, ``min_particles`` in bound particles;
        ``basis`` is ``"peak"`` (the object's best-resolved moment, the usual
        completeness criterion) or ``"final"`` (its last resolved snapshot).
        ``min_mass_fraction`` is a fraction of the host's z=0 M200c
        (recorded in the sidecar); the mass cut is the larger of it and
        ``min_mass``. ``inside_only`` additionally requires ``inside_now`` (inside the
        reference radius at the final snapshot). Unresolved satellites
        (no record) never pass. Raises if ``min_particles`` is requested but
        the tree carried no particle counts.
        """
        if basis not in ("peak", "final"):
            raise ValueError("basis must be 'peak' or 'final'.")
        mask = np.array([bool(r.series) for r in self.records], dtype=bool)
        if min_mass_fraction is not None:
            m200 = float(self.metadata.get("host_m200_final_msun", -1.0))
            if m200 <= 0:
                raise ValueError("this catalogue recorded no host M200.")
            min_mass = max(float(min_mass or 0.0), float(min_mass_fraction) * m200)
        if min_mass is not None:
            mask &= self.scalar(f"mass_{basis}") >= float(min_mass)
        if min_particles is not None:
            npart = self.scalar(f"npart_{basis}")
            if np.all(npart < 0):
                raise ValueError("this catalogue has no particle counts "
                                 "(tree lacked SubhaloLen).")
            mask &= npart >= int(min_particles)
        if inside_only:
            mask &= self.scalar("inside_now")
        return mask

    def scalar(self, name: str) -> np.ndarray:
        dtype, fill, _, _ = SCALARS[name]
        if name == "L_infall":
            out = np.full((len(self), 3), np.nan)
            for i, r in enumerate(self.records):
                out[i] = np.asarray(r.scalars[name], dtype=float) \
                    if np.ndim(r.scalars[name]) else np.nan
            return out
        return np.array([r.scalars[name] for r in self.records], dtype=dtype)

    # -- HDF5 --------------------------------------------------------
    def write_hdf5(self, path: str) -> None:
        with h5py.File(path, "w") as f:
            for k, v in self.metadata.items():
                f.attrs[k] = v
            f.attrs["n_satellites"] = len(self)
            g = f.create_group("Satellites")
            g.create_dataset("HaloID", data=self.ids)
            for name, (dtype, _, units, desc) in SCALARS.items():
                ds = g.create_dataset(name, data=self.scalar(name))
                ds.attrs["units"], ds.attrs["description"] = units, desc
            # events (flat)
            ev = [(i, e) for i, r in enumerate(self.records) for e in r.events]
            ge = f.create_group("Events")
            ge.create_dataset("sat_index", data=np.array(
                [i for i, _ in ev], dtype=np.int64))
            ge.create_dataset("kind", data=np.array(
                [e.kind for _, e in ev], dtype=np.int8))
            ge.attrs["kind_names"] = ", ".join(
                f"{k}={v}" for k, v in EVENT_NAMES.items())
            for fld, unit in (("t", "Gyr"), ("a", ""), ("r", "kpc"),
                              ("r_over_ref", "")):
                d = ge.create_dataset(fld, data=np.array(
                    [getattr(e, fld) for _, e in ev], dtype=float))
                d.attrs["units"] = unit
            # series (flat + offsets)
            gs = f.create_group("Series")
            lens = np.array([len(r.series.get("snapnum", [])) for r in
                             self.records], dtype=np.int64)
            gs.create_dataset("offset", data=np.concatenate([[0],
                                                             np.cumsum(lens)]))
            for fld in ("snapnum", "t", "a", "r", "vr", "vt", "L", "r_ref",
                        "mass", "extrapolated"):
                parts = [r.series[fld] for r in self.records
                         if "snapnum" in r.series]
                if parts:
                    gs.create_dataset(fld, data=np.concatenate(parts))

    @classmethod
    def read_hdf5(cls, path: str) -> "OrbitCatalogue":
        with h5py.File(path, "r") as f:
            md = {k: (v.decode() if isinstance(v, bytes) else
                      (v.item() if hasattr(v, "item") else v))
                  for k, v in f.attrs.items()}
            ids = f["Satellites/HaloID"][()]
            n = len(ids)
            # tolerate sidecars written before a scalar existed: fill value
            scal = {name: (f["Satellites"][name][()] if name in f["Satellites"]
                           else np.full((n, 3) if name == "L_infall" else n,
                                        SCALARS[name][1],
                                        dtype=SCALARS[name][0]))
                    for name in SCALARS}
            ev_idx = f["Events/sat_index"][()]
            ev = {k: f["Events"][k][()] for k in
                  ("kind", "t", "a", "r", "r_over_ref")}
            off = f["Series/offset"][()]
            ser = {k: f["Series"][k][()] for k in f["Series"] if k != "offset"}
        records = []
        for i in range(n):
            sel = np.flatnonzero(ev_idx == i)
            events = [Event(int(ev["kind"][j]), float(ev["t"][j]),
                            float(ev["a"][j]), float(ev["r"][j]),
                            float(ev["r_over_ref"][j])) for j in sel]
            S = {k: v[off[i]:off[i + 1]] for k, v in ser.items()} \
                if off[i + 1] > off[i] else {}
            sc = {k: (v[i] if v.ndim == 1 else v[i]) for k, v in scal.items()}
            records.append(OrbitRecord(int(ids[i]), S, events, sc))
        return cls(records, md)
