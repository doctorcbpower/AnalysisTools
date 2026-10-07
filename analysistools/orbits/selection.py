#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
analysistools.orbits.selection
------------------------------
Which subhaloes count as satellites of a host at a given snapshot.

Definition used for the Dorcha papers: a satellite is any subhalo whose
centre lies within the reference radius (default R200c, ``radius_factor``
adjustable) of the host centre, **regardless of which FOF group the tree puts
it in**, excluding the host's own central subhalo. FOF membership is *not*
the criterion: on the Dorcha host 25 per cent of z=0 FOF members lie
outside R200c (out to 1.8 R200c), and subhaloes of neighbouring FOF groups
can lie inside it.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

from .cosmology import FlatLCDM
from .kinematics import TrackUnits


@dataclass
class SatelliteSelection:
    subhalo_id: np.ndarray       # SubhaloNr of each selected subhalo
    index: np.ndarray            # row into the arrays passed in
    r: np.ndarray                # physical kpc from the host centre
    r_over_ref: np.ndarray
    in_host_fof: np.ndarray      # bool: also in the host's FOF group
    r_ref: float                 # physical kpc

    def __len__(self):
        return len(self.subhalo_id)


def select_satellites(pos, subhalo_id, host_pos, r_ref_phys: float, a: float,
                      units: TrackUnits, boxsize: Optional[float] = None,
                      group_nr=None, host_group: Optional[int] = None,
                      exclude_ids=(), mass=None, nparticles=None,
                      min_mass: Optional[float] = None,
                      min_particles: Optional[int] = None) -> SatelliteSelection:
    """Subhaloes of one snapshot within ``r_ref_phys`` (physical kpc) of
    ``host_pos``.

    ``pos``/``host_pos`` are in the tracks' own length units (converted via
    ``units``); ``boxsize`` is the comoving box side in those units.
    ``exclude_ids`` removes e.g. the host's central subhalo.

    Optional user-defined resolution thresholds, applied *at this snapshot*
    (bound values now, not peak -- use ``OrbitCatalogue.select(basis="peak")``
    for peak-based cuts): ``min_mass`` in Msun with ``mass`` in tree units
    (converted with ``units``), and/or ``min_particles`` with ``nparticles``.
    Both must hold if both are given. Supplying a threshold without its data
    raises, rather than silently selecting everything.
    """
    if min_mass is not None and mass is None:
        raise ValueError("min_mass given but no `mass` array.")
    if min_particles is not None and nparticles is None:
        raise ValueError("min_particles given but no `nparticles` array "
                         "(does the tree carry SubhaloLen?).")
    pos_c = units.comoving_kpc(np.asarray(pos, dtype=float),
                               np.full(len(pos), a))
    host_c = units.comoving_kpc(np.asarray(host_pos, dtype=float)[None],
                                np.array([a]))[0]
    dx = pos_c - host_c
    if boxsize:
        box_c = float(boxsize) * units.length_kpc * units._hfac()
        dx = (dx + 0.5 * box_c) % box_c - 0.5 * box_c
    r = np.linalg.norm(dx, axis=1) * a
    keep = r <= r_ref_phys
    ids = np.asarray(subhalo_id)
    if len(exclude_ids):
        keep &= ~np.isin(ids, list(exclude_ids))
    if min_mass is not None:
        keep &= units.mass(np.asarray(mass)) >= float(min_mass)
    if min_particles is not None:
        keep &= np.asarray(nparticles) >= int(min_particles)
    idx = np.flatnonzero(keep)
    in_fof = (np.asarray(group_nr)[idx] == host_group
              if group_nr is not None and host_group is not None
              else np.zeros(idx.size, dtype=bool))
    return SatelliteSelection(ids[idx], idx, r[idx], r[idx] / r_ref_phys,
                              in_fof, float(r_ref_phys))


def select_satellites_subfind(tree, host_track, cosmology: FlatLCDM,
                              units: TrackUnits, snapnum: int,
                              radius_key: str = "R200c_from_M200",
                              radius_factor: float = 1.0,
                              fixed_radius_kpc: Optional[float] = None,
                              min_mass: Optional[float] = None,
                              min_particles: Optional[int] = None,
                              min_mass_fraction: Optional[float] = None
                              ) -> SatelliteSelection:
    """``select_satellites`` for a ``MergerTreeTools`` SubFind-HBT tree: all
    subhaloes at ``snapnum`` within the reference radius of the host's
    central subhalo (the last point of ``host_track``)."""
    from .kinematics import host_reference
    d = tree.data
    rows = np.flatnonzero(np.asarray(d.SnapNum) == snapnum)
    k = int(np.flatnonzero(host_track.SnapNum == snapnum)[0])
    r_ref = float(host_reference(host_track, units, radius_key, radius_factor,
                                 fixed_radius_kpc, cosmology)[k])
    a = float(host_track.Time[k])
    host_id = int(host_track.halo_id)
    if min_mass_fraction is not None:
        m200 = float(units.mass(host_track.extra["GroupM200"][k]))
        min_mass = max(float(min_mass or 0.0), float(min_mass_fraction) * m200)
    return select_satellites(
        np.asarray(d.SubhaloPos)[rows], np.asarray(d.SubhaloNr)[rows],
        host_track.Pos[k], r_ref, a, units, tree.BoxSize,
        group_nr=np.asarray(d.GrpNr)[rows],
        host_group=int(host_track.extra["GroupNr"][k]),
        exclude_ids=(host_id,),
        mass=np.asarray(d.SubhaloMass)[rows],
        nparticles=(None if getattr(d, "SubhaloLen", None) is None
                    else np.asarray(d.SubhaloLen)[rows]),
        min_mass=min_mass, min_particles=min_particles)
