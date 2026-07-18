#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
merger_tree_types.py

Shared, format-agnostic types used across merger_tree_tools.py and the
treeio_*.py readers/walkers: the MergerTreeError exception, the
HaloTrack / OrbitAnalysis result containers, and small numerical helpers.

Keeping these here (rather than in merger_tree_tools.py) lets the
format-specific treeio_*.py modules build HaloTrack objects directly,
without importing back from merger_tree_tools.py and creating a circular
import -- mirroring the way halo_tools_standardise_names.py sits
underneath halo_tools.py's format readers.

Author: C. Power
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, Optional, Union

import numpy as np


class MergerTreeError(Exception):
    """Raised for tree-reading or tree-walking failures."""


# ---------------------------------------------------------------------
# Small numerical helpers
# ---------------------------------------------------------------------

def periodic_delta(dx: np.ndarray, boxsize: Optional[float]) -> np.ndarray:
    """Wrap a displacement (or array of displacements) into [-L/2, L/2)."""
    if boxsize is None or boxsize <= 0:
        return dx
    return (dx + 0.5 * boxsize) % boxsize - 0.5 * boxsize


# ---------------------------------------------------------------------
# HaloTrack: the format-agnostic result of walking a main branch
# ---------------------------------------------------------------------

@dataclass
class HaloTrack:
    """
    Time-ordered history of a single halo/subhalo, walking from the earliest
    snapshot in the tree (root) to the snapshot it was queried at.

    All arrays share the same length and ordering (ascending SnapNum / Time).
    """
    halo_id: Union[int, np.integer]
    query_snapnum: int
    treefileformat: str

    SnapNum: np.ndarray
    Redshift: np.ndarray
    Time: np.ndarray               # scale factor
    Mass: np.ndarray
    Pos: np.ndarray                # (N, 3)
    Vel: np.ndarray                # (N, 3)
    IsSubhalo: np.ndarray          # bool, True where the halo is a satellite
    HostID: np.ndarray             # id of the central/host at that snapshot
                                    # (equal to the halo's own id while central)
    extra: Dict[str, np.ndarray] = field(default_factory=dict)

    def __len__(self) -> int:
        return len(self.SnapNum)

    def __repr__(self) -> str:
        return (f"HaloTrack(id={self.halo_id}, format={self.treefileformat}, "
                f"nsnaps={len(self)}, "
                f"z_range=({self.Redshift.min():.2f}, {self.Redshift.max():.2f}))")

    def get(self, name: str) -> np.ndarray:
        """Fetch a quantity by name, checking core fields then `extra`."""
        if hasattr(self, name):
            return getattr(self, name)
        if name in self.extra:
            return self.extra[name]
        raise KeyError(f"'{name}' is not a field on this HaloTrack "
                        f"(core fields: SnapNum, Redshift, Time, Mass, Pos, "
                        f"Vel, IsSubhalo, HostID; extra: {list(self.extra)})")

    @property
    def is_currently_subhalo(self) -> bool:
        return bool(self.IsSubhalo[-1]) if len(self) else False

    def infall_snapshot(self) -> Optional[Dict[str, Any]]:
        """
        Identify the moment this halo first becomes a subhalo, i.e. the
        earliest-time transition from central (IsSubhalo == False) to
        satellite (IsSubhalo == True), walking forward in time.

        Returns
        -------
        dict with keys: index, snapshot, redshift, time, mass, host_id
        or None if the halo is never a subhalo anywhere along the track.
        """
        if not np.any(self.IsSubhalo):
            return None
        idx = int(np.argmax(self.IsSubhalo))  # first True, in ascending time
        return {
            "index": idx,
            "snapshot": int(self.SnapNum[idx]),
            "redshift": float(self.Redshift[idx]),
            "time": float(self.Time[idx]),
            "mass": float(self.Mass[idx]),
            "host_id": self.HostID[idx],
        }

    def to_dict(self) -> Dict[str, np.ndarray]:
        d = {
            "SnapNum": self.SnapNum, "Redshift": self.Redshift, "Time": self.Time,
            "Mass": self.Mass, "Pos": self.Pos, "Vel": self.Vel,
            "IsSubhalo": self.IsSubhalo, "HostID": self.HostID,
        }
        d.update(self.extra)
        return d


# ---------------------------------------------------------------------
# OrbitAnalysis: the format-agnostic result of analyse_orbit()
# ---------------------------------------------------------------------

@dataclass
class OrbitAnalysis:
    """
    A subhalo's orbit relative to a reference (host) track, over the
    snapshots the two tracks share, decomposed into radial/tangential
    components, with the reference's virial radius carried alongside for
    infall/crossing detection.

    All arrays are time-ordered (ascending SnapNum), matching HaloTrack.
    """
    halo_id: Union[int, np.integer]
    host_id: Union[int, np.integer]

    SnapNum: np.ndarray
    Redshift: np.ndarray
    Time: np.ndarray

    Distance: np.ndarray            # |pos - host_pos|, periodic-wrapped
    RelVel: np.ndarray              # |vel - host_vel|
    RadialVelocity: np.ndarray      # (vel - host_vel) . r_hat  (+ve = outgoing)
    TangentialVelocity: np.ndarray  # component of relative velocity perp. to r_hat
    Mass: np.ndarray                # subhalo mass at each snapshot
    Rvir: Optional[np.ndarray] = None   # host virial radius at each snapshot, if available

    def __len__(self) -> int:
        return len(self.SnapNum)

    def first_crossing(self) -> Optional[Dict[str, Any]]:
        """
        First time (walking forward from earliest snapshot) that Distance
        drops to or below Rvir -- i.e. the halo's first virial-radius
        crossing / infall into the host.

        Returns
        -------
        dict with keys: index, snapshot, redshift, time, mass, distance,
        rvir, radial_velocity, tangential_velocity, relative_speed
        or None if Rvir wasn't available, or the subhalo is never inside it.
        """
        if self.Rvir is None:
            return None
        inside = self.Distance <= self.Rvir
        if not np.any(inside):
            return None
        idx = int(np.argmax(inside))
        return {
            "index": idx,
            "snapshot": int(self.SnapNum[idx]),
            "redshift": float(self.Redshift[idx]),
            "time": float(self.Time[idx]),
            "mass": float(self.Mass[idx]),
            "distance": float(self.Distance[idx]),
            "rvir": float(self.Rvir[idx]),
            "radial_velocity": float(self.RadialVelocity[idx]),
            "tangential_velocity": float(self.TangentialVelocity[idx]),
            "relative_speed": float(self.RelVel[idx]),
        }
