#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
treeio_treefrog.py

Reader and main-branch walker for TreeFrog (VELOCIraptor) merger tree
files, for use with merger_tree_tools.MergerTreeTools.

Mirrors the shape of halo_tools's haloio_*.py readers: read_treefrog() is
a pure function that returns a small data container (TreeFrogTreeData)
rather than mutating an instance; walk_treefrog() takes that container as
an explicit argument.

Author: C. Power
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Tuple

import h5py
import numpy as np

from .merger_tree_types import HaloTrack, MergerTreeError


@dataclass
class TreeFrogTreeData:
    """Everything read from a TreeFrog tree file. TreeFrog stores data
    per-snapshot (one HDF5 group per "Snap_NNN"), so most fields here are
    dicts keyed by that group name; `lookup` makes this transparent to
    walk_treefrog()."""
    metadata: Dict[str, Any]
    HubbleParam: float

    snap_of_group: Dict[str, int]
    Time: Dict[int, float]                      # snapnum -> scale factor

    TreeProgenitor: Dict[str, np.ndarray]
    TreeProgenitorSnap: Dict[str, np.ndarray]
    TreeHalosID: Dict[str, np.ndarray]
    hostHaloID: Dict[str, np.ndarray]
    SubhaloMass: Dict[str, np.ndarray]
    SubhaloPos: Dict[str, np.ndarray]
    SubhaloVel: Dict[str, np.ndarray]
    SubhaloVmax: Dict[str, np.ndarray]
    SubhaloVmaxRad: Dict[str, np.ndarray]

    # (snapnum, halo id) -> (group key, local index within that group)
    lookup: Dict[Tuple[int, int], Tuple[str, int]]


def read_treefrog(filename: str, comoving: bool = False, logger=None) -> TreeFrogTreeData:
    """Read a TreeFrog tree file into a TreeFrogTreeData container."""
    with h5py.File(filename, "r") as f:
        header = dict(f["Header"].attrs.items())
        if "NSnaps" not in header:
            raise MergerTreeError(f"No trees found in '{filename}'.")
        if logger:
            logger.info(f"Found data for {header['NSnaps']:03d} snapshots in file.")

        HubbleParam = f["Header/Simulation"].attrs["h_val"]

        snap_of_group: Dict[str, int] = {}
        Time: Dict[int, float] = {}
        TreeProgenitor: Dict[str, np.ndarray] = {}
        TreeProgenitorSnap: Dict[str, np.ndarray] = {}
        TreeHalosID: Dict[str, np.ndarray] = {}
        hostHaloID: Dict[str, np.ndarray] = {}
        SubhaloMass: Dict[str, np.ndarray] = {}
        SubhaloPos: Dict[str, np.ndarray] = {}
        SubhaloVel: Dict[str, np.ndarray] = {}
        SubhaloVmax: Dict[str, np.ndarray] = {}
        SubhaloVmaxRad: Dict[str, np.ndarray] = {}

        for group in f.keys():
            if "Snap" not in group or f[group].attrs.get("NHalos", 0) <= 0:
                continue
            snapnum = int(f[group].attrs["Snapnum"])
            a = float(f[group].attrs["scalefactor"])
            snap_of_group[group] = snapnum
            Time[snapnum] = a

            TreeProgenitor[group] = f[group]["Progenitor"][()]
            TreeProgenitorSnap[group] = f[group]["ProgenitorSnap"][()]
            TreeHalosID[group] = f[group]["ID"][()]
            hostHaloID[group] = f[group]["hostHaloID"][()]

            mass = f[group]["Mass_tot"][()]
            pos = np.array([f[group]["Xcminpot"][()],
                             f[group]["Ycminpot"][()],
                             f[group]["Zcminpot"][()]]).T
            vel = np.array([f[group]["VXc"][()],
                             f[group]["VYc"][()],
                             f[group]["VZc"][()]]).T

            if comoving:
                pos = pos * HubbleParam / a
                mass = mass * HubbleParam

            SubhaloMass[group] = mass
            SubhaloPos[group] = pos
            SubhaloVel[group] = vel
            SubhaloVmax[group] = f[group]["Vmax"][()]
            SubhaloVmaxRad[group] = f[group]["Rmax"][()]

    # lookup: (snapnum, halo id) -> (group_key, local_index)
    lookup: Dict[Tuple[int, int], Tuple[str, int]] = {}
    for group, snapnum in snap_of_group.items():
        for local_idx, hid in enumerate(TreeHalosID[group]):
            lookup[(snapnum, int(hid))] = (group, local_idx)
    if logger:
        logger.info(f"Indexed {len(lookup):,} (snapshot, halo) tree nodes.")

    return TreeFrogTreeData(
        metadata=header, HubbleParam=HubbleParam,
        snap_of_group=snap_of_group, Time=Time,
        TreeProgenitor=TreeProgenitor, TreeProgenitorSnap=TreeProgenitorSnap,
        TreeHalosID=TreeHalosID, hostHaloID=hostHaloID,
        SubhaloMass=SubhaloMass, SubhaloPos=SubhaloPos, SubhaloVel=SubhaloVel,
        SubhaloVmax=SubhaloVmax, SubhaloVmaxRad=SubhaloVmaxRad,
        lookup=lookup,
    )


def walk_treefrog(data: TreeFrogTreeData, halo_id: int, snapnum: int,
                   max_length: int = 100_000) -> HaloTrack:
    """Walk a TreeFrog tree's main branch starting from (halo_id, snapnum)
    back to the root, returning a HaloTrack ordered earliest -> latest."""
    key = (int(snapnum), int(halo_id))
    if key not in data.lookup:
        raise MergerTreeError(f"(snapnum={snapnum}, id={halo_id}) not found in tree.")

    chain: list = []
    group, local_idx = data.lookup[key]
    for _ in range(max_length):
        chain.append((group, local_idx))
        prog_ids = data.TreeProgenitor[group][local_idx]
        prog_snaps = data.TreeProgenitorSnap[group][local_idx]
        prog_id = np.atleast_1d(prog_ids)
        prog_snap = np.atleast_1d(prog_snaps)
        # Root of the branch: TreeFrog conventions vary between a
        # self-loop (progenitor id == own id) and a -1 sentinel for "no
        # progenitor". Handle both.
        if prog_id.size == 0 or int(prog_id[0]) == -1 or \
                int(prog_id[0]) == int(data.TreeHalosID[group][local_idx]):
            break
        nxt_snap = int(prog_snap[0])
        nxt_group = f"Snap_{nxt_snap:03d}"
        if nxt_group not in data.TreeHalosID:
            break
        matches = np.where(data.TreeHalosID[nxt_group] == prog_id[0])[0]
        if matches.size == 0:
            break
        group, local_idx = nxt_group, int(matches[0])

    chain.reverse()  # earliest -> latest

    snap = np.array([data.snap_of_group[g] for g, _ in chain])
    mass = np.array([data.SubhaloMass[g][i] for g, i in chain])
    pos = np.array([data.SubhaloPos[g][i] for g, i in chain])
    vel = np.array([data.SubhaloVel[g][i] for g, i in chain])
    vmax = np.array([data.SubhaloVmax[g][i] for g, i in chain])
    vmaxrad = np.array([data.SubhaloVmaxRad[g][i] for g, i in chain])
    host = np.array([data.hostHaloID[g][i] for g, i in chain])
    own_id = np.array([data.TreeHalosID[g][i] for g, i in chain])
    is_subhalo = host != -1
    host_id = np.where(is_subhalo, host, own_id)
    time = np.array([data.Time[s] for s in snap])
    redshift = 1.0 / time - 1.0

    return HaloTrack(
        halo_id=halo_id, query_snapnum=snapnum, treefileformat="TreeFrog",
        SnapNum=snap, Redshift=redshift, Time=time,
        Mass=mass, Pos=pos, Vel=vel,
        IsSubhalo=is_subhalo, HostID=host_id,
        extra={"Vmax": vmax, "VmaxRad": vmaxrad, "ID": own_id},
    )
