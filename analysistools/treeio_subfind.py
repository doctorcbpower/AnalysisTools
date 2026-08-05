#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
treeio_subfind.py

Reader and main-branch walker for SubFind-HBT merger tree files
(the TreeHalos/* HDF5 layout), for use with
merger_tree_tools.MergerTreeTools.

Mirrors the shape of halo_tools's haloio_*.py readers: read_subfind_hbt()
is a pure function that returns a small data container (SubFindTreeData)
rather than mutating an instance; walk_subfind() and friends take that
container as an explicit argument.

Author: C. Power
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple

import h5py
import numpy as np

from .merger_tree_types import HaloTrack, MergerTreeError


@dataclass
class SubFindTreeData:
    """Everything read from a SubFind-HBT tree file, indexed and ready for
    walk_subfind() / resolve_group_first_sub() / get_group_members()."""
    metadata: Dict[str, Any]

    Redshift: np.ndarray
    Time: np.ndarray

    Omega0: float
    OmegaLambda: float
    OmegaBaryon: float
    HubbleParam: float
    BoxSize: float

    GrpNr: np.ndarray
    SubhaloNr: np.ndarray
    SnapNum: np.ndarray
    GrpM200: np.ndarray
    GrpR200: Optional[np.ndarray]
    SubhaloMass: np.ndarray
    SubhaloPos: np.ndarray
    SubhaloVel: np.ndarray
    SubhaloVelDisp: np.ndarray
    SubhaloVmax: np.ndarray
    SubhaloVmaxRad: np.ndarray

    TreeMainProgenitor: np.ndarray
    TreeFirstHaloInFOFgroup: np.ndarray
    TreeNextHaloInFOFgroup: np.ndarray

    # (snapnum, SubhaloNr) -> row index into the TreeHalos/* arrays
    lookup: Dict[Tuple[int, int], int]


def read_subfind_hbt(filename: str, logger=None) -> SubFindTreeData:
    """Read a SubFind-HBT tree file into a SubFindTreeData container.

    Returns raw values exactly as stored in the file -- comoving/little-h
    unit conversion is applied centrally by
    MergerTreeTools._apply_unit_conventions(), not here. Note SubFind-HBT's
    own native convention stores SubhaloPos/SubhaloVel *physical*, unlike
    the comoving convention group/subhalo catalogues normally use -- see
    MergerTreeTools.TREE_NATIVE_IS_COMOVING.
    """
    with h5py.File(filename, "r") as f:
        header = dict(f["Header"].attrs.items())
        if "Ntrees_Total" not in header:
            raise MergerTreeError(f"No trees found in '{filename}'.")
        if logger:
            logger.info(f"Found {header['Ntrees_Total']} merger trees in file.")

        Redshift = f["TreeTimes/Redshift"][()]
        Time = f["TreeTimes/Time"][()]

        Omega0 = f["Parameters"].attrs["Omega0"]
        OmegaLambda = f["Parameters"].attrs["OmegaLambda"]
        OmegaBaryon = f["Parameters"].attrs["OmegaBaryon"]
        HubbleParam = f["Parameters"].attrs["HubbleParam"]
        BoxSize = f["Parameters"].attrs["BoxSize"]

        GrpNr = f["TreeHalos/GroupNr"][()]
        SubhaloNr = f["TreeHalos/SubhaloNr"][()]
        SnapNum = f["TreeHalos/SnapNum"][()]
        GrpM200 = f["TreeHalos/Group_M_Crit200"][()]
        if "Group_R_Crit200" in f["TreeHalos"]:
            GrpR200 = f["TreeHalos/Group_R_Crit200"][()]
        else:
            if logger:
                logger.warning(
                    "TreeHalos/Group_R_Crit200 not found in this tree file -- "
                    "virial-radius overlays and analyse_orbit()'s crossing "
                    "detection won't be available unless you supply your own "
                    "radius array (e.g. from the halo catalogue) via a custom "
                    "radius_key in analyse_orbit().")
            GrpR200 = None
        SubhaloMass = f["TreeHalos/SubhaloMass"][()]
        SubhaloPos = f["TreeHalos/SubhaloPos"][()]
        SubhaloVel = f["TreeHalos/SubhaloVel"][()]
        SubhaloVelDisp = f["TreeHalos/SubhaloVelDisp"][()]
        SubhaloVmax = f["TreeHalos/SubhaloVmax"][()]
        SubhaloVmaxRad = f["TreeHalos/SubhaloVmaxRad"][()]

        TreeMainProgenitor = f["TreeHalos/TreeMainProgenitor"][()]
        TreeFirstHaloInFOFgroup = f["TreeHalos/TreeFirstHaloInFOFgroup"][()]
        TreeNextHaloInFOFgroup = f["TreeHalos/TreeNextHaloInFOFgroup"][()]

    # index lookup: (snapnum, subhalo id) -> array index. O(N) to build,
    # O(1) per lookup thereafter -- replaces repeated np.where chains.
    lookup = {(int(s), int(i)): idx
              for idx, (s, i) in enumerate(zip(SnapNum, SubhaloNr))}
    if logger:
        logger.info(f"Indexed {len(SubhaloNr):,} (snapshot, subhalo) tree nodes.")

    return SubFindTreeData(
        metadata=header, Redshift=Redshift, Time=Time,
        Omega0=Omega0, OmegaLambda=OmegaLambda, OmegaBaryon=OmegaBaryon,
        HubbleParam=HubbleParam, BoxSize=BoxSize,
        GrpNr=GrpNr, SubhaloNr=SubhaloNr, SnapNum=SnapNum,
        GrpM200=GrpM200, GrpR200=GrpR200,
        SubhaloMass=SubhaloMass, SubhaloPos=SubhaloPos, SubhaloVel=SubhaloVel,
        SubhaloVelDisp=SubhaloVelDisp, SubhaloVmax=SubhaloVmax, SubhaloVmaxRad=SubhaloVmaxRad,
        TreeMainProgenitor=TreeMainProgenitor,
        TreeFirstHaloInFOFgroup=TreeFirstHaloInFOFgroup,
        TreeNextHaloInFOFgroup=TreeNextHaloInFOFgroup,
        lookup=lookup,
    )


def resolve_group_first_sub(data: SubFindTreeData, group_id: int, snapnum: int) -> int:
    """Resolve a GroupNr to its central SubhaloNr at `snapnum`, using only
    data already in the tree (GrpNr + TreeFirstHaloInFOFgroup) -- no
    catalogue lookup required."""
    mask = (data.SnapNum == snapnum) & (data.GrpNr == group_id)
    candidates = np.where(mask)[0]
    if candidates.size == 0:
        raise MergerTreeError(f"No tree node found for GroupNr={group_id} at snapnum={snapnum}.")
    central = candidates[candidates == data.TreeFirstHaloInFOFgroup[candidates]]
    if central.size == 0:
        raise MergerTreeError(
            f"Could not find the central subhalo for GroupNr={group_id} "
            f"at snapnum={snapnum} (no node in the group satisfies "
            f"index == TreeFirstHaloInFOFgroup[index]).")
    return int(data.SubhaloNr[central[0]])


def get_group_members(data: SubFindTreeData, group_id: int, snapnum: int,
                       include_central: bool = False) -> List[int]:
    """List the SubhaloNr of every subhalo belonging to FOF group `group_id`
    at `snapnum`. If `include_central` is False (default), excludes the
    group's own central subhalo, i.e. returns only the satellites."""
    mask = (data.SnapNum == snapnum) & (data.GrpNr == group_id)
    idx = np.where(mask)[0]
    if idx.size == 0:
        raise MergerTreeError(f"No tree nodes found for GroupNr={group_id} at snapnum={snapnum}.")
    is_central = idx == data.TreeFirstHaloInFOFgroup[idx]
    keep = idx if include_central else idx[~is_central]
    return [int(x) for x in data.SubhaloNr[keep]]


def walk_subfind(data: SubFindTreeData, halo_id: int, snapnum: int,
                  max_length: int = 100_000) -> HaloTrack:
    """Walk a SubFind-HBT tree's main branch starting from (halo_id,
    snapnum) back to the root, returning a HaloTrack ordered earliest ->
    latest."""
    key = (int(snapnum), int(halo_id))
    if key not in data.lookup:
        raise MergerTreeError(f"(snapnum={snapnum}, id={halo_id}) not found in tree.")
    index = data.lookup[key]

    idx_list: List[int] = []
    idx = index
    for _ in range(max_length):
        idx_list.append(idx)
        idx = data.TreeMainProgenitor[idx]
        if idx == -1:
            break
    idx_list.reverse()  # earliest -> latest
    idx_arr = np.array(idx_list, dtype=np.int64)

    is_central = idx_arr == data.TreeFirstHaloInFOFgroup[idx_arr]
    host_index = data.TreeFirstHaloInFOFgroup[idx_arr]
    host_id = data.SubhaloNr[host_index]

    extra = {
        "GroupNr": data.GrpNr[idx_arr],
        "GroupM200": data.GrpM200[idx_arr],
        "VelDisp": data.SubhaloVelDisp[idx_arr],
        "Vmax": data.SubhaloVmax[idx_arr],
        "VmaxRad": data.SubhaloVmaxRad[idx_arr],
    }
    if data.GrpR200 is not None:
        extra["GroupR200"] = data.GrpR200[idx_arr]

    return HaloTrack(
        halo_id=halo_id, query_snapnum=snapnum, treefileformat="SubFind",
        SnapNum=data.SnapNum[idx_arr],
        Redshift=data.Redshift[data.SnapNum[idx_arr]],
        Time=data.Time[data.SnapNum[idx_arr]],
        Mass=data.SubhaloMass[idx_arr],
        Pos=data.SubhaloPos[idx_arr],
        Vel=data.SubhaloVel[idx_arr],
        IsSubhalo=~is_central,
        HostID=host_id,
        extra=extra,
    )
