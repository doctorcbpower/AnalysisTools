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


@dataclass
class TreeFrogWalkableData:
    """A TreeFrog *walkable tree* (Head/Tail pointer file). Unlike the full
    tree flavour, it stores only the tree topology per snapshot (groups
    under 'Snapshots/Snap_NNN'); halo properties (mass, position, ...)
    stay in the matching halo catalogues. walk_treefrog() fills properties
    from linked catalogues where available and NaN otherwise."""
    metadata: Dict[str, Any]

    snap_of_group: Dict[str, int]
    Time: Dict[int, float]                      # snapnum -> scale factor

    ID: Dict[str, np.ndarray]
    Tail: Dict[str, np.ndarray]                 # main progenitor (temporal ID)
    TailSnap: Dict[str, np.ndarray]
    Head: Dict[str, np.ndarray]                 # descendant (temporal ID)
    HeadSnap: Dict[str, np.ndarray]
    Num_progen: Dict[str, np.ndarray]

    # (snapnum, halo id) -> (group key, local index within that group)
    lookup: Dict[Tuple[int, int], Tuple[str, int]]
    # (head_snap, head_id) -> every (group, local_idx) whose Head/HeadSnap
    # points there -- the reverse of `lookup`'s forward Tail walk, and
    # what makes get_progenitors_treefrog() possible: this is the *only*
    # one of this codebase's tree readers with a forward (descendant)
    # pointer at all (see get_progenitors_treefrog's docstring), so full
    # (not just main-branch) progenitor lists are TreeFrog-walkable-tree
    # only, not available for the "full tree" TreeFrog flavour or
    # SubFind-HBT.
    descendant_lookup: Dict[Tuple[int, int], list]


def _read_treefrog_walkable(f, header, logger=None) -> TreeFrogWalkableData:
    """Read the walkable-tree layout (root groups Header + Snapshots)."""
    snap_of_group: Dict[str, int] = {}
    Time: Dict[int, float] = {}
    ID: Dict[str, np.ndarray] = {}
    Tail: Dict[str, np.ndarray] = {}
    TailSnap: Dict[str, np.ndarray] = {}
    Head: Dict[str, np.ndarray] = {}
    HeadSnap: Dict[str, np.ndarray] = {}
    Num_progen: Dict[str, np.ndarray] = {}

    snaps = f["Snapshots"]
    for group in snaps.keys():
        g = snaps[group]
        if int(g.attrs.get("NHalos", 0)) <= 0:
            continue
        snapnum = int(g.attrs["Snapnum"])
        snap_of_group[group] = snapnum
        Time[snapnum] = float(g.attrs.get("scalefactor", 1.0))
        ID[group] = g["ID"][()]
        Tail[group] = g["Tail"][()]
        TailSnap[group] = g["TailSnap"][()]
        Head[group] = g["Head"][()]
        HeadSnap[group] = g["HeadSnap"][()]
        Num_progen[group] = g["Num_progen"][()] if "Num_progen" in g \
            else np.zeros(len(ID[group]), dtype=np.uint32)

    lookup: Dict[Tuple[int, int], Tuple[str, int]] = {}
    descendant_lookup: Dict[Tuple[int, int], list] = {}
    for group, snapnum in snap_of_group.items():
        for local_idx, hid in enumerate(ID[group]):
            lookup[(snapnum, int(hid))] = (group, local_idx)
            desc_key = (int(HeadSnap[group][local_idx]),
                       int(Head[group][local_idx]))
            descendant_lookup.setdefault(desc_key, []).append(
                (group, local_idx))
    if logger:
        logger.info(f"Walkable tree: indexed {len(lookup):,} "
                    f"(snapshot, halo) tree nodes.")

    return TreeFrogWalkableData(
        metadata=header, snap_of_group=snap_of_group, Time=Time,
        ID=ID, Tail=Tail, TailSnap=TailSnap, Head=Head, HeadSnap=HeadSnap,
        Num_progen=Num_progen, lookup=lookup,
        descendant_lookup=descendant_lookup,
    )


def read_treefrog(filename: str, logger=None):
    """Read a TreeFrog tree file.

    Handles both flavours:
    * full tree (per-snapshot root groups with Progenitor + properties)
      -> TreeFrogTreeData;
    * walkable tree (topology-only groups under 'Snapshots/')
      -> TreeFrogWalkableData.

    Returns raw values exactly as stored in the file -- comoving/little-h
    unit conversion is applied centrally by
    MergerTreeTools._apply_unit_conventions(), not here.
    """
    with h5py.File(filename, "r") as f:
        header = dict(f["Header"].attrs.items())
        if "Snapshots" in f:
            return _read_treefrog_walkable(f, header, logger=logger)
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


def _walkable_local_index(hid: int, snapnum: int, tid_val: int) -> int:
    """TreeFrog temporal-ID convention: ID = snapnum*tid_val + (index+1)."""
    return int(hid - snapnum * tid_val) - 1


def _walk_treefrog_walkable(data: TreeFrogWalkableData, halo_id: int,
                             snapnum: int, max_length: int = 100_000,
                             halo_tools_obj=None) -> HaloTrack:
    """Walk a walkable tree's main branch via Tail pointers.

    Halo properties are not stored in the walkable tree itself; they are
    filled from linked halo catalogue(s) where available and NaN otherwise.

    Parameters
    ----------
    halo_tools_obj : optional
        Either a single HaloTools-like object (with .snapnum and
        .standardised_halos), or a dict {snapnum: HaloTools-like}.
    """
    key = (int(snapnum), int(halo_id))
    if key not in data.lookup:
        raise MergerTreeError(f"(snapnum={snapnum}, id={halo_id}) not found in tree.")

    tid_val = int(data.metadata.get("Temporal_halo_id_value", 1_000_000_000_000))

    # normalise catalogue link to {snapnum: standardised table}
    tables: Dict[int, Dict[str, np.ndarray]] = {}
    if halo_tools_obj is not None:
        objs = halo_tools_obj if isinstance(halo_tools_obj, dict) \
            else {getattr(halo_tools_obj, "snapnum", None): halo_tools_obj}
        for snum, obj in objs.items():
            table = getattr(obj, "standardised_halos", None)
            if snum is not None and table is not None:
                tables[int(snum)] = table

    chain: list = []
    group, local_idx = data.lookup[key]
    for _ in range(max_length):
        chain.append((group, local_idx))
        own_id = int(data.ID[group][local_idx])
        tail = int(data.Tail[group][local_idx])
        tail_snap = int(data.TailSnap[group][local_idx])
        # Root: Tail is a self-loop (walkable-tree convention) or missing.
        if tail == own_id or tail <= 0:
            break
        nxt = (tail_snap, tail)
        if nxt not in data.lookup:
            break
        group, local_idx = data.lookup[nxt]

    chain.reverse()  # earliest -> latest

    n = len(chain)
    snap = np.array([data.snap_of_group[g] for g, _ in chain])
    own_id = np.array([data.ID[g][i] for g, i in chain], dtype=np.int64)
    time = np.array([data.Time[s] for s in snap])
    redshift = 1.0 / time - 1.0

    mass = np.full(n, np.nan)
    pos = np.full((n, 3), np.nan)
    vel = np.full((n, 3), np.nan)
    is_subhalo = np.zeros(n, dtype=bool)
    host_id = own_id.copy()

    for j, (s, hid) in enumerate(zip(snap, own_id)):
        table = tables.get(int(s))
        if table is None:
            continue
        cat_ids = table.get("halo_id")
        row = None
        if cat_ids is not None:
            matches = np.where(np.asarray(cat_ids) == hid)[0]
            if matches.size:
                row = int(matches[0])
        if row is None:
            # temporal-ID convention fallback: row = ID - snap*tid_val - 1
            cand = _walkable_local_index(int(hid), int(s), tid_val)
            ncat = len(next(iter(table.values()))) if table else 0
            if 0 <= cand < ncat:
                row = cand
        if row is None:
            continue
        if table.get("mass") is not None:
            mass[j] = table["mass"][row]
        if table.get("pos") is not None:
            pos[j] = table["pos"][row]
        if table.get("vel") is not None:
            vel[j] = table["vel"][row]
        hh = table.get("hostHaloID")
        if hh is not None and hh[row] != -1:
            is_subhalo[j] = True
            host_id[j] = hh[row]

    return HaloTrack(
        halo_id=halo_id, query_snapnum=snapnum, treefileformat="TreeFrog",
        SnapNum=snap, Redshift=redshift, Time=time,
        Mass=mass, Pos=pos, Vel=vel,
        IsSubhalo=is_subhalo, HostID=host_id,
        extra={"ID": own_id,
               "Num_progen": np.array([data.Num_progen[g][i]
                                       for g, i in chain])},
    )


def get_progenitors_treefrog(data: TreeFrogWalkableData, halo_id: int,
                             snapnum: int) -> list:
    """Every node the tree builder considers a progenitor of (halo_id,
    snapnum) -- not just the main-branch one ``Tail``/``walk_treefrog``
    follows. Uses ``descendant_lookup`` (every node's forward Head/
    HeadSnap pointer, inverted), so this returns the *actual* set of
    merging progenitors, cross-checked against (and normally identical
    in count to) the tree builder's own ``Num_progen`` for this node.

    A node at the most recent snapshot in the tree self-loops (Head ==
    its own ID, since it has no descendant yet), which would otherwise
    make a node its own "progenitor" in the reverse lookup -- filtered
    out explicitly.

    Returns
    -------
    list of dict, each ``{"halo_id", "snapnum", "is_main"}`` -- "is_main"
    is True for the one entry matching this node's own Tail/TailSnap
    (the branch ``walk_treefrog`` itself would follow), False for every
    other (i.e. merging-in) progenitor. Empty list for a root (no
    progenitor at all).
    """
    key = (int(snapnum), int(halo_id))
    if key not in data.lookup:
        raise MergerTreeError(
            f"(snapnum={snapnum}, id={halo_id}) not found in tree.")
    group, local_idx = data.lookup[key]
    own_tail = int(data.Tail[group][local_idx])
    own_tail_snap = int(data.TailSnap[group][local_idx])

    progenitors = []
    for prog_group, prog_idx in data.descendant_lookup.get(key, []):
        prog_id = int(data.ID[prog_group][prog_idx])
        prog_snap = data.snap_of_group[prog_group]
        if (prog_id, prog_snap) == (halo_id, snapnum):
            continue  # the final-snapshot self-loop, not a real progenitor
        progenitors.append({
            "halo_id": prog_id, "snapnum": prog_snap,
            "is_main": (prog_id, prog_snap) == (own_tail, own_tail_snap),
        })
    return progenitors


def walk_treefrog(data, halo_id: int, snapnum: int,
                   max_length: int = 100_000, halo_tools_obj=None) -> HaloTrack:
    """Walk a TreeFrog tree's main branch starting from (halo_id, snapnum)
    back to the root, returning a HaloTrack ordered earliest -> latest.

    Dispatches on the container type: TreeFrogTreeData (full tree, has
    its own properties) or TreeFrogWalkableData (topology only; properties
    filled from `halo_tools_obj` catalogues where available)."""
    if isinstance(data, TreeFrogWalkableData):
        return _walk_treefrog_walkable(data, halo_id, snapnum, max_length,
                                        halo_tools_obj=halo_tools_obj)
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
