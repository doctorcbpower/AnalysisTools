#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
treeio_ahf.py

Reader and main-branch walker for AHF MergerTree link files
(.AHF_mtree-style), for use with merger_tree_tools.MergerTreeTools.

AHF MergerTree files only encode descendant<->progenitor ID links, not
halo properties -- masses/positions/velocities must come from the halo
catalogue itself via a linked HaloTools instance, passed in to
walk_ahf() as `halo_tools_obj` (make sure its catalogue for each relevant
snapshot has been loaded before calling MergerTreeTools.get_track()).

Mirrors the shape of halo_tools's haloio_*.py readers: read_ahf_mergertree()
is a pure function that returns a small data container (AHFTreeData)
rather than mutating an instance; walk_ahf() takes that container (plus
the linked halo_tools_obj and a time/scale-factor array) as explicit
arguments.

Author: C. Power
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from .merger_tree_types import HaloTrack, MergerTreeError

# AHF's own convention (used both in .AHF_mtree link files and, commonly,
# baked into the "HaloID" column of .AHF_halos itself) is
# ID = SnapNum * multiplier + local_halo_index. This is the default;
# merger_tree_tools.py's AHF_SNAPNUM_ID_MULTIPLIER is the place to
# override it project-wide.
DEFAULT_AHF_SNAPNUM_ID_MULTIPLIER = 10 ** 12


@dataclass
class AHFTreeData:
    """Everything read from an AHF MergerTree link file: just the
    progenitor<->descendant ID links, plus the ID-decoding convention used
    to recover (snapnum, local_id) from a single packed AHF ID."""
    metadata: Dict[str, Any]

    desc_ids: np.ndarray
    prog_len: np.ndarray
    prog_ids: np.ndarray
    desc_index: np.ndarray

    prog_to_desc: Dict[int, int]              # progenitor id -> descendant id
    snapnum_id_multiplier: Optional[int]


def read_ahf_mergertree(filename: str,
                         snapnum_id_multiplier: Optional[int] = DEFAULT_AHF_SNAPNUM_ID_MULTIPLIER,
                         logger=None) -> AHFTreeData:
    """
    Read an AHF MergerTree link file into an AHFTreeData container.

    No comoving/little_h conversion here (unlike the other tree readers) --
    AHF MergerTree link files carry no physical (position/mass/velocity)
    properties of their own to convert.
    """
    with open(filename) as fh:
        next(fh); next(fh); next(fh)
        nlines, ndata = 0, 0
        for line in fh:
            if len(line.split()) > 1:
                ndata += 1
            nlines += 1

    mtreedata = np.zeros([nlines, 2], dtype=np.int64)
    with open(filename) as fh:
        next(fh); next(fh); next(fh)
        i = 0
        for line in fh:
            if "END" in line:
                break
            mtreedata[i] = line.split()
            i += 1

    if logger:
        logger.info(f"Number of halo links: {ndata}; number of halo IDs: {nlines}")

    desc_index = np.where(mtreedata[:, 0] != mtreedata[:, 1])[0]
    desc_ids = mtreedata[desc_index][:, 0]
    prog_len = mtreedata[desc_index][:, 1]
    prog_ids = mtreedata[np.where(mtreedata[:, 0] == mtreedata[:, 1])[0]][:, 0]

    # progenitor id -> descendant id, one step at a time
    prog_to_desc: Dict[int, int] = {}
    for k, d in enumerate(desc_ids):
        start = desc_index[k] - k
        finish = start + prog_len[k]
        for pid in prog_ids[start:finish]:
            prog_to_desc.setdefault(int(pid), int(d))

    metadata = {"n_links": ndata, "n_ids": nlines}

    return AHFTreeData(
        metadata=metadata, desc_ids=desc_ids, prog_len=prog_len, prog_ids=prog_ids,
        desc_index=desc_index, prog_to_desc=prog_to_desc,
        snapnum_id_multiplier=snapnum_id_multiplier,
    )


def ahf_decode(data: AHFTreeData, ahf_id: int) -> Tuple[int, int]:
    """
    Decode a packed AHF ID into (snapnum, local_id).

    ASSUMPTION: AHF halo IDs pack the snapshot number as
    ID = SnapNum * data.snapnum_id_multiplier + local_index. Set
    `snapnum_id_multiplier` at read_ahf_mergertree() time to match your
    AHF ID convention if it differs, or to None if your AHF IDs already
    come as plain (snapnum, id) pairs (in which case this function can't
    be used -- decode them yourself before calling get_track()).
    """
    if data.snapnum_id_multiplier is None:
        raise MergerTreeError(
            "snapnum_id_multiplier is None; set AHF_SNAPNUM_ID_MULTIPLIER "
            "in merger_tree_tools.py to match your AHF ID convention, or "
            "supply (snapnum, id) pairs directly.")
    snapnum = ahf_id // data.snapnum_id_multiplier
    local_id = ahf_id % data.snapnum_id_multiplier
    return int(snapnum), int(local_id)


def ahf_lookup_properties(halo_tools_obj, snapnum: int, local_id: int) -> Dict[str, Any]:
    """
    Look up mass/pos/vel/host status for a single AHF halo at a given
    snapshot via the linked HaloTools catalogue.

    Uses the *raw* (non-standardised) field names from haloio_ahf.py --
    "HaloID", "Mass", "Pos", "Vel", "HostHaloID" -- since
    halo_tools_standardise_names.py maps "id" -> "ID" for AHF, which
    doesn't match what read_ahf() actually produces ("HaloID"). AHF's
    convention is HostHaloID == -1 for a central halo, and the ID of the
    host halo otherwise, so subhalo status is read directly rather than
    inferred.

    NOTE: this expects `halo_tools_obj` to already have the catalogue for
    `snapnum` loaded (i.e. you've called read_catalogue for that
    snapshot, with standardise=False). This module does not manage
    per-snapshot catalogue loading/caching for you.
    """
    if halo_tools_obj is None:
        raise MergerTreeError(
            "AHF MergerTree has no property data of its own; pass a "
            "linked halo_tools=HaloTools(...) instance with the relevant "
            "snapshots already read (or re-readable) so properties can be "
            "looked up per-snapshot.")
    halos = halo_tools_obj.halos
    if halos is None:
        raise MergerTreeError(
            f"No catalogue currently loaded on the linked HaloTools "
            f"instance for snapshot {snapnum}. Load it (read_catalogue "
            f"with standardise=False) before calling get_track() for AHF "
            f"trees.")
    if "HaloID" not in halos:
        raise MergerTreeError(
            "Linked catalogue has no 'HaloID' field -- pass a raw "
            "(non-standardised) AHF catalogue, or adjust the field names "
            "used in ahf_lookup_properties() to match your setup.")
    idx = np.where(halos["HaloID"] == local_id)[0]
    if idx.size == 0:
        raise MergerTreeError(f"Halo id {local_id} not found at snapshot {snapnum}.")
    i = int(idx[0])
    host_field = halos.get("HostHaloID", np.full(len(halos["HaloID"]), -1))
    is_sub = bool(host_field[i] != -1)
    host_id = host_field[i] if is_sub else halos["HaloID"][i]
    return {
        "mass": halos["Mass"][i], "pos": halos["Pos"][i], "vel": halos["Vel"][i],
        "is_subhalo": is_sub, "host_id": host_id,
    }


def ahf_time_for_snaps(time_array: Optional[np.ndarray], snaps: np.ndarray) -> np.ndarray:
    """Look up scale factors for an array of snapshot numbers. AHF
    MergerTree link files carry no time table of their own, so this has
    to be supplied by the caller (see walk_ahf's `time_array` argument)."""
    if isinstance(time_array, np.ndarray):
        return time_array[snaps]
    raise MergerTreeError(
        "No scale-factor table available for AHF trees; supply one via "
        "mt.Time (e.g. from the linked HaloTools metadata) before calling "
        "get_track().")


def walk_ahf(data: AHFTreeData, halo_id: int, snapnum: int, max_length: int,
             halo_tools_obj=None, time_array: Optional[np.ndarray] = None) -> HaloTrack:
    """Walk an AHF MergerTree link chain starting from `halo_id` (at
    `snapnum`, used only for the returned HaloTrack.query_snapnum) forward
    through progenitor->descendant links, returning a HaloTrack ordered
    earliest -> latest."""
    chain_ids: List[int] = [int(halo_id)]
    cur = int(halo_id)
    for _ in range(max_length):
        desc = data.prog_to_desc.get(cur)
        if desc is None or desc == cur:
            break
        chain_ids.append(desc)
        cur = desc
    # chain_ids currently walks progenitor -> descendant, i.e. earliest ->
    # latest already, matching the convention used elsewhere.

    snaps, ids_local = zip(*(ahf_decode(data, i) for i in chain_ids))
    snaps = np.array(snaps)
    mass_l, pos_l, vel_l, host_l, sub_flag_l = [], [], [], [], []
    for s, local_id in zip(snaps, ids_local):
        props = ahf_lookup_properties(halo_tools_obj, int(s), int(local_id))
        mass_l.append(props["mass"])
        pos_l.append(props["pos"])
        vel_l.append(props["vel"])
        host_l.append(props["host_id"])
        sub_flag_l.append(props["is_subhalo"])

    time = ahf_time_for_snaps(time_array, snaps)
    redshift = 1.0 / time - 1.0

    return HaloTrack(
        halo_id=halo_id, query_snapnum=snapnum, treefileformat="MergerTree",
        SnapNum=snaps, Redshift=redshift, Time=time,
        Mass=np.array(mass_l), Pos=np.array(pos_l), Vel=np.array(vel_l),
        IsSubhalo=np.array(sub_flag_l, dtype=bool), HostID=np.array(host_l),
        extra={},
    )
