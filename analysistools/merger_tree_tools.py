#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
merger_tree_tools.py
Refactored: 2026-07-03

Unified interface for reading and analysing halo merger trees produced by
SubFind-HBT MergerTree, TreeFrog (VELOCIraptor), and AHF MergerTree, designed
to work seamlessly alongside halo_tools.HaloTools.

Design summary
--------------
Every supported tree format is parsed once, by a dedicated reader in its
own module (treeio_subfind.py, treeio_treefrog.py, treeio_ahf.py), into a
small data container plus a fast (SnapNum, ID) -> index lookup table. From
there, a matching format-specific `walk_*()` function walks the main
branch and returns a `HaloTrack`: a small container of time-ordered numpy
arrays (mass, position, velocity, subhalo status, host id, ...). This
class just dispatches to those readers/walkers -- all downstream analysis
and plotting works off `HaloTrack` objects and is completely
format-agnostic.

Module layout
-------------
merger_tree_types.py   HaloTrack, OrbitAnalysis, MergerTreeError, periodic_delta
treeio_subfind.py       SubFind-HBT reader (read_subfind_hbt) + walker (walk_subfind)
treeio_treefrog.py      TreeFrog/VELOCIraptor reader (read_treefrog) + walker (walk_treefrog)
treeio_ahf.py           AHF MergerTree reader (read_ahf_mergertree) + walker (walk_ahf)
merger_tree_tools.py    MergerTreeTools (this file): dispatch + analysis + plotting

Typical usage
-------------
>>> ht = HaloTools(comoving_units=True)
>>> ht.read_catalogue("groups_010.hdf5", fileformat="SubFind", standardise=True)
>>>
>>> mt = MergerTreeTools("trees.hdf5", treefileformat="SubFind")
>>> track = mt.from_halo(ht, index=42, object_type="Subhalo")
>>> mt.plot_track(track)
>>>
>>> host_track = mt.get_track(host_id, snapnum)
>>> mt.plot_relative(track, host_track)
>>>
>>> event = track.infall_snapshot()
>>> print(event["snapshot"], event["mass"])

Author: C. Power
"""
from __future__ import annotations

import os
import logging
from typing import Optional, Dict, Any, List, Tuple, Union

import numpy as np
import matplotlib.pyplot as plt

from .merger_tree_types import MergerTreeError, HaloTrack, OrbitAnalysis, periodic_delta
from .treeio_subfind import (
    read_subfind_hbt, walk_subfind, resolve_group_first_sub,
    get_group_members as _subfind_get_group_members,
)
from .treeio_treefrog import read_treefrog, walk_treefrog
from .treeio_ahf import read_ahf_mergertree, walk_ahf

# ---------------------------------------------------------------------
# ASSUMPTIONS / CONFIGURATION
# ---------------------------------------------------------------------
# These are the only places you should need to touch if your
# halo_tools_standardise_names.py schema, or your AHF MergerTree ID/SnapNum
# encoding, differs from what's assumed below. I don't have visibility into
# either of those files, so these are best-effort defaults.

# Field name expected in ht.halos / ht.subhalos *after* standardise_names()
# has been applied, used by MergerTreeTools.from_halo(). This matches
# STANDARD_KEYS in halo_tools_standardise_names.py (lowercase).
STD_ID_FIELD = "halo_id"

# Sentinel meaning "this table has no stored ID field; the ID is the row
# index" -- matches ROW_INDEX_SENTINEL in halo_tools_standardise_names.py.
# SubFind subhalos are numbered this way (no SubhaloNr dataset; subhalo N is
# row N of the Subhalo/* arrays), which is also how the SubFind-HBT tree
# file's TreeHalos/SubhaloNr values were assigned in the first place -- so
# this only holds for a catalogue read from a single, unsplit file.
ROW_INDEX_SENTINEL = "__ROW_INDEX__"

# Native (pre-standardisation) ID field name per (treefileformat, object_type),
# used as a fallback / for AHF property look-ups, since it's what the raw
# reader output in haloio_*.py actually contains today. Adjust if your local
# copies of the haloio_*.py readers differ.
RAW_ID_FIELD = {
    ("SubFind", "Group"): "GroupNr",
    ("SubFind", "Subhalo"): ROW_INDEX_SENTINEL,
    ("MergerTree", "Group"): "HaloID",     # haloio_ahf.py
    ("MergerTree", "Subhalo"): "HaloID",
}

# AHF MergerTree encodes a single integer ID per halo per snapshot. AHF's own
# convention (used both in the .AHF_mtree link files and, commonly, baked
# into the "HaloID" column of .AHF_halos itself) is
# ID = SnapNum * 10**12 + local_halo_index. Set this to match your files;
# set to None if your AHF IDs already come as plain (snapnum, id) pairs.
AHF_SNAPNUM_ID_MULTIPLIER: Optional[int] = 10**12


# ---------------------------------------------------------------------
# Main user-facing class
# ---------------------------------------------------------------------

class MergerTreeTools:
    """
    Unified high-level interface to load a merger tree and extract/plot
    HaloTrack objects for individual halos or subhalos.

    Parameters
    ----------
    treefilename : str
        Path to the tree file.
    treefileformat : str
        One of "SubFind", "TreeFrog", "MergerTree" (AHF).
    comoving_units : bool
        If True, positions/masses are converted to comoving units on read.
    halo_tools : HaloTools, optional
        A linked catalogue reader, enabling `from_halo()` and (for AHF)
        per-snapshot property lookups.
    """

    SUPPORTED_FORMATS = ("SubFind", "TreeFrog", "MergerTree")

    def __init__(
        self,
        treefilename: str,
        treefileformat: str,
        comoving_units: bool = False,
        halo_tools: Optional["object"] = None,
        **kwargs,
    ):
        if treefileformat not in self.SUPPORTED_FORMATS:
            raise ValueError(f"Unsupported tree format '{treefileformat}'. "
                              f"Must be one of {self.SUPPORTED_FORMATS}.")

        self.treefilename = treefilename
        self.treefileformat = treefileformat
        self.comoving_units = comoving_units
        self.halo_tools = halo_tools

        self.metadata: Dict[str, Any] = {}
        self.data = None            # format-specific data container, set by read_tree()
        self.BoxSize: Optional[float] = None
        self._loaded = False

        self.logger = logging.getLogger(__name__)
        if not self.logger.handlers:
            handler = logging.StreamHandler()
            handler.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))
            self.logger.addHandler(handler)
            self.logger.propagate = False
        self.logger.setLevel(kwargs.get("loglevel", logging.INFO))

        self.read_tree()

    # ------------------------------------------------------------------
    # Reading (dispatch to treeio_*.py)
    # ------------------------------------------------------------------

    def read_tree(self) -> None:
        """Dispatch to the format-specific reader. Safe to call once; the
        constructor already calls this, so you shouldn't normally need to."""
        if self._loaded:
            return
        if self.treefileformat == "SubFind":
            self.data = read_subfind_hbt(
                self.treefilename, comoving=self.comoving_units, logger=self.logger)
        elif self.treefileformat == "TreeFrog":
            self.data = read_treefrog(
                self.treefilename, comoving=self.comoving_units, logger=self.logger)
        elif self.treefileformat == "MergerTree":
            self.data = read_ahf_mergertree(
                self.treefilename, snapnum_id_multiplier=AHF_SNAPNUM_ID_MULTIPLIER,
                logger=self.logger)
        else:
            raise MergerTreeError(f"Unhandled format '{self.treefileformat}'")

        self.metadata = self.data.metadata
        self.BoxSize = getattr(self.data, "BoxSize", None)  # only SubFind-HBT carries this
        self._loaded = True

    # ------------------------------------------------------------------
    # Main-branch walk -> HaloTrack (dispatch to treeio_*.py)
    # ------------------------------------------------------------------

    def get_track(self, halo_id: int, snapnum: int, object_type: str = "Subhalo",
                   max_length: int = 100_000) -> HaloTrack:
        """
        Walk the main branch of the tree starting from (halo_id, snapnum)
        back to the root, and return it as a HaloTrack ordered from the
        earliest snapshot to `snapnum`.

        Parameters
        ----------
        halo_id : int
            Halo/subhalo ID as it appears in the halo catalogue at `snapnum`
            (for SubFind/TreeFrog), or the AHF MergerTree ID (for AHF).
        snapnum : int
            Snapshot number at which `halo_id` is defined.
        object_type : {"Subhalo", "Group"}
            If "Group", `halo_id` (a GroupNr) is resolved to its central
            subhalo *using the tree file itself* -- the halo catalogue
            doesn't carry a first-subhalo pointer, but the tree's
            TreeFirstHaloInFOFgroup does. Only meaningful for SubFind.
        """
        if object_type not in ("Subhalo", "Group"):
            raise ValueError("object_type must be 'Subhalo' or 'Group'")

        if object_type == "Group":
            halo_id = self._resolve_group_first_sub(halo_id, snapnum)

        if self.treefileformat == "SubFind":
            return walk_subfind(self.data, halo_id, snapnum, max_length)
        elif self.treefileformat == "TreeFrog":
            return walk_treefrog(self.data, halo_id, snapnum, max_length)
        elif self.treefileformat == "MergerTree":
            return walk_ahf(self.data, halo_id, snapnum, max_length,
                             halo_tools_obj=self.halo_tools,
                             time_array=getattr(self, "Time", None))
        raise MergerTreeError(f"Unhandled format '{self.treefileformat}'")

    def _resolve_group_first_sub(self, group_id: int, snapnum: int) -> int:
        if self.treefileformat != "SubFind":
            raise MergerTreeError(
                f"object_type='Group' resolution is only implemented for "
                f"SubFind-HBT trees, not '{self.treefileformat}'.")
        return resolve_group_first_sub(self.data, group_id, snapnum)

    def get_group_members(self, group_id: int, snapnum: int,
                           include_central: bool = False) -> List[int]:
        """
        List the SubhaloNr of every subhalo belonging to FOF group `group_id`
        at `snapnum` (SubFind-HBT only -- TreeFrog/AHF track subhalo
        membership via hostHaloID per-object rather than a group table, so
        there's no single group-membership query for them here).

        Parameters
        ----------
        include_central : bool
            If False (default), excludes the group's own central subhalo,
            i.e. returns only the satellites.
        """
        if self.treefileformat != "SubFind":
            raise MergerTreeError(
                f"get_group_members() is only implemented for SubFind-HBT "
                f"trees, not '{self.treefileformat}'.")
        return _subfind_get_group_members(self.data, group_id, snapnum,
                                           include_central=include_central)

    def get_group_subhalo_tracks(self, group_id: int, snapnum: int,
                                  include_central: bool = False) -> Dict[int, HaloTrack]:
        """
        Convenience wrapper around get_group_members() + get_track(): returns
        {subhalo_id: HaloTrack} for every (satellite) member of the group.
        """
        members = self.get_group_members(group_id, snapnum, include_central=include_central)
        return {sub_id: self.get_track(sub_id, snapnum, object_type="Subhalo") for sub_id in members}

    # ------------------------------------------------------------------
    # Integration with HaloTools
    # ------------------------------------------------------------------

    def from_halo(self, halo_tools_obj, index: int, object_type: str = "Subhalo",
                   snapnum: Optional[int] = None, id_field: Optional[str] = None,
                   standardised: bool = True) -> HaloTrack:
        """
        Convenience wrapper: given a HaloTools instance and the row index of
        a halo/subhalo of interest in its currently-loaded catalogue,
        extract the ID and return the corresponding HaloTrack.

        Parameters
        ----------
        snapnum : int, optional
            Snapshot the catalogue was read at. If not given, taken from
            `halo_tools_obj.snapnum` (set automatically if you passed
            `snapnum=` to `read_catalogue`). One of the two is required.
        id_field : str, optional
            Field name to read the halo ID from. Defaults to the
            standardised key "halo_id" if `standardised=True`, otherwise to
            RAW_ID_FIELD[(treefileformat, object_type)]. Override this if
            neither matches your catalogue.
        standardised : bool
            Whether `halo_tools_obj.halos`/`.subhalos` were produced with
            `read_catalogue(..., standardise=True)`. SubFind Subhalo tables
            have no stored ID field (subhalo N is just row N of the
            Subhalo/* arrays), so this is handled either way -- via the
            standardised "halo_id" column (populated as arange(n) by
            halo_tools_standardise_names.py) or, with standardised=False,
            directly as the row index.
        """
        if standardised:
            table = halo_tools_obj.standardised_subhalos if object_type == "Subhalo" \
                else halo_tools_obj.standardised_halos
        else:
            table = halo_tools_obj.subhalos if object_type == "Subhalo" else halo_tools_obj.halos
        if table is None:
            raise MergerTreeError(
                f"No {'standardised ' if standardised else ''}{object_type.lower()} "
                f"catalogue loaded on halo_tools_obj. "
                f"{'Call read_catalogue(..., standardise=True) first.' if standardised else ''}")

        if snapnum is None:
            snapnum = getattr(halo_tools_obj, "snapnum", None)
        if snapnum is None:
            raise MergerTreeError(
                "Could not determine the snapshot number. Either pass "
                "snapnum=... explicitly, or call "
                "halo_tools_obj.read_catalogue(..., snapnum=...) so it's "
                "recorded automatically.")

        if id_field is None:
            id_field = STD_ID_FIELD if standardised else \
                RAW_ID_FIELD.get((self.treefileformat, object_type))

        if id_field == ROW_INDEX_SENTINEL:
            halo_id = index
        else:
            if id_field is None or id_field not in table or table[id_field] is None:
                raise MergerTreeError(
                    f"Could not find a usable ID field ('{id_field}') in the "
                    f"{object_type.lower()} catalogue (fields: {list(table)}). "
                    f"Try standardised={not standardised}, or pass "
                    f"id_field=... explicitly.")
            halo_id = table[id_field][index]

        return self.get_track(int(halo_id), int(snapnum), object_type=object_type)

    # ------------------------------------------------------------------
    # Analysis (format-agnostic: works purely off HaloTrack objects)
    # ------------------------------------------------------------------

    def find_infall(self, track: HaloTrack) -> Optional[Dict[str, Any]]:
        """Thin wrapper around HaloTrack.infall_snapshot() for symmetry
        with the rest of the API."""
        return track.infall_snapshot()

    def analyse_orbit(self, track: HaloTrack, reference_track: HaloTrack,
                       radius_key: str = "GroupR200",
                       boxsize: Optional[float] = None) -> OrbitAnalysis:
        """
        Compute `track`'s orbit relative to `reference_track` (typically the
        host halo's main progenitor), decomposed into radial and tangential
        velocity components, over the snapshots the two tracks share.

        Parameters
        ----------
        radius_key : str
            Key into `reference_track.extra` holding the reference's virial
            (or other infall) radius at each snapshot. Defaults to
            "GroupR200", which SubFind-HBT tracks carry automatically. If
            absent, the returned OrbitAnalysis.Rvir is None and
            first_crossing() will return None -- pass your own radius_key
            (and make sure it's on the reference track's `extra`) for other
            tree formats.
        boxsize : float, optional
            Defaults to self.BoxSize (available for SubFind-HBT trees) for
            periodic wrapping of the separation vector.
        """
        if boxsize is None:
            boxsize = getattr(self, "BoxSize", None)

        snaps_common, i_self, i_ref = np.intersect1d(
            track.SnapNum, reference_track.SnapNum, return_indices=True)
        if snaps_common.size == 0:
            raise MergerTreeError("Track and reference track share no common snapshots.")

        dpos = periodic_delta(track.Pos[i_self] - reference_track.Pos[i_ref], boxsize)
        dvel = track.Vel[i_self] - reference_track.Vel[i_ref]

        distance = np.linalg.norm(dpos, axis=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            r_hat = np.where(distance[:, None] > 0, dpos / distance[:, None], 0.0)
        radial_vel = np.sum(dvel * r_hat, axis=1)
        tangential_vec = dvel - radial_vel[:, None] * r_hat
        tangential_vel = np.linalg.norm(tangential_vec, axis=1)
        rel_speed = np.linalg.norm(dvel, axis=1)

        rvir = reference_track.extra.get(radius_key)
        rvir_common = rvir[i_ref] if rvir is not None else None

        return OrbitAnalysis(
            halo_id=track.halo_id, host_id=reference_track.halo_id,
            SnapNum=snaps_common, Redshift=track.Redshift[i_self], Time=track.Time[i_self],
            Distance=distance, RelVel=rel_speed,
            RadialVelocity=radial_vel, TangentialVelocity=tangential_vel,
            Mass=track.Mass[i_self], Rvir=rvir_common,
        )

    # ------------------------------------------------------------------
    # Plotting
    # ------------------------------------------------------------------

    def plot_track(self, track: HaloTrack, quantities: Tuple[str, ...] = ("Mass", "Speed", "Position"),
                    use_redshift: bool = False, mark_infall: bool = True,
                    figsize: Tuple[float, float] = (7, 8)) -> "plt.Figure":
        """
        Plot a HaloTrack in isolation: mass, speed, and position components
        vs. time (or redshift), stacked as subplots. `quantities` can
        include "Mass", "Speed", "Position", "Vmax", "VmaxRad", or any
        key in track.extra.
        """
        x = 1.0 / track.Time - 1.0 if use_redshift else track.Time
        xlabel = "Redshift" if use_redshift else "Expansion Factor"

        fig, axes = plt.subplots(len(quantities), 1, figsize=figsize, sharex=True)
        if len(quantities) == 1:
            axes = [axes]

        infall = track.infall_snapshot() if mark_infall else None
        x_infall = (1.0 / infall["time"] - 1.0) if (infall and use_redshift) else \
                   (infall["time"] if infall else None)

        for ax, q in zip(axes, quantities):
            if q == "Mass":
                ax.plot(x, track.Mass, color="black")
                ax.set_yscale("log")
                ax.set_ylabel("Mass")
            elif q == "Speed":
                speed = np.linalg.norm(track.Vel, axis=1)
                ax.plot(x, speed, color="darkred")
                ax.set_ylabel("|Velocity|")
            elif q == "Position":
                for i, label in enumerate(("x", "y", "z")):
                    ax.plot(x, track.Pos[:, i], label=label)
                ax.set_ylabel("Position")
                ax.legend(fontsize=8, ncol=3)
            else:
                ax.plot(x, track.get(q), color="steelblue")
                ax.set_ylabel(q)

            if x_infall is not None:
                ax.axvline(x_infall, color="grey", ls=":", lw=1.5)

        axes[-1].set_xlabel(xlabel)
        if not use_redshift:
            axes[-1].set_xlim(0, 1)
        fig.suptitle(f"Halo {track.halo_id} ({track.treefileformat})")
        fig.tight_layout()
        return fig

    def plot_relative(self, track: HaloTrack, reference_track: HaloTrack,
                       boxsize: Optional[float] = None, use_redshift: bool = False,
                       mark_infall: bool = True,
                       figsize: Tuple[float, float] = (7, 6)) -> "plt.Figure":
        """
        Plot a HaloTrack's separation and relative velocity from a reference
        (e.g. host/central) track, over the snapshots the two tracks share.
        Positions are periodically wrapped using `boxsize` (defaults to
        self.BoxSize if available, from a SubFind-HBT tree).
        """
        if boxsize is None:
            boxsize = getattr(self, "BoxSize", None)

        snaps_common, i_self, i_ref = np.intersect1d(
            track.SnapNum, reference_track.SnapNum, return_indices=True)
        if snaps_common.size == 0:
            raise MergerTreeError("Track and reference track share no common snapshots.")

        dpos = track.Pos[i_self] - reference_track.Pos[i_ref]
        dpos = periodic_delta(dpos, boxsize)
        dvel = track.Vel[i_self] - reference_track.Vel[i_ref]

        distance = np.linalg.norm(dpos, axis=1)
        rel_speed = np.linalg.norm(dvel, axis=1)
        time = track.Time[i_self]
        x = 1.0 / time - 1.0 if use_redshift else time
        xlabel = "Redshift" if use_redshift else "Expansion Factor"

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize, sharex=True)
        ax1.plot(x, distance, color="black")
        ax1.set_ylabel("Separation")
        ax2.plot(x, rel_speed, color="darkred")
        ax2.set_ylabel("|Relative Velocity|")
        ax2.set_xlabel(xlabel)

        if mark_infall:
            infall = track.infall_snapshot()
            if infall is not None:
                x_infall = (1.0 / infall["time"] - 1.0) if use_redshift else infall["time"]
                for ax in (ax1, ax2):
                    ax.axvline(x_infall, color="grey", ls=":", lw=1.5)
                ax1.annotate(f"infall  z={infall['redshift']:.2f}\n"
                             f"M={infall['mass']:.3g}",
                             xy=(x_infall, distance[np.searchsorted(x if not use_redshift else -x,
                                                                     x_infall if not use_redshift else -x_infall)]
                                 if False else distance.max()),
                             xytext=(5, 0), textcoords="offset points", fontsize=8, color="grey")

        if not use_redshift:
            ax2.set_xlim(0, 1)
        fig.suptitle(f"Halo {track.halo_id} relative to {reference_track.halo_id}")
        fig.tight_layout()
        return fig

    # ------------------------------------------------------------------

    def plot_orbits(self, orbits: Dict[Any, "OrbitAnalysis"], reference_track: HaloTrack,
                     use_redshift: bool = False, normalize_by_rvir: bool = False,
                     mark_crossing: bool = True, figsize: Tuple[float, float] = (7, 5)) -> "plt.Figure":
        """
        Overlay several subhalos' orbits (distance from host vs. time) on one
        axes, with the host's virial radius drawn as a reference curve and
        each subhalo's first virial-radius crossing marked.

        Parameters
        ----------
        orbits : dict
            {subhalo_id: OrbitAnalysis}, e.g. from
            {sid: mt.analyse_orbit(track, host_track) for sid, track in
             mt.get_group_subhalo_tracks(group_id, snapnum).items()}
        reference_track : HaloTrack
            The host track the orbits were computed relative to (used here
            only to draw Rvir(t); pass the same one used in analyse_orbit).
        normalize_by_rvir : bool
            If True, plot Distance / Rvir(t) instead of raw Distance (so the
            host's virial radius becomes the horizontal line y=1).
        """
        x_ref = 1.0 / reference_track.Time - 1.0 if use_redshift else reference_track.Time
        xlabel = "Redshift" if use_redshift else "Expansion Factor"

        fig, ax = plt.subplots(figsize=figsize)
        rvir = reference_track.extra.get("GroupR200")

        if not normalize_by_rvir and rvir is not None:
            ax.plot(x_ref, rvir, color="black", lw=2, label=r"Host $R_{200}$", zorder=5)
        elif normalize_by_rvir:
            ax.axhline(1.0, color="black", lw=2, label=r"Host $R_{200}$", zorder=5)

        cmap = plt.get_cmap("tab10")
        for i, (sub_id, orb) in enumerate(orbits.items()):
            x = 1.0 / orb.Time - 1.0 if use_redshift else orb.Time
            y = orb.Distance / orb.Rvir if (normalize_by_rvir and orb.Rvir is not None) else orb.Distance
            color = cmap(i % 10)
            ax.plot(x, y, color=color, lw=1.2, label=f"Subhalo {sub_id}")

            if mark_crossing:
                crossing = orb.first_crossing()
                if crossing is not None:
                    x_c = (1.0 / crossing["time"] - 1.0) if use_redshift else crossing["time"]
                    y_c = 1.0 if normalize_by_rvir else crossing["distance"]
                    ax.scatter([x_c], [y_c], color=color, edgecolor="black", zorder=6, s=40)

        ax.set_xlabel(xlabel)
        ax.set_ylabel(r"Distance / $R_{200}$" if normalize_by_rvir else "Distance from host")
        if not use_redshift:
            ax.set_xlim(0, 1)
        ax.legend(fontsize=8, ncol=2)
        fig.suptitle(f"Subhalo orbits about host {reference_track.halo_id}")
        fig.tight_layout()
        return fig

    def summary(self) -> None:
        print("Merger Tree Summary")
        print(f"  Format: {self.treefileformat}")
        print(f"  File:   {os.path.basename(self.treefilename)}")
        for k, v in self.metadata.items():
            print(f"  {k:20s}: {v}")
