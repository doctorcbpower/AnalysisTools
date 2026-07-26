#!/usr/bin/env python3
"""
analysistools.api.trees
-----------------------
MergerTree: unified access to merger trees, and TrackDataset: the Dataset
view of a single halo's main-branch history (a HaloTrack).

A tree file is not itself a flat table, so MergerTree is a thin handle
whose queries return TrackDataset objects; those follow the full Dataset
contract (fields, selection, common names: mass/pos/vel/redshift), so a
track can be fed to anything that accepts a Dataset.
"""
from __future__ import annotations

from typing import Dict, Optional, Union

import numpy as np

from ..merger_tree_tools import MergerTreeTools
from ..merger_tree_types import HaloTrack
from .dataset import Dataset
from .halos import HaloCatalogue


class TrackDataset(Dataset):
    """Dataset view of one HaloTrack (rows = snapshots, earliest first)."""

    kind = "tree"

    def __init__(self, track: HaloTrack, label: Optional[str] = None):
        super().__init__(path="", fileformat=track.treefileformat,
                         label=label or f"track:{track.halo_id}")
        self.track = track                     # underlying HaloTrack
        self._columns = {k: v for k, v in track.to_dict().items()
                         if isinstance(v, np.ndarray)}
        self.meta.update({
            "halo_id": track.halo_id,
            "query_snapnum": track.query_snapnum,
            "redshift": float(track.Redshift[-1]) if len(track) else None,
            "is_currently_subhalo": track.is_currently_subhalo,
        })
        self._loaded = True

    def _load(self) -> None:      # data injected in __init__
        pass

    def infall(self) -> Optional[dict]:
        """First central->satellite transition (HaloTrack.infall_snapshot)."""
        return self.track.infall_snapshot()


def _to_halo_tools(obj):
    """Accept a HaloCatalogue adapter or a raw HaloTools; return HaloTools."""
    if isinstance(obj, HaloCatalogue):
        obj.preload()                # ensures standardised tables exist
        return obj.backend
    return obj


class MergerTree:
    """
    Unified handle on a merger tree file.

    Parameters
    ----------
    path : str
        Tree file path.
    fileformat : str
        "TreeFrog" (full or walkable flavour), "SubFind", "MergerTree" (AHF).
    halos : optional
        Linked halo catalogue(s) used to resolve IDs and fill track
        properties: a HaloCatalogue adapter, a HaloTools instance, or a
        dict {snapnum: either} (walkable TreeFrog trees store topology
        only, so properties come from these catalogues).
    comoving : bool, optional
    label : str, optional

    Examples
    --------
    >>> tree = at.load("VELOCIraptor.walkabletree.hdf5", halos=cat)
    >>> tr = tree.track(halo_id, snapnum=31)      # -> TrackDataset
    >>> tr["mass"], tr["redshift"]
    >>> tr.infall()
    >>> tr2 = tree.from_halo(cat, index=0)        # via catalogue row
    """

    kind = "tree"

    def __init__(self, path: str, fileformat: str = "TreeFrog",
                 halos=None, comoving: bool = False,
                 label: Optional[str] = None, **backend_kwargs):
        self.path = path
        self.fileformat = fileformat
        self.label = label or path.rsplit("/", 1)[-1]
        self._comoving = comoving
        self._backend_kwargs = backend_kwargs
        self._halos = halos
        self._backend: Optional[MergerTreeTools] = None

    # ------------------------------------------------------------------

    def _halo_link(self):
        if self._halos is None:
            return None
        if isinstance(self._halos, dict):
            return {int(k): _to_halo_tools(v) for k, v in self._halos.items()}
        return _to_halo_tools(self._halos)

    @property
    def backend(self) -> MergerTreeTools:
        """The underlying MergerTreeTools (reads the file on first use)."""
        if self._backend is None:
            kwargs = dict(self._backend_kwargs)
            kwargs.setdefault("loglevel", 30)
            self._backend = MergerTreeTools(
                self.path, treefileformat=self.fileformat,
                comoving_units=self._comoving,
                halo_tools=self._halo_link(), **kwargs)
        return self._backend

    def preload(self) -> "MergerTree":
        self.backend
        return self

    # ------------------------------------------------------------------
    # Queries -> TrackDataset
    # ------------------------------------------------------------------

    def track(self, halo_id: int, snapnum: int,
              object_type: str = "Subhalo") -> TrackDataset:
        """Main-branch history of (halo_id, snapnum) as a TrackDataset."""
        tr = self.backend.get_track(int(halo_id), int(snapnum),
                                    object_type=object_type)
        return TrackDataset(tr, label=f"{self.label}:{halo_id}")

    def from_halo(self, catalogue, index: int,
                  object_type: str = "Group",
                  snapnum: Optional[int] = None, **kwargs) -> TrackDataset:
        """
        Track of the halo at row `index` of a catalogue (HaloCatalogue
        adapter or HaloTools instance).

        Note: for VELOCIraptor/TreeFrog, catalogue IDs are per-snapshot
        (1..N) while walkable-tree IDs are temporal
        (snapnum*Temporal_halo_id_value + N); the temporal conversion is
        applied automatically when needed.
        """
        ht = _to_halo_tools(catalogue)
        if snapnum is None:
            snapnum = getattr(ht, "snapnum", None)
            if snapnum is None and isinstance(catalogue, HaloCatalogue):
                snapnum = catalogue.meta.get("snapnum")
        if snapnum is None:
            raise ValueError("from_halo(): pass snapnum=, or load the "
                             "catalogue with snapnum= so it is recorded.")

        table = ht.standardised_halos if object_type == "Group" \
            else ht.standardised_subhalos
        halo_id = int(table["halo_id"][index])

        # temporal-ID conversion for walkable TreeFrog trees
        data = self.backend.data
        lookup = getattr(data, "lookup", None)
        if lookup is not None and (int(snapnum), halo_id) not in lookup:
            tid_val = int(getattr(data, "metadata", {}).get(
                "Temporal_halo_id_value", 0) or 0)
            if tid_val:
                cand = int(snapnum) * tid_val + halo_id
                if (int(snapnum), cand) in lookup:
                    halo_id = cand

        tr = self.backend.get_track(halo_id, int(snapnum),
                                    object_type="Subhalo", **kwargs)
        return TrackDataset(tr, label=f"{self.label}:{halo_id}")

    # ------------------------------------------------------------------

    def find_infall(self, track) -> Optional[dict]:
        tr = track.track if isinstance(track, TrackDataset) else track
        return self.backend.find_infall(tr)

    def summary(self) -> None:
        b = self.backend
        print(f"MergerTree ({self.fileformat})")
        print(f"  label: {self.label}")
        print(f"  path:  {self.path}")
        lookup = getattr(b.data, "lookup", None)
        if lookup is not None:
            print(f"  nodes: {len(lookup):,}")
        if b.metadata:
            for k in ("NSnaps", "Total_number_of_halos", "Name", "Title"):
                if k in b.metadata:
                    print(f"  {k}: {b.metadata[k]}")

    def __repr__(self) -> str:
        return (f"<MergerTree '{self.label}' format={self.fileformat}"
                f"{' loaded' if self._backend is not None else ''}>")
