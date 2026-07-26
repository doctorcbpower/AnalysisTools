#!/usr/bin/env python3
"""
analysistools.api
-----------------
The unified interface (DEVELOPMENT.md section 3): ``load()`` opens any
supported data product and returns a Dataset with common syntax.

Phase 1 covers snapshots and halo catalogues. Merger trees, SHARK galaxy
catalogues, the Simulation umbrella, and the plotting layer follow in
Phases 2-4.
"""
from __future__ import annotations

import os
from typing import Optional

from .dataset import Dataset, FIELD_ALIASES
from .snapshot import SnapshotDataset
from .halos import HaloCatalogue

__all__ = ["load", "Dataset", "SnapshotDataset", "HaloCatalogue",
           "FIELD_ALIASES"]

#: kinds load() understands today; "tree", "galaxies", "field" arrive in
#: later phases (see DEVELOPMENT.md section 5).
_SUPPORTED_KINDS = ("snapshot", "halos")


def _sniff_hdf5(path: str):
    """Return (kind, fileformat, convention) guessed from HDF5 contents."""
    import h5py
    with h5py.File(path, "r") as f:
        keys = set(f.keys())

        # --- snapshot: PartTypeN groups ---
        if any(k.startswith("PartType") for k in keys):
            if "Cosmology" in keys:                     # SWIFT
                return "snapshot", "HDF5", "SWIFT"
            if "Parameters" in keys and "Omega0" in f["Parameters"].attrs:
                return "snapshot", "HDF5", "GADGET4"    # GADGET-4 / AREPO
            return "snapshot", "HDF5", "GADGET2/3"

        # --- SHARK galaxies (Phase 2) ---
        if "galaxies" in keys:
            return "galaxies", "SHARK", None

        # --- SubFind group catalogue ---
        if "Group" in keys or "Subhalo" in keys:
            return "halos", "SUBFIND", None

        # --- SWIFT FOF ---
        if "Groups" in keys and "MetaData" not in keys:
            grp = f["Groups"]
            if "Centres" in grp or "Masses" in grp:
                return "halos", "SWIFT_FOF", None

        # --- VELOCIraptor: grouped or flat .properties layout ---
        if "Groups" in keys and "MetaData" in keys:
            return "halos", "VELOCIraptor", None
        if {"Xc", "Yc", "Zc"} <= keys or "Total_num_of_groups" in keys:
            return "halos", "VELOCIraptor", None

    return None, None, None


def load(path: str, kind: Optional[str] = None,
         fileformat: Optional[str] = None,
         label: Optional[str] = None, **kwargs) -> Dataset:
    """
    Open a simulation data product with a single call.

    The data kind and on-disk format are sniffed from the file where
    possible; pass ``kind=`` / ``fileformat=`` to override. No data is read
    until the first field access (lazy).

    Parameters
    ----------
    path : str
        File to open.
    kind : {"snapshot", "halos"}, optional
        Data kind ("tree"/"galaxies"/"field" arrive in later phases).
    fileformat : str, optional
        Backend format ("HDF5", "SNAP2", "SUBFIND", "VELOCIraptor", ...).
    label : str, optional
        Name for plot legends.
    **kwargs :
        Passed to the adapter (e.g. convention=, comoving=, extra_blocks=).

    Examples
    --------
    >>> import analysistools as at
    >>> snap  = at.load("data/snap_0031.hdf5")               # sniffed
    >>> halos = at.load("data/halos/snap_0031.VELOCIraptor.properties.0")
    >>> snap.dm["pos"], halos["mass"]
    """
    if not os.path.exists(path):
        raise FileNotFoundError(path)

    sniffed_kind = sniffed_fmt = sniffed_conv = None
    try:
        import h5py
        if h5py.is_hdf5(path):
            sniffed_kind, sniffed_fmt, sniffed_conv = _sniff_hdf5(path)
    except (ImportError, OSError):
        pass

    # Non-HDF5 fallbacks by filename convention
    if sniffed_kind is None:
        name = os.path.basename(path)
        if "AHF_halos" in name:
            sniffed_kind, sniffed_fmt = "halos", "AHF"

    kind = kind or sniffed_kind
    fileformat = fileformat or sniffed_fmt

    if kind is None:
        raise ValueError(
            f"Could not determine the data kind of '{path}'. "
            f"Pass kind= explicitly (one of {_SUPPORTED_KINDS}).")

    if kind == "snapshot":
        kwargs.setdefault("convention", sniffed_conv)
        return SnapshotDataset(path, fileformat=fileformat or "HDF5",
                               label=label, **kwargs)
    if kind == "halos":
        if fileformat is None:
            raise ValueError("load(): fileformat is required for halo "
                             "catalogues that cannot be sniffed.")
        return HaloCatalogue(path, fileformat=fileformat, label=label,
                             **kwargs)
    if kind in ("tree", "galaxies", "field"):
        raise NotImplementedError(
            f"kind='{kind}' is planned for a later phase (see "
            f"DEVELOPMENT.md section 5); use the existing "
            f"MergerTreeTools / SharkModel / FDMFieldTools APIs meanwhile.")
    raise ValueError(f"Unknown kind '{kind}' (use one of {_SUPPORTED_KINDS}).")
