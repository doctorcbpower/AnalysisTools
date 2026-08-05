#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
haloio_velociraptor.py

Reader for VELOCIraptor-format halo catalogues (HDF5).
"""

import os
import h5py
import numpy as np
from typing import Dict, Any, Tuple


def _read_siminfo(filename: str) -> Dict[str, Any]:
    """Best-effort parse of a VELOCIraptor '<base>.siminfo' sidecar file
    (one 'key : value # dtype' pair per line, e.g. 'h_val : 0.703000 #
    float64'), if one exists alongside `filename`. Not every VELOCIraptor
    run produces this file, so an absent one is not an error -- just
    returns {}."""
    base = filename.split(".properties")[0]
    siminfo_path = base + ".siminfo"
    if not os.path.exists(siminfo_path):
        return {}

    info: Dict[str, Any] = {}
    with open(siminfo_path) as f:
        for line in f:
            line = line.split("#", 1)[0].strip()
            if not line or ":" not in line:
                continue
            key, _, value = line.partition(":")
            key, value = key.strip(), value.strip()
            try:
                info[key] = float(value)
            except ValueError:
                info[key] = value
    return info


def read_velociraptor(filename: str) -> Tuple[Dict[str, Any], Dict[str, np.ndarray], Dict[str, np.ndarray]]:
    """Read a VELOCIraptor-format halo catalogue (HDF5).

    Returns raw values exactly as stored in the file -- comoving/little-h
    unit conversion is applied centrally by HaloTools.standardise_names(),
    not here (see its comoving/little_h parameters). HubbleParam isn't part
    of the main HDF5 file's own attrs; if a '<base>.siminfo' sidecar file is
    present alongside `filename` (common but not universal -- depends on how
    VELOCIraptor was run), its 'h_val' is used. Without one, metadata has no
    HubbleParam and the little_h conversion silently can't do anything for
    this catalogue; the comoving (scale-factor) conversion still works via
    ScaleFactor either way.

    Parameters
    ----------
    filename : str
        Path to the halo catalogue file.

    Returns
    -------
    metadata : dict
    halos : dict[str, np.ndarray]
    subhalos : dict[str, np.ndarray]
    """
    with h5py.File(filename, "r") as f:
        if "Groups" in f and "MetaData" in f:
            # Grouped layout (MetaData/Groups[/SubGroups])
            header = dict(f["MetaData"].attrs)

            metadata = {
                "BoxSize": header.get("BoxSizeComovingKpch"),
                "ScaleFactor": header.get("ScaleFactor", 1.0),
                "TotNgroups": header.get("Num_of_groups", 0),
            }

            halos = {
                key: f["Groups/" + key][()] for key in f["Groups"].keys()
                if isinstance(f["Groups/" + key], h5py.Dataset)
            }

            subhalos = {}
            if "SubGroups" in f:
                subhalos = {
                    key: f["SubGroups/" + key][()] for key in f["SubGroups"].keys()
                    if isinstance(f["SubGroups/" + key], h5py.Dataset)
                }
        else:
            # Flat layout: standard VELOCIraptor .properties output --
            # per-halo datasets at file root, metadata in root attributes.
            # Field halos and subhalos share one table; subhalos are the
            # rows with hostHaloID != -1.
            header = dict(f.attrs)

            ntot = f["Total_num_of_groups"][0] if "Total_num_of_groups" in f \
                else (f["Num_of_groups"][0] if "Num_of_groups" in f else 0)
            metadata = {
                "BoxSize": header.get("Period"),
                "ScaleFactor": header.get("Time", 1.0),
                "TotNgroups": ntot,
            }

            # Per-halo datasets only (skip bookkeeping scalars and non-datasets)
            skip = {"Num_of_groups", "Total_num_of_groups", "File_id",
                    "Num_of_files", "Configuration", "SimulationInfo",
                    "UnitInfo"}
            table = {
                key: f[key][()] for key in f.keys()
                if key not in skip and isinstance(f[key], h5py.Dataset)
            }
            # Guard against stray scalar datasets
            n = len(table["ID"]) if "ID" in table else max(
                (len(v) for v in table.values() if hasattr(v, "__len__")),
                default=0)
            table = {k: v for k, v in table.items()
                     if hasattr(v, "__len__") and len(v) == n}

            halos = table
            subhalos = {}
            if "hostHaloID" in table:
                is_sub = table["hostHaloID"] != -1
                if is_sub.any():
                    subhalos = {k: v[is_sub] for k, v in table.items()}

    siminfo = _read_siminfo(filename)
    if siminfo.get("h_val") is not None:
        metadata["HubbleParam"] = siminfo["h_val"]

    return metadata, halos, subhalos

