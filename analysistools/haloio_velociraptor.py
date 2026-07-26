#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
haloio_velociraptor.py

Reader for VELOCIraptor-format halo catalogues (HDF5).
"""

import h5py
import numpy as np
from typing import Dict, Any, Tuple


def read_velociraptor(filename: str, comoving: bool = False) -> Tuple[Dict[str, Any], Dict[str, np.ndarray], Dict[str, np.ndarray]]:
    """Read a VELOCIraptor-format halo catalogue (HDF5).

    Parameters
    ----------
    filename : str
        Path to the halo catalogue file.
    comoving : bool, optional
        Convert positions to comoving coordinates.

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

    if comoving:
        scale = metadata.get("ScaleFactor", 1.0)
        for arrname in ["Xc", "Yc", "Zc"]:
            if arrname in halos:
                halos[arrname] *= scale
            if arrname in subhalos:
                subhalos[arrname] *= scale

    return metadata, halos, subhalos

