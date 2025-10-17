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

    if comoving:
        scale = metadata.get("ScaleFactor", 1.0)
        for arrname in ["Xc", "Yc", "Zc"]:
            if arrname in halos:
                halos[arrname] *= scale
            if arrname in subhalos:
                subhalos[arrname] *= scale

    return metadata, halos, subhalos

