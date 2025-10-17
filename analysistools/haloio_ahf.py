#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
haloio_ahf.py

Reader for AHF-format halo catalogues (text-based).
"""

import os
import numpy as np
from typing import Dict, Any, Tuple


def read_ahf(filename: str, comoving: bool = False) -> Tuple[Dict[str, Any], Dict[str, np.ndarray], Dict[str, np.ndarray]]:
    """Read an AHF halo catalogue.

    Parameters
    ----------
    filename : str
        Path to the `.AHF_halos` file.
    comoving : bool, optional
        Convert positions to comoving coordinates.

    Returns
    -------
    metadata : dict
    halos : dict[str, np.ndarray]
    subhalos : dict[str, np.ndarray]
    """
    metadata: Dict[str, Any] = {}
    halos: Dict[str, np.ndarray] = {}
    subhalos: Dict[str, np.ndarray] = {}

    # Read main halos
    if filename.endswith(".AHF_halos"):
        data = np.genfromtxt(filename, comments="#")
        halos["HaloID"] = data[:, 0].astype(np.int64)
        halos["HostHaloID"] = data[:, 1].astype(np.int64)
        halos["NumSubStruct"] = data[:, 2].astype(np.int64)
        halos["Mass"] = data[:, 3].astype(np.float64)
        halos["NumPart"] = data[:, 4].astype(np.int64)
        halos["Pos"] = data[:, 5:8].astype(np.int32)
        halos["Vel"] = data[:, 8:11].astype(np.int32)
        halos["Rvir"] = data[:, 11].astype(np.float64)
        halos["com_offset"] = data[:,15].astype(np.float64)
        halos["Vmax"] = data[:,16].astype(np.float64)
        metadata["Nhalos"] = len(data)

    # Read optional substructure file
    subfile = filename.replace(".AHF_halos", ".AHF_substructure")
    if os.path.exists(subfile):
        subdata = np.genfromtxt(subfile, comments="#")
        subhalos["ParentHaloID"] = subdata[:, 0].astype(np.int64)
        subhalos["SubhaloMass"] = subdata[:, 1]
        subhalos["Pos"] = subdata[:, 5:8]
        subhalos["Vel"] = subdata[:, 8:11]
        subhalos["Rvir"] = subdata[:, 11]
        metadata["Nsubhalos"] = len(subdata)

    if comoving and "Pos" in halos:
        halos["Pos"] *= 1.0  # placeholder for unit conversion

    return metadata, halos, subhalos

