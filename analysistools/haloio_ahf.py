#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
haloio_ahf.py

Reader for AHF-format halo catalogues (text-based).
"""

import os
import numpy as np
from typing import Dict, Any, Tuple


def read_ahf(filename: str) -> Tuple[Dict[str, Any], Dict[str, np.ndarray], Dict[str, np.ndarray]]:
    """Read an AHF halo catalogue.

    Returns raw values exactly as stored in the file -- comoving/little-h
    unit conversion is applied centrally by HaloTools.standardise_names(),
    not here (see its comoving/little_h parameters). Note this reader
    doesn't currently parse HubbleParam/redshift from the file, so that
    conversion can't do anything for AHF yet regardless of what's requested.

    Parameters
    ----------
    filename : str
        Path to the `.AHF_halos` file.

    Returns
    -------
    metadata : dict
    halos : dict[str, np.ndarray]
    subhalos : dict[str, np.ndarray]
    """
    if not filename.endswith(".AHF_halos"):
        raise ValueError(
            f"read_ahf() expects a '.AHF_halos' file, got '{filename}'. "
            "(The '.AHF_substructure' sibling file, if present, is picked "
            "up automatically.)"
        )

    metadata: Dict[str, Any] = {}
    halos: Dict[str, np.ndarray] = {}
    subhalos: Dict[str, np.ndarray] = {}

    # Read main halos. ndmin=2 forces a 2-D array even for a single-row
    # file (genfromtxt otherwise collapses it to 1-D, breaking the [:, col]
    # indexing below).
    data = np.genfromtxt(filename, comments="#", ndmin=2)
    halos["HaloID"] = data[:, 0].astype(np.int64)
    halos["HostHaloID"] = data[:, 1].astype(np.int64)
    halos["NumSubStruct"] = data[:, 2].astype(np.int64)
    halos["Mass"] = data[:, 3].astype(np.float64)
    halos["NumPart"] = data[:, 4].astype(np.int64)
    halos["Pos"] = data[:, 5:8].astype(np.float64)
    halos["Vel"] = data[:, 8:11].astype(np.float64)
    halos["Rvir"] = data[:, 11].astype(np.float64)
    halos["com_offset"] = data[:,15].astype(np.float64)
    halos["Vmax"] = data[:,16].astype(np.float64)
    metadata["Nhalos"] = len(data)

    # Read optional substructure file
    subfile = filename.replace(".AHF_halos", ".AHF_substructure")
    if os.path.exists(subfile):
        subdata = np.genfromtxt(subfile, comments="#", ndmin=2)
        subhalos["ParentHaloID"] = subdata[:, 0].astype(np.int64)
        subhalos["SubhaloMass"] = subdata[:, 1]
        subhalos["Pos"] = subdata[:, 5:8]
        subhalos["Vel"] = subdata[:, 8:11]
        subhalos["Rvir"] = subdata[:, 11]
        metadata["Nsubhalos"] = len(subdata)

    return metadata, halos, subhalos

