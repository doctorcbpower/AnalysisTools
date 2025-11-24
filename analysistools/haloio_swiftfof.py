#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
haloio_swiftfof.py

Reader for SWIFT-format FOF catalogues (HDF5).
"""

import h5py
import numpy as np
from typing import Dict, Any, Tuple


def read_swiftfof(filename: str, comoving: bool = False) -> Tuple[Dict[str, Any], Dict[str, np.ndarray], Dict[str, np.ndarray]]:
    """Read a SWIFT-format FOF catalogue (HDF5).

    Parameters
    ----------
    filename : str
        Path to the SWIFT FOF catalogue file.
    comoving : bool, optional
        If True, convert positions to comoving coordinates.

    Returns
    -------
    metadata : dict
        Header and parameter attributes.
    halos : dict[str, np.ndarray]
        Group-level properties.
    subhalos : dict[str, np.ndarray]
        Subhalo-level properties (if present).
    """
    with h5py.File(filename, "r") as f:
        header = dict(f["Header"].attrs)
        cosmo = dict(f["Cosmology"].attrs)

        metadata = {
            "BoxSize": header.get("BoxSize"),
            "NumFiles": header.get("NumFilesPerSnapshot", 1),
            "HubbleParam": cosmo.get("h"),
            "TotNgroups": header.get("NumGroups_Total"),
        }

        halos = {key: f["Groups/" + key][()] for key in f["Groups"].keys() if isinstance(f["Groups/" + key], h5py.Dataset)}
        subhalos = {}
        
    # Optional comoving conversion
    if comoving and "GroupPos" in halos and metadata.get("HubbleParam"):
        hubble = metadata["HubbleParam"]
        halos["GroupPos"] /= hubble

    return metadata, halos, subhalos

