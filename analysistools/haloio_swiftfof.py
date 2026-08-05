#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
haloio_swiftfof.py

Reader for SWIFT-format FOF catalogues (HDF5).
"""

import h5py
import numpy as np
from typing import Dict, Any, Tuple


def read_swiftfof(filename: str) -> Tuple[Dict[str, Any], Dict[str, np.ndarray], Dict[str, np.ndarray]]:
    """Read a SWIFT-format FOF catalogue (HDF5).

    Returns raw values exactly as stored in the file -- comoving/little-h
    unit conversion is applied centrally by HaloTools.standardise_names(),
    not here (see its comoving/little_h parameters). Note SWIFT's own native
    convention already excludes h (lengths in Mpc, not Mpc/h).

    Parameters
    ----------
    filename : str
        Path to the SWIFT FOF catalogue file.

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
            "OmegaDM": cosmo.get("Omega_cdm"),
            "OmegaB": cosmo.get("Omega_b"),
            "OmegaLambda": cosmo.get("Omega_lambda"),
            "Redshift": cosmo.get("Redshift"), 
            "TotNgroups": header.get("NumGroups_Total"),
        }

        halos = {key: f["Groups/" + key][()] for key in f["Groups"].keys() if isinstance(f["Groups/" + key], h5py.Dataset)}
        subhalos = {}
        
    return metadata, halos, subhalos

