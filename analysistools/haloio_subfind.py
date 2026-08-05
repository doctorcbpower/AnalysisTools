#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
haloio_subfind.py

Reader for SubFind-format halo catalogues (HDF5).
"""

import h5py
import numpy as np
from typing import Dict, Any, Tuple


def read_subfind(filename: str) -> Tuple[Dict[str, Any], Dict[str, np.ndarray], Dict[str, np.ndarray]]:
    """Read a SubFind-format halo catalogue (HDF5).

    Returns raw values exactly as stored in the file -- comoving/little-h
    unit conversion is applied centrally by HaloTools.standardise_names(),
    not here (see its comoving/little_h parameters).

    Parameters
    ----------
    filename : str
        Path to the SubFind halo catalogue file.

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
        params = dict(f["Parameters"].attrs)

        metadata = {
            "BoxSize": header.get("BoxSize"),
            "NumFiles": header.get("NumFiles", 1),
            "HubbleParam": params.get("HubbleParam"),
            "OmegaDM": params.get("Omega_cdm"),
            "OmegaB": params.get("Omega_b"),
            "OmegaLambda": params.get("OmegaLambda"),
            "Redshift": params.get("Redshift"), 
            "TotNgroups": header.get("Ngroups_Total"),
            "TotNsubgroups": header.get("Nsubgroups_Total", header.get("Nsubhalos_Total")),
        }

        halos = {key: f["Group/" + key][()] for key in f["Group"].keys() if isinstance(f["Group/" + key], h5py.Dataset)}

        subhalos = {}
        if "Subhalo" in f:
            subhalos = {
                key: f["Subhalo/" + key][()] for key in f["Subhalo"].keys()
                if isinstance(f["Subhalo/" + key], h5py.Dataset)
            }

    return metadata, halos, subhalos

