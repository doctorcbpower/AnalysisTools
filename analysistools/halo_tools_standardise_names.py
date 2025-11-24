#!/usr/bin/env python3
"""
halo_tools_unifier.py
Version: 2025-11-12

Provides schema harmonisation across different halo catalogue formats
(AHF, SUBFIND, VELOCIraptor, SWIFTFOF, etc.).

Each reader returns native fields with their own conventions; this module
maps them into a common "standardised" schema.

Author: C. Power
"""

import numpy as np
import logging
from typing import Dict, Optional

STANDARD_KEYS = [
    "mass",     # total halo mass (M200, Mvir, etc.)
    "pos",      # 3D position [x,y,z]
    "vel",      # 3D velocity [vx,vy,vz]
    "radius",   # virial or FOF radius
    "halo_id",       # halo ID
    "num_part",       # number of particles
]

CATALOGUE_MAPPINGS: Dict[str, Dict[str, str]] = {
    "SUBFIND": {
        "mass": "Group_M_Crit200",
        "pos": "GroupPos",
        "vel": "GroupVel",
        "radius": "Group_R_Crit200",
        "id": "GroupNr",
        "np": "GroupLen",
    },
    "AHF": {
        "mass": "Mvir",
        "pos": "Xc",   # can be split from Xc, Yc, Zc if stored separately
        "vel": "Vcm",  # likewise Vx, Vy, Vz may need stacking
        "radius": "Rvir",
        "id": "ID",
        "np": "npart",
    },
    "VELOCIraptor": {
        "mass": "Mass_200crit",
        "pos": "CentreOfPotential",
        "vel": "Velocity",
        "radius": "R_200crit",
        "id": "ID",
        "np": "npart",
    },
    "SWIFT_FOF": {
        "mass": "Masses",
        "pos": "Centres",
        "vel": "",
        "radius": "Radii",
        "id": "GroupIDs",
        "np": "Sizes",
    },
}


def standardise_catalogue_names(halos: Dict[str, np.ndarray],
                                halotype: str,
                                logger: Optional[logging.Logger] = None) -> Dict[str, np.ndarray]:
    """
    Convert halo dictionaries to common naming conventions.

    Parameters
    ----------
    halos : dict
        Raw halo data as returned by a reader (e.g. read_subfind).
    halotype : str
        One of "SUBFIND", "AHF", "VELOCIraptor", "SWIFT_FOF".
    logger : logging.Logger, optional
        Logger for messages.

    Returns
    -------
    dict
        Dictionary containing keys from STANDARD_KEYS.
    """
    if halos is None:
        raise ValueError("No halo data provided for normalisation.")

    halotype = halotype.upper()
    mapping = CATALOGUE_MAPPINGS.get(halotype)
    if mapping is None:
        raise ValueError(f"No normalisation mapping defined for {halotype}.")

    if logger:
        logger.info(f"Normalising catalogue fields for {halotype}")

    standardised = {}

    for stdkey in STANDARD_KEYS:
        nativekey = mapping.get(stdkey)
        if nativekey is None:
            continue
        if nativekey in halos:
            standardised[stdkey] = halos[nativekey]
        else:
            if logger:
                logger.warning(f"Missing field '{nativekey}' for {halotype}")
            standardised[stdkey] = None

    # Optional: convert units or reformat arrays
    if halotype == "AHF":
        # Example: assemble coordinates if stored separately
        if "Xc" in halos and "Yc" in halos and "Zc" in halos:
            standardised["pos"] = np.vstack((halos["Xc"], halos["Yc"], halos["Zc"])).T
        if "Vx" in halos and "Vy" in halos and "Vz" in halos:
            standardised["vel"] = np.vstack((halos["Vx"], halos["Vy"], halos["Vz"])).T

    # Example unit scaling (to comoving Mpc/h, km/s)
    # You can extend this to detect internal units from metadata
    # and rescale consistently here.

    return standardised
