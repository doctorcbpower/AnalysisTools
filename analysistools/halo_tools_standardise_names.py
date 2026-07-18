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
from typing import Dict, Optional, Tuple

STANDARD_KEYS = [
    "mass",     # total halo mass (M200, Mvir, etc.)
    "pos",      # 3D position [x,y,z]
    "vel",      # 3D velocity [vx,vy,vz]
    "radius",   # virial or FOF radius
    "halo_id",       # halo ID
    "num_part",       # number of particles
]

# Sentinel used in CATALOGUE_MAPPINGS in place of a native field name, for
# cases where a catalogue has no explicit stored ID and the ID is instead
# just the halo's position in the array (SubFind's subhalo numbering
# convention: subhalo N is row N of the Subhalo/* datasets).
ROW_INDEX_SENTINEL = "__ROW_INDEX__"

# Mappings are keyed by (halotype, object_type). "Group"-level and
# "Subhalo"-level tables commonly use different native field names (e.g.
# SubFind's GroupPos vs SubhaloPos), so they need separate entries -- using
# one mapping for both, as this file previously did, silently produces
# all-None fields for whichever table the mapping doesn't actually match.
CATALOGUE_MAPPINGS: Dict[Tuple[str, str], Dict[str, str]] = {
    ("SUBFIND", "Group"): {
        "mass": "Group_M_Crit200",
        "pos": "GroupPos",
        "vel": "GroupVel",
        "radius": "Group_R_Crit200",
        "halo_id": ROW_INDEX_SENTINEL,   # no stored SubhaloNr; ID is the row index
        "num_part": "GroupLen",
    },
    ("SUBFIND", "Subhalo"): {
        "mass": "SubhaloMass",
        "pos": "SubhaloPos",
        "vel": "SubhaloVel",
        "radius": "SubhaloHalfmassRad",
        "halo_id": ROW_INDEX_SENTINEL,   # no stored SubhaloNr; ID is the row index
        "num_part": "SubhaloLen",
    },
    ("AHF", "Group"): {
        "mass": "Mvir",
        "pos": "Xc",   # can be split from Xc, Yc, Zc if stored separately
        "vel": "Vcm",  # likewise Vx, Vy, Vz may need stacking
        "radius": "Rvir",
        "halo_id": "ID",
        "num_part": "npart",
    },
    ("VELOCIraptor", "Group"): {
        "mass": "Mass_200crit",
        "pos": "CentreOfPotential",
        "vel": "Velocity",
        "radius": "R_200crit",
        "halo_id": "ID",
        "num_part": "npart",
    },
    ("VELOCIraptor", "Subhalo"): {
        "mass": "Mass_200crit",
        "pos": "CentreOfPotential",
        "vel": "Velocity",
        "radius": "R_200crit",
        "halo_id": "ID",
        "num_part": "npart",
    },
    ("SWIFT_FOF", "Group"): {
        "mass": "Masses",
        "pos": "Centres",
        "vel": "",
        "radius": "Radii",
        "halo_id": "GroupIDs",
        "num_part": "Sizes",
    },
}


def standardise_catalogue_names(halos: Dict[str, np.ndarray],
                                halotype: str,
                                object_type: str = "Group",
                                logger: Optional[logging.Logger] = None,
                                keep_extra_fields: bool = True) -> Dict[str, np.ndarray]:
    """
    Convert halo dictionaries to common naming conventions.

    Parameters
    ----------
    halos : dict
        Raw halo data as returned by a reader (e.g. read_subfind), for
        *either* the Group table or the Subhalo table -- pass the matching
        `object_type` so the right native field names are used.
    halotype : str
        One of "SUBFIND", "AHF", "VELOCIraptor", "SWIFT_FOF".
    object_type : {"Group", "Subhalo"}
        Which table `halos` holds. Group- and Subhalo-level tables often use
        different native field names (SubFind's GroupPos vs SubhaloPos), so
        this matters -- passing the wrong one will silently produce None (or
        wrong) fields rather than raising, since a missing field looks the
        same either way.
    logger : logging.Logger, optional
        Logger for messages.
    keep_extra_fields : bool, optional
        If True (default), every native field that isn't consumed by
        CATALOGUE_MAPPINGS is passed through into the returned dict under
        its original native name, so format-specific fields (e.g.
        SubFind's GroupFirstSub, AHF's n2, VELOCIraptor's Vmax, ...) aren't
        silently dropped. Set False to get back only STANDARD_KEYS, which
        was the previous behaviour.

    Returns
    -------
    dict
        Dictionary containing the STANDARD_KEYS entries, plus (by default)
        any remaining native fields under their original names.
    """
    if halos is None:
        raise ValueError("No halo data provided for normalisation.")

    halotype = halotype.upper()
    mapping = CATALOGUE_MAPPINGS.get((halotype, object_type))
    if mapping is None:
        raise ValueError(
            f"No normalisation mapping defined for ({halotype}, {object_type}). "
            f"Available: {list(CATALOGUE_MAPPINGS.keys())}")

    if logger:
        logger.info(f"Normalising catalogue fields for {halotype} ({object_type})")

    standardised = {}
    n = len(next(iter(halos.values()))) if halos else 0

    # Native keys that CATALOGUE_MAPPINGS consumes for this (halotype,
    # object_type) -- these are excluded from the pass-through of "extra"
    # fields below, since they're already represented under their standard
    # name.
    mapped_native_keys = {v for v in mapping.values() if v != ROW_INDEX_SENTINEL}

    for stdkey in STANDARD_KEYS:
        nativekey = mapping.get(stdkey)
        if nativekey is None:
            continue
        if nativekey == ROW_INDEX_SENTINEL:
            # This table has no explicit stored ID field; the ID convention
            # is the row index (e.g. SubFind subhalo numbering). Only valid
            # for a catalogue read from a single file with no ID offset.
            standardised[stdkey] = np.arange(n, dtype=np.int64)
        elif nativekey in halos:
            standardised[stdkey] = halos[nativekey]
        else:
            if logger:
                logger.warning(f"Missing field '{nativekey}' for {halotype} ({object_type})")
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

    # Pass through every remaining native field so format-specific data
    # (e.g. GroupFirstSub, SubhaloCM, Vmax, ...) isn't lost just because it
    # has no entry in CATALOGUE_MAPPINGS.
    if keep_extra_fields:
        for nativekey, value in halos.items():
            if nativekey in mapped_native_keys:
                continue
            if nativekey in standardised:
                # A native field name happens to collide with one of
                # STANDARD_KEYS (e.g. a catalogue that already has a field
                # literally called "mass"). Don't clobber the standardised
                # entry silently.
                if logger:
                    logger.warning(
                        f"Native field '{nativekey}' collides with a standard "
                        f"key name for {halotype} ({object_type}); keeping the "
                        f"standardised value and skipping the raw pass-through.")
                continue
            standardised[nativekey] = value

    return standardised
