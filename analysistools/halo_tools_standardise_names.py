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
        # Native names as actually produced by haloio_ahf.read_ahf() (which
        # combines Xc/Yc/Zc and VXc/VYc/VZc into single Pos/Vel arrays
        # itself) -- this previously named fields ("Mvir", "Xc", "Vcm", "ID",
        # "npart") that read_ahf() never produces, so mass/pos/vel/halo_id/
        # num_part all silently came back None for every AHF catalogue.
        "mass": "Mass",
        "pos": "Pos",
        "vel": "Vel",
        "radius": "Rvir",
        "halo_id": "HaloID",
        "num_part": "NumPart",
    },
    ("AHF", "Subhalo"): {
        # Native names from haloio_ahf.read_ahf()'s .AHF_substructure
        # parsing. No distinct subhalo ID or particle-count column is read,
        # hence ROW_INDEX_SENTINEL and no "num_part" entry (omitted, rather
        # than pointing at a field that doesn't exist).
        "mass": "SubhaloMass",
        "pos": "Pos",
        "vel": "Vel",
        "radius": "Rvir",
        "halo_id": ROW_INDEX_SENTINEL,
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

    # Case-insensitive halotype lookup ("VELOCIraptor" == "VELOCIRAPTOR")
    halotype = halotype.upper()
    _by_upper = {(ht.upper(), ot): m for (ht, ot), m in CATALOGUE_MAPPINGS.items()}
    mapping = _by_upper.get((halotype, object_type))
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
            # Don't warn for pos/vel when per-component columns are present:
            # the stacking block below will assemble them.
            stackable = (
                (stdkey == "pos" and all(k in halos for k in ("Xc", "Yc", "Zc")))
                or (stdkey == "vel" and (
                    all(k in halos for k in ("Vx", "Vy", "Vz"))
                    or all(k in halos for k in ("VXc", "VYc", "VZc")))))
            if logger and not stackable:
                logger.warning(f"Missing field '{nativekey}' for {halotype} ({object_type})")
            standardised[stdkey] = None

    # Assemble (N, 3) pos/vel from per-component columns where the mapped
    # field is absent or 1-D. Covers AHF (Xc/Yc/Zc, Vx/Vy/Vz) and flat
    # VELOCIraptor .properties files (Xc/Yc/Zc, VXc/VYc/VZc).
    def _needs_stacking(key):
        arr = standardised.get(key)
        return arr is None or (hasattr(arr, "ndim") and arr.ndim == 1)

    if _needs_stacking("pos") and all(k in halos for k in ("Xc", "Yc", "Zc")):
        standardised["pos"] = np.vstack((halos["Xc"], halos["Yc"], halos["Zc"])).T
    for vx, vy, vz in (("Vx", "Vy", "Vz"), ("VXc", "VYc", "VZc")):
        if _needs_stacking("vel") and all(k in halos for k in (vx, vy, vz)):
            standardised["vel"] = np.vstack((halos[vx], halos[vy], halos[vz])).T
            break

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
