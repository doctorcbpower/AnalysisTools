#!/usr/bin/env python3
"""
halo_tools.py
Refactored: 2025-10-16

Unified interface for reading and analysing halo catalogues
from AHF, SubFind, VELOCIraptor, and SWIFT outputs.

Author: C. Power
"""

import os
import logging
from typing import Optional, Dict, Any, Tuple, Union

import h5py
import numpy as np

from .haloio_subfind import read_subfind
from .haloio_swiftfof import read_swiftfof
from .haloio_ahf import read_ahf
from .haloio_velociraptor import read_velociraptor
from .halo_tools_standardise_names import standardise_catalogue_names

FORMAT_READERS = {
    "SUBFIND": read_subfind,
    "AHF": read_ahf,
    "VELOCIraptor": read_velociraptor,
    "SWIFT_FOF": read_swiftfof,
#    "SWIFT_HBT": read_swifthbt,
}

# ---------------------------------------------------------------------
# Main user-facing class
# ---------------------------------------------------------------------

class HaloTools:
    """
    Unified high-level interface to load and inspect halo catalogues.

    Parameters
    ----------
    comoving_units : bool, optional
        See individual reader docstrings (divides positions by HubbleParam).
    centre_on_subhalo : bool, optional
        SUBFIND only. If True, standardise_names() replaces each group's
        'pos' (and 'vel', if present) with its primary subhalo's
        SubhaloPos/SubhaloVel (looked up via GroupFirstSub) instead of
        GroupPos/GroupVel -- the FOF group centre can be offset from the
        actual halo centre by substructure, so the primary subhalo's
        position is usually the better proxy. Groups with no identified
        subhalo (GroupFirstSub invalid) keep GroupPos/GroupVel. Adds a
        'centred_on_subhalo' bool column recording which groups were
        actually re-centred. SubhaloVel is SubFind's own subhalo bulk
        velocity, not a self-consistent recomputation from bound
        particles -- treat the velocity substitution as an approximation.

    Examples
    --------
    >>> ht = HaloTools(comoving_units=True, centre_on_subhalo=True)
    >>> halos = ht.read_catalogue(filename="groups_010.hdf5",fileformat="SubFind",
    ...                            standardise=True)
    >>> ht.summary()
    """

    def __init__(
        self,
        comoving_units: bool = False,
        centre_on_subhalo: bool = False,
        **kwargs,
    ):
        self.fmtmap = {
            "1": "SUBFIND",
            "2": "AHF",
            "3": "VELOCIraptor",
            "4": "SWIFT_FOF",
            "5": "SWIFT_HBT",
            1: "SUBFIND",
            2: "AHF",
            3: "VELOCIraptor",
            4: "SWIFT_FOF",
            5: "SWIFT_HBT",
        }

        self.comoving_units = comoving_units
        self.centre_on_subhalo = centre_on_subhalo
        self.metadata: Dict[str, Any] = {}
        self.halos: Optional[Dict[str, np.ndarray]] = None
        self.subhalos: Optional[Dict[str, np.ndarray]] = None

        self.usehalocatonly = kwargs.get("usehalocatonly", False)
        self.usesubstructure_file = kwargs.get("usesubstructure_file", False)

        # Logging setup
        self.logger = logging.getLogger(__name__)
        if not self.logger.handlers:
            handler = logging.StreamHandler()
            handler.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))
            self.logger.addHandler(handler)
            self.logger.propagate = False
        self.logger.setLevel(kwargs.get("loglevel", logging.INFO))

    def read_catalogue(self, filename: str, fileformat: Union[str, int], standardise: bool = False,
                        snapnum: Optional[int] = None) -> Dict[str, np.ndarray]:
        """
        Read a halo catalogue from disk.
        
        Parameters
        ----------
        filename : str
            Path to the halo catalogue file.
        fileformat : str or int
            Format of the halo catalogue (e.g. "SUBFIND", "AHF", "VELOCIraptor", "SWIFT_FOF").
        standardise : bool, optional
            If True, normalise field names to common schema.
        snapnum : int, optional
            Snapshot number this catalogue corresponds to. Not used internally
            by HaloTools itself, but recorded on `self.snapnum` so that other
            tools (e.g. MergerTreeTools.from_halo) can look a halo up in a
            merger tree without you having to pass the snapshot number
            separately every time.
        
        Returns
        -------
        halos : dict[str, np.ndarray]
            Dictionary of halo properties.
        subhalos : dict[str, np.ndarray]
            Dictionary of subhalo properties (if present).                      
        metadata : dict
            Dictionary of metadata from the catalogue header.            
        """
        
        self.filename = filename
        self.snapnum = snapnum
        # Normalise: integer codes via fmtmap, then case-insensitive match
        # against FORMAT_READERS keys (so "VELOCIraptor"/"velociraptor"/3 all work).
        fmt = self.fmtmap.get(fileformat, self.fmtmap.get(str(fileformat), str(fileformat)))
        canonical = {k.upper(): k for k in FORMAT_READERS}
        self.halocatfileformat = canonical.get(str(fmt).upper(), str(fmt).upper())
        if self.halocatfileformat not in FORMAT_READERS:
            raise ValueError(f"Unsupported halo catalogue format: {fileformat}")

        reader = FORMAT_READERS[self.halocatfileformat]

        self.logger.info(f"Reading halo catalogue '{filename}' using {self.halocatfileformat}")
        self.metadata, self.halos, self.subhalos = reader(filename, comoving=self.comoving_units)

        nh = len(next(iter(self.halos.values()))) if self.halos else 0
        self.logger.info(f"Loaded {nh:,} halos from {self.halocatfileformat} file.")
        
        # Optional post-processing
        if standardise:
            self.standardise_names()
            return self.metadata, self.standardised_halos, self.standardised_subhalos

        return self.metadata, self.halos, self.subhalos

    def standardise_names(self):
        """Normalise raw halo and subhalo data to the common schema.

        Stores results on `self.standardised_halos` and
        `self.standardised_subhalos`. Group- and Subhalo-level tables are
        standardised separately (with object_type="Group"/"Subhalo") since
        they generally use different native field names.
        """
        if self.halos is None:
            self.logger.warning("No halo data loaded.")
            self.standardised_halos = None
        else:
            self.standardised_halos = standardise_catalogue_names(
                self.halos, self.halocatfileformat, object_type="Group", logger=self.logger)

        if self.subhalos:
            self.standardised_subhalos = standardise_catalogue_names(
                self.subhalos, self.halocatfileformat, object_type="Subhalo", logger=self.logger)
        else:
            self.standardised_subhalos = None

        if (self.centre_on_subhalo and self.halocatfileformat == "SUBFIND"
                and self.standardised_halos is not None
                and self.standardised_subhalos is not None):
            self._centre_on_first_subhalo()

    def _centre_on_first_subhalo(self) -> None:
        """Re-centre self.standardised_halos on each group's primary
        subhalo (GroupFirstSub -> SubhaloPos/SubhaloVel). See
        centre_on_subhalo in the class docstring."""
        halos = self.standardised_halos
        subs = self.standardised_subhalos

        first_sub = halos.get("GroupFirstSub")
        if first_sub is None:
            self.logger.warning(
                "centre_on_subhalo requested but 'GroupFirstSub' is not "
                "present in this SUBFIND catalogue; leaving GroupPos/GroupVel "
                "as-is.")
            return

        # Upcast before comparing so both sentinel conventions are caught:
        # signed -1, or (for an unsigned dtype, which can't hold -1) a huge
        # positive fill value -- both fail one of the two bounds below.
        first_sub = np.asarray(first_sub).astype(np.int64)
        n_sub = len(subs["pos"]) if subs.get("pos") is not None else 0
        valid = (first_sub >= 0) & (first_sub < n_sub)

        did_pos = did_vel = False
        if halos.get("pos") is not None and subs.get("pos") is not None:
            new_pos = halos["pos"].copy()
            new_pos[valid] = subs["pos"][first_sub[valid]]
            halos["pos"] = new_pos
            did_pos = True

        if halos.get("vel") is not None and subs.get("vel") is not None:
            new_vel = halos["vel"].copy()
            new_vel[valid] = subs["vel"][first_sub[valid]]
            halos["vel"] = new_vel
            did_vel = True

        halos["centred_on_subhalo"] = valid
        self.logger.info(
            f"Re-centred {int(valid.sum()):,}/{len(valid):,} groups on their "
            f"primary subhalo's "
            f"{'position and velocity' if did_vel else 'position' if did_pos else '(nothing -- no pos/vel to replace)'}"
            f"; the rest kept GroupPos/GroupVel (no valid GroupFirstSub).")

    def summary(self) -> None:
        """Print a concise summary of loaded halo catalogue."""
        if self.halos is None:
            print("No halo data loaded.")
            return

        nh = len(next(iter(self.halos.values()))) if self.halos else 0
        ns = len(next(iter(self.subhalos.values()))) if self.subhalos else 0
        fmt = self.halocatfileformat

        print(f"Halo Catalogue Summary")
        print(f"  Format: {fmt}")
        print(f"  File:   {os.path.basename(self.filename) if self.filename else 'N/A'}")
        print(f"  Halos:  {nh:,}")
        print(f"  Subhalos: {ns:,}")
        if self.metadata:
            for k, v in self.metadata.items():
                print(f"  {k:20s}: {v}")
