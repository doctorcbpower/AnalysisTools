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

    Examples
    --------
    >>> ht = HaloTools(comoving_units=True)
    >>> halos = ht.read_catalogue(halocat_filename="groups_010.hdf5",halocat_fileformat="SubFind")
    >>> ht.summary()
    """

    def __init__(
        self,
        comoving_units: bool = False,
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
        # self.halocatfileformat = fmtmap.get(str(halocatfileformat).upper(), str(halocatfileformat).upper())
        # if self.halocatfileformat not in FORMAT_READERS:
        #     raise ValueError(f"Unsupported halo catalogue format: {halocatfileformat}")

        self.comoving_units = comoving_units
        self.filename: Optional[str] = None
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

    def read_catalogue(self, filename: str, fileformat: Union[str, int],) -> Dict[str, np.ndarray]:
        """Read a halo catalogue from disk."""
        self.filename = filename
        self.halocatfileformat = self.fmtmap.get(str(fileformat).upper(), str(fileformat).upper())
        if self.halocatfileformat not in FORMAT_READERS:
            raise ValueError(f"Unsupported halo catalogue format: {fileformat}")

        reader = FORMAT_READERS[self.halocatfileformat]

        self.logger.info(f"Reading halo catalogue '{filename}' using {self.halocatfileformat}")
        self.metadata, self.halos, self.subhalos = reader(filename, comoving=self.comoving_units)

        nh = len(next(iter(self.halos.values()))) if self.halos else 0
        self.logger.info(f"Loaded {nh:,} halos from {self.halocatfileformat} file.")
        return self.metadata, self.halos, self.subhalos

    def standardise_names(self):
        """Normalise raw halo data to the common schema."""
        if self.halos is None:
            self.logger.warning("No halo data loaded.")
            return

        self.standardised_names = standardise_catalogue_names(self.halos, self.halocatfileformat, self.logger)

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

