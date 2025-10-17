#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
halo_tools.py
Refactored: 2025-10-16

Unified interface for reading and analysing halo catalogues
from AHF, SubFind, and VELOCIraptor outputs.

Author: C. Power
Refactor: ChatGPT (GPT-5)
"""

import os
import logging
from typing import Optional, Dict, Any, Tuple, Union

import h5py
import numpy as np

from .haloio_subfind import read_subfind
from .haloio_ahf import read_ahf
from .haloio_velociraptor import read_velociraptor

FORMAT_READERS = {
    "SUBFIND": read_subfind,
    "AHF": read_ahf,
    "VELOCIrapTOR": read_velociraptor,
}

# ---------------------------------------------------------------------
# Main user-facing class
# ---------------------------------------------------------------------

class HaloTools:
    """
    Unified high-level interface to load and inspect halo catalogues.

    Examples
    --------
    >>> ht = HaloTools("SubFind", comoving_units=True)
    >>> halos = ht.read("groups_010.hdf5")
    >>> ht.summary()
    """

    def __init__(
        self,
        halocatfileformat: Union[str, int],
        comoving_units: bool = False,
        **kwargs,
    ):
        fmtmap = {
            "1": "SUBFIND",
            "2": "AHF",
            "3": "VELOCIrapTOR",
            1: "SUBFIND",
            2: "AHF",
            3: "VELOCIrapTOR",
        }
        self.halocatfileformat = fmtmap.get(str(halocatfileformat).upper(), str(halocatfileformat).upper())
        if self.halocatfileformat not in FORMAT_READERS:
            raise ValueError(f"Unsupported halo catalogue format: {halocatfileformat}")

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

    # -----------------------------------------------------------------
    # Core I/O methods
    # -----------------------------------------------------------------

    def read(self, filename: str) -> Dict[str, np.ndarray]:
        """Read a halo catalogue from disk."""
        self.filename = filename
        reader = FORMAT_READERS[self.halocatfileformat]

        self.logger.info(f"Reading halo catalogue '{filename}' using {self.halocatfileformat}")
        self.metadata, self.halos, self.subhalos = reader(filename, comoving=self.comoving_units)

        nh = len(next(iter(self.halos.values()))) if self.halos else 0
        self.logger.info(f"Loaded {nh:,} halos from {self.halocatfileformat} file.")
        return self.metadata, self.halos, self.subhalos

    # -----------------------------------------------------------------
    # Utility methods
    # -----------------------------------------------------------------

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

    def get_property(self, key: str, subhalos: bool = False) -> Optional[np.ndarray]:
        """Retrieve a property array for halos or subhalos."""
        data = self.subhalos if subhalos else self.halos
        if data is None:
            self.logger.warning("No data loaded.")
            return None
        return data.get(key)

    # -----------------------------------------------------------------
    # Optional: convenience computation
    # -----------------------------------------------------------------

    def compute_virial_masses(self) -> Optional[np.ndarray]:
        """Example derived quantity computation."""
        if self.halos is None or "Group_R_Crit200" not in self.halos:
            self.logger.warning("No virial radii found.")
            return None

        r = self.halos["Group_R_Crit200"]
        rho_crit = 2.775e11  # h^2 Msun/Mpc^3
        mass = 4 / 3 * np.pi * (r**3) * 200 * rho_crit
        self.halos["M200"] = mass
        return mass
