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

# Whether each format's *raw*, on-disk length/mass values already include
# the little-h factor (Mpc/h, 1e10 Msun/h) or have it factored out (Mpc,
# Msun) -- a property of the simulation code the catalogue was built from,
# not of AnalysisTools. SWIFT stores h-free units by convention; SubFind (run
# directly on GADGET/Arepo output) stores h-scaled units. AHF and
# VELOCIraptor are conventionally run against GADGET/Arepo-family output too
# (hence True here), but either can in principle be run against a SWIFT
# snapshot instead -- if that's your case, this default is wrong for your
# catalogue and little_h will silently apply the wrong conversion, so verify
# against your catalogue's own documented units before trusting it.
NATIVE_INCLUDES_LITTLE_H = {
    "SUBFIND": True,
    "AHF": True,
    "VELOCIraptor": True,
    "SWIFT_FOF": False,
}

# ---------------------------------------------------------------------
# Main user-facing class
# ---------------------------------------------------------------------

class HaloTools:
    """
    Unified high-level interface to load and inspect halo catalogues.

    Two independent unit axes, both applied in standardise_names() (raw
    read_catalogue() output without standardise=True is untouched, exactly
    as stored in the file):

    - comoving vs physical (scale factor 'a'): catalogues are always stored
      comoving on disk, so comoving=True (default) is a no-op; comoving=False
      multiplies pos/boxsize by 'a' to convert to physical coordinates.
    - little-h (whether Mpc/h or Mpc, 1e10 Msun/h or Msun): *not* the same
      axis as the above, and the two must not be conflated -- a value can be
      comoving-and-h-scaled, physical-and-h-free, or any other combination.
      Whether little_h=True/False requires actually dividing/multiplying by h
      depends on the *catalogue format's own native convention*: SWIFT
      already stores h-free values, GADGET/Arepo-family codes (SubFind, and
      conventionally AHF/VELOCIraptor) store h-scaled values -- see
      NATIVE_INCLUDES_LITTLE_H. Comparing a SWIFT-origin catalogue against a
      SubFind/Arepo one without accounting for this is a common source of
      factor-of-h mismatches.

      SubFind's convention is fixed (it's GADGET/Arepo's own group finder).
      AHF and VELOCIraptor's are *not* fixed -- both are standalone halo
      finders that inherit whatever convention their input snapshot used,
      so NATIVE_INCLUDES_LITTLE_H's True default for them is only a common-
      case guess, not a guarantee. Confirmed by example: a VELOCIraptor
      catalogue in this repo's own test data was run against a SWIFT
      snapshot and its Period is SWIFT's raw h-free BoxSize passed straight
      through, not VELOCIraptor's own h-scaled convention -- the True
      default would silently double-strip h for that catalogue. Pass
      native_includes_h= explicitly for AHF/VELOCIraptor rather than
      trusting the default, unless you've verified it against your
      catalogue's own units documentation.

    Parameters
    ----------
    comoving : bool, optional
        If True (default), ensure pos/boxsize are comoving (a no-op, since
        that's how these catalogues are stored). If False, convert to
        physical coordinates (multiply by the scale factor).
    little_h : bool, optional
        If True, ensure pos/mass/boxsize are in little-h units (Mpc/h,
        1e10 Msun/h). If False (default), ensure h is divided out (Mpc,
        Msun) -- physical units most analysis code expects. Whether this
        actually changes anything depends on the catalogue format's native
        convention (NATIVE_INCLUDES_LITTLE_H); a no-op is correctly applied
        for SWIFT_FOF even when little_h=False, since SWIFT is already
        h-free. No-ops (with a warning) if the catalogue's HubbleParam isn't
        available (currently: AHF; VELOCIraptor only if a '.siminfo' sidecar file is present).
    native_includes_h : bool, optional
        Override NATIVE_INCLUDES_LITTLE_H's per-format guess for whether
        this catalogue's *raw* pos/mass/BoxSize already include h. Default
        None (use the table). Strongly recommended for AHF/VELOCIraptor,
        whose native convention actually depends on the snapshot they were
        run against, not a fixed property of the format -- see above.
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
    >>> ht = HaloTools(comoving=True, little_h=False, centre_on_subhalo=True)
    >>> halos = ht.read_catalogue(filename="groups_010.hdf5",fileformat="SubFind",
    ...                            standardise=True)
    >>> ht.summary()
    """

    def __init__(
        self,
        comoving: bool = True,
        little_h: bool = False,
        native_includes_h: Optional[bool] = None,
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

        self.comoving = comoving
        self.little_h = little_h
        self.native_includes_h = native_includes_h
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
        self.metadata, self.halos, self.subhalos = reader(filename)

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

        self._apply_unit_conventions()

    def _get_scale_factor(self) -> Optional[float]:
        """Best-effort scale factor 'a' from self.metadata (ScaleFactor
        directly, or derived from Redshift), or None if unavailable."""
        meta = self.metadata or {}
        if meta.get("ScaleFactor") is not None:
            return float(meta["ScaleFactor"])
        z = meta.get("Redshift")
        if z is not None:
            return 1.0 / (1.0 + float(z))
        return None

    def _get_hubble_param(self) -> Optional[float]:
        """Best-effort HubbleParam from self.metadata, or None if this
        format's reader doesn't provide one (currently: AHF; VELOCIraptor only if a '.siminfo' sidecar is present)."""
        meta = self.metadata or {}
        for key in ("HubbleParam", "h", "hubble"):
            if meta.get(key) is not None:
                return float(meta[key])
        return None

    def _apply_unit_conventions(self) -> None:
        """Apply the comoving (scale-factor) and little_h conversions to
        self.standardised_halos/self.standardised_subhalos['pos'/'mass'] and
        self.metadata['BoxSize'], so all three end up in the same units --
        see the comoving/little_h parameters in the class docstring.

        Both axes reduce to a single combined scalar factor per quantity
        (length_factor for pos/BoxSize, mass_factor for mass), applied once,
        out-of-place, to whichever standardised tables/fields exist.
        """
        length_factor = 1.0
        mass_factor = 1.0

        a = self._get_scale_factor()
        if not self.comoving:
            if a is not None:
                length_factor *= a
            else:
                self.logger.warning(
                    "comoving=False requested but no scale factor could be "
                    "determined for this %s catalogue; pos/BoxSize left "
                    "comoving.", self.halocatfileformat)

        native_includes_h = (
            self.native_includes_h if self.native_includes_h is not None
            else NATIVE_INCLUDES_LITTLE_H.get(self.halocatfileformat, True))
        if self.little_h != native_includes_h:
            h = self._get_hubble_param()
            if h:
                # native h-included -> h-free: divide. native h-free ->
                # h-included: multiply. (Mpc/h = Mpc * h, per the standard
                # little-h convention -- see docs/unified_interface.md.)
                h_factor = (1.0 / h) if native_includes_h else h
                length_factor *= h_factor
                mass_factor *= h_factor
            else:
                self.logger.warning(
                    "little_h=%s requested but no HubbleParam is available "
                    "for this %s catalogue; pos/mass/BoxSize left as-is.",
                    self.little_h, self.halocatfileformat)

        if length_factor == 1.0 and mass_factor == 1.0:
            return

        for table in (self.standardised_halos, self.standardised_subhalos):
            if table is None:
                continue
            if table.get("pos") is not None and length_factor != 1.0:
                table["pos"] = table["pos"] * length_factor
            if table.get("mass") is not None and mass_factor != 1.0:
                table["mass"] = table["mass"] * mass_factor

        if self.metadata.get("BoxSize") is not None and length_factor != 1.0:
            # Convert from the raw, on-disk BoxSize (stashed on first call)
            # rather than the already-converted value, so re-running
            # standardise_names() on the same HaloTools instance doesn't
            # compound the factor.
            raw_boxsize = self.metadata.setdefault(
                "_raw_boxsize", self.metadata["BoxSize"])
            self.metadata["BoxSize"] = raw_boxsize * length_factor

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
