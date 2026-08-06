#!/usr/bin/env python3
"""
analysistools.api.halos
-----------------------
HaloCatalogue: the Dataset adapter over HaloTools.

The main object is the Group (host halo) table; ``cat.subhalos`` returns a
second HaloCatalogue over the Subhalo table where the format provides one.
Fields use the standardised schema (mass, pos, vel, radius, halo_id,
num_part) with all native fields passed through.
"""
from __future__ import annotations

from typing import Any, Dict, Optional

import numpy as np

from ..halo_tools import HaloTools
from .dataset import Dataset


class HaloCatalogue(Dataset):
    """
    Lazy Dataset view of a halo catalogue.

    Parameters
    ----------
    path : str
        Catalogue path.
    fileformat : str or int
        "SUBFIND", "AHF", "VELOCIraptor", "SWIFT_FOF" (or codes 1-4).
    comoving : bool, optional
        If True (default), pos/boxsize are comoving (a no-op, since that's
        the native storage convention). If False, converted to physical
        coordinates. This is the scale-factor axis, independent of little_h
        below -- see HaloTools' docstring for why the two must not be
        conflated (a common source of factor-of-h mismatches between
        catalogue formats).
    little_h : bool, optional
        If True, pos/mass/boxsize are in little-h units (Mpc/h, 1e10
        Msol/h). If False (default), h is divided out -- but whether that's
        actually a no-op or a real conversion depends on the catalogue
        format's native convention: SWIFT_FOF is already h-free; SUBFIND
        (and conventionally AHF/VELOCIraptor) are not. See
        HaloTools.NATIVE_INCLUDES_LITTLE_H.
    snapnum : int, optional
        Recorded for merger-tree lookups (MergerTreeTools.from_halo).
    label : str, optional
        Name for plot legends.
    **backend_kwargs :
        Any further HaloTools constructor options, e.g. centre_on_subhalo=True
        (SUBFIND only: use each group's primary subhalo -- GroupFirstSub ->
        SubhaloPos/SubhaloVel -- as its 'pos'/'vel' instead of GroupPos/
        GroupVel), or native_includes_h=False/True to override
        NATIVE_INCLUDES_LITTLE_H's per-format guess -- strongly recommended
        for AHF/VELOCIraptor catalogues, whose native little-h convention
        depends on the snapshot they were run against, not a fixed property
        of the format (see HaloTools' docstring).

    Examples
    --------
    >>> cat = HaloCatalogue("data/halos/snap_0031.VELOCIraptor.properties.0",
    ...                     fileformat="VELOCIraptor")
    >>> cat["mass"], cat["pos"]
    >>> massive = cat.select(mass=(1e3, None))
    >>> cat.subhalos["mass"]
    >>> sf = HaloCatalogue("groups_010.hdf5", fileformat="SUBFIND",
    ...                    centre_on_subhalo=True)
    >>> sf["centred_on_subhalo"]   # False where a group had no bound subhalo
    """

    kind = "halos"

    def __init__(self, path: str, fileformat, comoving: bool = True,
                 little_h: bool = False,
                 snapnum: Optional[int] = None,
                 label: Optional[str] = None, **backend_kwargs):
        super().__init__(path=path, fileformat=str(fileformat), label=label)
        backend_kwargs.setdefault("loglevel", 30)
        self._backend = HaloTools(comoving=comoving, little_h=little_h,
                                  **backend_kwargs)
        self._snapnum = snapnum
        self._sub_columns: Optional[Dict[str, np.ndarray]] = None

    # ------------------------------------------------------------------

    def _load(self) -> None:
        meta, halos, subhalos = self._backend.read_catalogue(
            self.path, fileformat=self.fileformat, standardise=True,
            snapnum=self._snapnum)

        self._columns = {k: v for k, v in halos.items()
                         if isinstance(v, np.ndarray)}
        self._sub_columns = ({k: v for k, v in subhalos.items()
                              if isinstance(v, np.ndarray)}
                             if subhalos else None)

        a = float(meta.get("ScaleFactor", 1.0) or 1.0)
        comoving = self._backend.comoving
        little_h = self._backend.little_h
        h_suffix = ", h-included" if little_h else ", h-free"
        native_meta = {k: v for k, v in meta.items() if not k.startswith("_")}

        self.meta.update({
            "scale_factor": a,
            "redshift": 1.0 / a - 1.0 if a > 0 else None,
            "boxsize": (float(meta["BoxSize"])
                        if meta.get("BoxSize") is not None else None),
            "h0": meta.get("HubbleParam"),
            "comoving": comoving,
            "little_h": little_h,
            "n_groups": int(meta.get("TotNgroups", 0) or 0),
            "snapnum": self._snapnum,
            "native_meta": native_meta,
            "units": {"length": "catalogue-native" + h_suffix,
                      "mass": "catalogue-native" + h_suffix,
                      "velocity": "km/s"},
        })
        # normalise the canonical format string chosen by HaloTools
        self.fileformat = self._backend.halocatfileformat

    # ------------------------------------------------------------------

    @property
    def subhalos(self) -> "HaloCatalogue":
        """Subhalo table as its own HaloCatalogue (raises if absent)."""
        self._ensure_loaded()
        if not self._sub_columns:
            raise AttributeError(
                f"No subhalo table in this {self.fileformat} catalogue.")
        sub = object.__new__(HaloCatalogue)
        sub.__dict__.update(self.__dict__)
        sub._columns = self._sub_columns
        sub._sub_columns = None
        sub._index = None
        sub.label = f"{self.label}:subhalos"
        return sub

    @property
    def has_subhalos(self) -> bool:
        """
        Whether this catalogue has a subhalo table.

        Returns
        -------
        bool
        """
        self._ensure_loaded()
        return bool(self._sub_columns)

    @property
    def backend(self) -> HaloTools:
        """The underlying HaloTools instance (full legacy API)."""
        return self._backend

    def summary(self) -> None:
        """
        Print the base dataset summary plus the subhalo count, if present.
        """
        super().summary()
        if self.has_subhalos:
            n = len(next(iter(self._sub_columns.values())))
            print(f"  subhalos: {n:,} (access via .subhalos)")
