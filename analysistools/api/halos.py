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
        Return positions in comoving coordinates (HaloTools comoving_units).
    snapnum : int, optional
        Recorded for merger-tree lookups (MergerTreeTools.from_halo).
    label : str, optional
        Name for plot legends.

    Examples
    --------
    >>> cat = HaloCatalogue("data/halos/snap_0031.VELOCIraptor.properties.0",
    ...                     fileformat="VELOCIraptor")
    >>> cat["mass"], cat["pos"]
    >>> massive = cat.select(mass=(1e3, None))
    >>> cat.subhalos["mass"]
    """

    kind = "halos"

    def __init__(self, path: str, fileformat, comoving: bool = True,
                 snapnum: Optional[int] = None,
                 label: Optional[str] = None, **backend_kwargs):
        super().__init__(path=path, fileformat=str(fileformat), label=label)
        backend_kwargs.setdefault("loglevel", 30)
        self._backend = HaloTools(comoving_units=comoving, **backend_kwargs)
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
        self.meta.update({
            "scale_factor": a,
            "redshift": 1.0 / a - 1.0 if a > 0 else None,
            "boxsize": (float(meta["BoxSize"])
                        if meta.get("BoxSize") is not None else None),
            "n_groups": int(meta.get("TotNgroups", 0) or 0),
            "snapnum": self._snapnum,
            "native_meta": dict(meta),
            "units": {"length": "catalogue-native", "mass": "catalogue-native",
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
        self._ensure_loaded()
        return bool(self._sub_columns)

    @property
    def backend(self) -> HaloTools:
        """The underlying HaloTools instance (full legacy API)."""
        return self._backend

    def summary(self) -> None:
        super().summary()
        if self.has_subhalos:
            n = len(next(iter(self._sub_columns.values())))
            print(f"  subhalos: {n:,} (access via .subhalos)")
