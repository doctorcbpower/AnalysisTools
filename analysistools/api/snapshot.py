#!/usr/bin/env python3
"""
analysistools.api.snapshot
--------------------------
SnapshotDataset: the Dataset adapter over SnapshotTools/SnapshotData.

Per-species access (design decision, DEVELOPMENT.md section 7.2):
``ds["pos"]`` returns all particles; ``ds.dm``, ``ds.gas``, ``ds.star``,
``ds.bh`` return per-species Dataset views, so ``ds.dm.pos`` and
``ds.dm["pos"]`` both work. ``select(species="dm")`` is the programmatic
equivalent.
"""
from __future__ import annotations

from typing import Any, Dict, Optional, Sequence

import numpy as np

from ..snapshot_tools import SnapshotTools
from .dataset import Dataset

#: species name -> SnapshotTools constructor attribute holding the type code
SPECIES_TYPE_ATTR = {"gas": "gas_type", "dm": "dm_type",
                     "star": "star_type", "bh": "bh_type"}


class SnapshotDataset(Dataset):
    """
    Lazy Dataset view of a particle snapshot (HDF5 or GADGET binary).

    Parameters
    ----------
    path : str
        Snapshot path (multi-file snapshots detected as in SnapshotTools).
    fileformat : str, optional
        "HDF5" (default), "SNAP1", "SNAP2" -- passed to SnapshotTools.
    convention : str, optional
        "SWIFT", "GADGET4", "AREPO", "GADGET2/3". The load() factory sniffs
        this from the file; default here is SnapshotTools' default.
    label : str, optional
        Name for plot legends.
    **backend_kwargs :
        Any further SnapshotTools constructor options (extra_blocks, ...).

    Examples
    --------
    >>> snap = SnapshotDataset("data/snap_0031.hdf5", convention="SWIFT")
    >>> snap["pos"].shape          # all particles
    >>> snap.dm["pos"].shape       # dark matter only (same as snap.dm.pos)
    >>> sub = snap.select(centre=[70, 70, 70], size=5.0)
    """

    kind = "snapshot"

    def __init__(self, path: str, fileformat: str = "HDF5",
                 convention: Optional[str] = None,
                 label: Optional[str] = None, **backend_kwargs):
        super().__init__(path=path, fileformat=fileformat, label=label)
        backend_kwargs.setdefault("loglevel", 30)   # WARNING: quiet library
        self._backend = SnapshotTools(snapfileformat=fileformat,
                                      **backend_kwargs)
        self._convention = convention
        self._data = None                  # SnapshotData, set by _load()

    # ------------------------------------------------------------------

    def _load(self) -> None:
        self._data = self._backend.read_snapshot(self.path,
                                                 convention=self._convention)
        d = self._data

        for name in ("pos", "vel", "pids", "mass", "ptype", "potential",
                     "rho", "u", "sfr", "age", "groupid", "hsml"):
            arr = getattr(d, name, None)
            if isinstance(arr, np.ndarray) and arr.ndim >= 1 and len(arr):
                self._columns[name] = arr

        def _scalar(x):
            arr = np.atleast_1d(np.asarray(x, dtype=float))
            return float(arr[0])

        a = _scalar(getattr(d, "scale_factor", 1.0))
        self.meta.update({
            "scale_factor": a,
            "redshift": 1.0 / a - 1.0 if a > 0 else None,
            "boxsize": _scalar(getattr(d, "box_size", 0.0)) or None,
            "h0": _scalar(getattr(d, "hubble_param", 0.0)) or None,
            "omega_0": _scalar(getattr(d, "omega_0", 0.0)) or None,
            "omega_lambda": _scalar(getattr(d, "omega_lambda", 0.0)) or None,
            "num_part_total": np.asarray(getattr(d, "num_part_total", [])),
            "convention": self._convention or self._backend.convention,
            # plain strings (no astropy) -- GADGET-family defaults
            "units": {"length": "code (kpc-family, see backend unit_* attrs)",
                      "mass": "code (1e10 Msol-family)",
                      "velocity": "km/s"},
        })

    # ------------------------------------------------------------------
    # Species views
    # ------------------------------------------------------------------

    def _species_view(self, species: str) -> "SnapshotDataset":
        if species not in SPECIES_TYPE_ATTR:
            raise ValueError(f"Unknown species '{species}' "
                             f"(use one of {sorted(SPECIES_TYPE_ATTR)}).")
        self._ensure_loaded()
        code = getattr(self._backend, SPECIES_TYPE_ATTR[species])
        ptype = self._resolve("ptype")
        if self._index is not None:
            ptype = ptype[self._index]
        rows = np.flatnonzero(ptype == code)
        view = self._make_view(rows)
        view.label = f"{self.label}:{species}"
        return view

    @property
    def gas(self) -> "SnapshotDataset":
        return self._species_view("gas")

    @property
    def dm(self) -> "SnapshotDataset":
        return self._species_view("dm")

    @property
    def star(self) -> "SnapshotDataset":
        return self._species_view("star")

    @property
    def bh(self) -> "SnapshotDataset":
        return self._species_view("bh")

    def select(self, mask=None, *, species: Optional[str] = None,
               **kwargs) -> "SnapshotDataset":
        base = self._species_view(species) if species else self
        return Dataset.select(base, mask, **kwargs)

    # ------------------------------------------------------------------

    @property
    def backend(self) -> SnapshotTools:
        """The underlying SnapshotTools instance (full legacy API)."""
        return self._backend

    @property
    def data(self):
        """The underlying SnapshotData instance (post-read), or None."""
        return self._data

    def summary(self) -> None:
        super().summary()
        npt = self.meta.get("num_part_total")
        if npt is not None and len(npt):
            counts = {s: int(np.sum(npt[getattr(self._backend,
                                                SPECIES_TYPE_ATTR[s])]))
                      for s in SPECIES_TYPE_ATTR}
            present = ", ".join(f"{k}={v:,}" for k, v in counts.items() if v)
            print(f"  species: {present or 'none flagged'}")
