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

#: header/config scalars SnapshotData may carry after read_snapshot(), beyond
#: the handful promoted to top-level meta keys -- mirrors the attribute list
#: SnapshotTools._transfer_attributes_to_writer uses for the same purpose.
#: Collected into meta["native_meta"], parity with HaloCatalogue's native_meta.
_HEADER_ATTRS = (
    "box_size", "dimension", "scale_factor", "redshift", "time",
    "omega_dm", "omega_bar", "omega_b", "omega_0", "omega_lambda",
    "hubble_param", "h",
    "num_files", "num_part_total", "num_part_this_file", "num_type",
    "num_part_type", "mass_table",
    "flag_cooling", "flag_stellar_age", "flag_sfr", "flag_metals",
    "flag_feedback", "flag_double_precision",
    "ispotential", "isgroupid",
    "name_of_mass_block", "name_of_u_block",
    "unit_length_in_cgs", "unit_mass_in_cgs", "unit_velocity_in_cgs",
    "unit_time_in_cgs", "unit_density_in_cgs", "unit_sfr_in_cgs",
    "convention",
)


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
    comoving : bool, optional
        If True (default, matching HaloCatalogue's default), 'pos', 'mass',
        and 'boxsize' are divided by the snapshot's HubbleParam -- i.e. the
        little-h factor is stripped, exactly as HaloCatalogue(comoving=True)
        does for GroupPos/SubhaloPos. Pass comoving=False to get values
        exactly as stored in the file (typically h^-1 code units for
        GADGET/Arepo/SWIFT). Set both adapters the same way when
        cross-matching particles against halo positions.
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
                 comoving: bool = True,
                 label: Optional[str] = None, **backend_kwargs):
        super().__init__(path=path, fileformat=fileformat, label=label)
        backend_kwargs.setdefault("loglevel", 30)   # WARNING: quiet library
        self._backend = SnapshotTools(snapfileformat=fileformat,
                                      **backend_kwargs)
        self._convention = convention
        self._comoving = comoving
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
        h0 = _scalar(getattr(d, "hubble_param", 0.0)) or None
        boxsize = _scalar(getattr(d, "box_size", 0.0)) or None

        # strip the little-h factor, matching HaloCatalogue(comoving=True)'s
        # GroupPos/SubhaloPos conversion -- out-of-place, so self._backend's
        # own arrays (reachable via .backend/.data) are left untouched.
        strip_h = bool(self._comoving and h0)
        if strip_h:
            for key in ("pos", "mass"):
                if key in self._columns:
                    self._columns[key] = self._columns[key] / h0
            boxsize = boxsize / h0 if boxsize is not None else None

        length_unit = "code (kpc-family, see backend unit_* attrs)"
        mass_unit = "code (1e10 Msol-family)"
        if h0:
            suffix = ", h-free" if strip_h else ", per h"
            length_unit += suffix
            mass_unit += suffix

        def _unwrap(x):
            # HDF5 attrs often arrive as 0-d/1-element arrays; reduce those
            # to plain Python scalars but leave real arrays (mass_table,
            # num_part_total, ...) and non-numeric values (strings, bools)
            # alone -- unlike _scalar(), this never assumes a numeric dtype.
            if isinstance(x, np.ndarray) and (x.ndim == 0 or x.size == 1):
                return x.item()
            return x

        # raw (un-h-stripped) header/config scalars, for fields that don't
        # have a standardised meta key -- parity with HaloCatalogue's
        # meta["native_meta"].
        native_meta = {
            attr: _unwrap(getattr(d, attr))
            for attr in _HEADER_ATTRS if getattr(d, attr, None) is not None
        }

        self.meta.update({
            "scale_factor": a,
            "redshift": 1.0 / a - 1.0 if a > 0 else None,
            "boxsize": boxsize,
            "h0": h0,
            "comoving": strip_h,
            "omega_0": _scalar(getattr(d, "omega_0", 0.0)) or None,
            "omega_lambda": _scalar(getattr(d, "omega_lambda", 0.0)) or None,
            "num_part_total": np.asarray(getattr(d, "num_part_total", [])),
            "convention": self._convention or self._backend.convention,
            "native_meta": native_meta,
            # plain strings (no astropy) -- GADGET-family defaults
            "units": {"length": length_unit,
                      "mass": mass_unit,
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
