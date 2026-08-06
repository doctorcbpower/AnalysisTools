#!/usr/bin/env python3
"""
analysistools.api.dataset
-------------------------
The common Dataset interface: one field vocabulary and one access pattern
across snapshots, halo catalogues, merger trees, and galaxy catalogues.

See DEVELOPMENT.md sections 3.1-3.2. Phase 1 implements the base class and
the field-alias registry; adapters live in the sibling modules.

Design decisions (DEVELOPMENT.md section 7):
* plain-string units in meta["units"], no astropy;
* per-species views via attributes (ds.dm, ds.gas, ...) on snapshots;
* select()/__len__ are optional per kind (mesh data may not support them).

Author: C. Power (design), implemented 2026-07.
"""
from __future__ import annotations

import logging
from typing import Any, Callable, Dict, List, Optional, Sequence, Union

import numpy as np

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Field-alias registry
# ---------------------------------------------------------------------------
# Per-kind aliases mapping standardised names -> native column names (tried
# in order) or a callable(dataset) for derived fields. The registry extends
# STANDARD_KEYS from halo_tools_standardise_names.py across all data kinds.

FIELD_ALIASES: Dict[str, Dict[str, Union[Sequence[str], Callable]]] = {
    "snapshot": {
        "pos":  ("pos",),
        "vel":  ("vel",),
        "mass": ("mass",),
        "id":   ("pids",),
        "pids": ("pids",),          # muscle memory
        "type": ("ptype",),
        "ptype": ("ptype",),
        "potential": ("potential",),
        "density": ("rho",),
        "internal_energy": ("u",),
        "sfr":  ("sfr",),
        "groupid": ("groupid",),
    },
    "halos": {
        "pos":    ("pos",),
        "vel":    ("vel",),
        "mass":   ("mass",),
        "id":     ("halo_id",),
        "halo_id": ("halo_id",),
        "radius": ("radius",),
        "num_part": ("num_part",),
        "vmax":   ("Vmax", "vmax_subhalo", "SubhaloVmax"),
    },
    "tree": {
        "pos":  ("Pos",),
        "vel":  ("Vel",),
        "mass": ("Mass",),
        "redshift": ("Redshift",),
        "snapnum":  ("SnapNum",),
        "is_subhalo": ("IsSubhalo",),
        "host_id":  ("HostID",),
    },
    "galaxies": {
        # SHARK catalogues. Callables build derived/stacked fields from
        # native columns via ds._resolve() (see note in Dataset._resolve).
        "pos": lambda ds: np.stack([ds._resolve("position_x"),
                                    ds._resolve("position_y"),
                                    ds._resolve("position_z")], axis=1),
        "vel": lambda ds: np.stack([ds._resolve("velocity_x"),
                                    ds._resolve("velocity_y"),
                                    ds._resolve("velocity_z")], axis=1),
        "mass": lambda ds: (ds._resolve("mstars_disk")
                            + ds._resolve("mstars_bulge")),
        "mstar": lambda ds: (ds._resolve("mstars_disk")
                             + ds._resolve("mstars_bulge")),
        "mgas": lambda ds: (ds._resolve("mgas_disk")
                            + ds._resolve("mgas_bulge")),
        "sfr": lambda ds: (ds._resolve("sfr_disk")
                           + ds._resolve("sfr_burst")),
        "id": ("id_galaxy",),
        "type": ("type",),
        "mhalo": ("mvir_hosthalo",),
        "msubhalo": ("mvir_subhalo",),
        "mbh": ("m_bh", "mbh"),
        "radius": ("rstar_disk",),
    },
    "satellites": {
        # Master science catalogues (analysistools.catalogue,
        # DEVELOPMENT.md Phase 6). Field names are the catalogue's own
        # (already standardised at build time -- see
        # docs/dorcha_master_catalogue_design.md section 2.2), so aliases
        # here just bridge the cross-kind vocabulary shared with halos/
        # galaxies for code that treats a catalogue like any other Dataset
        # (e.g. api/plotting.py's mass_function).
        "id":     ("SatelliteID",),
        "mass":   ("StellarMass",),
        "halo_id": ("HostHaloID",),
        "radius": ("HalfLightRadius",),
    },
}


class Dataset:
    """
    Uniform, mapping-style access to one table of simulation data.

    Concrete adapters (SnapshotDataset, HaloCatalogue, ...) populate
    ``_columns`` (native name -> ndarray) lazily and fill ``meta``. Field
    lookup order for ``ds[name]``:

    1. exact native column name;
    2. standardised alias from FIELD_ALIASES[kind];
    3. KeyError listing what is available.

    Selections (``ds.select(...)``) return lightweight views: same class of
    object, same columns, with a row-index applied on access.
    """

    kind: str = "dataset"

    def __init__(self, path: str = "", fileformat: str = "",
                 label: Optional[str] = None):
        self.path = path
        self.fileformat = fileformat
        self.label = label or (path.rsplit("/", 1)[-1] if path else self.kind)
        self.meta: Dict[str, Any] = {"units": {}}
        self._columns: Dict[str, np.ndarray] = {}
        self._index: Optional[np.ndarray] = None   # None = full table
        self._loaded = False

    # ------------------------------------------------------------------
    # Loading (lazy) -- adapters override _load()
    # ------------------------------------------------------------------

    def _load(self) -> None:
        """Read data from disk into self._columns / self.meta."""
        raise NotImplementedError

    def _ensure_loaded(self) -> None:
        if not self._loaded:
            self._load()
            self._loaded = True

    def preload(self) -> "Dataset":
        """Force the read now (otherwise it happens on first field access)."""
        self._ensure_loaded()
        return self

    # ------------------------------------------------------------------
    # Field access
    # ------------------------------------------------------------------

    def _fetch_native(self, field: str) -> Optional[np.ndarray]:
        """Hook for adapters with an on-demand backend (e.g. a SharkModel):
        fetch one native column by name, or return None if unavailable.
        Fetched columns are cached in self._columns."""
        return None

    def _resolve(self, field: str) -> np.ndarray:
        self._ensure_loaded()
        if field in self._columns:
            return self._columns[field]
        alias = FIELD_ALIASES.get(self.kind, {}).get(field)
        if alias is not None:
            if callable(alias):
                # NOTE: callables must build from ds._resolve() (unindexed
                # arrays); __getitem__ applies any view index afterwards.
                return alias(self)
            for native in alias:
                if native in self._columns:
                    return self._columns[native]
        # on-demand backend fetch (field itself, then its aliases)
        candidates = [field] + (list(alias) if alias is not None
                                and not callable(alias) else [])
        for native in candidates:
            arr = self._fetch_native(native)
            if arr is not None:
                self._columns[native] = arr
                return arr
        raise KeyError(
            f"'{field}' not found in this {self.kind} dataset. "
            f"Native columns: {sorted(self._columns)}")

    def __getitem__(self, field: str) -> np.ndarray:
        arr = self._resolve(field)
        if arr is None:
            raise KeyError(f"Field '{field}' is present but empty (None) in "
                           f"this {self.kind} dataset.")
        if self._index is not None and isinstance(arr, np.ndarray) \
                and arr.ndim >= 1:
            return arr[self._index]
        return arr

    def get(self, field: str, default=None):
        """
        Like `self[field]` but returns `default` instead of raising when
        the field is missing.

        Parameters
        ----------
        field : str
        default : optional
            Value returned if `field` cannot be resolved. Default None.
        """
        try:
            return self[field]
        except KeyError:
            return default

    def __contains__(self, field: str) -> bool:
        try:
            self._resolve(field)
            return True
        except KeyError:
            return False

    def __getattr__(self, name: str):
        # Attribute sugar: ds.pos == ds["pos"] (only reached when normal
        # attribute lookup fails; never shadows real attributes/methods).
        if name.startswith("_"):
            raise AttributeError(name)
        try:
            return self[name]
        except KeyError:
            raise AttributeError(
                f"'{type(self).__name__}' has no attribute or field '{name}'")

    @property
    def fields(self) -> List[str]:
        """Standardised aliases that resolve, plus all native column names."""
        self._ensure_loaded()
        names = set(self._columns)
        for std in FIELD_ALIASES.get(self.kind, {}):
            try:
                self._resolve(std)
                names.add(std)
            except KeyError:
                pass
        return sorted(names)

    def __len__(self) -> int:
        self._ensure_loaded()
        if self._index is not None:
            return len(self._index)
        for v in self._columns.values():
            if isinstance(v, np.ndarray) and v.ndim >= 1:
                return len(v)
        return 0

    # ------------------------------------------------------------------
    # Selection
    # ------------------------------------------------------------------

    def select(self, mask: Optional[np.ndarray] = None, *,
               centre: Optional[Sequence[float]] = None,
               size: Optional[float] = None,
               geometry: str = "spherical",
               periodic: Optional[bool] = None,
               **cuts) -> "Dataset":
        """
        Return a view of this dataset restricted to a subset of rows.

        Parameters
        ----------
        mask : bool or int ndarray, optional
            Explicit row mask/index (in this view's frame).
        centre, size, geometry, periodic :
            Geometric cut on ``self["pos"]``: 'spherical' (size = radius) or
            'cubic' (size = half-width). Periodic wrapping uses
            meta["boxsize"] and defaults to on when a boxsize is known.
        **cuts :
            Attribute cuts as field=(lo, hi); either bound may be None.
            e.g. ``ds.select(mass=(1e12, None))``.

        Returns
        -------
        Dataset view (chainable; anything accepting a Dataset accepts it).
        """
        self._ensure_loaded()
        n = len(self)
        keep = np.ones(n, dtype=bool)

        if mask is not None:
            mask = np.asarray(mask)
            if mask.dtype == bool:
                keep &= mask
            else:                     # integer indices
                m = np.zeros(n, dtype=bool)
                m[mask] = True
                keep &= m

        if centre is not None:
            if size is None:
                raise ValueError("select(): 'size' is required with 'centre'.")
            pos = self["pos"]
            delta = pos - np.asarray(centre, dtype=float)
            boxsize = self.meta.get("boxsize")
            if periodic is None:
                periodic = boxsize is not None
            if periodic and boxsize:
                delta = (delta + 0.5 * boxsize) % boxsize - 0.5 * boxsize
            if geometry == "spherical":
                keep &= np.sum(delta**2, axis=1) < size**2
            elif geometry in ("cubic", "cuboidal"):
                keep &= np.all(np.abs(delta) < size, axis=1)
            else:
                raise ValueError(f"Unknown geometry '{geometry}' "
                                 f"(use 'spherical' or 'cubic').")

        for field, bounds in cuts.items():
            arr = self[field]
            try:
                lo, hi = bounds
            except TypeError:
                raise ValueError(
                    f"select(): cut on '{field}' must be a (lo, hi) pair; "
                    f"got {bounds!r}. Use None for an open bound.")
            if lo is not None:
                keep &= arr >= lo
            if hi is not None:
                keep &= arr < hi

        return self._make_view(np.flatnonzero(keep))

    def _make_view(self, rows: np.ndarray) -> "Dataset":
        """Build a view sharing columns/meta, with `rows` (frame of this
        view) composed onto any existing index."""
        view = object.__new__(type(self))
        view.__dict__.update(self.__dict__)
        if self._index is not None:
            view._index = self._index[rows]
        else:
            view._index = rows
        return view

    @property
    def is_view(self) -> bool:
        """
        Whether this dataset is a row-indexed view onto another dataset's
        columns (e.g. as produced by `select`).

        Returns
        -------
        bool
        """
        return self._index is not None

    # ------------------------------------------------------------------
    # Reporting
    # ------------------------------------------------------------------

    def summary(self) -> None:
        """Print a concise human-readable description."""
        self._ensure_loaded()
        print(f"{type(self).__name__} ({self.kind})"
              + (" [view]" if self.is_view else ""))
        print(f"  label:  {self.label}")
        if self.path:
            print(f"  path:   {self.path}")
        if self.fileformat:
            print(f"  format: {self.fileformat}")
        print(f"  rows:   {len(self):,}")
        for key in ("redshift", "boxsize", "h0", "scale_factor"):
            if self.meta.get(key) is not None:
                print(f"  {key}: {self.meta[key]}")
        print(f"  fields: {', '.join(self.fields)}")

    def __repr__(self) -> str:
        z = self.meta.get("redshift")
        ztxt = f", z={z:.3g}" if isinstance(z, (int, float)) else ""
        return (f"<{type(self).__name__} '{self.label}' kind={self.kind} "
                f"n={len(self) if self._loaded else '?'}{ztxt}"
                f"{' view' if self.is_view else ''}>")
