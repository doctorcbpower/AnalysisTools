#!/usr/bin/env python3
"""
analysistools.halo_model
------------------------
HaloModel: lazy, cached access to halo catalogues across multiple snapshots,
mirroring the SharkModel interface from the shark package.

Supports SUBFIND, AHF, VELOCIraptor, and SWIFT_FOF formats via the existing
HaloTools readers.

Typical usage
-------------
    from analysistools.halo_model import HaloModel

    catalogues = {
        0.0: "./catalogues/snap_200/fof_output.hdf5",
        1.0: "./catalogues/snap_151/fof_output.hdf5",
        2.0: "./catalogues/snap_110/fof_output.hdf5",
    }

    hm = HaloModel(
        catalogues  = catalogues,
        fileformat  = "SWIFT_FOF",
        label       = "CDM",
        colour      = "#e41a1c",
        comoving    = True,
        standardise = True,
    )

    # Access is lazy — file is only read on first request
    mvir = hm.get("Group_M_Crit200", redshift=0.0)

    # Or use derived-quantity methods
    pos  = hm.position(redshift=0.0)
    mvir = hm.mvir(redshift=0.0)

    # Works identically with a list of models for comparison plots
    models = [
        HaloModel(catalogues_cdm, "SWIFT_FOF", label="CDM",        colour="#e41a1c"),
        HaloModel(catalogues_wdm, "SWIFT_FOF", label="WDM 0.5 keV", colour="#377eb8"),
    ]
"""

from __future__ import annotations

import os
import logging
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np

from .halo_tools import HaloTools, FORMAT_READERS

# ---------------------------------------------------------------------------
# Redshift table helper
# ---------------------------------------------------------------------------

class _SnapshotTable:
    """Maps redshift floats to catalogue file paths.

    Accepts either:
    * a dict  {redshift: filepath}  — exact lookup with nearest-neighbour fallback
    * a callable  redshift -> filepath  — full user control
    """

    def __init__(self, source: Union[Dict[float, str], Callable[[float], str]]):
        if callable(source):
            self._fn = source
            self._redshifts = None
            self._paths = None
        else:
            zs = sorted(source.keys())
            self._redshifts = np.array(zs)
            self._paths = [source[z] for z in zs]
            self._fn = None

    def __getitem__(self, redshift: float) -> str:
        if self._fn is not None:
            return self._fn(redshift)
        idx = int(np.argmin(np.abs(self._redshifts - redshift)))
        return self._paths[idx]

    @property
    def available_redshifts(self) -> Optional[np.ndarray]:
        return self._redshifts


# ---------------------------------------------------------------------------
# HaloModel
# ---------------------------------------------------------------------------

class HaloModel:
    """Lazy, cached interface to halo catalogues across multiple snapshots.

    Mirrors the SharkModel API so models of both types can be passed to the
    same Analysis and Plotter classes.

    Parameters
    ----------
    catalogues : dict {float: str} or callable float -> str
        Mapping from redshift to catalogue file path.  For a dict the nearest
        available redshift is used; for a callable full control is yours.
    fileformat : str or int
        Halo catalogue format — one of "SUBFIND", "AHF", "VELOCIraptor",
        "SWIFT_FOF" (or the integer codes 1–4).
    label : str, optional
        Human-readable name for plot legends.
    colour : str, optional
        Matplotlib colour string.  Assigned by Plotter if omitted.
    comoving : bool, optional
        If True, positions are returned in comoving coordinates.
    standardise : bool, optional
        If True, field names are normalised to the common schema via
        ``standardise_catalogue_names`` before caching.
    loglevel : int, optional
        Logging level (default logging.WARNING — less verbose than HaloTools
        default so notebooks stay quiet).

    Examples
    --------
    >>> hm = HaloModel(
    ...     {0.0: "./snap200/fof.hdf5", 1.0: "./snap151/fof.hdf5"},
    ...     fileformat="SWIFT_FOF", label="CDM", comoving=True,
    ... )
    >>> mvir = hm.get("Group_M_Crit200", redshift=0.0)
    >>> pos  = hm.position(redshift=0.0)
    """

    def __init__(
        self,
        catalogues:  Union[Dict[float, str], Callable[[float], str]],
        fileformat:  Union[str, int],
        label:       Optional[str]  = None,
        colour:      Optional[str]  = None,
        comoving:    bool           = True,
        standardise: bool           = False,
        loglevel:    int            = logging.WARNING,
    ):
        self._table      = _SnapshotTable(catalogues)
        self.fileformat  = fileformat
        self.label       = label or "HaloModel"
        self.colour      = colour
        self.comoving    = comoving
        self.standardise = standardise

        # Internal HaloTools instance — reused for every load
        self._ht = HaloTools(comoving_units=comoving, loglevel=loglevel)

        # Cache: redshift_key -> {"halos": dict, "subhalos": dict, "meta": dict}
        self._cache: Dict[float, Dict[str, Any]] = {}

    # ------------------------------------------------------------------
    # Core access — mirrors SharkModel.get / get_meta
    # ------------------------------------------------------------------

    def get(self, field: str, redshift: float) -> np.ndarray:
        """Return a halo field array at the snapshot nearest *redshift*.

        Searches halos first, then subhalos.  The result is cached so
        repeated calls are free.

        Parameters
        ----------
        field : str
            Field name as it appears in the catalogue (raw or standardised,
            depending on how the model was constructed).
        redshift : float
            Target redshift; matched to the nearest available snapshot.

        Returns
        -------
        arr : ndarray
        """
        snap = self._load(redshift)
        for store in ("halos", "subhalos"):
            data = snap[store]
            if data and field in data:
                return data[field]
        available = sorted(
            list(snap["halos"].keys()) + list(snap["subhalos"].keys())
        )
        raise KeyError(
            f"Field '{field}' not found in catalogue at z≈{redshift:.2f}. "
            f"Available: {available}"
        )

    def get_halos(self, redshift: float) -> Dict[str, np.ndarray]:
        """Return the full halo dict at the snapshot nearest *redshift*."""
        return self._load(redshift)["halos"]

    def get_subhalos(self, redshift: float) -> Dict[str, np.ndarray]:
        """Return the full subhalo dict at the snapshot nearest *redshift*."""
        return self._load(redshift)["subhalos"]

    def get_meta(self, redshift: float) -> Dict[str, Any]:
        """Return catalogue metadata (header + cosmology) for the snapshot nearest *redshift*.

        Always includes ``h0`` and ``redshift`` keys for compatibility with
        SharkModel.get_meta().
        """
        return self._load(redshift)["meta"]

    def available_fields(self, redshift: float) -> List[str]:
        """Return all field names available at the snapshot nearest *redshift*."""
        snap = self._load(redshift)
        return sorted(list(snap["halos"].keys()) + list(snap["subhalos"].keys()))

    def available_redshifts(self) -> Optional[np.ndarray]:
        """Return the redshifts for which catalogue files are registered.

        Returns None if the catalogue mapping was provided as a callable.
        """
        return self._table.available_redshifts

    # ------------------------------------------------------------------
    # Derived quantities — mirrors SharkModel derived methods
    # ------------------------------------------------------------------

    def mvir(self, redshift: float) -> np.ndarray:
        """Halo virial mass.  Tries common field names across formats."""
        return self._first_available(
            redshift,
            ["Group_M_Crit200", "Mass_200crit", "Mvir", "m200c",
             "Group_M_Mean200", "Mass_tot", "Masses"],
            label="virial mass",
        )

    def rvir(self, redshift: float) -> np.ndarray:
        """Halo virial radius."""
        return self._first_available(
            redshift,
            ["Group_R_Crit200", "R_200crit", "Rvir", "r200c",
             "Group_R_Mean200", "Radii"],
            label="virial radius",
        )

    def vmax(self, redshift: float) -> np.ndarray:
        """Maximum circular velocity."""
        return self._first_available(
            redshift,
            ["VMax", "vmax", "SubhaloVmax", "Vmax"],
            label="Vmax",
        )

    def position(self, redshift: float) -> np.ndarray:
        """Halo positions as array of shape (n, 3).

        Tries 3-column array fields first, then assembles from x/y/z scalars.
        """
        snap = self._load(redshift)
        halos = snap["halos"]

        for key in ["GroupPos", "Xcminpot", "Xc", "position", "Centres"]:
            if key in halos:
                arr = halos[key]
                if arr.ndim == 2 and arr.shape[1] == 3:
                    return arr
                if arr.ndim == 1:
                    # Some formats store a flat (3N,) array
                    return arr.reshape(-1, 3)

        # Try assembling from x/y/z component fields
        for x, y, z in [("Xc", "Yc", "Zc"), ("x", "y", "z"),
                         ("pos_x", "pos_y", "pos_z")]:
            if x in halos and y in halos and z in halos:
                return np.column_stack([halos[x], halos[y], halos[z]])

        raise KeyError(
            f"No position field found in catalogue at z≈{redshift:.2f}. "
            f"Available: {self.available_fields(redshift)}"
        )

    def velocity(self, redshift: float) -> np.ndarray:
        """Halo bulk velocities as array of shape (n, 3)."""
        snap = self._load(redshift)
        halos = snap["halos"]

        for key in ["GroupVel", "Vel", "velocity"]:
            if key in halos:
                arr = halos[key]
                if arr.ndim == 2 and arr.shape[1] == 3:
                    return arr
                if arr.ndim == 1:
                    return arr.reshape(-1, 3)

        for vx, vy, vz in [("VXc", "VYc", "VZc"), ("vel_x", "vel_y", "vel_z")]:
            if vx in halos and vy in halos and vz in halos:
                return np.column_stack([halos[vx], halos[vy], halos[vz]])

        raise KeyError(
            f"No velocity field found in catalogue at z≈{redshift:.2f}."
        )

    def npart(self, redshift: float) -> np.ndarray:
        """Number of particles per halo."""
        return self._first_available(
            redshift,
            ["GroupLen", "npart", "Npart", "NumPart", "n_part"],
            label="particle count",
        )

    def concentration(self, redshift: float) -> np.ndarray:
        """NFW concentration parameter."""
        return self._first_available(
            redshift,
            ["cNFW", "cnfw", "concentration", "Concentration"],
            label="concentration",
        )

    def h0(self, redshift: float) -> float:
        """Hubble parameter h for the snapshot nearest *redshift*."""
        meta = self.get_meta(redshift)
        for key in ("h0", "HubbleParam", "h"):
            if key in meta and meta[key] is not None:
                return float(meta[key])
        raise KeyError(f"Hubble parameter not found in metadata: {meta}")

    def volume(self, redshift: float) -> float:
        """Box volume in Mpc^3 inferred from metadata BoxSize and h0."""
        meta = self.get_meta(redshift)
        box = meta.get("BoxSize")
        if box is None:
            raise KeyError("BoxSize not found in catalogue metadata.")
        h = self.h0(redshift)
        # BoxSize is typically in Mpc/h — convert to Mpc^3
        return float(box / h) ** 3

    def is_loaded(self, redshift: float) -> bool:
        """True if the catalogue for the snapshot nearest *redshift* is cached."""
        z_key = self._z_key(redshift)
        return z_key in self._cache

    def clear_cache(self) -> None:
        """Release all cached arrays, freeing memory."""
        self._cache.clear()

    def preload(self, redshifts: Sequence[float]) -> "HaloModel":
        """Eagerly load catalogue files for all *redshifts* into the cache.

        Returns self for chaining.
        """
        for z in redshifts:
            self._load(z)
        return self

    # ------------------------------------------------------------------
    # Repr / summary
    # ------------------------------------------------------------------

    def __repr__(self) -> str:
        n = len(self._cache)
        return (
            f"HaloModel(label={self.label!r}, "
            f"format={self.fileformat!r}, "
            f"snapshots_cached={n})"
        )

    def summary(self, redshift: float) -> None:
        """Print a concise summary for the snapshot nearest *redshift*."""
        snap  = self._load(redshift)
        meta  = snap["meta"]
        halos = snap["halos"]
        subs  = snap["subhalos"]
        nh = len(next(iter(halos.values()))) if halos else 0
        ns = len(next(iter(subs.values()))) if subs else 0

        print(f"HaloModel: {self.label}")
        print(f"  Format   : {self.fileformat}")
        print(f"  z ≈      : {redshift:.2f}")
        print(f"  Halos    : {nh:,}")
        print(f"  Subhalos : {ns:,}")
        if meta:
            for k, v in meta.items():
                print(f"  {k:20s}: {v}")

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _z_key(self, redshift: float) -> float:
        """Canonical redshift key: the z stored in the snapshot table."""
        zs = self._table.available_redshifts
        if zs is not None:
            return float(zs[int(np.argmin(np.abs(zs - redshift)))])
        return float(redshift)

    def _load(self, redshift: float) -> Dict[str, Any]:
        """Return cached snapshot data, loading from disk if necessary."""
        z_key = self._z_key(redshift)
        if z_key in self._cache:
            return self._cache[z_key]

        filepath = self._table[redshift]
        meta, halos, subhalos = self._ht.read_catalogue(
            filename    = filepath,
            fileformat  = self.fileformat,
            standardise = self.standardise,
        )

        # Normalise metadata to always include h0 and redshift keys
        meta = dict(meta)
        for h_key in ("HubbleParam", "h", "hubble"):
            if h_key in meta and "h0" not in meta:
                meta["h0"] = meta[h_key]
        if "redshift" not in meta and "Redshift" in meta:
            meta["redshift"] = meta["Redshift"]

        self._cache[z_key] = {
            "halos":    halos    or {},
            "subhalos": subhalos or {},
            "meta":     meta,
        }
        return self._cache[z_key]

    def _first_available(
        self,
        redshift: float,
        candidates: List[str],
        label: str = "field",
    ) -> np.ndarray:
        """Return the first candidate field name that exists in the catalogue."""
        snap = self._load(redshift)
        for store in ("halos", "subhalos"):
            data = snap[store]
            for name in candidates:
                if data and name in data:
                    return data[name]
        raise KeyError(
            f"No {label} field found at z≈{redshift:.2f}. "
            f"Tried: {candidates}. "
            f"Available: {self.available_fields(redshift)}"
        )
