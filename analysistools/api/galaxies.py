#!/usr/bin/env python3
"""
analysistools.api.galaxies
--------------------------
GalaxyCatalogue: the Dataset adapter for SHARK galaxy catalogues.

Two construction routes:

* file-backed -- one ``galaxies.hdf5`` file:
      gals = at.load(".../199/0/galaxies.hdf5")
* model-backed -- a SharkModel frozen at one redshift (fields fetched
  lazily through the model's cache):
      gals = GalaxyCatalogue.from_model(model, redshift=0.0)
      gals = model.at(0.0)                      # same thing

Standardised fields: pos, vel (stacked from *_x/y/z), mass/mstar
(mstars_disk+bulge), mgas, sfr, id, type, mhalo, msubhalo, mbh, radius --
plus every native SHARK column.
"""
from __future__ import annotations

from typing import Any, Dict, Optional

import numpy as np

from .dataset import Dataset


class GalaxyCatalogue(Dataset):
    """Dataset over one epoch of a SHARK galaxy catalogue."""

    kind = "galaxies"

    def __init__(self, path: str, label: Optional[str] = None):
        super().__init__(path=path, fileformat="SHARK", label=label)
        self._model = None
        self._model_z: Optional[float] = None

    # ------------------------------------------------------------------
    # Construction from a SharkModel (multi-epoch -> frozen view)
    # ------------------------------------------------------------------

    @classmethod
    def from_model(cls, model, redshift: float,
                   label: Optional[str] = None) -> "GalaxyCatalogue":
        """Freeze a SharkModel at one redshift as a GalaxyCatalogue."""
        self = cls(path=getattr(model, "model_dir", ""),
                   label=label or f"{getattr(model, 'label', 'shark')}"
                                  f"@z={redshift:g}")
        self._model = model
        self._model_z = float(redshift)
        return self

    # ------------------------------------------------------------------

    def _load(self) -> None:
        if self._model is not None:
            self._load_from_model()
        else:
            self._load_from_file()

    def _load_from_file(self) -> None:
        import h5py
        with h5py.File(self.path, "r") as f:
            gals = f["galaxies"]
            self._columns = {k: gals[k][()] for k in gals.keys()
                             if isinstance(gals[k], h5py.Dataset)}
            h0 = float(f["cosmology/h"][()]) if "cosmology/h" in f else None
            vol = (float(f["run_info/effective_volume"][()])
                   if "run_info/effective_volume" in f else None)
            z = (float(f["run_info/redshift"][()])
                 if "run_info/redshift" in f else None)
        self.meta.update({
            "redshift": z, "h0": h0, "volume": vol,
            "units": {"length": "cMpc/h", "mass": "Msun/h",
                      "velocity": "km/s", "sfr": "Msun/Gyr/h"},
        })

    def _load_from_model(self) -> None:
        z = self._model_z
        # establish table length with one cheap field
        self._columns["id_galaxy"] = np.asarray(
            self._model.get("id_galaxy", z))
        try:
            meta = self._model.get_meta(z)
        except Exception:
            meta = {}
        self.meta.update({
            "redshift": z,
            "h0": meta.get("h0"),
            "volume": meta.get("volume"),
            "units": {"length": "cMpc/h", "mass": "Msun/h",
                      "velocity": "km/s", "sfr": "Msun/Gyr/h"},
        })

    # ------------------------------------------------------------------
    # Lazy per-field fetch (model-backed only)
    # ------------------------------------------------------------------

    def _fetch_native(self, field: str) -> Optional[np.ndarray]:
        if self._model is None:
            return None
        # try the name as-is (logical SharkModel name), then reverse-map a
        # native HDF5 dataset name to its logical key (e.g. sfr_burst ->
        # sfr_bulge in GALAXY_FIELDS).
        names = [field]
        try:
            from ..shark.model import GALAXY_FIELDS
            names += [logical for logical, (_, dset) in GALAXY_FIELDS.items()
                      if dset == field and logical != field]
        except ImportError:
            pass
        for name in names:
            try:
                return np.asarray(self._model.get(name, self._model_z))
            except (KeyError, ValueError):
                continue
            except OSError:
                return None
        return None

    @property
    def fields(self):
        """
        Available field names: the base dataset's fields plus, for
        model-backed catalogues, all logical SHARK galaxy field names.

        Returns
        -------
        list of str, sorted.
        """
        base = set(super().fields)
        if self._model is not None:
            try:
                from ..shark.model import GALAXY_FIELDS
                base |= set(GALAXY_FIELDS)
            except ImportError:
                pass
        return sorted(base)

    @property
    def model(self):
        """The underlying SharkModel (model-backed catalogues), or None."""
        return self._model
