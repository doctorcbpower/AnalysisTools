#!/usr/bin/env python3
"""
analysistools.api.catalogue
-----------------------------
CatalogueDataset: the Dataset adapter for a master science catalogue
built by ``analysistools.catalogue`` (DEVELOPMENT.md Phase 6,
docs/master_catalogue.md).

Unlike the other adapters, the on-disk layout is nested groups
(``Satellites/HaloProperties/Mpeak``, ``Satellites/GalaxyProperties/
StellarMass``, ...) sharing one canonical row order keyed on
``SatelliteID`` (see docs/dorcha_master_catalogue_design.md section 2.1).
``_load`` flattens every leaf dataset under ``Satellites/`` into one
column namespace on first access, by its dataset name (last path
component) -- so ``cat["Mpeak"]``, ``cat["StellarMass"]`` work exactly
like any other Dataset, and every existing ``api/plotting.py`` function
(mass_function, overlay_points, profile) works unchanged on a built
catalogue. Per-field units/description/provenance/is_derived HDF5
attributes are exposed via :meth:`field_info`.

This module only reads; ``analysistools.catalogue.pipeline.WriteStage`` is
the only code permitted to write a catalogue file.
"""
from __future__ import annotations

from typing import Any, Dict, Optional

import logging

import numpy as np

from .dataset import Dataset

logger = logging.getLogger(__name__)


class CatalogueDataset(Dataset):
    """Dataset over one released master catalogue file."""

    kind = "satellites"

    def __init__(self, path: str, label: Optional[str] = None):
        super().__init__(path=path, fileformat="HDF5-catalogue", label=label)
        self._field_attrs: Dict[str, Dict[str, Any]] = {}

    def _load(self) -> None:
        import h5py

        columns: Dict[str, np.ndarray] = {}
        attrs: Dict[str, Dict[str, Any]] = {}

        def visit(name: str, obj: "h5py.HLObject") -> None:
            if not isinstance(obj, h5py.Dataset):
                return
            if not name.startswith("Satellites/"):
                return
            field_name = name.rsplit("/", 1)[-1]
            if field_name in columns:
                logger.warning(
                    "Duplicate field name '%s' at multiple paths in %s; "
                    "keeping the first one seen (%s).",
                    field_name, self.path, name)
                return
            columns[field_name] = obj[()]
            attrs[field_name] = dict(obj.attrs)

        with h5py.File(self.path, "r") as f:
            if "Satellites" not in f:
                raise ValueError(
                    f"'{self.path}' has no 'Satellites' group -- is this a "
                    f"master catalogue file? (analysistools.catalogue"
                    f".pipeline.WriteStage output expected.)")
            f.visititems(visit)

            meta_grp = f.get("Metadata")
            self.meta.update({
                "catalogue_version": (meta_grp.attrs.get("catalogue_version")
                                     if meta_grp is not None else None),
                "schema_version": (meta_grp.attrs.get("schema_version")
                                   if meta_grp is not None else None),
                "n_satellites": len(next(iter(columns.values())))
                               if columns else 0,
            })
            if "Cosmology" in f:
                cosmo = f["Cosmology"]
                self.meta["h0"] = float(cosmo.attrs.get("H0", 0)) / 100.0 \
                    if "H0" in cosmo.attrs else None
            self.meta["units"] = {
                name: a.get("units") for name, a in attrs.items()
                if a.get("units")
            }

        self._columns = columns
        self._field_attrs = attrs

    def field_info(self, field: str) -> Dict[str, Any]:
        """Units/description/provenance/is_derived for one field, as
        stored in its HDF5 attributes at write time."""
        self._ensure_loaded()
        try:
            return self._field_attrs[field]
        except KeyError:
            raise KeyError(
                f"No attribute metadata for '{field}' (native columns: "
                f"{sorted(self._field_attrs)})")
