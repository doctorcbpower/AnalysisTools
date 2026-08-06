#!/usr/bin/env python3
"""
analysistools.api.epoch
-----------------------
EpochModel: the formal interface shared by HaloModel and SharkModel
(lazy multi-redshift access: ``get(field, redshift)``, ``get_meta``,
``preload``, ``clear_cache``, ``label``/``colour``), plus ``at(model, z)``
-- the bridge from the multi-epoch world to single-epoch Datasets.

``HaloModel.at(z)`` and ``SharkModel.at(z)`` delegate here, so everything
downstream (selection, plotting, cross-matching) only ever needs to
understand Dataset.
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any, Dict, Optional

import numpy as np

from .dataset import Dataset
from .galaxies import GalaxyCatalogue


class EpochModel(ABC):
    """Lazy, cached, multi-redshift data source.

    HaloModel and SharkModel are registered as virtual subclasses; use
    ``isinstance(obj, EpochModel)`` to accept either.
    """

    @abstractmethod
    def get(self, field: str, redshift: float) -> np.ndarray:
        """
        Fetch a field's values at the given redshift, loading/caching that
        epoch's data as needed.

        Parameters
        ----------
        field : str
        redshift : float

        Returns
        -------
        array
        """
        ...

    @abstractmethod
    def get_meta(self, redshift: float) -> Dict[str, Any]:
        """
        Fetch this model's metadata dict at the given redshift.

        Parameters
        ----------
        redshift : float

        Returns
        -------
        dict
        """
        ...

    def at(self, redshift: float) -> Dataset:
        """Freeze this model at one redshift as a Dataset view."""
        return at(self, redshift)


class EpochHaloView(Dataset):
    """HaloModel frozen at one redshift: a halos-kind Dataset."""

    kind = "halos"

    def __init__(self, model, redshift: float, label: Optional[str] = None):
        super().__init__(path="", fileformat=str(model.fileformat),
                         label=label or f"{model.label}@z={redshift:g}")
        self._model = model
        self._z = float(redshift)

    def _load(self) -> None:
        halos = self._model.get_halos(self._z)
        self._columns = {k: v for k, v in halos.items()
                         if isinstance(v, np.ndarray)}
        meta = {}
        try:
            meta = dict(self._model.get_meta(self._z))
        except Exception:
            pass
        self.meta.update({
            "redshift": self._z,
            "boxsize": meta.get("BoxSize", meta.get("boxsize")),
            "native_meta": meta,
        })
        for name, key in (("h0", "h0"), ("volume", "volume")):
            try:
                self.meta[key] = getattr(self._model, name)(self._z)
            except Exception:
                pass

    @property
    def subhalos(self) -> "EpochHaloView":
        """
        Subhalo table at this view's redshift, as its own EpochHaloView.

        Returns
        -------
        EpochHaloView

        Raises
        ------
        AttributeError
            If no subhalo table exists at this redshift.
        """
        subs = self._model.get_subhalos(self._z)
        if not subs:
            raise AttributeError("No subhalo table at this redshift.")
        view = EpochHaloView(self._model, self._z,
                             label=f"{self.label}:subhalos")
        view._columns = {k: v for k, v in subs.items()
                         if isinstance(v, np.ndarray)}
        view._loaded = True
        view.meta.update(self.meta if self._loaded else {"redshift": self._z})
        return view

    @property
    def model(self):
        """The underlying HaloModel."""
        return self._model


def at(model, redshift: float) -> Dataset:
    """Freeze any EpochModel at one redshift as a Dataset.

    HaloModel -> EpochHaloView (kind "halos");
    SharkModel -> GalaxyCatalogue (kind "galaxies").
    """
    # duck-typed dispatch: SharkModel has GALAXY_FIELDS-style get + model_dir
    if hasattr(model, "get_halos"):
        return EpochHaloView(model, redshift)
    if hasattr(model, "model_dir"):
        return GalaxyCatalogue.from_model(model, redshift)
    raise TypeError(f"Don't know how to freeze {type(model).__name__} "
                    f"at a redshift; expected HaloModel or SharkModel.")


def _register_models() -> None:
    """Register HaloModel/SharkModel as virtual EpochModel subclasses and
    give them .at(z) (additive; existing APIs untouched)."""
    from ..halo_model import HaloModel
    EpochModel.register(HaloModel)
    if not hasattr(HaloModel, "at"):
        HaloModel.at = lambda self, redshift: at(self, redshift)
        HaloModel.at.__doc__ = at.__doc__

    try:
        from ..shark.model import SharkModel
    except ImportError:      # optional deps missing
        return
    EpochModel.register(SharkModel)
    if not hasattr(SharkModel, "at"):
        SharkModel.at = lambda self, redshift: at(self, redshift)
        SharkModel.at.__doc__ = at.__doc__


_register_models()
