"""
shark.photometry.cosmology
==========================
Cosmology helpers for the photometry pipeline.

SharkModel already reads H0/OmegaM/OmegaB from the HDF5 run_info group and
stores them in _meta.  This module builds an astropy cosmology object from
those stored values so the pipeline can compute ages without re-reading HDF5.

Usage
-----
>>> cosmo = SharkCosmology.from_model(model, redshift=0.1)
>>> tage  = cosmo.age_at_z(0.1)
"""

from __future__ import annotations

import numpy as np
from astropy.cosmology import FlatLambdaCDM
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ..model import SharkModel


class SharkCosmology:
    """
    Thin astropy wrapper initialised from a SharkModel instance.

    Parameters
    ----------
    H0 : float   Hubble constant in km/s/Mpc  (little-h * 100)
    Om0 : float  Total matter density parameter
    Ob0 : float  Baryon density parameter

    Prefer the factory method ``SharkCosmology.from_model`` over direct
    construction.
    """

    def __init__(self, H0: float, Om0: float, Ob0: float):
        self.H0  = H0
        self.Om0 = Om0
        self.Ob0 = Ob0
        self.cosmo = FlatLambdaCDM(H0=H0, Om0=Om0, Ob0=Ob0)

    @classmethod
    def from_model(cls, model: "SharkModel", redshift: float) -> "SharkCosmology":
        """
        Build a SharkCosmology from the cosmological parameters stored
        inside a SharkModel.

        SharkModel._meta stores h0 (little-h) and vol but not OmegaM/OmegaB
        directly.  We read OmegaM and OmegaB from the HDF5 cosmology group
        of the first subvolume file at the relevant snapshot — this is a
        single cheap read, cached implicitly by the OS after the first call.

        Parameters
        ----------
        model : SharkModel
        redshift : float
            Any redshift whose snapshot has already been (or will be) loaded.
        """
        import h5py
        import os

        snapshot    = int(model.redshift_table[redshift])
        first_subvol = sorted(model.subvols)[0]
        fname = os.path.join(
            model.model_dir, str(snapshot), str(first_subvol), "galaxies.hdf5"
        )
        with h5py.File(fname, "r") as f:
            h0      = float(f["cosmology/h"][()])
            omega_m = float(f["cosmology/OmegaM"][()])
            omega_b = float(f["cosmology/OmegaB"][()])

        return cls(H0=h0 * 100.0, Om0=omega_m, Ob0=omega_b)

    def age_at_z(self, z: float) -> float:
        """Age of the Universe [Gyr] at redshift z."""
        return float(self.cosmo.age(z).value)

    def lookback_at_z(self, z: float) -> float:
        """Lookback time [Gyr] to redshift z."""
        return float(self.cosmo.lookback_time(z).value)
