"""
shark.photometry.photometry
============================
Top-level pipeline: SharkModel → absolute magnitudes.

PhotometryPipeline accepts a SharkModel instance (which already owns all
HDF5 reading, caching, and subvolume handling) and drives the SPS convolution
via SPSEngine.

The SFH time axis used for FSPS comes from ``lbt_mean`` in
star_formation_histories.hdf5 (accessed via SharkModel.get_sfh_meta),
converted to age-of-Universe using the simulation cosmology.

Example
-------
>>> from shark.model import SharkModel
>>> from shark.common import _redshift_table, parse_subvolumes
>>> from shark.photometry import PhotometryPipeline
>>> import numpy as np

>>> rt    = _redshift_table("redshift_list.txt")
>>> sv    = parse_subvolumes("0-63")
>>> model = SharkModel("./CDM/base_model", rt, sv, label="CDM")

>>> pipe  = PhotometryPipeline(model, z_obs=0.1, bands=["v"])
>>> M_V   = pipe.abs_mag_V(gal_indices=np.arange(1000))
>>> ML_V  = pipe.mass_to_light(gal_indices=np.arange(1000))
"""

from __future__ import annotations

import numpy as np
from typing import Optional

from .io import galaxy_data, sfh_ages_from_model
from .sps import SPSEngine

from ..model import SharkModel


class PhotometryPipeline:
    """
    End-to-end pipeline from SharkModel SFH output to absolute magnitudes.

    Parameters
    ----------
    model : SharkModel
        Fully initialised SharkModel instance.  The pipeline reads SFH
        fields on first use (lazy via SharkModel's cache).
    z_obs : float
        Redshift at which to evaluate magnitudes.  Should correspond to a
        snapshot in the Shark run.
    bands : list of str, optional
        FSPS filter names.  Default ['v'] (Buser V-band).
        Call SPSEngine().available_bands() to list all options.
    imf_type : int
        FSPS IMF flag.  0=Salpeter, 1=Chabrier (default), 2=Kroupa.
        Must match the IMF assumed in your Shark run.
    add_dust : bool
        Apply Calzetti (2000) attenuation.  Default False (intrinsic).
    add_neb_emission : bool
        Include nebular line and continuum emission.  Default False.
    progress : bool
        Print progress every 1000 galaxies.  Default True.

    Notes
    -----
    The isochrone/spectral library (MIST, Padova, etc.) is fixed at
    python-fsps install time and cannot be chosen here — see
    ``SPSEngine.libraries`` after construction to check which one your
    installed python-fsps actually uses.
    """

    def __init__(
        self,
        model: SharkModel,
        z_obs: float,
        bands: Optional[list[str]] = None,
        imf_type: int = 1,
        add_dust: bool = False,
        add_neb_emission: bool = False,
        progress: bool = True,
    ):
        self.model    = model
        self.z_obs    = z_obs
        self.bands    = bands or ["v"]
        self.progress = progress

        self.engine = SPSEngine(
            imf_type=imf_type,
            add_dust=add_dust,
            add_neb_emission=add_neb_emission,
            bands=self.bands,
        )

        print(f"PhotometryPipeline initialized with bands={self.bands}, "
              f"IMF type={imf_type}, add_dust={add_dust}, "
              f"add_neb_emission={add_neb_emission}.")

        # SFH bin age array — derived from lbt_mean in the SFH file
        # Shape (n_sfh_bins,), ascending age of Universe in Gyr
        self._sfh_ages = sfh_ages_from_model(model, z_obs)
        self._tage     = model.age_at_z(z_obs)

        self._validate_sfh_ages()

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    @property
    def n_galaxies(self) -> int:
        """Total galaxies in the SharkModel catalogue at z_obs."""
        return len(self.model.get("mstars_disk", self.z_obs))

    def abs_mag_V(self, gal_indices: Optional[np.ndarray] = None) -> np.ndarray:
        """
        Absolute V-band magnitude for selected galaxies.

        Parameters
        ----------
        gal_indices : array-like, optional
            Indices into the SharkModel catalogue.  If None, runs all galaxies.

        Returns
        -------
        M_V : ndarray, shape (n_selected,)
            Absolute AB V-band magnitude.  NaN for galaxies with no stars.
        """
        if "v" not in self.bands:
            raise ValueError(
                "V band not in self.bands. Reinitialise with bands=['v', ...]"
            )
        mags  = self.abs_mags(gal_indices)
        v_idx = self.bands.index("v")
        return mags[:, v_idx]

    def abs_mags(self, gal_indices: Optional[np.ndarray] = None) -> np.ndarray:
        """
        Absolute magnitudes in all requested bands.

        Parameters
        ----------
        gal_indices : array-like, optional
            Galaxy indices.  If None, runs all galaxies.

        Returns
        -------
        mags : ndarray, shape (n_selected, n_bands)
            Absolute AB magnitudes.  Column order matches self.bands.
        """
        idx  = self._resolve_indices(gal_indices)
        data = galaxy_data(self.model, self.z_obs, idx)

        mags = self.engine.mag_galaxies(
            sfr_disk    = data["sfr_disk"],
            sfr_bulge   = data["sfr_bulge"],
            Z_disk      = data["Z_disk"],
            Z_bulge     = data["Z_bulge"],
            snap_ages   = self._sfh_ages,
            mstar_disk  = data["mstar_disk"],
            mstar_bulge = data["mstar_bulge"],
            tage        = self._tage,
            z_obs       = self.z_obs,
            progress    = self.progress,
        )
        return mags

    def mass_to_light(self, gal_indices: Optional[np.ndarray] = None) -> np.ndarray:
        """
        V-band stellar mass-to-light ratio (M/L)_V in solar units.

        Returns
        -------
        ML_V : ndarray, shape (n_selected,)
        """
        idx  = self._resolve_indices(gal_indices)
        M_V  = self.abs_mag_V(idx)

        data = galaxy_data(self.model, self.z_obs, idx)
        mstar_total = data["mstar_disk"] + data["mstar_bulge"]

        # Sun's absolute V-band magnitude (AB)
        M_V_sun = 4.83
        L_V = 10.0 ** (-0.4 * (M_V - M_V_sun))

        with np.errstate(divide="ignore", invalid="ignore"):
            ML_V = np.where(L_V > 0, mstar_total / L_V, np.nan)

        return ML_V

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _resolve_indices(self, gal_indices) -> np.ndarray:
        if gal_indices is None:
            return np.arange(self.n_galaxies)
        return np.atleast_1d(np.asarray(gal_indices, dtype=int))

    def _validate_sfh_ages(self) -> None:
        """Sanity-check that the SFH age array is well-formed for FSPS."""
        ages = self._sfh_ages
        if len(ages) == 0:
            raise ValueError("SFH age array is empty.")
        if not np.all(np.diff(ages) > 0):
            raise ValueError(
                "sfh_ages must be strictly increasing (oldest bin first). "
                "Check sfh_ages_from_model in io.py."
            )
        if ages[0] <= 0:
            raise ValueError(
                f"First SFH age is {ages[0]:.4f} Gyr — must be > 0 for FSPS."
            )
