"""
shark.photometry.sps
====================
python-fsps wrapper for convolving Shark SFHs with SSP templates.

Design principles
-----------------
- One StellarPopulation object is instantiated per SPSEngine and reused
  across all galaxies via set_tabular_sfh, avoiding the expensive C
  initialisation overhead.
- Disk and bulge are treated as independent SSPs then combined in
  luminosity space.
- Normalisation: FSPS set_tabular_sfh works in units of Msun/yr and
  normalises internally to 1 Msun *formed*. We recover the correct
  total luminosity by scaling by (mstar_surviving / sp.stellar_mass),
  where sp.stellar_mass is FSPS's surviving fraction for the same SFH.
- Metallicity is passed as Z(t) per snapshot, linearly interpolated by
  FSPS (zcontinuous=1).

Supported IMFs (imf_type)
-------------------------
  0 : Salpeter (1955)
  1 : Chabrier (2003)  ← Shark default
  2 : Kroupa (2001)
"""

from __future__ import annotations

import numpy as np
from typing import Sequence
import warnings


# FSPS SSP metallicity grid bounds (Padova isochrones)
Z_MIN = 1.0e-4
Z_MAX = 0.03

# Minimum SFR to pass to FSPS — avoids numerical noise in empty snapshots
SFR_FLOOR = 1.0e-20  # Msun/yr


class SPSEngine:
    """
    Reusable FSPS engine for tabulated SFH convolution.

    Parameters
    ----------
    imf_type : int
        FSPS IMF flag. Default 1 (Chabrier).
    add_dust : bool
        Apply Calzetti (2000) attenuation. Default False (intrinsic).
    add_neb_emission : bool
        Include nebular line and continuum emission. Useful at high-z.
    bands : list of str
        FSPS filter names to compute. Default ['v'] (Buser V).
        Run fsps.list_filters() to see all available.

    Notes
    -----
    The isochrone/spectral library (e.g. MIST, Padova) is fixed at
    python-fsps *install* time via compiler flags — it cannot be chosen
    here at runtime. Check ``self.libraries`` (= ``sp.libraries``) to see
    which library your installed python-fsps is actually using; if you
    need a different one, reinstall python-fsps with the relevant
    ``FFLAGS`` (see python-fsps installation docs).
    """

    def __init__(
        self,
        imf_type: int = 1,
        add_dust: bool = False,
        add_neb_emission: bool = False,
        bands: list[str] = None,
    ):
        self.bands = bands or ["v"]
        self._imf_type = imf_type

        try:
            import fsps
        except ImportError as e:
            raise ImportError(
                "python-fsps is required to use SPSEngine. "
                "Install with: pip install fsps"
            ) from e
        self._fsps = fsps

        self._sp = fsps.StellarPopulation(
            compute_vega_mags=False,      # AB magnitudes throughout
            zcontinuous=3,                # interpolate in log Z
            sfh=3,                        # tabulated SFH input mode
            imf_type=imf_type,
            add_dust_emission=add_dust,
            add_neb_emission=add_neb_emission,
            dust_type=2 if add_dust else 0,   # Calzetti if dust on
            vactoair_flag=False,
            smooth_velocity=True,
        )

        # Whatever isochrone/spectral library this install was compiled
        # with — informational only, not settable here.
        self.libraries = self._sp.libraries

    # ------------------------------------------------------------------
    # Single-galaxy interface
    # ------------------------------------------------------------------

    def mag_one_component(
        self,
        sfr: np.ndarray,
        Z_hist: np.ndarray,
        snap_ages: np.ndarray,
        mstar_surviving: float,
        tage: float,
        z_obs: float = 0.0,
    ) -> np.ndarray:
        """
        Compute absolute AB magnitudes for one galaxy component (disk or bulge).

        Parameters
        ----------
        sfr : ndarray, shape (n_snaps,)
            SFR at each snapshot in Msun/yr.
        Z_hist : ndarray, shape (n_snaps,)
            Metallicity Z (not log Z) at each snapshot.
        snap_ages : ndarray, shape (n_snaps,)
            Age of the Universe [Gyr] at each snapshot, strictly ascending.
        mstar_surviving : float
            Surviving stellar mass [Msun] at z_obs from Shark.
        tage : float
            Age of the Universe [Gyr] at z_obs.
        z_obs : float
            Observed redshift (for K-correction in get_mags). Default 0.

        Returns
        -------
        mags : ndarray, shape (n_bands,)
            Absolute AB magnitudes in self.bands.
            Returns array of NaN if component has negligible mass.
        """
        # Guard: skip components with no stars
        if mstar_surviving <= 0 or np.sum(sfr) <= 0:
            return np.full(len(self.bands), np.nan)

        sfr_clean  = np.clip(sfr,   SFR_FLOOR, None).astype(np.float64)
        Z_clean    = np.clip(Z_hist, Z_MIN, Z_MAX).astype(np.float64)
        ages_clean = np.clip(snap_ages, 1e-3, None).astype(np.float64)

        self._sp.set_tabular_sfh(age=ages_clean, sfr=sfr_clean, Z=Z_clean)

        # FSPS normalises to 1 Msun formed; get_mags returns absolute mag
        mags_1msun = self._sp.get_mags(
            tage=tage,
            redshift=z_obs,
            bands=self.bands,
        )

        # FSPS surviving fraction for this SFH
        fsps_surviving = self._sp.stellar_mass  # Msun surviving per Msun formed

        if fsps_surviving <= 0:
            return np.full(len(self.bands), np.nan)

        # Scale: mstar_surviving [Msun] / fsps_surviving = Msun formed
        log_scale = np.log10(mstar_surviving / fsps_surviving)
        mags = mags_1msun - 2.5 * log_scale

        return mags

    def mag_galaxy(
        self,
        sfr_disk: np.ndarray,
        sfr_bulge: np.ndarray,
        Z_disk: np.ndarray,
        Z_bulge: np.ndarray,
        snap_ages: np.ndarray,
        mstar_disk: float,
        mstar_bulge: float,
        tage: float,
        z_obs: float = 0.0,
    ) -> np.ndarray:
        """
        Compute absolute magnitudes for a full galaxy (disk + bulge combined).

        Disk and bulge SEDs are computed independently and combined in
        linear flux space before converting back to magnitudes.

        Parameters
        ----------
        sfr_disk, sfr_bulge : ndarray (n_snaps,)
        Z_disk, Z_bulge : ndarray (n_snaps,)
        snap_ages : ndarray (n_snaps,)
        mstar_disk, mstar_bulge : float  [Msun]
        tage : float  [Gyr]
        z_obs : float

        Returns
        -------
        mags : ndarray (n_bands,)  absolute AB magnitudes
        """
        mags_disk  = self.mag_one_component(
            sfr_disk,  Z_disk,  snap_ages, mstar_disk,  tage, z_obs
        )
        mags_bulge = self.mag_one_component(
            sfr_bulge, Z_bulge, snap_ages, mstar_bulge, tage, z_obs
        )

        return _combine_mags(mags_disk, mags_bulge)

    # ------------------------------------------------------------------
    # Batch interface — loop over galaxy array, reusing sp object
    # ------------------------------------------------------------------

    def mag_galaxies(
        self,
        sfr_disk: np.ndarray,
        sfr_bulge: np.ndarray,
        Z_disk: np.ndarray,
        Z_bulge: np.ndarray,
        snap_ages: np.ndarray,
        mstar_disk: np.ndarray,
        mstar_bulge: np.ndarray,
        tage: float,
        z_obs: float = 0.0,
        progress: bool = True,
    ) -> np.ndarray:
        """
        Compute absolute magnitudes for an array of galaxies.

        Parameters
        ----------
        sfr_disk, sfr_bulge : ndarray (n_gal, n_snaps)
        Z_disk, Z_bulge     : ndarray (n_gal, n_snaps)
        snap_ages           : ndarray (n_snaps,)
        mstar_disk, mstar_bulge : ndarray (n_gal,)
        tage                : float [Gyr]
        z_obs               : float
        progress            : bool  print progress every 1000 galaxies

        Returns
        -------
        mags : ndarray (n_gal, n_bands)  absolute AB magnitudes
        """
        n_gal   = sfr_disk.shape[0]
        n_bands = len(self.bands)
        mags    = np.full((n_gal, n_bands), np.nan)

        for i in range(n_gal):
            if progress and i % 1000 == 0:
                print(f"  SPSEngine: galaxy {i}/{n_gal}", flush=True)
            try:
                mags[i] = self.mag_galaxy(
                    sfr_disk[i],  sfr_bulge[i],
                    Z_disk[i],    Z_bulge[i],
                    snap_ages,
                    float(mstar_disk[i]), float(mstar_bulge[i]),
                    tage, z_obs,
                )
            except Exception as exc:
                warnings.warn(f"Galaxy {i} failed: {exc}")
                # Leave as NaN and continue

        return mags

    # ------------------------------------------------------------------
    # Utility
    # ------------------------------------------------------------------

    @property
    def imf_type(self) -> int:
        return self._imf_type

    def available_bands(self) -> list[str]:
        """List all FSPS filter names available for get_mags."""
        return self._fsps.list_filters()


# ------------------------------------------------------------------
# Module-level helper
# ------------------------------------------------------------------

def _combine_mags(mags_a: np.ndarray, mags_b: np.ndarray) -> np.ndarray:
    """
    Combine two magnitude arrays by summing luminosities.

    NaN components are treated as zero luminosity (absent component).
    If both are NaN, returns NaN.
    """
    result = np.full_like(mags_a, np.nan)
    for i in range(len(mags_a)):
        a_nan = np.isnan(mags_a[i])
        b_nan = np.isnan(mags_b[i])
        if a_nan and b_nan:
            result[i] = np.nan
        elif a_nan:
            result[i] = mags_b[i]
        elif b_nan:
            result[i] = mags_a[i]
        else:
            L_a = 10.0 ** (-0.4 * mags_a[i])
            L_b = 10.0 ** (-0.4 * mags_b[i])
            result[i] = -2.5 * np.log10(L_a + L_b)
    return result
