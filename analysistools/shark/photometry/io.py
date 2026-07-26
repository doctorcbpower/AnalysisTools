"""
shark.photometry.io
===================
Connects SharkModel and the photometry pipeline.

SharkModel (shark.model) is the authoritative source for all HDF5 reading —
it already handles subvolume concatenation, caching, and field discovery for
both galaxies.hdf5 and star_formation_histories.hdf5 via common.read_sfh.

This module provides one function, ``galaxy_data``, which pulls the arrays
the SPSEngine needs from a SharkModel instance at a given redshift.  There
is no direct HDF5 access here.

SFH time axis
-------------
star_formation_histories.hdf5 contains ``lbt_mean`` — the mean lookback time
[Gyr] of each SFH bin.  For FSPS set_tabular_sfh we need the *age of the
Universe* at each bin, which is:

    snap_age = age_of_universe_at_z_obs - lbt_mean

This is computed in ``sfh_ages_from_model``, using ``SharkModel.age_at_z``
for the cosmology (no separate cosmology object needed).

Metallicity
-----------
The SFH file stores metal *mass* per bin (sfh_metals_disk / sfh_metals_bulge)
and SFR per bin.  Metallicity Z(t) is derived as mz_bin / (sfr_bin * delta_t),
i.e. metal mass formed per bin divided by stellar mass formed per bin.
SharkModel._sfh_metallicity performs this calculation; we call it via the
Z_disk_history / Z_bulge_history convenience methods.
"""

from __future__ import annotations

import numpy as np
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ..model import SharkModel

# Z bounds matching FSPS SSP grid (Padova / MIST)
Z_MIN = 1.0e-4
Z_MAX = 0.03


def sfh_ages_from_model(model: "SharkModel", redshift: float) -> np.ndarray:
    """
    Return the age of the Universe [Gyr] at each SFH bin midpoint.

    FSPS set_tabular_sfh requires ages in ascending order starting from
    the earliest bin.  Shark's lbt_mean is in *descending* order (most
    recent bin last has the smallest lookback time), so we reverse it.

    Parameters
    ----------
    model : SharkModel
    redshift : float
        The observation redshift — used to get the age of the Universe
        at z_obs from model cosmology.

    Returns
    -------
    ages : ndarray, shape (n_sfh_bins,)
        Age of the Universe [Gyr] at each SFH bin, strictly ascending.
    """
    meta   = model.get_sfh_meta(redshift)
    lbt    = meta["lbt_mean"]                      # Gyr, Shark ordering

    t_obs  = model.age_at_z(redshift)               # Gyr, age of Universe at z_obs

    # Convert lookback time → age of Universe
    # lbt_mean is stored high-lbt first (oldest bins first in Shark convention)
    ages = t_obs - lbt
    ages = np.clip(ages, 1e-3, None)

    # Ensure strictly ascending (oldest = smallest age first)
    if ages[0] > ages[-1]:
        ages = ages[::-1]

    return ages.astype(np.float64)


def galaxy_data(
    model: "SharkModel",
    redshift: float,
    gal_idx: np.ndarray,
) -> dict:
    """
    Extract SFH, metallicity, and stellar mass arrays for selected galaxies.

    Reads from the SharkModel cache (triggering HDF5 reads on first call).
    All arrays are indexed along axis 0 by gal_idx; the SFH bin axis is
    already in ascending-time order as required by FSPS.

    Parameters
    ----------
    model : SharkModel
    redshift : float
    gal_idx : ndarray of int
        Galaxy indices into the full SharkModel catalogue at this redshift.

    Returns
    -------
    dict with keys:
        sfr_disk        (n_gal, n_sfh_bins)  Msun/yr
        sfr_bulge       (n_gal, n_sfh_bins)  Msun/yr
        Z_disk          (n_gal, n_sfh_bins)  dimensionless, clipped to FSPS range
        Z_bulge         (n_gal, n_sfh_bins)
        mstar_disk      (n_gal,)  Msun surviving at z_obs
        mstar_bulge     (n_gal,)  Msun surviving at z_obs
    """
    idx = np.atleast_1d(gal_idx)

    # SFH arrays — shape (n_gal_total, n_sfh_bins) from SharkModel
    sfr_disk  = model.sfh_disk(redshift)[idx]
    sfr_bulge = model.sfh_bulge(redshift)[idx]
    Z_disk    = model.Z_disk_history(redshift)[idx]
    Z_bulge   = model.Z_bulge_history(redshift)[idx]

    # Ensure SFH bin axis is ascending in time (oldest bin first)
    # SharkModel._load_sfh_snapshot does not guarantee ordering; check here.
    sfr_disk, sfr_bulge, Z_disk, Z_bulge = _ensure_ascending_sfh(
        sfr_disk, sfr_bulge, Z_disk, Z_bulge
    )

    # Surviving stellar mass at z_obs — from galaxies.hdf5 via SharkModel
    mstar_disk  = model.get("mstars_disk",  redshift)[idx]
    mstar_bulge = model.get("mstars_bulge", redshift)[idx]

    # Sanitise: clip negatives, replace non-finite values
    out = {
        "sfr_disk":   _sanitise(sfr_disk),
        "sfr_bulge":  _sanitise(sfr_bulge),
        "Z_disk":     np.clip(_sanitise(Z_disk,  fill=0.02), Z_MIN, Z_MAX),
        "Z_bulge":    np.clip(_sanitise(Z_bulge, fill=0.02), Z_MIN, Z_MAX),
        "mstar_disk":  np.clip(mstar_disk,  0.0, None),
        "mstar_bulge": np.clip(mstar_bulge, 0.0, None),
    }
    return out


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _ensure_ascending_sfh(*arrays: np.ndarray):
    """
    Flip the SFH bin axis (axis=1) if SFR is higher in the first bin
    than the last for the median galaxy — a sign the bins run
    newest-to-oldest rather than oldest-to-newest.
    """
    sfr = arrays[0]
    if sfr.ndim < 2 or sfr.shape[1] < 2:
        return arrays
    if np.median(sfr[:, 0]) > np.median(sfr[:, -1]):
        return tuple(a[:, ::-1] for a in arrays)
    return arrays


def _sanitise(arr: np.ndarray, fill: float = 0.0) -> np.ndarray:
    """Replace non-finite values with *fill* and clip negatives to zero."""
    arr = np.where(np.isfinite(arr), arr, fill)
    return np.clip(arr, 0.0, None)
