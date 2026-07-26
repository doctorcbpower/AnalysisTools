"""
shark.analysis
--------------
Analysis: computes binned statistics from one or more SharkModel instances.

Design principles
~~~~~~~~~~~~~~~~~
* Analysis wraps a list of SharkModels and knows about binning and statistics.
* It never touches matplotlib — results are plain numpy arrays.
* Each ``compute_*`` method returns a dict of arrays keyed by model label,
  ready to be passed directly to HaloMFPlotter or inspected in a notebook.
* Bin configurations are plain dicts so they are trivially serialisable.
"""

from __future__ import annotations

from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np

from .model import SharkModel

# ---------------------------------------------------------------------------
# Binning helpers
# ---------------------------------------------------------------------------

def make_bins(low: float, upp: float, dbin: float) -> Tuple[np.ndarray, np.ndarray]:
    """Return ``(bin_edges, bin_centres)`` for a uniform log10-mass grid."""
    edges   = np.arange(low, upp, dbin)
    centres = edges + dbin / 2.0
    return edges, centres


def number_density(
    values:   np.ndarray,
    weights:  Optional[np.ndarray],
    edges:    np.ndarray,
    vol:      float,
    dbin:     float,
) -> np.ndarray:
    """Compute dn/dlog10(x) [Mpc^-3] from *values* in log10 space.

    Parameters
    ----------
    values : ndarray
        Log10-transformed quantity (e.g. log10 mass) for each galaxy.
    weights : ndarray or None
        Per-galaxy weights.  Pass None for unweighted counts.
    edges : ndarray
        Bin edges in log10 space.
    vol : float
        Comoving survey volume in Mpc^3.
    dbin : float
        Bin width (= spacing between edges).

    Returns
    -------
    phi : ndarray
        Log10 number density per bin.  NaN where the bin is empty.
    """
    H, _ = np.histogram(values, bins=np.append(edges, edges[-1] + dbin),
                        weights=weights)
    with np.errstate(divide="ignore", invalid="ignore"):
        phi = np.where(H > 0, np.log10(H / (vol * dbin)), np.nan)
    return phi


# ---------------------------------------------------------------------------
# Analysis
# ---------------------------------------------------------------------------

class Analysis:
    """Compute binned statistics from a collection of SharkModel instances.

    Parameters
    ----------
    models : list of SharkModel
        Models to analyse.  All are processed identically.

    Examples
    --------
    >>> from shark import SharkModel, Analysis
    >>> an = Analysis([cdm_model, wdm_model])
    >>> results = an.compute_halo_mf(redshifts=[0.0, 1.0])
    >>> results["CDM"]["z=0.0"]["phi"]   # log10 number density array
    """

    def __init__(self, models: List[SharkModel]):
        if not models:
            raise ValueError("models list must not be empty.")
        self.models = models

    # ------------------------------------------------------------------
    # Mass functions
    # ------------------------------------------------------------------

    def compute_halo_mf(
        self,
        redshifts:  Sequence[float],
        low:  float = 4.0,
        upp:  float = 15.0,
        dbin: float = 0.3,
        centrals_only: bool = True,
    ) -> Dict[str, Dict]:
        """Halo mass function  dn/dlog10(M_halo) [Mpc^-3].

        Parameters
        ----------
        redshifts : sequence of float
        low, upp, dbin : float
            Bin configuration in log10(M_halo / M_sun).
        centrals_only : bool
            If True (default) use only central galaxies (type <= 0)
            for the host halo MF.  Set False for the subhalo MF.

        Returns
        -------
        results : dict
            ``results[label][z_key]`` contains:

            ``phi``     — log10 number density array, shape (n_bins,)
            ``x``       — bin centres, shape (n_bins,)
            ``edges``   — bin edges
            ``plotz``   — bool, True if volume was valid
        """
        edges, centres = make_bins(low, upp, dbin)
        results = {}

        for model in self.models:
            model_results = {}
            for z in redshifts:
                z_key  = f"z={z:.2f}"
                meta   = model.get_meta(z)
                h0     = meta["h0"]
                vol    = meta["vol"]
                plotz  = vol > 0.0

                phi = np.full(len(edges), np.nan)
                if plotz:
                    typeg = model.galaxy_type(z)
                    ms    = model.mstars(z)
                    mh    = model.mhalo(z) if centrals_only else model.msubhalo(z)

                    sel = (ms > 1e5)
                    if centrals_only:
                        sel &= (typeg <= 0)
                    else:
                        sel &= (typeg <= 1)

                    log_mh = np.where(
                        sel & (mh > 0),
                        np.log10(mh) - np.log10(h0),
                        0.0,
                    )
                    phi = number_density(log_mh[sel & (mh > 0)],
                                         None, edges, vol, dbin)

                model_results[z_key] = dict(
                    phi=phi, x=centres, edges=edges, plotz=plotz, z=z
                )
            results[model.label] = model_results

        return results

    def compute_stellar_mf(
        self,
        redshifts:  Sequence[float],
        low:  float = 7.0,
        upp:  float = 12.5,
        dbin: float = 0.2,
    ) -> Dict[str, Dict]:
        """Stellar mass function  dn/dlog10(M_star) [Mpc^-3].

        Parameters
        ----------
        redshifts : sequence of float
        low, upp, dbin : float
            Bin configuration in log10(M_star / M_sun).
        """
        edges, centres = make_bins(low, upp, dbin)
        results = {}

        for model in self.models:
            model_results = {}
            for z in redshifts:
                z_key = f"z={z:.2f}"
                meta  = model.get_meta(z)
                h0, vol = meta["h0"], meta["vol"]
                plotz   = vol > 0.0

                phi = np.full(len(edges), np.nan)
                if plotz:
                    ms  = model.mstars(z)
                    sel = ms > 1e5
                    log_ms = np.log10(ms[sel]) - np.log10(h0)
                    phi = number_density(log_ms, None, edges, vol, dbin)

                model_results[z_key] = dict(
                    phi=phi, x=centres, edges=edges, plotz=plotz, z=z
                )
            results[model.label] = model_results

        return results

    # ------------------------------------------------------------------
    # Scaling relations  (return scatter arrays, not binned densities)
    # ------------------------------------------------------------------

    def compute_sfr_main_sequence(
        self,
        redshifts:      Sequence[float],
        ms_min:         float = 1e7,
        ssfr_min:       float = 1e-12,
        centrals_only:  bool  = False,
    ) -> Dict[str, Dict]:
        """Star formation rate vs stellar mass (main sequence).

        Returns per-galaxy log10(M_star) and log10(SFR) arrays — the caller
        controls how to bin or plot them.

        Parameters
        ----------
        redshifts : sequence of float
        ms_min : float
            Minimum stellar mass cut  [M_sun / h].
        ssfr_min : float
            Minimum sSFR to be considered star-forming  [yr^-1].
        centrals_only : bool
            Restrict to central galaxies only.
        """
        results = {}

        for model in self.models:
            model_results = {}
            for z in redshifts:
                z_key = f"z={z:.2f}"
                meta  = model.get_meta(z)
                h0    = meta["h0"]

                ms    = model.mstars(z)
                sfr   = model.sfr(z)
                ssfr  = model.ssfr(z)
                typeg = model.galaxy_type(z)

                sel = (ms > ms_min) & (ssfr > ssfr_min)
                if centrals_only:
                    sel &= (typeg <= 0)

                log_ms  = np.log10(ms[sel])  - np.log10(h0)
                log_sfr = np.log10(sfr[sel]) - np.log10(h0)

                model_results[z_key] = dict(
                    log_mstar=log_ms, log_sfr=log_sfr, z=z
                )
            results[model.label] = model_results

        return results

    def compute_size_mass(
        self,
        redshifts:      Sequence[float],
        ms_min:         float = 1e7,
        centrals_only:  bool  = False,
    ) -> Dict[str, Dict]:
        """Half-stellar-mass radius vs stellar mass.

        Returns per-galaxy log10(M_star) and log10(r_star) arrays.

        Parameters
        ----------
        redshifts : sequence of float
        ms_min : float
            Minimum stellar mass cut  [M_sun / h].
        centrals_only : bool
            Restrict to central galaxies only.
        """
        results = {}

        for model in self.models:
            model_results = {}
            for z in redshifts:
                z_key = f"z={z:.2f}"
                meta  = model.get_meta(z)
                h0    = meta["h0"]

                ms    = model.mstars(z)
                r     = model.rstar(z)
                typeg = model.galaxy_type(z)

                sel = (ms > ms_min) & (r > 0)
                if centrals_only:
                    sel &= (typeg <= 0)

                log_ms = np.log10(ms[sel]) - np.log10(h0)
                log_r  = np.log10(r[sel])  - np.log10(h0)

                model_results[z_key] = dict(
                    log_mstar=log_ms, log_rstar=log_r, z=z
                )
            results[model.label] = model_results

        return results

    def compute_bh_bulge(
        self,
        redshifts:      Sequence[float],
        mbulge_min:     float = 1e7,
        mbh_min:        float = 1e4,
        centrals_only:  bool  = True,
    ) -> Dict[str, Dict]:
        """Black hole mass vs bulge stellar mass.

        Returns per-galaxy log10(M_bulge) and log10(M_BH) arrays.

        Parameters
        ----------
        redshifts : sequence of float
        mbulge_min, mbh_min : float
            Minimum mass cuts  [M_sun / h].
        centrals_only : bool
            Restrict to central galaxies (default True — satellites rarely
            host significantly grown BHs).
        """
        results = {}

        for model in self.models:
            model_results = {}
            for z in redshifts:
                z_key = f"z={z:.2f}"
                meta  = model.get_meta(z)
                h0    = meta["h0"]

                mbulge = model.mbulge(z)
                mbh    = model.mbh(z)
                typeg  = model.galaxy_type(z)

                sel = (mbulge > mbulge_min) & (mbh > mbh_min)
                if centrals_only:
                    sel &= (typeg <= 0)

                log_mbulge = np.log10(mbulge[sel]) - np.log10(h0)
                log_mbh    = np.log10(mbh[sel])    - np.log10(h0)

                model_results[z_key] = dict(
                    log_mbulge=log_mbulge, log_mbh=log_mbh, z=z
                )
            results[model.label] = model_results

        return results

    # ------------------------------------------------------------------
    # Generic / custom statistic
    # ------------------------------------------------------------------

    def compute_custom(
        self,
        redshifts:  Sequence[float],
        x_func:     Callable[[SharkModel, float], np.ndarray],
        y_func:     Callable[[SharkModel, float], np.ndarray],
        sel_func:   Optional[Callable[[SharkModel, float], np.ndarray]] = None,
        x_log:      bool = True,
        y_log:      bool = True,
    ) -> Dict[str, Dict]:
        """Compute a fully custom x vs y relation.

        Parameters
        ----------
        redshifts : sequence of float
        x_func : callable(model, z) -> ndarray
            Returns the x-axis quantity (before optional log10).
        y_func : callable(model, z) -> ndarray
            Returns the y-axis quantity (before optional log10).
        sel_func : callable(model, z) -> bool ndarray, optional
            Returns a selection mask.  Defaults to all galaxies.
        x_log, y_log : bool
            If True (default), take log10 of the respective axis.

        Examples
        --------
        Compute cold gas fraction vs stellar mass::

            results = analysis.compute_custom(
                redshifts=[0.0],
                x_func=lambda m, z: m.mstars(z),
                y_func=lambda m, z: m.mgas(z) / (m.mstars(z) + m.mgas(z)),
                sel_func=lambda m, z: m.mstars(z) > 1e8,
                x_log=True,
                y_log=False,
            )
        """
        results = {}

        for model in self.models:
            model_results = {}
            for z in redshifts:
                z_key = f"z={z:.2f}"
                meta  = model.get_meta(z)
                h0    = meta["h0"]

                x_raw = x_func(model, z)
                y_raw = y_func(model, z)

                if sel_func is not None:
                    sel = sel_func(model, z)
                    x_raw = x_raw[sel]
                    y_raw = y_raw[sel]

                with np.errstate(divide="ignore", invalid="ignore"):
                    x_out = np.log10(x_raw) - np.log10(h0) if x_log else x_raw
                    y_out = np.log10(y_raw) - np.log10(h0) if y_log else y_raw

                model_results[z_key] = dict(x=x_out, y=y_out, z=z)
            results[model.label] = model_results

        return results
