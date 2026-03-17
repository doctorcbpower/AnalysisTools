"""
shark.plots
-----------
Plotter: turns Analysis result dicts into publication-quality figures.

Design principles
~~~~~~~~~~~~~~~~~
* Plotter knows nothing about HDF5 or SHARK internals.
* Every ``plot_*`` method accepts the dict returned by the corresponding
  ``Analysis.compute_*`` method, plus optional style overrides.
* Returns the list of saved file paths so callers can chain or inspect.
* A ``plot_custom`` method handles any x/y result dict from
  ``Analysis.compute_custom``.
"""

from __future__ import annotations

import os
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from . import common

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------

DEFAULT_COLOURS = [
    "#e41a1c",  # red
    "#377eb8",  # blue
    "#4daf4a",  # green
    "#ff7f00",  # orange
    "#984ea3",  # purple
    "#a65628",  # brown
    "#f781bf",  # pink
    "#999999",  # grey
]

# Axis configuration bundles for each built-in plot type.
_AXIS_DEFAULTS: Dict[str, dict] = {
    "halo_mf": dict(
        xmin=6.0, xmax=15.0, ymin=-10.0, ymax=-1.0,
        xtit=r"$\log_{10}(M_{\rm halo}\,/\,M_\odot)$",
        ytit=r"$\log_{10}(dn/d\log_{10}M\,[\rm Mpc^{-3}\,dex^{-1}])$",
        locators=(0.1, 1, 0.1),
    ),
    "stellar_mf": dict(
        xmin=7.0, xmax=12.5, ymin=-7.0, ymax=-1.0,
        xtit=r"$\log_{10}(M_\star\,/\,M_\odot)$",
        ytit=r"$\log_{10}(\phi\,[\rm Mpc^{-3}\,dex^{-1}])$",
        locators=(0.1, 1, 0.1),
    ),
    "sfr_ms": dict(
        xmin=7.0, xmax=12.0, ymin=-3.0, ymax=3.0,
        xtit=r"$\log_{10}(M_\star\,/\,M_\odot)$",
        ytit=r"$\log_{10}(\rm SFR\,[M_\odot\,yr^{-1}])$",
        locators=(0.1, 1, 0.1),
    ),
    "size_mass": dict(
        xmin=7.0, xmax=12.0, ymin=-1.5, ymax=1.5,
        xtit=r"$\log_{10}(M_\star\,/\,M_\odot)$",
        ytit=r"$\log_{10}(r_{1/2}\,[\rm kpc])$",
        locators=(0.1, 1, 0.1),
    ),
    "bh_bulge": dict(
        xmin=7.0, xmax=12.5, ymin=4.0, ymax=11.0,
        xtit=r"$\log_{10}(M_{\rm bulge}\,/\,M_\odot)$",
        ytit=r"$\log_{10}(M_{\rm BH}\,/\,M_\odot)$",
        locators=(0.1, 1, 0.1),
    ),
}


# ---------------------------------------------------------------------------
# Plotter
# ---------------------------------------------------------------------------

class Plotter:
    """Render Analysis results to figures and save them to *outdir*.

    Parameters
    ----------
    outdir : str
        Directory for saved figures.  Created if it does not exist.
    colours : list of str, optional
        Colour cycle for models.  Defaults to ``DEFAULT_COLOURS``.
    figsize : tuple of float, optional
        Default figure size in inches.  Default: ``(7, 7)``.
    scatter_alpha : float, optional
        Alpha for scatter plot points.  Default: ``0.3``.
    scatter_size : float, optional
        Point size for scatter plots.  Default: ``0.5``.

    Examples
    --------
    >>> from shark import SharkModel, Analysis
    >>> from shark.plots import Plotter
    >>> an = Analysis([cdm, wdm])
    >>> p  = Plotter("./plots")
    >>> p.plot_halo_mf(an.compute_halo_mf([0.0, 1.0]))
    >>> p.plot_stellar_mf(an.compute_stellar_mf([0.0, 1.0]))
    """

    def __init__(
        self,
        outdir:        str,
        colours:       Optional[List[str]] = None,
        figsize:       Tuple[float, float] = (7, 7),
        scatter_alpha: float = 0.3,
        scatter_size:  float = 0.5,
    ):
        self.outdir        = outdir
        self.colours       = colours or list(DEFAULT_COLOURS)
        self.figsize       = figsize
        self.scatter_alpha = scatter_alpha
        self.scatter_size  = scatter_size
        os.makedirs(self.outdir, exist_ok=True)

    # ------------------------------------------------------------------
    # Built-in plot types
    # ------------------------------------------------------------------

    def plot_halo_mf(
        self,
        results:  Dict[str, Dict],
        redshifts: Optional[Sequence[float]] = None,
        axis_cfg: Optional[dict] = None,
        fname_prefix: str = "halomf",
    ) -> List[str]:
        """Plot halo mass functions from ``Analysis.compute_halo_mf`` results.

        Parameters
        ----------
        results : dict
            As returned by ``Analysis.compute_halo_mf``.
        redshifts : sequence of float, optional
            Subset of redshifts to plot.  Defaults to all available.
        axis_cfg : dict, optional
            Override axis defaults (xmin, xmax, ymin, ymax, xtit, ytit, locators).
        fname_prefix : str
            Prefix for output filenames.

        Returns
        -------
        saved_paths : list of str
        """
        cfg = {**_AXIS_DEFAULTS["halo_mf"], **(axis_cfg or {})}
        return self._plot_phi(results, redshifts, cfg, fname_prefix,
                              x_key="x", y_key="phi")

    def plot_stellar_mf(
        self,
        results:   Dict[str, Dict],
        redshifts: Optional[Sequence[float]] = None,
        axis_cfg:  Optional[dict] = None,
        fname_prefix: str = "stellarmf",
    ) -> List[str]:
        """Plot stellar mass functions from ``Analysis.compute_stellar_mf`` results."""
        cfg = {**_AXIS_DEFAULTS["stellar_mf"], **(axis_cfg or {})}
        return self._plot_phi(results, redshifts, cfg, fname_prefix,
                              x_key="x", y_key="phi")

    def plot_sfr_main_sequence(
        self,
        results:   Dict[str, Dict],
        redshifts: Optional[Sequence[float]] = None,
        axis_cfg:  Optional[dict] = None,
        fname_prefix: str = "sfr_ms",
        use_median: bool = True,
        bin_width:  float = 0.2,
    ) -> List[str]:
        """Plot the SFR main sequence from ``Analysis.compute_sfr_main_sequence``.

        Parameters
        ----------
        use_median : bool
            If True (default), plot the running median ± 16th/84th percentiles.
            If False, plot individual galaxies as a scatter.
        bin_width : float
            Bin width in log10(M_star) for the running median.
        """
        cfg = {**_AXIS_DEFAULTS["sfr_ms"], **(axis_cfg or {})}
        return self._plot_scatter_or_median(
            results, redshifts, cfg, fname_prefix,
            x_key="log_mstar", y_key="log_sfr",
            use_median=use_median, bin_width=bin_width,
        )

    def plot_size_mass(
        self,
        results:   Dict[str, Dict],
        redshifts: Optional[Sequence[float]] = None,
        axis_cfg:  Optional[dict] = None,
        fname_prefix: str = "sizemass",
        use_median: bool = True,
        bin_width:  float = 0.2,
    ) -> List[str]:
        """Plot the size-mass relation from ``Analysis.compute_size_mass``."""
        cfg = {**_AXIS_DEFAULTS["size_mass"], **(axis_cfg or {})}
        return self._plot_scatter_or_median(
            results, redshifts, cfg, fname_prefix,
            x_key="log_mstar", y_key="log_rstar",
            use_median=use_median, bin_width=bin_width,
        )

    def plot_bh_bulge(
        self,
        results:   Dict[str, Dict],
        redshifts: Optional[Sequence[float]] = None,
        axis_cfg:  Optional[dict] = None,
        fname_prefix: str = "bh_bulge",
        use_median: bool = True,
        bin_width:  float = 0.2,
    ) -> List[str]:
        """Plot BH-bulge mass relation from ``Analysis.compute_bh_bulge``."""
        cfg = {**_AXIS_DEFAULTS["bh_bulge"], **(axis_cfg or {})}
        return self._plot_scatter_or_median(
            results, redshifts, cfg, fname_prefix,
            x_key="log_mbulge", y_key="log_mbh",
            use_median=use_median, bin_width=bin_width,
        )

    def plot_custom(
        self,
        results:      Dict[str, Dict],
        redshifts:    Optional[Sequence[float]] = None,
        axis_cfg:     Optional[dict] = None,
        fname_prefix: str = "custom",
        use_median:   bool = False,
        bin_width:    float = 0.2,
    ) -> List[str]:
        """Plot any custom x/y result from ``Analysis.compute_custom``.

        Parameters
        ----------
        results : dict
            As returned by ``Analysis.compute_custom``.
            Expected keys per z-slice: ``x``, ``y``.
        axis_cfg : dict
            Must supply at least ``xtit`` and ``ytit``.
            Axis limits default to data range if not given.
        use_median : bool
            Plot running median instead of scatter.
        """
        # For custom plots, axis limits default to data range
        all_x = np.concatenate([
            z_data["x"]
            for m_data in results.values()
            for z_data in m_data.values()
            if len(z_data.get("x", [])) > 0
        ])
        all_y = np.concatenate([
            z_data["y"]
            for m_data in results.values()
            for z_data in m_data.values()
            if len(z_data.get("y", [])) > 0
        ])

        defaults = dict(
            xmin=float(np.nanmin(all_x)) if len(all_x) else 0.0,
            xmax=float(np.nanmax(all_x)) if len(all_x) else 1.0,
            ymin=float(np.nanmin(all_y)) if len(all_y) else 0.0,
            ymax=float(np.nanmax(all_y)) if len(all_y) else 1.0,
            xtit="x", ytit="y",
            locators=(0.1, 1, 0.1),
        )
        cfg = {**defaults, **(axis_cfg or {})}
        return self._plot_scatter_or_median(
            results, redshifts, cfg, fname_prefix,
            x_key="x", y_key="y",
            use_median=use_median, bin_width=bin_width,
        )

    # ------------------------------------------------------------------
    # Private rendering engine
    # ------------------------------------------------------------------

    def _plot_phi(
        self,
        results:      Dict[str, Dict],
        redshifts:    Optional[Sequence[float]],
        cfg:          dict,
        fname_prefix: str,
        x_key:        str,
        y_key:        str,
    ) -> List[str]:
        """Render one number-density (phi) figure per redshift."""
        plt    = common.load_matplotlib()
        z_keys = self._resolve_z_keys(results, redshifts)
        saved  = []

        for z_key in z_keys:
            fig, ax = plt.subplots(figsize=self.figsize)
            common.prepare_ax(
                ax,
                cfg["xmin"], cfg["xmax"],
                cfg["ymin"], cfg["ymax"],
                cfg["xtit"], cfg["ytit"],
                locators=cfg["locators"],
            )

            z_val = self._z_val_from_key(results, z_key)
            ax.text(
                cfg["xmax"] - 0.02 * (cfg["xmax"] - cfg["xmin"]),
                cfg["ymax"] - 0.05 * (cfg["ymax"] - cfg["ymin"]),
                f"z = {z_val:.1f}", ha="right", va="top",
            )

            any_plotted = False
            for i, (label, m_data) in enumerate(results.items()):
                if z_key not in m_data:
                    continue
                z_data = m_data[z_key]
                if not z_data.get("plotz", True):
                    continue
                x   = z_data[x_key]
                y   = z_data[y_key]
                col = self.colours[i % len(self.colours)]
                valid = np.isfinite(y)
                if valid.any():
                    ax.plot(x[valid], y[valid],
                            color=col, linewidth=1.8, label=label)
                    any_plotted = True

            if not any_plotted:
                plt.close(fig)
                continue

            self._add_legend(ax)
            path = self._save(plt, fig, f"{fname_prefix}_{z_key}.pdf")
            saved.append(path)

        return saved

    def _plot_scatter_or_median(
        self,
        results:      Dict[str, Dict],
        redshifts:    Optional[Sequence[float]],
        cfg:          dict,
        fname_prefix: str,
        x_key:        str,
        y_key:        str,
        use_median:   bool,
        bin_width:    float,
    ) -> List[str]:
        """Render one scatter / running-median figure per redshift."""
        plt    = common.load_matplotlib()
        z_keys = self._resolve_z_keys(results, redshifts)
        saved  = []

        for z_key in z_keys:
            fig, ax = plt.subplots(figsize=self.figsize)
            common.prepare_ax(
                ax,
                cfg["xmin"], cfg["xmax"],
                cfg["ymin"], cfg["ymax"],
                cfg["xtit"], cfg["ytit"],
                locators=cfg["locators"],
            )

            z_val = self._z_val_from_key(results, z_key)
            ax.text(
                cfg["xmax"] - 0.02 * (cfg["xmax"] - cfg["xmin"]),
                cfg["ymax"] - 0.05 * (cfg["ymax"] - cfg["ymin"]),
                f"z = {z_val:.1f}", ha="right", va="top",
            )

            any_plotted = False
            for i, (label, m_data) in enumerate(results.items()):
                if z_key not in m_data:
                    continue
                z_data = m_data[z_key]
                x = z_data[x_key]
                y = z_data[y_key]
                col = self.colours[i % len(self.colours)]

                if len(x) == 0:
                    continue

                if use_median:
                    xm, ym, ylo, yhi = self._running_median(
                        x, y, cfg["xmin"], cfg["xmax"], bin_width
                    )
                    valid = np.isfinite(ym)
                    if valid.any():
                        ax.plot(xm[valid], ym[valid],
                                color=col, linewidth=1.8, label=label)
                        ax.fill_between(
                            xm[valid], ylo[valid], yhi[valid],
                            color=col, alpha=0.2,
                        )
                        any_plotted = True
                else:
                    finite = np.isfinite(x) & np.isfinite(y)
                    if finite.any():
                        ax.scatter(
                            x[finite], y[finite],
                            color=col, s=self.scatter_size,
                            alpha=self.scatter_alpha, label=label,
                            rasterized=True,
                        )
                        any_plotted = True

            if not any_plotted:
                plt.close(fig)
                continue

            self._add_legend(ax)
            path = self._save(plt, fig, f"{fname_prefix}_{z_key}.pdf")
            saved.append(path)

        return saved

    # ------------------------------------------------------------------
    # Utilities
    # ------------------------------------------------------------------

    @staticmethod
    def _running_median(
        x: np.ndarray,
        y: np.ndarray,
        xmin: float,
        xmax: float,
        dbin: float,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Compute running median and 16th/84th percentile envelope."""
        edges   = np.arange(xmin, xmax + dbin, dbin)
        centres = edges[:-1] + dbin / 2.0
        med = np.full(len(centres), np.nan)
        lo  = np.full(len(centres), np.nan)
        hi  = np.full(len(centres), np.nan)

        for j, (e0, e1) in enumerate(zip(edges[:-1], edges[1:])):
            sel = (x >= e0) & (x < e1) & np.isfinite(y)
            if sel.sum() > 5:
                med[j] = np.median(y[sel])
                lo[j]  = np.percentile(y[sel], 16)
                hi[j]  = np.percentile(y[sel], 84)

        return centres, med, lo, hi

    @staticmethod
    def _resolve_z_keys(
        results:   Dict[str, Dict],
        redshifts: Optional[Sequence[float]],
    ) -> List[str]:
        """Return the list of z-key strings to iterate over."""
        # Gather all z_keys present in any model
        all_keys = []
        for m_data in results.values():
            for k in m_data:
                if k not in all_keys:
                    all_keys.append(k)
        if redshifts is None:
            return all_keys
        wanted = {f"z={z:.2f}" for z in redshifts}
        return [k for k in all_keys if k in wanted]

    @staticmethod
    def _z_val_from_key(results: Dict[str, Dict], z_key: str) -> float:
        """Extract the stored z float from the first model that has this z_key."""
        for m_data in results.values():
            if z_key in m_data:
                return m_data[z_key].get("z", float(z_key.split("=")[1]))
        return 0.0

    def _add_legend(self, ax) -> None:
        handles, labels = ax.get_legend_handles_labels()
        if handles:
            ax.legend(handles, labels, loc="lower left",
                      fontsize="small", framealpha=0.7)

    def _save(self, plt, fig, fname: str) -> str:
        common.savefig(self.outdir, fig, fname)
        plt.close(fig)
        path = os.path.join(self.outdir, fname)
        print(f"  [plotter] saved {path}")
        return path
