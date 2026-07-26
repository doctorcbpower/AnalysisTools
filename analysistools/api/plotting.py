#!/usr/bin/env python3
"""
analysistools.api.plotting
--------------------------
Dataset-aware plotting: thin wrappers that accept Dataset objects (or
selection views) directly, extract the arrays, and delegate the numerics
to GriddingTools / ProfileTools. Because every argument is a Dataset,
combining sources (snapshot + halos + galaxies on one axes) is just
several calls against the same ax.

    from analysistools.api import plotting as plt2

    ax = plt2.density_map(snap.dm, centre=c, size=5.0)
    plt2.overlay_points(ax, halos.select(mass=(1e3, None)), marker="o")
    plt2.overlay_points(ax, gals, colour_by="sfr")

    prof = plt2.profile(snap.dm, centre=c, kind="density",
                        rmin=0.05, rmax=1.0)
    plt2.plot_profile(prof, label=snap.label)

    plt2.mass_function([halos_cdm, halos_wdm])
"""
from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple, Union

import numpy as np

from ..gridding_tools import GriddingTools
from ..profile_tools import ProfileTools
from .dataset import Dataset

_PROJECTIONS = {"xy": (0, 1), "xz": (0, 2), "yz": (1, 2),
                "yx": (1, 0), "zx": (2, 0), "zy": (2, 1)}


def _axes_of(projection: str) -> Tuple[int, int]:
    try:
        return _PROJECTIONS[projection.lower()]
    except KeyError:
        raise ValueError(f"Unknown projection '{projection}' "
                         f"(use one of {sorted(_PROJECTIONS)}).")


def _get_ax(ax=None, figsize=(7, 6)):
    import matplotlib.pyplot as plt
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    return ax


# ---------------------------------------------------------------------------
# Gridded density maps
# ---------------------------------------------------------------------------

def density_map(ds: Dataset, *, projection: str = "xy",
                grid_size: Union[int, Tuple[int, int]] = 512,
                method: str = "CIC", weight: str = "mass",
                centre: Optional[Sequence[float]] = None,
                size: Optional[float] = None,
                species: Optional[str] = None,
                log: bool = True, smooth: Optional[float] = None,
                cmap: str = "inferno", ax=None, colorbar: bool = True,
                **imshow_kwargs):
    """
    Projected surface-density map of a Dataset (snapshot, halos, galaxies).

    Parameters
    ----------
    ds : Dataset
        Any dataset with "pos" (+ `weight` field). Selection views work.
    projection : "xy", "xz", "yz", ...
    grid_size : int or (Nx, Ny)
    method : "NGP", "CIC", or "Gaussian" (GriddingTools.smooth_to_grid)
    weight : field summed onto the grid (default "mass"; use e.g. "sfr")
    centre, size :
        Optional region: select |delta| < size around centre first (cubic,
        periodic) and set the plot extent to centre +/- size in the
        projected plane.
    species : str, optional
        Snapshot species shortcut (ds.select(species=...)).
    log : plot log10 of the surface density.
    smooth : optional Gaussian smoothing of the final grid (pixels).

    Returns
    -------
    ax : matplotlib Axes (grid available as ax._at_grid)
    """
    if species:
        ds = ds.select(species=species)
    if centre is not None:
        if size is None:
            raise ValueError("density_map(): size= is required with "
                             "centre=.")
        ds = ds.select(centre=centre, size=size, geometry="cubic")

    i, j = _axes_of(projection)
    pos = np.asarray(ds["pos"], dtype=float)
    w = np.asarray(ds[weight], dtype=float) if weight in ds \
        else np.ones(len(ds))

    if centre is not None:
        # wrap into the frame centred on `centre` so the region is
        # contiguous even across the periodic boundary
        box = ds.meta.get("boxsize")
        delta = pos - np.asarray(centre, dtype=float)
        if box:
            delta = (delta + 0.5 * box) % box - 0.5 * box
        p2 = delta[:, (i, j)]
        limits = (-size, size, -size, size)
    else:
        p2 = pos[:, (i, j)]
        box = ds.meta.get("boxsize")
        if box:
            limits = (0.0, float(box), 0.0, float(box))
        else:
            limits = (p2[:, 0].min(), p2[:, 0].max(),
                      p2[:, 1].min(), p2[:, 1].max())

    if isinstance(grid_size, int):
        grid_size = (grid_size, grid_size)

    gt = GriddingTools()
    grid = gt.smooth_to_grid(p2, w, tuple(grid_size), limits,
                             method=method, filter_sigma=smooth)

    # mass per cell -> surface density
    dx = (limits[1] - limits[0]) / grid_size[0]
    dy = (limits[3] - limits[2]) / grid_size[1]
    grid = grid / (dx * dy)

    img = np.log10(grid, out=np.full_like(grid, np.nan),
                   where=grid > 0) if log else grid

    ax = _get_ax(ax)
    im = ax.imshow(img.T, origin="lower",
                   extent=(limits[0], limits[1], limits[2], limits[3]),
                   cmap=cmap, **imshow_kwargs)
    names = "xyz"
    unit = ds.meta.get("units", {}).get("length", "")
    suffix = f" [{unit}]" if unit and "code" not in unit else ""
    ax.set_xlabel(f"{names[i]}{suffix}")
    ax.set_ylabel(f"{names[j]}{suffix}")
    if colorbar:
        import matplotlib.pyplot as plt
        label = (r"$\log_{10}\,\Sigma$" if log else r"$\Sigma$") \
            + (f" ({weight}/area)" if weight != "mass" else "")
        plt.colorbar(im, ax=ax, label=label)
    if ds.label:
        ax.set_title(ds.label)
    ax._at_grid = grid          # numerics stay reachable
    return ax


def overlay_points(ax, ds: Dataset, *, projection: str = "xy",
                   centre: Optional[Sequence[float]] = None,
                   colour_by: Optional[str] = None,
                   log_colour: bool = True, s: float = 30,
                   **scatter_kwargs):
    """
    Scatter a Dataset's positions on existing axes (e.g. halos or SHARK
    galaxies over a density_map). Pass the same `centre` you gave
    density_map so coordinates land in the same (wrapped) frame.
    """
    i, j = _axes_of(projection)
    pos = np.asarray(ds["pos"], dtype=float)
    if centre is not None:
        box = ds.meta.get("boxsize")
        delta = pos - np.asarray(centre, dtype=float)
        if box:
            delta = (delta + 0.5 * box) % box - 0.5 * box
        p2 = delta[:, (i, j)]
    else:
        p2 = pos[:, (i, j)]

    kwargs = dict(edgecolors="white", facecolors="none", linewidths=1.2)
    if colour_by is not None:
        c = np.asarray(ds[colour_by], dtype=float)
        if log_colour:
            c = np.log10(np.clip(c, np.finfo(float).tiny, None))
        kwargs = dict(c=c)
    kwargs.update(scatter_kwargs)
    sc = ax.scatter(p2[:, 0], p2[:, 1], s=s,
                    label=ds.label or None, **kwargs)
    return sc


# ---------------------------------------------------------------------------
# Radial profiles
# ---------------------------------------------------------------------------

#: kind -> (ProfileTools method, needed fields, y key, y label)
_PROFILE_KINDS = {
    "density": ("volume_density", ("pos", "mass"), "density",
                r"$\rho(r)$"),
    "surface_density": ("surface_density", ("pos", "mass"), "density",
                        r"$\Sigma(R)$"),
    "velocity_dispersion": ("velocity_dispersion", ("pos", "vel"), "sigma",
                            r"$\sigma(r)$"),
    "vertical_velocity_dispersion": ("vertical_velocity_dispersion",
                                     ("pos", "vel"), "sigma",
                                     r"$\sigma_z(R)$"),
}


def profile(ds: Dataset, centre: Sequence[float], *, kind: str = "density",
            rmin: float, rmax: float, nbins: int = 40,
            species: Optional[str] = None) -> Dict:
    """
    Radial profile of a Dataset around `centre` via ProfileTools.

    kind : "density" (3D), "surface_density" (cylindrical),
           "velocity_dispersion", "vertical_velocity_dispersion"

    Returns the ProfileTools result dict, annotated with "kind", "ykey",
    "ylabel", and "label" so plot_profile() needs no further context.
    """
    if kind not in _PROFILE_KINDS:
        raise ValueError(f"Unknown profile kind '{kind}' "
                         f"(use one of {sorted(_PROFILE_KINDS)}).")
    if species:
        ds = ds.select(species=species)

    method, fields, ykey, ylabel = _PROFILE_KINDS[kind]
    pt = ProfileTools(numbins=nbins)
    args = [np.asarray(ds[f], dtype=float) for f in fields]
    result = getattr(pt, method)(*args, np.asarray(centre, dtype=float),
                                 rmin, rmax)
    result.update({"kind": kind, "ykey": ykey, "ylabel": ylabel,
                   "label": ds.label})
    return result


def plot_profile(result: Union[Dict, List[Dict]], *, ax=None,
                 logx: bool = True, logy: bool = True,
                 xlabel: str = "r", ylabel: Optional[str] = None,
                 label: Optional[str] = None, **plot_kwargs):
    """Plot one or more profile() results on shared axes."""
    ax = _get_ax(ax, figsize=(6, 5))
    results = result if isinstance(result, (list, tuple)) else [result]
    for res in results:
        ykey = res.get("ykey") or ("density" if "density" in res
                                   else "sigma")
        ax.plot(res["r"], res[ykey],
                label=label or res.get("label"), **plot_kwargs)
    if logx:
        ax.set_xscale("log")
    if logy:
        ax.set_yscale("log")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel or results[0].get("ylabel") or "")
    if any(res.get("label") for res in results) or label:
        ax.legend()
    return ax


# ---------------------------------------------------------------------------
# Mass functions
# ---------------------------------------------------------------------------

def mass_function(datasets: Union[Dataset, Sequence[Dataset]], *,
                  mass_field: str = "mass",
                  volume: Union[str, float, None] = "auto",
                  dlog10m: float = 0.2,
                  bins: Optional[np.ndarray] = None,
                  cumulative: bool = False, ax=None,
                  **plot_kwargs) -> Dict[str, Dict]:
    """
    (Halo or stellar) mass function dn/dlog10 M for one or more Datasets
    -- any mix of halo catalogues and galaxy catalogues.

    volume : "auto" (meta["volume"], else boxsize^3), a number, or None
        (raw counts per bin, no volume normalisation).

    Returns {label: {"log10m", "phi", "counts", "volume"}}; the plot lands
    on `ax` (created if needed, reachable via results["_ax"]).
    """
    if isinstance(datasets, Dataset):
        datasets = [datasets]

    ax = _get_ax(ax, figsize=(6, 5))
    results: Dict[str, Dict] = {}

    for ds in datasets:
        m = np.asarray(ds[mass_field], dtype=float)
        m = m[np.isfinite(m) & (m > 0)]
        logm = np.log10(m)

        if bins is None:
            lo = np.floor(logm.min() / dlog10m) * dlog10m
            hi = np.ceil(logm.max() / dlog10m) * dlog10m
            edges = np.arange(lo, hi + 0.5 * dlog10m, dlog10m)
        else:
            edges = np.asarray(bins, dtype=float)

        counts, edges = np.histogram(logm, bins=edges)
        centres = 0.5 * (edges[1:] + edges[:-1])
        if cumulative:
            counts = np.cumsum(counts[::-1])[::-1]

        vol = None
        if volume == "auto":
            vol = ds.meta.get("volume")
            if not vol and ds.meta.get("boxsize"):
                vol = float(ds.meta["boxsize"]) ** 3
        elif volume is not None:
            vol = float(volume)

        phi = counts.astype(float)
        if vol:
            phi = phi / vol
            if not cumulative:
                phi = phi / np.diff(edges)

        keep = counts > 0
        ax.plot(centres[keep], phi[keep], drawstyle="steps-mid",
                label=ds.label or None,
                color=getattr(ds, "colour", None), **plot_kwargs)
        results[ds.label or f"dataset{len(results)}"] = {
            "log10m": centres, "phi": phi, "counts": counts, "volume": vol,
        }

    ax.set_yscale("log")
    ax.set_xlabel(r"$\log_{10} M$")
    if cumulative:
        ax.set_ylabel(r"$n(>M)$" + (r" [$V^{-1}$]" if vol else ""))
    else:
        ax.set_ylabel(r"$\mathrm{d}n/\mathrm{d}\log_{10}M$"
                      + (r" [$V^{-1}$]" if vol else ""))
    if len(results) > 1 or any(ds.label for ds in datasets):
        ax.legend()
    results["_ax"] = ax
    return results
