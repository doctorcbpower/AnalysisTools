from scipy.spatial import cKDTree
from scipy.ndimage import gaussian_filter
from scipy.spatial import Voronoi
from matplotlib.collections import LineCollection

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

import plotly.graph_objects as go
import plotly.express as px

import numpy as np

from scipy.ndimage import gaussian_filter
from scipy.spatial import cKDTree

from numba import njit


# ---------------------------------------------------------------------------
# Jitted assignment kernels. 2D and 3D are separate functions because numba
# type-checks every branch against the actual grid dimensionality, so a
# single kernel with an `if dim == 3` cannot compile for both.
# ---------------------------------------------------------------------------

@njit
def _ngp_assign_2d(grid, coords, values, Nx, Ny):
    for p in range(coords.shape[0]):
        i = int(np.floor(coords[p, 0]))
        j = int(np.floor(coords[p, 1]))
        grid[i % Nx, j % Ny] += values[p]


@njit
def _ngp_assign_3d(grid, coords, values, Nx, Ny, Nz):
    for p in range(coords.shape[0]):
        i = int(np.floor(coords[p, 0]))
        j = int(np.floor(coords[p, 1]))
        k = int(np.floor(coords[p, 2]))
        grid[i % Nx, j % Ny, k % Nz] += values[p]


@njit
def _cic_assign_2d(grid, coords, values, Nx, Ny):
    for p in range(coords.shape[0]):
        x, y = coords[p, 0], coords[p, 1]
        i, j = int(np.floor(x)), int(np.floor(y))
        fx, fy = x - i, y - j
        grid[i % Nx, j % Ny] += values[p] * (1-fx)*(1-fy)
        grid[(i+1) % Nx, j % Ny] += values[p] * fx*(1-fy)
        grid[i % Nx, (j+1) % Ny] += values[p] * (1-fx)*fy
        grid[(i+1) % Nx, (j+1) % Ny] += values[p] * fx*fy


@njit
def _cic_assign_3d(grid, coords, values, Nx, Ny, Nz):
    for p in range(coords.shape[0]):
        x, y, z = coords[p, 0], coords[p, 1], coords[p, 2]
        i, j, k = int(np.floor(x)), int(np.floor(y)), int(np.floor(z))
        fx, fy, fz = x - i, y - j, z - k

        grid[i % Nx, j % Ny, k % Nz] += values[p] * (1-fx)*(1-fy)*(1-fz)
        grid[(i+1) % Nx, j % Ny, k % Nz] += values[p] * fx*(1-fy)*(1-fz)
        grid[i % Nx, (j+1) % Ny, k % Nz] += values[p] * (1-fx)*fy*(1-fz)
        grid[i % Nx, j % Ny, (k+1) % Nz] += values[p] * (1-fx)*(1-fy)*fz
        grid[(i+1) % Nx, j % Ny, (k+1) % Nz] += values[p] * fx*(1-fy)*fz
        grid[i % Nx, (j+1) % Ny, (k+1) % Nz] += values[p] * (1-fx)*fy*fz
        grid[(i+1) % Nx, (j+1) % Ny, k % Nz] += values[p] * fx*fy*(1-fz)
        grid[(i+1) % Nx, (j+1) % Ny, (k+1) % Nz] += values[p] * fx*fy*fz


@njit
def _cubic_spline_kernel(q):
    """Normalized 3D cubic spline (Monaghan & Lattanzio 1985), compact
    support out to q=2. Caller supplies the 1/(pi h^3) normalization by
    construction of the deposit weights (see _sph_assign_3d) -- this
    returns the dimensionless shape function only."""
    if q < 1.0:
        return 1.0 - 1.5 * q * q + 0.75 * q * q * q
    elif q < 2.0:
        t = 2.0 - q
        return 0.25 * t * t * t
    else:
        return 0.0


@njit
def _sph_assign_3d(grid, pos, values, h, xmin, ymin, zmin, dx, dy, dz, Nx, Ny, Nz):
    """
    Per-particle SPH-kernel deposition onto a 3D grid, using each
    particle's own smoothing length `h` (compact-support cubic spline,
    support radius 2h) rather than one global smoothing scale applied to
    every particle alike.

    Works entirely in physical coordinates (not grid-index units), so it
    is correct for anisotropic grid spacing (dx != dy != dz, e.g. a
    thinner z-extent than x/y). Non-periodic: kernel support is clipped
    to the grid, not wrapped (unlike _ngp_assign_3d/_cic_assign_3d), since
    this targets isolated (non-cosmological) galaxy examples.

    Each particle's deposited weights are renormalized to sum to exactly
    `values[p]` (computed in a first pass over its own affected cells,
    then applied in a second pass) -- this is what keeps a particle's
    contribution mass-conserving regardless of whether its own h happens
    to be small or large relative to the grid spacing, which a single
    global smoothing scale cannot do (see the CIC+Gaussian gap/blob
    artifact this method replaces). A particle whose kernel support
    doesn't overlap any cell centre (h << spacing) falls back to
    nearest-grid-point deposition so its value is never silently dropped.
    """
    n = pos.shape[0]
    for p in range(n):
        hp = h[p]
        x, y, z = pos[p, 0], pos[p, 1], pos[p, 2]

        if hp <= 0.0:
            i = int(np.floor((x - xmin) / dx))
            j = int(np.floor((y - ymin) / dy))
            k = int(np.floor((z - zmin) / dz))
            if 0 <= i < Nx and 0 <= j < Ny and 0 <= k < Nz:
                grid[i, j, k] += values[p]
            continue

        rmax = 2.0 * hp
        imin = max(int(np.floor((x - xmin - rmax) / dx)), 0)
        imax = min(int(np.ceil((x - xmin + rmax) / dx)), Nx - 1)
        jmin = max(int(np.floor((y - ymin - rmax) / dy)), 0)
        jmax = min(int(np.ceil((y - ymin + rmax) / dy)), Ny - 1)
        kmin = max(int(np.floor((z - zmin - rmax) / dz)), 0)
        kmax = min(int(np.ceil((z - zmin + rmax) / dz)), Nz - 1)

        if imin > imax or jmin > jmax or kmin > kmax:
            continue

        wsum = 0.0
        for i in range(imin, imax + 1):
            cx = xmin + (i + 0.5) * dx
            ddx = cx - x
            for j in range(jmin, jmax + 1):
                cy = ymin + (j + 0.5) * dy
                ddy = cy - y
                for k in range(kmin, kmax + 1):
                    cz = zmin + (k + 0.5) * dz
                    ddz = cz - z
                    r = np.sqrt(ddx * ddx + ddy * ddy + ddz * ddz)
                    q = r / hp
                    if q < 2.0:
                        wsum += _cubic_spline_kernel(q)

        if wsum <= 0.0:
            i = int(np.floor((x - xmin) / dx))
            j = int(np.floor((y - ymin) / dy))
            k = int(np.floor((z - zmin) / dz))
            if 0 <= i < Nx and 0 <= j < Ny and 0 <= k < Nz:
                grid[i, j, k] += values[p]
            continue

        v = values[p]
        for i in range(imin, imax + 1):
            cx = xmin + (i + 0.5) * dx
            ddx = cx - x
            for j in range(jmin, jmax + 1):
                cy = ymin + (j + 0.5) * dy
                ddy = cy - y
                for k in range(kmin, kmax + 1):
                    cz = zmin + (k + 0.5) * dz
                    ddz = cz - z
                    r = np.sqrt(ddx * ddx + ddy * ddy + ddz * ddz)
                    q = r / hp
                    if q < 2.0:
                        w = _cubic_spline_kernel(q)
                        grid[i, j, k] += v * (w / wsum)


class GriddingTools:
    def __init__(self):
        pass

    @staticmethod
    def ngp_assign(grid, coords, values, grid_size):
        """
        Nearest-Grid-Point (NGP) assignment (2D or 3D, periodic wrap).
        """
        gs = np.asarray(grid_size, dtype=np.int64)
        if grid.ndim == 3:
            _ngp_assign_3d(grid, coords, values, gs[0], gs[1], gs[2])
        else:
            _ngp_assign_2d(grid, coords, values, gs[0], gs[1])

    @staticmethod
    def cic_assign(grid, coords, values, grid_size):
        """
        Cloud-In-Cell (CIC) assignment (2D or 3D, periodic wrap).
        """
        gs = np.asarray(grid_size, dtype=np.int64)
        if grid.ndim == 3:
            _cic_assign_3d(grid, coords, values, gs[0], gs[1], gs[2])
        else:
            _cic_assign_2d(grid, coords, values, gs[0], gs[1])

    @staticmethod
    def sph_assign(grid, positions, values, smoothing_lengths, grid_limits, grid_size):
        """
        Per-particle SPH-kernel (cubic spline) deposition onto a 3D grid,
        using each particle's own smoothing length rather than a single
        global smoothing scale. See _sph_assign_3d for the deposition
        details (mass-conserving per particle, non-periodic, correct for
        anisotropic grid spacing). 3D only -- unlike ngp_assign/cic_assign,
        there is no 2D variant (no current caller needs one).
        """
        if grid.ndim != 3:
            raise ValueError("sph_assign only supports 3D grids")
        dx = (grid_limits[1] - grid_limits[0]) / grid_size[0]
        dy = (grid_limits[3] - grid_limits[2]) / grid_size[1]
        dz = (grid_limits[5] - grid_limits[4]) / grid_size[2]
        _sph_assign_3d(grid, positions, values, smoothing_lengths,
                        grid_limits[0], grid_limits[2], grid_limits[4],
                        dx, dy, dz, grid_size[0], grid_size[1], grid_size[2])

    def smooth_to_grid(self, positions, values, grid_size, grid_limits,
                       method="NGP", sigma=1.0, filter_sigma=None,
                       smoothing_lengths=None):
        """
        Assign particle values to a 2D or 3D grid.

        Args:
            positions (ndarray): (N,2) or (N,3) array of input coordinates.
            values (ndarray): (N,) array of particle values.
            grid_size (tuple): Grid shape (Nx, Ny[, Nz]).
            grid_limits (tuple): (xmin, xmax, ymin, ymax[, zmin, zmax]).
            method (str): "NGP", "CIC", "Gaussian", or "SPH".
            sigma (float): Gaussian width in grid cells (for 'Gaussian' method).
            filter_sigma (float): Optional Gaussian smoothing of final grid.
            smoothing_lengths (ndarray): (N,) array of per-particle physical
                smoothing lengths, same units as `positions`/`grid_limits`.
                Required for method="SPH" -- e.g. each cell's own effective
                radius (3*Volume/(4*pi))**(1/3) for a Voronoi/moving-mesh
                code, where a single global Gaussian sigma either
                over-smooths dense regions or under-samples (leaves
                particle-scale gaps in) sparse ones.

        Returns:
            ndarray: Grid of assigned values.
        """
        dim = len(grid_size)
        grid = np.zeros(grid_size, dtype=float)
        # numba-jitted assigners need an ndarray (tuple indexing is
        # branch-checked statically and fails for 2D input)
        grid_size = np.asarray(grid_size, dtype=np.int64)

        if method.upper() == "SPH":
            if smoothing_lengths is None:
                raise ValueError("method='SPH' requires smoothing_lengths")
            if dim != 3:
                raise ValueError("method='SPH' only supports 3D grids")
            self.sph_assign(grid, positions, values, smoothing_lengths,
                             grid_limits, grid_size)
            if filter_sigma is not None:
                grid = gaussian_filter(grid, sigma=filter_sigma)
            return grid

        # Grid spacing
        spacing = [(grid_limits[2*i+1] - grid_limits[2*i]) / grid_size[i]
                   for i in range(dim)]

        # Normalized coordinates in [0, grid_size)
        coords = [(positions[:,i] - grid_limits[2*i]) / spacing[i] for i in range(dim)]
        coords = np.stack(coords, axis=1)

        if method.upper() == "NGP":
            self.ngp_assign(grid, coords, values, grid_size)
        elif method.upper() == "CIC":
            self.cic_assign(grid, coords, values, grid_size)
        elif method.upper() == "GAUSSIAN":
            self.cic_assign(grid, coords, values, grid_size)
            grid = gaussian_filter(grid, sigma=sigma)

        # Optional additional smoothing
        if filter_sigma is not None:
            grid = gaussian_filter(grid, sigma=filter_sigma)

        return grid
            
    @staticmethod
    def _collapse(grid_3d, grid_limits, slice_axis='z', slice_index=None,
                  slice_width=None, slice_average=True, mode='projection',
                  projection='mean', normalize_by_area=False):
        """
        Shared collapse logic for plot_3d_slice/plot_3d_projections/
        plot_ratio_projections: reduces a 3D grid to a 2D (data, extent,
        xlabel, ylabel, title_str) along one axis, either a thin slice or a
        full projection (mean/sum/max).

        normalize_by_area (only meaningful with projection='sum'): divides
        the summed data by the physical area of one pixel in the two
        surviving axes -- turns a per-column *total* (e.g. total mass in a
        column, code mass units) into a per-unit-area *density* (mass per
        unit area, in [input mass units] / [input grid_limits units]^2).
        Callers wanting a specific physical unit (e.g. Msun/pc^2) should
        pre-scale the values passed into smooth_to_grid accordingly (see
        GriddingTools docstring / calling code) -- this method has no
        notion of physical units itself, only of the numbers it's given.
        """
        nx, ny, nz = grid_3d.shape

        if mode == 'slice':
            if slice_axis == 'z':
                n = nz
                if slice_index is None:
                    slice_index = n // 2
                if slice_width is None:
                    slice_width = 1
                start_idx = max(slice_index - slice_width // 2, 0)
                end_idx = min(slice_index + slice_width // 2 + 1, n)
                sub = grid_3d[:, :, start_idx:end_idx]
                data = sub.mean(axis=2) if slice_average else sub.sum(axis=2)
                extent = [grid_limits[0], grid_limits[1], grid_limits[2], grid_limits[3]]
                xlabel, ylabel, title_str = 'X', 'Y', f'XY slice (Z={slice_index})'

            elif slice_axis == 'y':
                n = ny
                if slice_index is None:
                    slice_index = n // 2
                if slice_width is None:
                    slice_width = 1
                start_idx = max(slice_index - slice_width // 2, 0)
                end_idx = min(slice_index + slice_width // 2 + 1, n)
                sub = grid_3d[:, start_idx:end_idx, :]
                data = sub.mean(axis=1) if slice_average else sub.sum(axis=1)
                extent = [grid_limits[0], grid_limits[1], grid_limits[4], grid_limits[5]]
                xlabel, ylabel, title_str = 'X', 'Z', f'XZ slice (Y={slice_index})'

            elif slice_axis == 'x':
                n = nx
                if slice_index is None:
                    slice_index = n // 2
                if slice_width is None:
                    slice_width = 1
                start_idx = max(slice_index - slice_width // 2, 0)
                end_idx = min(slice_index + slice_width // 2 + 1, n)
                sub = grid_3d[start_idx:end_idx, :, :]
                data = sub.mean(axis=0) if slice_average else sub.sum(axis=0)
                extent = [grid_limits[2], grid_limits[3], grid_limits[4], grid_limits[5]]
                xlabel, ylabel, title_str = 'Y', 'Z', f'YZ slice (X={slice_index})'

            else:
                raise ValueError("slice_axis must be 'x', 'y', or 'z'")

        elif mode == 'projection':
            collapse_axis = {'z': 2, 'y': 1, 'x': 0}.get(slice_axis)
            if collapse_axis is None:
                raise ValueError("slice_axis must be 'x', 'y', or 'z'")

            if projection == 'mean':
                data = grid_3d.mean(axis=collapse_axis)
            elif projection == 'sum':
                data = grid_3d.sum(axis=collapse_axis)
            elif projection == 'max':
                data = grid_3d.max(axis=collapse_axis)
            else:
                raise ValueError("projection must be 'mean', 'sum', or 'max'")

            if slice_axis == 'z':
                extent = [grid_limits[0], grid_limits[1], grid_limits[2], grid_limits[3]]
                xlabel, ylabel, title_str = 'X', 'Y', f'XY projection ({projection} along Z)'
            elif slice_axis == 'y':
                extent = [grid_limits[0], grid_limits[1], grid_limits[4], grid_limits[5]]
                xlabel, ylabel, title_str = 'X', 'Z', f'XZ projection ({projection} along Y)'
            else:
                extent = [grid_limits[2], grid_limits[3], grid_limits[4], grid_limits[5]]
                xlabel, ylabel, title_str = 'Y', 'Z', f'YZ projection ({projection} along X)'

            if normalize_by_area:
                pixel_area = ((extent[1] - extent[0]) / data.shape[0]) * \
                             ((extent[3] - extent[2]) / data.shape[1])
                data = data / pixel_area

        else:
            raise ValueError("mode must be 'slice' or 'projection'")

        return data, extent, xlabel, ylabel, title_str

    @staticmethod
    def _plot_panel(ax, data, extent, xlabel, ylabel, title_str, cmap='viridis',
                     vmin=None, vmax=None, log=True):
        """Shared imshow + colorbar for a single (already-collapsed) 2D panel."""
        with np.errstate(divide='ignore', invalid='ignore'):
            plotted = np.log10(data.T) if log else data.T
        im = ax.imshow(plotted, origin='lower', extent=extent, cmap=cmap,
                        aspect='auto', vmin=vmin, vmax=vmax)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title_str)
        return im

    @staticmethod
    def _shared_clim(data_list, log=True):
        """Finite-value min/max across several panels' data, for a colour
        scale shared across all of them by default (rather than each panel
        auto-scaling independently, which makes the panels incomparable)."""
        with np.errstate(divide='ignore', invalid='ignore'):
            all_vals = np.concatenate([
                (np.log10(d) if log else d).ravel() for d in data_list
            ])
        finite = all_vals[np.isfinite(all_vals)]
        if finite.size == 0:
            return None, None
        return float(finite.min()), float(finite.max())

    def plot_3d_slice(self, grid_3d, grid_limits,
                 slice_axis='z', slice_index=None,
                 slice_width=None, slice_average=True,
                 mode='slice', projection='mean',
                 title="3D Grid Slice", cmap='viridis', figsize=(12, 4),
                 vmin=None, vmax=None, normalize_by_area=False,
                 cbar_label=None):
        """
        Visualize 3D grid by showing orthogonal slices or projections.

        Parameters
        ----------
        grid_3d : ndarray
            3D array of data.
        grid_limits : list
            [xmin, xmax, ymin, ymax, zmin, zmax]
        slice_axis : str, optional
            Axis along which to slice/project: 'x', 'y', or 'z'
        slice_index : int, optional
            Index at which to slice (used if mode='slice'). Defaults to centre.
        slice_width : int, optional
            Thickness of slice (used if mode='slice'). Defaults to 1.
        slice_average : bool, optional
            Sum or average over cells in slice (used if mode='slice'). Defaults to average.
        mode : str, optional
            'slice' = single slice at slice_index
            'projection' = collapse along slice_axis
        projection : str, optional
            If mode='projection', how to collapse: 'mean', 'sum', 'max'
        vmin, vmax : float, optional
            Fixed log10-scale colour limits. Default: matplotlib auto-scales
            to this panel's own data range.
        normalize_by_area : bool, optional
            With mode='projection', projection='sum': divide the summed
            column total by pixel area, turning it into a density. See
            _collapse's docstring. Default: False (raw sum).
        cbar_label : str, optional
            Colorbar label (e.g. units).
        """
        data, extent, xlabel, ylabel, title_str = self._collapse(
            grid_3d, grid_limits, slice_axis=slice_axis, slice_index=slice_index,
            slice_width=slice_width, slice_average=slice_average, mode=mode,
            projection=projection, normalize_by_area=normalize_by_area,
        )

        fig, ax = plt.subplots(figsize=figsize)
        im = self._plot_panel(ax, data, extent, xlabel, ylabel, title_str,
                               cmap=cmap, vmin=vmin, vmax=vmax)
        fig.suptitle(title, fontsize=14)
        cbar = plt.colorbar(im, ax=ax)
        if cbar_label is not None:
            cbar.set_label(cbar_label)
        plt.tight_layout()
        return fig, ax

    def plot_3d_projections(self,grid_3d, grid_limits,
                            mode='projection', projection='sum',
                            cmap='viridis', figsize=(12, 4), title=None,
                            slice_axis='z', slice_index=None,
                            slice_width=None, slice_average=True,
                            vmin=None, vmax=None, normalize_by_area=False,
                            cbar_label=None):
        """
        Plot a row of 3 orthogonal projections (XY, XZ, YZ) of the same
        3D grid.

        Parameters
        ----------
        grid_3d : ndarray
            3D array of data values.
        grid_limits : list
            [xmin, xmax, ymin, ymax, zmin, zmax]
        projection : str, optional
            How to collapse along the chosen axis ('mean', 'sum', 'max').
        cmap : str, optional
            Colormap to use for imshow.
        figsize : tuple, optional
            Figure size (width, height).
        title : str, optional
            Title for the figure.
        slice_index : int, optional
            Index at which to slice (used if mode='slice'). Defaults to centre.
        slice_width : int, optional
            Thickness of slice (used if mode='slice'). Defaults to 1.
        slice_average : bool, optional
            Sum or average over cells in slice (used if mode='slice'). Defaults to average.
        vmin, vmax : float, optional
            Fixed log10-scale colour limits shared by all three panels.
            Default: computed from the combined finite data range across
            all three panels, so they share one scale by default (unlike
            each panel auto-scaling independently, which makes side-by-side
            panels visually incomparable).
        normalize_by_area : bool, optional
            With projection='sum': divide each panel's summed column total
            by its own pixel area, turning a raw sum into a density. See
            GriddingTools._collapse's docstring for units caveats.
        cbar_label : str, optional
            Shared colorbar label (e.g. units) for all three panels.

        Returns
        -------
        fig : matplotlib.figure.Figure
        axes : list of matplotlib.axes.Axes
        """
        fig, axes = plt.subplots(1, 3, figsize=figsize)

        axis_labels = ['z', 'y', 'x']
        panels = []
        for axis in axis_labels:
            idx = slice_index
            if mode == 'slice' and isinstance(slice_index, dict):
                idx = slice_index.get(axis, None)
            panels.append(self._collapse(
                grid_3d, grid_limits, slice_axis=axis, slice_index=idx,
                slice_width=slice_width, slice_average=slice_average, mode=mode,
                projection=projection, normalize_by_area=normalize_by_area,
            ))

        if vmin is None or vmax is None:
            auto_vmin, auto_vmax = self._shared_clim([p[0] for p in panels])
            vmin = auto_vmin if vmin is None else vmin
            vmax = auto_vmax if vmax is None else vmax

        im = None
        for ax, lbl, (data, extent, xlabel, ylabel, _) in zip(axes, ['XY', 'XZ', 'YZ'], panels):
            im = self._plot_panel(ax, data, extent, xlabel, ylabel, lbl,
                                   cmap=cmap, vmin=vmin, vmax=vmax)
            cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
            if cbar_label is not None:
                cbar.set_label(cbar_label)

        if title:
            fig.suptitle(title, fontsize=14)

        plt.tight_layout()
        return fig, axes

    def plot_ratio_projections(self, numerator_3d, denominator_3d, grid_limits,
                               projection='sum', cmap='viridis', figsize=(12, 4),
                               title=None, vmin=None, vmax=None,
                               min_denominator=0.0, cbar_label=None):
        """
        Plot a row of 3 orthogonal projections of the ratio of two 3D
        grids built on the same grid_size/grid_limits (e.g. dust mass /
        gas mass, for a dust-to-gas-ratio map) -- collapses each grid with
        `projection` independently, then divides, so numerator and
        denominator are each properly summed (or averaged) before the
        ratio is taken, not divided cell-by-cell before projection.

        Any mass-unit/area-normalization factor cancels in the ratio, so
        (unlike plot_3d_projections) there's no normalize_by_area option
        here -- the two grids just need to share units with each other,
        not represent any particular physical unit.

        Parameters
        ----------
        numerator_3d, denominator_3d : ndarray
            3D grids of the same shape (e.g. from two smooth_to_grid calls
            with the same grid_size/grid_limits).
        min_denominator : float, optional
            Columns whose collapsed denominator is <= this are masked to
            NaN (not plotted) rather than producing a divide-by-zero
            artifact. Default: 0.0 (mask exact zeros only).
        (see plot_3d_projections for the remaining parameters)

        Returns
        -------
        fig : matplotlib.figure.Figure
        axes : list of matplotlib.axes.Axes
        """
        if numerator_3d.shape != denominator_3d.shape:
            raise ValueError("numerator_3d and denominator_3d must have the same shape")

        fig, axes = plt.subplots(1, 3, figsize=figsize)

        axis_labels = ['z', 'y', 'x']
        panels = []
        for axis in axis_labels:
            num, extent, xlabel, ylabel, _ = self._collapse(
                numerator_3d, grid_limits, slice_axis=axis, mode='projection',
                projection=projection,
            )
            den, _, _, _, _ = self._collapse(
                denominator_3d, grid_limits, slice_axis=axis, mode='projection',
                projection=projection,
            )
            with np.errstate(divide='ignore', invalid='ignore'):
                ratio = np.where(den > min_denominator, num / den, np.nan)
            panels.append((ratio, extent, xlabel, ylabel))

        if vmin is None or vmax is None:
            auto_vmin, auto_vmax = self._shared_clim([p[0] for p in panels])
            vmin = auto_vmin if vmin is None else vmin
            vmax = auto_vmax if vmax is None else vmax

        for ax, lbl, (data, extent, xlabel, ylabel) in zip(axes, ['XY', 'XZ', 'YZ'], panels):
            im = self._plot_panel(ax, data, extent, xlabel, ylabel, lbl,
                                   cmap=cmap, vmin=vmin, vmax=vmax)
            cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
            if cbar_label is not None:
                cbar.set_label(cbar_label)

        if title:
            fig.suptitle(title, fontsize=14)

        plt.tight_layout()
        return fig, axes
        
    def slice_nearest_neighbour(self, positions, field,
                                 width, height,
                                 Xpixels, Ypixels,
                                 slice_axis='z', slice_value=0.0,
                                 max_distance=None):
        """
        Produce a 2D image by nearest-neighbour lookup along a slice plane.

        Parameters
        ----------
        positions : ndarray (N, 3)
            Particle positions, already centred on the region of interest.
        field : ndarray (N,)
            Field values (e.g. density) per particle.
        width, height : float
            Physical extent of the output image.
        Xpixels, Ypixels : int
            Output image resolution.
        slice_axis : str
            Axis normal to the slice plane: 'x', 'y', or 'z'.
        slice_value : float
            Position along slice_axis to query (default 0, i.e. centred).
        max_distance: float
            Maximum distance a particle can be from the plane
        Returns
        -------
        image : ndarray (Ypixels, Xpixels)
        extent : list [xmin, xmax, ymin, ymax]
        """
        axis_map = {'x': 0, 'y': 1, 'z': 2}
        normal   = axis_map[slice_axis]
        plane    = [i for i in range(3) if i != normal]

        tree = cKDTree(positions)

        x = np.linspace(-width/2,  width/2,  Xpixels)
        y = np.linspace(-height/2, height/2, Ypixels)
        X, Y = np.meshgrid(x, y)

        # Build query points in 3D
        query = np.zeros((Xpixels * Ypixels, 3))
        query[:, plane[0]] = X.ravel()
        query[:, plane[1]] = Y.ravel()
        query[:, normal]   = slice_value

        dist, idx = tree.query(query, k=1)
        field_mapped = field[idx].astype(float)

        # Mask pixels where nearest neighbour was too far from the slice plane
        if max_distance is not None:
            field_mapped[dist > max_distance] = np.nan

        image  = field[idx].reshape(Ypixels, Xpixels)
        extent = [-width/2, width/2, -height/2, height/2]

        return image, extent

    def plot_voronoi_slice(self, positions,
                           width, height,
                           slice_axis='z', slice_value=0.0,
                           max_distance=None,
                           figsize=(8, 8), title=None):
        """
        Draw the 2D Voronoi tessellation of particles within `max_distance`
        of a slice plane, using mirrored points to bound edge cells cleanly.

        Parameters
        ----------
        positions : (N, 3) array
        width, height : float
            Extent of the plotting box along the two in-plane axes.
        slice_axis : {'x', 'y', 'z'}, optional
            Axis normal to the slice plane. Default 'z'.
        slice_value : float, optional
            Coordinate of the slice plane along `slice_axis`. Default 0.0.
        max_distance : float, optional
            Half-thickness of the slab selected around the slice plane.
            Defaults to an estimate of the mean particle spacing in-plane.
        figsize : tuple, optional
            Matplotlib figure size. Default (8, 8).
        title : str, optional

        Returns
        -------
        (fig, ax) : matplotlib Figure and Axes with the tessellation drawn.
        """
        axis_map = {'x': 0, 'y': 1, 'z': 2}
        normal   = axis_map[slice_axis]
        plane_ax = [i for i in range(3) if i != normal]

        # Estimate default before mask is applied
        if max_distance is None:
            n_total = np.sum(
                (positions[:, plane_ax[0]] > -width/2)  &
                (positions[:, plane_ax[0]] <  width/2)  &
                (positions[:, plane_ax[1]] > -height/2) &
                (positions[:, plane_ax[1]] <  height/2)
            )
            max_distance = np.sqrt((width * height) / n_total) if n_total > 0 else width * 0.01

        mask = (
            (np.abs(positions[:, normal] - slice_value) < max_distance) &
            (positions[:, plane_ax[0]] > -width/2)  &
            (positions[:, plane_ax[0]] <  width/2)  &
            (positions[:, plane_ax[1]] > -height/2) &
            (positions[:, plane_ax[1]] <  height/2)
        )

        pts_2d = positions[mask][:, plane_ax]

        print(f"Plotting Voronoi for {mask.sum()} particles, "
              f"mean spacing ~ {max_distance:.4f}")

        # Mirror points to cleanly bound edge cells — no infinite ridges
        pad = max(width, height) * 2
        for px, py in [( pad, 0), (-pad, 0), (0,  pad), (0, -pad),
                       ( pad,  pad), (-pad,  pad),
                       ( pad, -pad), (-pad, -pad)]:
            pass  # built below
        mirrors = np.array([
            [ pad,  0], [-pad,  0], [0,  pad], [0, -pad],
            [ pad,  pad], [-pad,  pad], [ pad, -pad], [-pad, -pad]
        ])
        pts_ext = np.vstack([pts_2d, mirrors])

        vor = Voronoi(pts_ext)

        # Collect only finite ridge segments, clipped to plot box
        xmin, xmax = -width/2,  width/2
        ymin, ymax = -height/2, height/2

        segments = []
        for ridge in vor.ridge_vertices:
            if -1 in ridge:
                continue  # skip infinite ridges entirely
            p1, p2 = vor.vertices[ridge[0]], vor.vertices[ridge[1]]
            # Only keep segments whose midpoint is inside the plot box
            mid = (p1 + p2) / 2
            if xmin <= mid[0] <= xmax and ymin <= mid[1] <= ymax:
                segments.append([p1, p2])

        fig, ax = plt.subplots(figsize=figsize)
        lc = LineCollection(segments, colors='k', linewidths=0.4)
        ax.add_collection(lc)

        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)

        axis_labels = ['x (pc)', 'y (pc)', 'z (pc)']
        ax.set_xlabel(axis_labels[plane_ax[0]])
        ax.set_ylabel(axis_labels[plane_ax[1]])
        ax.set_aspect('equal')

        if title:
            ax.set_title(title)

        plt.tight_layout()
        return fig, ax
