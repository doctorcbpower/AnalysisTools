from scipy.spatial import cKDTree
from scipy.ndimage import gaussian_filter

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

import plotly.graph_objects as go
import plotly.express as px

import numpy as np

from scipy.ndimage import gaussian_filter
from scipy.spatial import cKDTree

from numba import njit

@njit
def ngp_assign(grid, coords, values, grid_size):
    """
    Nearest-Grid-Point (NGP) assignment.
    """
    n_particles, dim = coords.shape
    # Ensure grid_size elements are int64
    grid_size = grid_size.astype(np.int64)

    for p in range(n_particles):
        idx = np.empty(dim, dtype=np.int64)
        fx = np.empty(dim, dtype=np.float64)
        
        for d in range(dim):
            x = coords[p, d]
            i = int(np.floor(x))
            idx[d] = i
        
        if dim == 3:
            i, j, k = idx
            grid[i, j, k] += values[p]
        else:
            # Implement 2D case if needed
            raise NotImplementedError("Only 3D CIC implemented.")


@njit
def cic_assign(grid, coords, values, grid_size):
    """
    Cloud-In-Cell (CIC) assignment.
    """
    n_particles, dim = coords.shape
    # Ensure grid_size elements are int64
    grid_size = grid_size.astype(np.int64)

    for p in range(n_particles):
        idx = np.empty(dim, dtype=np.int64)
        fx = np.empty(dim, dtype=np.float64)
        
        for d in range(dim):
            x = coords[p, d]
            i = int(np.floor(x))
            f = x - i
            idx[d] = i
            fx[d] = f
        
        if dim == 3:
            i, j, k = idx
            fx_i, fx_j, fx_k = fx
            
            # Compute weights
            w000 = (1-fx_i)*(1-fx_j)*(1-fx_k)
            w100 = fx_i*(1-fx_j)*(1-fx_k)
            w010 = (1-fx_i)*fx_j*(1-fx_k)
            w001 = (1-fx_i)*(1-fx_j)*fx_k
            w101 = fx_i*(1-fx_j)*fx_k
            w011 = (1-fx_i)*fx_j*fx_k
            w110 = fx_i*fx_j*(1-fx_k)
            w111 = fx_i*fx_j*fx_k
            
            # Periodic boundaries
            Nx, Ny, Nz = grid_size
            grid[i % Nx, j % Ny, k % Nz] += values[p] * w000
            grid[(i+1) % Nx, j % Ny, k % Nz] += values[p] * w100
            grid[i % Nx, (j+1) % Ny, k % Nz] += values[p] * w010
            grid[i % Nx, j % Ny, (k+1) % Nz] += values[p] * w001
            grid[(i+1) % Nx, j % Ny, (k+1) % Nz] += values[p] * w101
            grid[i % Nx, (j+1) % Ny, (k+1) % Nz] += values[p] * w011
            grid[(i+1) % Nx, (j+1) % Ny, k % Nz] += values[p] * w110
            grid[(i+1) % Nx, (j+1) % Ny, (k+1) % Nz] += values[p] * w111

        else:
            # Implement 2D case if needed
            raise NotImplementedError("Only 3D CIC implemented.")

class GriddingTools:
    def __init__(self):
        pass

    def smooth_to_grid(self, positions, values, grid_size, grid_limits,
                       method="NGP", sigma=1.0, filter_sigma=None):
        """
        Assign particle values to a 2D or 3D grid.

        Args:
            positions (ndarray): (N,2) or (N,3) array of input coordinates.
            values (ndarray): (N,) array of particle values.
            grid_size (tuple): Grid shape (Nx, Ny[, Nz]).
            grid_limits (tuple): (xmin, xmax, ymin, ymax[, zmin, zmax]).
            method (str): "NGP", "CIC", or "Gaussian".
            sigma (float): Gaussian width (for 'Gaussian' method).
            filter_sigma (float): Optional Gaussian smoothing of final grid.

        Returns:
            ndarray: Grid of assigned values.
        """
        dim = len(grid_size)
        grid = np.zeros(grid_size, dtype=float)

        # Grid spacing
        spacing = [(grid_limits[2*i+1] - grid_limits[2*i]) / grid_size[i] 
                   for i in range(dim)]

        # Normalized coordinates in [0, grid_size)
        coords = [(positions[:,i] - grid_limits[2*i]) / spacing[i] for i in range(dim)]
        coords = np.stack(coords, axis=1)

        if method.upper() == "NGP":
            ngp_assign(grid, coords, values, grid_size)
        elif method.upper() == "CIC":
            cic_assign(grid, coords, values, grid_size)
        elif method.upper() == "GAUSSIAN":
            cic_assign(grid, coords, values, grid_size)
            grid = gaussian_filter(grid, sigma=sigma)

        # Optional additional smoothing
        if filter_sigma is not None:
            grid = gaussian_filter(grid, sigma=filter_sigma)

        return grid
            
    def plot_3d_slice(self, grid_3d, grid_limits,
                 slice_axis='z', slice_index=None,
                 slice_width=None, slice_average=True,
                 mode='slice', projection='mean',
                 title="3D Grid Slice", cmap='viridis', figsize=(12, 4)):
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
        """
        nx, ny, nz = grid_3d.shape

        # Handle mode
        if mode == 'slice':
            if slice_axis == 'z':
                if slice_index is None:
                    slice_index = nz // 2
    
                if slice_width is None:
                    slice_width = 1  # e.g., 1 cell thick; can make this a function argument

                # Compute start/end indices safely
                start_idx = max(slice_index - slice_width//2, 0)
                end_idx = min(slice_index + slice_width//2 + 1, nz)

                # Extract thin slice (average or sum over the thickness)
                if slice_average is True:
                    data = grid_3d[:, :, start_idx:end_idx].mean(axis=2)
                else:
                    data = grid_3d[:, :, start_idx:end_idx].sum(axis=2)
                    
                extent = [grid_limits[0], grid_limits[1], grid_limits[2], grid_limits[3]]
                xlabel, ylabel, title_str = 'X', 'Y', f'XY slice (Z={slice_index})'

            elif slice_axis == 'y':
                if slice_index is None:
                    slice_index = ny // 2

                if slice_width is None:
                    slice_width = 1  # e.g., 1 cells thick; can make this a function argument

                # Compute start/end indices safely
                start_idx = max(slice_index - slice_width//2, 0)
                end_idx = min(slice_index + slice_width//2 + 1, ny)

                # Extract thin slice (average or sum over the thickness)
                if slice_average is True:
                    data = grid_3d[:, start_idx:end_idx, :].mean(axis=1)
                else:
                    data = grid_3d[:, start_idx:end_idx, :].sum(axis=1)
                    
                extent = [grid_limits[0], grid_limits[1], grid_limits[4], grid_limits[5]]
                xlabel, ylabel, title_str = 'X', 'Z', f'XZ slice (Y={slice_index})'

            elif slice_axis == 'x':
                if slice_index is None:
                    slice_index = nx // 2

                if slice_width is None:
                    slice_width = 1  # e.g., 3 cells thick; can make this a function argument

                # Compute start/end indices safely
                start_idx = max(slice_index - slice_width//2, 0)
                end_idx = min(slice_index + slice_width//2 + 1, nx)

                # Extract thin slice (average or sum over the thickness)
                if slice_average is True:
                    data = grid_3d[start_idx:end_idx, :, :].mean(axis=0)
                else:
                    data = grid_3d[start_idx:end_idx, :, :].sum(axis=0)
                    
                extent = [grid_limits[2], grid_limits[3], grid_limits[4], grid_limits[5]]
                xlabel, ylabel, title_str = 'Y', 'Z', f'YZ slice (X={slice_index})'

            else:
                raise ValueError("slice_axis must be 'x', 'y', or 'z'")

        elif mode == 'projection':
            if slice_axis == 'z':
                if projection == 'mean':
                    data = grid_3d.mean(axis=2)
                elif projection == 'sum':
                    data = grid_3d.sum(axis=2)
                elif projection == 'max':
                    data = grid_3d.max(axis=2)
                extent = [grid_limits[0], grid_limits[1], grid_limits[2], grid_limits[3]]
                xlabel, ylabel, title_str = 'X', 'Y', f'XY projection ({projection} along Z)'

            elif slice_axis == 'y':
                if projection == 'mean':
                    data = grid_3d.mean(axis=1)
                elif projection == 'sum':
                    data = grid_3d.sum(axis=1)
                elif projection == 'max':
                    data = grid_3d.max(axis=1)
                extent = [grid_limits[0], grid_limits[1], grid_limits[4], grid_limits[5]]
                xlabel, ylabel, title_str = 'X', 'Z', f'XZ projection ({projection} along Y)'

            elif slice_axis == 'x':
                if projection == 'mean':
                    data = grid_3d.mean(axis=0)
                elif projection == 'sum':
                    data = grid_3d.sum(axis=0)
                elif projection == 'max':
                    data = grid_3d.max(axis=0)
                extent = [grid_limits[2], grid_limits[3], grid_limits[4], grid_limits[5]]
                xlabel, ylabel, title_str = 'Y', 'Z', f'YZ projection ({projection} along X)'

            else:
                raise ValueError("slice_axis must be 'x', 'y', or 'z'")

        else:
            raise ValueError("mode must be 'slice' or 'projection'")

        # Plot
        fig, ax = plt.subplots(figsize=figsize)
        im = ax.imshow(np.log10(data.T), origin='lower', extent=extent, cmap=cmap, aspect='auto')
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title_str)
        fig.suptitle(title, fontsize=14)
        plt.colorbar(im, ax=ax)
        plt.tight_layout()
        return fig, ax

    def plot_3d_projections(self,grid_3d, grid_limits,
                            mode='projection', projection='sum',
                            cmap='viridis', figsize=(12, 4), title=None,
                            slice_axis='z', slice_index=None,
                            slice_width=None, slice_average=True,):
        """
        Wrapper to plot a row of 3 orthogonal projections (XY, XZ, YZ)
        by calling plot_3d_slice for each axis.

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

        Returns
        -------
        fig : matplotlib.figure.Figure
        axes : list of matplotlib.axes.Axes
        """
        fig, axes = plt.subplots(1, 3, figsize=figsize)

        for ax, axis, lbl in zip(
            axes,
            ['z', 'y', 'x'],
            ['XY', 'XZ', 'YZ']
        ):
            # Handle slice_index per axis
            idx = None
            if mode == 'slice' and slice_index is not None:
                if isinstance(slice_index, dict):
                    idx = slice_index.get(axis, None)
                else:
                    idx = slice_index
        
            # Call your existing slice function in projection mode
            _, single_ax = self.plot_3d_slice(
                grid_3d, grid_limits,
                slice_axis=axis,
                mode=mode,
                projection=projection,
                cmap=cmap,
                figsize=(5, 5)  # ignored, since we reuse fig/ax
            )

            # Transfer the image from the returned ax to our chosen subplot
            im = single_ax.images[0]
            ax.imshow(im.get_array(), origin='lower',
                      extent=im.get_extent(),
                      cmap=im.get_cmap() if hasattr(im, 'get_cmap') else cmap,
                      aspect='auto')  # set aspect manually
            ax.set_title(lbl)
            ax.set_xlabel(single_ax.get_xlabel())
            ax.set_ylabel(single_ax.get_ylabel())
            plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
            plt.close(single_ax.figure)  # avoid extra figures popping up

        if title:
            fig.suptitle(title, fontsize=14)

        plt.tight_layout()
        return fig, axes

