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
        # numba-jitted assigners need an ndarray (tuple indexing is
        # branch-checked statically and fails for 2D input)
        grid_size = np.asarray(grid_size, dtype=np.int64)
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
