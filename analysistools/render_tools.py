"""
Rendering hook for vortrace (https://github.com/gusbeane/vortrace), which performs
exact projections through Voronoi meshes. This is a natural fit for Arepo output:
cell positions + densities already define the tessellation vortrace integrates
through, so no NGP/CIC/Gaussian grid deposition (see GriddingTools) is needed.

vortrace is an optional dependency -- install with:
    pip install analysistools[rendering]
or directly from its repo:
    pip install git+https://github.com/gusbeane/vortrace

vortrace.ProjectionCloud.grid_projection always integrates along the *third*
column of the (N, 3) `pos` array passed to it, over `bounds`, onto an image
plane spanned by the first two columns, over `extent`. To project along a
different line of sight, reorder the columns of `pos` before constructing the
cloud -- see `VortraceRenderer.project_orthogonal`.
"""
import numpy as np
import matplotlib.pyplot as plt


class VortraceRenderer:
    """Thin wrapper around vortrace.ProjectionCloud for AnalysisTools snapshot data.

    Parameters
    ----------
    pos : ndarray (N, 3)
        Cell/particle coordinates, e.g. SnapshotTools.pos for a given ptype.
    rho : ndarray (N,) or (N, 4)
        Density values, or RGBA values per cell for volume rendering
        (reduction='volume'), e.g. SnapshotTools.rho.
    """

    def __init__(self, pos, rho):
        try:
            import vortrace as vt
        except ImportError as exc:
            raise ImportError(
                "VortraceRenderer requires the 'vortrace' package. Install with "
                "`pip install analysistools[rendering]` or "
                "`pip install git+https://github.com/gusbeane/vortrace`."
            ) from exc

        self.pos = np.asarray(pos)
        self.rho = np.asarray(rho)
        self._cloud = vt.ProjectionCloud(self.pos, self.rho)

    def project(self, extent, npix, bounds, reduction='sum'):
        """
        Integrate through the Voronoi mesh along the cloud's third coordinate.

        Parameters
        ----------
        extent : list [xmin, xmax, ymin, ymax]
            Image-plane region, in the same (absolute, not centre-relative)
            units as `pos`.
        npix : int
            Output resolution (npix x npix).
        bounds : list [zmin, zmax]
            Integration bounds along the line of sight.
        reduction : str, optional
            'sum' for a standard projection (returns (npix, npix)), or
            'volume' for RGBA volume rendering (returns (npix, npix, 3));
            requires `rho` to have been passed as an (N, 4) RGBA array.

        Returns
        -------
        ndarray
        """
        return self._cloud.grid_projection(extent, npix, bounds, reduction=reduction)

    def project_orthogonal(self, axis, half_extent, npix, half_depth, centre=None,
                            reduction='sum'):
        """
        Project along one of the three coordinate axes, centred on `centre`.

        Reorders `pos`/`rho` so the requested axis becomes the line of sight
        (vortrace always integrates along the third column), then delegates
        to a fresh ProjectionCloud -- vortrace has no API to change the
        integration axis on an existing cloud.

        Parameters
        ----------
        axis : {'x', 'y', 'z'}
            Axis to integrate along (i.e. the line of sight).
        half_extent : float
            Half-width of the square image plane, centred on `centre`.
        npix : int
            Output resolution (npix x npix).
        half_depth : float
            Half-depth of the integration range along `axis`, centred on `centre`.
        centre : ndarray (3,), optional
            Centre of the projected region. Defaults to the origin.
        reduction : str, optional
            See `project`.

        Returns
        -------
        image : ndarray
        extent : list [xmin, xmax, ymin, ymax] of the two in-plane axes.
        """
        import vortrace as vt

        if centre is None:
            centre = np.zeros(3)
        centre = np.asarray(centre)

        axis_order = {'z': (0, 1, 2), 'y': (0, 2, 1), 'x': (1, 2, 0)}
        if axis not in axis_order:
            raise ValueError("axis must be 'x', 'y', or 'z'")
        perm = axis_order[axis]

        pos_reordered = self.pos[:, perm]
        centre_reordered = centre[list(perm)]

        extent = [
            centre_reordered[0] - half_extent, centre_reordered[0] + half_extent,
            centre_reordered[1] - half_extent, centre_reordered[1] + half_extent,
        ]
        bounds = [centre_reordered[2] - half_depth, centre_reordered[2] + half_depth]

        cloud = vt.ProjectionCloud(pos_reordered, self.rho)
        image = cloud.grid_projection(extent, npix, bounds, reduction=reduction)
        return image, extent

    def plot_projection(self, extent, npix, bounds, reduction='sum', log=True,
                         cmap='viridis', figsize=(6, 6), xlabel='X', ylabel='Y',
                         title=None):
        """Project and render as a single imshow figure. Returns (fig, ax)."""
        data = self.project(extent, npix, bounds, reduction=reduction)

        fig, ax = plt.subplots(figsize=figsize)
        if reduction == 'volume':
            im = ax.imshow(np.transpose(data, (1, 0, 2)), origin='lower',
                            extent=extent, aspect='auto')
        else:
            with np.errstate(divide="ignore"):
                plot_data = np.log10(data.T) if log else data.T
            im = ax.imshow(plot_data, origin='lower', extent=extent, cmap=cmap,
                            aspect='auto')
            plt.colorbar(im, ax=ax)

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        if title:
            ax.set_title(title)
        plt.tight_layout()
        return fig, ax

    def plot_projections(self, half_extent, npix, half_depth, centre=None,
                          reduction='sum', log=True, cmap='viridis',
                          figsize=(15, 5), title=None):
        """
        Three orthogonal projections (XY, XZ, YZ) through `centre`, mirroring
        GriddingTools.plot_3d_projections. Returns (fig, axes).
        """
        fig, axes = plt.subplots(1, 3, figsize=figsize)
        axis_labels = {'z': ('X', 'Y'), 'y': ('X', 'Z'), 'x': ('Y', 'Z')}

        for ax, axis, lbl in zip(axes, ['z', 'y', 'x'], ['XY', 'XZ', 'YZ']):
            image, extent = self.project_orthogonal(
                axis, half_extent, npix, half_depth, centre=centre,
                reduction=reduction,
            )

            if reduction == 'volume':
                plot_data = np.transpose(image, (1, 0, 2))
            else:
                with np.errstate(divide="ignore"):
                    plot_data = np.log10(image.T) if log else image.T

            im = ax.imshow(plot_data, origin='lower', extent=extent,
                            cmap=cmap if reduction != 'volume' else None,
                            aspect='auto')
            ax.set_title(lbl)
            xlabel, ylabel = axis_labels[axis]
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            if reduction != 'volume':
                plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

        if title:
            fig.suptitle(title, fontsize=14)

        plt.tight_layout()
        return fig, axes
