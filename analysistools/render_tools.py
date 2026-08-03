"""
Rendering hooks for two external, publication-quality renderers, kept in one
module because they're both alternatives to GriddingTools' NGP/CIC grid
deposition -- but they have very different shapes:

- VortraceRenderer wraps vortrace (https://github.com/gusbeane/vortrace), an
  in-process Python library that performs exact projections through Voronoi
  meshes. A natural fit for Arepo output: cell positions + densities already
  define the tessellation vortrace integrates through.

- ImbasRenderer wraps imbas_renderer
  (https://github.com/doctorcbpower/imbas_renderer), a standalone, OpenMP/MPI-
  parallel C executable (`render_image`) that reads the snapshot itself and is
  driven by a YAML config with CLI overrides, rather than numpy arrays passed
  in-process. It supports SPH-kernel and CIC projection for GADGET4/SWIFT/Arepo
  snapshots, and -- since it's MPI-parallel -- is the renderer of choice for
  large snapshots on HPC systems. This wrapper's job is to build that
  YAML/CLI invocation (and, optionally, a SLURM sbatch script) and shell out
  to it; it does not talk to the executable in-process.

Both are optional dependencies -- install with:
    pip install analysistools[rendering]
vortrace has no PyPI release (installed from GitHub via pyproject.toml).
imbas_renderer is a compiled binary and is not pip-installable at all; build
it separately (see its README) and either put `render_image` on PATH or pass
`executable=<path>` to ImbasRenderer.

vortrace.ProjectionCloud.grid_projection always integrates along the *third*
column of the (N, 3) `pos` array passed to it, over `bounds`, onto an image
plane spanned by the first two columns, over `extent`. To project along a
different line of sight, reorder the columns of `pos` before constructing the
cloud -- see `VortraceRenderer.project_orthogonal`.
"""
import shlex
import subprocess

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


class ImbasRenderer:
    """Wrapper around the imbas_renderer `render_image` executable.

    Unlike VortraceRenderer this never touches particle arrays directly --
    imbas_renderer reads the snapshot itself, so this class's job is to
    assemble the YAML config it expects, build the command line (optionally
    under an MPI/SLURM launcher), and either run it directly or write a
    batch script for later submission on an HPC system.

    Parameters
    ----------
    executable : str, optional
        Path to the `render_image` binary. Defaults to 'render_image' (must
        then be on PATH).
    launcher : list of str, optional
        Command prefix used to invoke the executable under MPI, e.g.
        ['mpirun', '-np', '32'] or ['srun']. Defaults to [] (run directly,
        OpenMP-only).
    """

    def __init__(self, executable="render_image", launcher=None):
        self.executable = executable
        self.launcher = list(launcher) if launcher else []

    @staticmethod
    def build_config(input, output, xc, yc, zc, lbox, is_hdf5=True, **extra):
        """
        Assemble the YAML config dict for a render.

        Parameters
        ----------
        input : str
            Snapshot path.
        output : str
            Output prefix; render_image writes `<output>.NNNN.png`.
        xc, yc, zc, lbox : float
            View centre and box side length, in snapshot units.
        is_hdf5 : bool, optional
            False for GADGET binary input.
        **extra
            Passed through verbatim as additional YAML keys (e.g.
            dark_matter, scene, colormap, opacity_function, rotate, nframes,
            ...) -- see the imbas_renderer README for the full option set.

        Returns
        -------
        dict
        """
        config = {
            "input": input,
            "output": output,
            "isHDF5": bool(is_hdf5),
            "xc": xc, "yc": yc, "zc": zc,
            "lbox": lbox,
        }
        config.update(extra)
        return config

    @staticmethod
    def write_config(config, path):
        """Write `config` (from build_config) to a YAML file at `path`."""
        import yaml
        with open(path, "w") as f:
            yaml.safe_dump(config, f, sort_keys=False)
        return path

    def command(self, config_path, cli_overrides=None):
        """
        Build the full argv list, without running it -- e.g. for embedding
        in an sbatch script rather than running interactively.

        Parameters
        ----------
        config_path : str
            Path to a YAML file written by `write_config`.
        cli_overrides : dict, optional
            Extra `-key value` flags appended after `-config`, taking
            precedence over the YAML file per imbas_renderer's three-tier
            config precedence.

        Returns
        -------
        list of str
        """
        cmd = list(self.launcher) + [self.executable, "-config", str(config_path)]
        for key, value in (cli_overrides or {}).items():
            cmd += [f"-{key}", str(value)]
        return cmd

    def render(self, config, config_path=None, cli_overrides=None,
               check=True, **subprocess_kwargs):
        """
        Write `config` to YAML and run render_image directly (optionally
        under `self.launcher`). Suitable for a quick test render on a login
        node or inside an interactive allocation; for a full HPC batch run,
        use `write_sbatch_script` instead so the job goes through the
        scheduler.

        Returns
        -------
        subprocess.CompletedProcess
        """
        if config_path is None:
            config_path = f"{config.get('output', 'imbas_render')}.config.yaml"
        self.write_config(config, config_path)

        cmd = self.command(config_path, cli_overrides=cli_overrides)
        return subprocess.run(cmd, check=check, **subprocess_kwargs)

    def write_sbatch_script(self, config_path, script_path, cli_overrides=None,
                             job_name="imbas_render", nodes=1, ntasks=32,
                             time="01:00:00", partition=None, account=None,
                             modules=None, extra_directives=None):
        """
        Write a SLURM batch script that runs render_image under
        `self.launcher` (typically ['srun'] or ['mpirun', '-np', str(ntasks)])
        against `config_path`. Does not submit it -- run `sbatch script_path`
        yourself once you've checked it over.

        Parameters
        ----------
        config_path : str
            Path to a YAML file written by `write_config`.
        script_path : str
            Where to write the batch script.
        cli_overrides : dict, optional
            See `command`.
        job_name, nodes, ntasks, time, partition, account : various
            Standard SLURM directives.
        modules : list of str, optional
            Modules to `module load` before running (e.g. ['openmpi/4.1']).
        extra_directives : list of str, optional
            Additional raw `#SBATCH ...` lines.

        Returns
        -------
        str
            `script_path`, for convenience.
        """
        lines = [
            "#!/bin/bash",
            f"#SBATCH --job-name={job_name}",
            f"#SBATCH --nodes={nodes}",
            f"#SBATCH --ntasks={ntasks}",
            f"#SBATCH --time={time}",
        ]
        if partition:
            lines.append(f"#SBATCH --partition={partition}")
        if account:
            lines.append(f"#SBATCH --account={account}")
        for directive in (extra_directives or []):
            lines.append(f"#SBATCH {directive}")
        lines.append("")

        for module in (modules or []):
            lines.append(f"module load {module}")
        if modules:
            lines.append("")

        cmd = self.command(config_path, cli_overrides=cli_overrides)
        lines.append(" ".join(shlex.quote(c) for c in cmd))
        lines.append("")

        with open(script_path, "w") as f:
            f.write("\n".join(lines))
        return script_path

    @staticmethod
    def output_paths(output_prefix, nframes=None):
        """
        Predict the PNG path(s) render_image will produce, following its
        `<prefix>.NNNN.png` numbering.

        Parameters
        ----------
        output_prefix : str
        nframes : int, optional
            If given, returns a list of `nframes` paths (0000, 0001, ...).
            Otherwise returns the single-frame path (0000).

        Returns
        -------
        str or list of str
        """
        if nframes is None:
            return f"{output_prefix}.0000.png"
        return [f"{output_prefix}.{i:04d}.png" for i in range(nframes)]
