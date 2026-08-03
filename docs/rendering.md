# Gridding and Rendering

Three ways to turn particle/cell data into an image, in increasing order of fidelity (and dependency weight):

| Tool | Method | In-process? | Best for |
|------|--------|-------------|----------|
| `GriddingTools` | NGP/CIC/Gaussian deposition onto a Cartesian grid | Yes | Quick, dependency-free diagnostics |
| `VortraceRenderer` | Exact projection through the Arepo Voronoi mesh (`vortrace`) | Yes | Publication-quality Arepo figures, no grid choice needed |
| `ImbasRenderer` | SPH-kernel/CIC projection via the standalone `imbas_renderer` executable | No — shells out | Large snapshots on HPC (MPI-parallel), animations |

## Gridding and Mesh Tools

```python
from analysistools.gridding_tools import GriddingTools

gt = GriddingTools()
grid = gt.smooth_to_grid(positions, values, grid_size, grid_limits)
gt.plot_3d_slice(grid, grid_limits)
```

`smooth_to_grid` supports `method="NGP"`, `"CIC"`, or `"Gaussian"` (CIC + `scipy.ndimage.gaussian_filter`); `method="NGP"` is unsmoothed and will look pixelated on a log colour scale — prefer `"Gaussian"` for anything you're going to look at closely.

## Rendering with vortrace

`VortraceRenderer` wraps [vortrace](https://github.com/gusbeane/vortrace)'s `ProjectionCloud`, which integrates exactly through the Arepo Voronoi mesh rather than depositing onto a Cartesian grid -- no NGP/CIC/Gaussian smoothing choice needed. Requires the `rendering` extra (`pip install -e ".[rendering]"`).

```python
from analysistools import SnapshotTools, VortraceRenderer

snap = SnapshotTools(snapfileformat="HDF5", convention="Arepo")
snap.read_snapshot("snap_010.hdf5")
snap.LoadParticlesByType("gas")

renderer = VortraceRenderer(snap.gas.pos, snap.gas.density)

# single projection along z
fig, ax = renderer.plot_projection(
    extent=[40, 60, 40, 60], npix=512, bounds=[45, 55],
)

# three orthogonal projections through a centre, mirroring GriddingTools.plot_3d_projections
fig, axes = renderer.plot_projections(
    half_extent=10.0, npix=512, half_depth=2.5, centre=[50, 50, 50],
)
```

`vortrace` always integrates along the third column of the `pos` array passed to `ProjectionCloud`; `project_orthogonal`/`plot_projections` handle reordering columns for the other two axes.

## Rendering with imbas_renderer (HPC)

`ImbasRenderer` wraps [imbas_renderer](https://github.com/doctorcbpower/imbas_renderer)'s `render_image`, a standalone OpenMP/MPI-parallel C executable, not a Python library -- it reads the snapshot itself, driven by a YAML config with CLI overrides. This wrapper builds that config/command line; it doesn't pass particle arrays in-process the way `VortraceRenderer` does. Because it's MPI-parallel, it's the renderer to reach for on large snapshots on an HPC cluster.

Quick interactive render (login node or inside an interactive allocation):

```python
from analysistools import ImbasRenderer

renderer = ImbasRenderer(launcher=["mpirun", "-np", "8"])

config = renderer.build_config(
    input="snap_122.hdf5", output="frames/frame",
    xc=56.18, yc=43.35, zc=49.13, lbox=0.2,
    dark_matter=True, scene="cluster",
)
renderer.render(config, config_path="frame.config.yaml")
# -> frames/frame.0000.png
```

For a full batch run, generate a SLURM script instead of running directly, then submit it yourself:

```python
renderer = ImbasRenderer(launcher=["srun"])
renderer.write_config(config, "frame.config.yaml")
renderer.write_sbatch_script(
    "frame.config.yaml", "render.sbatch",
    job_name="disc_render", nodes=2, ntasks=64, time="02:00:00",
    partition="work", modules=["openmpi/4.1"],
)
```

```bash
sbatch render.sbatch
```

`executable` defaults to `render_image` on `PATH`; pass `executable="/path/to/render_image"` if it isn't. `output_paths()` predicts the `<prefix>.NNNN.png` filename(s) render_image will produce, for locating results afterwards.
