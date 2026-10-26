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

`smooth_to_grid` supports `method="NGP"`, `"CIC"`, `"Gaussian"` (CIC + `scipy.ndimage.gaussian_filter`) and `"SPH"` (3D grids only; needs per-particle `smoothing_lengths`, e.g. each Voronoi cell's own effective radius `(3V/4π)^(1/3)`). `grid_size` is `(Nx, Ny)` or `(Nx, Ny, Nz)` and `grid_limits` is `(xmin, xmax, ymin, ymax[, zmin, zmax])`; 2D grids take `(N, 2)` positions, so you can deposit a projection directly. `method="NGP"` is unsmoothed and will look pixelated on a log colour scale — prefer `"Gaussian"` for anything you're going to look at closely. Grid arrays are indexed `[x, y]`, so transpose (`.T`) before `imshow(..., origin="lower")`.

## Rendering with vortrace

`VortraceRenderer` wraps [vortrace](https://github.com/gusbeane/vortrace)'s `ProjectionCloud`, which integrates exactly through the Arepo Voronoi mesh rather than depositing onto a Cartesian grid -- no NGP/CIC/Gaussian smoothing choice needed. Requires the `rendering` extra (`pip install -e ".[rendering]"`).

```python
from analysistools import SnapshotTools, VortraceRenderer

snap = SnapshotTools(snapfileformat="HDF5", convention="AREPO")
data = snap.read_snapshot("snap_010.hdf5")     # returns SnapshotData
data.LoadParticlesByType("gas")

renderer = VortraceRenderer(data.gas.pos, data.gas.density)

# single projection along z
fig, ax = renderer.plot_projection(
    extent=[40, 60, 40, 60], npix=512, bounds=[45, 55],
)

# three orthogonal projections through a centre, mirroring GriddingTools.plot_3d_projections
fig, axes = renderer.plot_projections(
    half_extent=10.0, npix=512, half_depth=2.5, centre=[50, 50, 50],
)
```

`vortrace` always integrates along the third column of the `pos` array passed to `ProjectionCloud`; `project_orthogonal`/`plot_projections` handle reordering columns for the other two axes. `project_orthogonal(axis, ...)` returns `(image, extent)`, where `axis` is the line of sight; `image` is indexed `[x, y]` of the in-plane axes (`'z'`: x,y; `'y'`: x,z; `'x'`: y,z), so plot `image.T` with `origin="lower"`.

A projection of `rho` is a **column density**. Because vortrace integrates whatever per-cell value you hand it, a density-weighted quantity is the ratio of two projections (see the galaxy recipe below).

## Recipe: face-on and edge-on views of an isolated galaxy

Gas via vortrace (density-weighted density, ∫ρ² dz / ∫ρ dz), stars via the internal gridding (stellar surface density). Stars are particles, not mesh cells, so they go through `GriddingTools`, not vortrace.

```python
import numpy as np
import matplotlib.pyplot as plt
from analysistools import SnapshotTools, VortraceRenderer
from analysistools.gridding_tools import GriddingTools

snap = SnapshotTools(snapfileformat="HDF5", convention="AREPO")
data = snap.read_snapshot("snap_010.hdf5")
data.LoadParticlesByType("gas")
data.LoadParticlesByType("star")
gas, star = data.gas, data.star

# centre on the stellar centre of mass
centre = np.average(star.pos, weights=star.mass, axis=0)
gpos, spos = gas.pos - centre, star.pos - centre

# optional: rotate so the gas angular momentum (inside r < 5, code units) is along +z
sel = np.linalg.norm(gpos, axis=1) < 5.0
L = np.sum(np.cross(gpos[sel], gas.vel[sel]) * gas.mass[sel, None], axis=0)
zhat = L / np.linalg.norm(L)
xhat = np.cross([0, 1, 0], zhat); xhat /= np.linalg.norm(xhat)
R = np.vstack([xhat, np.cross(zhat, xhat), zhat])      # rows = new axes
gpos, spos = gpos @ R.T, spos @ R.T

half, npix = 15.0, 512          # same units as the snapshot positions

# gas: density-weighted density, face-on (line of sight z) and edge-on (y)
num = VortraceRenderer(gpos, gas.density**2)
den = VortraceRenderer(gpos, gas.density)
fig, axes = plt.subplots(2, 2, figsize=(10, 10))
for ax, axis, lbl in zip(axes[0], ["z", "y"], ["face-on", "edge-on"]):
    n, ext = num.project_orthogonal(axis, half, npix, half_depth=half)
    d, _   = den.project_orthogonal(axis, half, npix, half_depth=half)
    ax.imshow(np.log10((n / d).T), origin="lower", extent=ext, cmap="magma")
    ax.set_title(f"gas, {lbl}")

# stars: mass surface density on a 2D grid
gt = GriddingTools()
lim, dA = (-half, half, -half, half), (2 * half / npix) ** 2
for ax, cols, lbl in zip(axes[1], ([0, 1], [0, 2]), ["face-on", "edge-on"]):
    m = gt.smooth_to_grid(spos[:, cols], star.mass, (npix, npix), lim,
                          method="Gaussian", sigma=1.5)
    with np.errstate(divide="ignore"):
        ax.imshow(np.log10(m.T / dA), origin="lower", extent=lim, cmap="gray")
    ax.set_title(f"stars, {lbl}")
plt.tight_layout()
```

If the disc already lies in the xy plane (typical for an isolated-galaxy IC), skip the rotation. To avoid vortrace for gas, deposit `gas.mass * gas.density` and `gas.mass` with `smooth_to_grid` and divide; that gives a mass-weighted density per pixel.

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
