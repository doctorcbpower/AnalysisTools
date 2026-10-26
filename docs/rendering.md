# Gridding and Rendering

Ways to turn particle/cell data into an image, in increasing order of fidelity (and dependency weight):

| Tool | Method | In-process? | Best for |
|------|--------|-------------|----------|
| `GriddingTools` | NGP/CIC/Gaussian deposition onto a Cartesian grid | Yes | Quick, dependency-free diagnostics |
| `VortraceRenderer` | Exact projection through the Arepo Voronoi mesh (`vortrace`) | Yes | Publication-quality Arepo figures, no grid choice needed |
| `analysistools.projection` | Backend-neutral SPH projection (numba or cosmokit deposit), weighting, periodic boxes, orientation | Yes | Quantitative maps from particles/cells; see below |
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

## Shared projection layer (`analysistools.projection`)

One backend-neutral set of routines for SPH-smoothed projections, so that figures from different deposit engines are directly comparable. It handles particle selection, periodic boxes, orientation, smoothing lengths and weighting; a *backend* only has to supply a k-nearest-neighbour search, an SPH density, and a 3D deposit summed along the line of sight.

```python
from analysistools import projection as pj

H   = pj.smoothing_lengths(pos, k=32)                    # support radius = distance to k-th neighbour
rho = pj.sph_density(pos, mass, k=32)                    # local density for particles that carry none (stars)
R   = pj.orientation_from_angular_momentum(pos - c, vel, mass, r_max=5.0)

m = pj.project(pos, mass, rotation=R, view="face", centre=c, half_extent=15.0, npix=512)
m.plot(cmap="magma")                                     # Map.image is [ix, iy]; plot() handles the transpose
```

| Backend | Selected by | Notes |
|---|---|---|
| `numba` (default) | `backend="numba"` | scipy k-d tree + `GriddingTools.sph_assign`; no extra dependency |
| `cosmokit` | `backend="cosmokit"` | C/OpenMP k-d tree and SPH deposit; needs `cosmokit` on `PYTHONPATH`. Isotropic cells (the `nz` hint is ignored) and a padded cubic mesh, so memory is ~`nmesh³` doubles; it raises with advice if that exceeds 4 GiB. Typically 1.5-2x faster. |

`backend="auto"` is `numba` unless `ANALYSISTOOLS_PROJECTION_BACKEND` is set. Any object with `knn_h`, `sph_density` and `deposit_projected` can be passed as `backend`.

Key points:

- **Weighting.** `project(pos, weights)` is the surface density Σw/dA. With `weight_by=f` it returns Σ w f / Σ w along the line of sight. Since a cell's mass is m = ρV, `weights=mass, weight_by=rho` is exactly the ρ-weighted density ∫ρ² dz / ∫ρ dz.
- **Periodic boxes.** Pass `boxsize=L`. Positions are wrapped about `centre`, and periodic images are added when the view region (plus kernel margin) is larger than half the box, so views that straddle or exceed the box edge are correct.
- **Edges.** Particles just outside the window whose kernels reach in are included (`drop_outside` deposit): mass that falls off the grid is lost, not piled onto edge pixels.
- **Smoothing.** `H` is the distance to the k-th neighbour (k=32 by default) for both backends, so they agree exactly; the cubic spline has support `H`. Pass `h=` to use your own, e.g. a cell's own radius.
- **`sph_density` bias.** `include_self=True` (default, as in SPH and cosmokit) adds the particle's own kernel and biases a Poisson-sampled field high (~20% for k=48); use `include_self=False` for an unbiased estimate.
- From cosmokit: `import cosmokit.projection as cp; cp.project(...)` delegates to the same code with `backend="cosmokit"` (AnalysisTools must be importable).

### Recipe: face-on and edge-on views of an isolated galaxy

```python
import numpy as np
import matplotlib.pyplot as plt
from analysistools import SnapshotTools, projection as pj

snap = SnapshotTools(snapfileformat="HDF5", convention="AREPO")
data = snap.read_snapshot("snap_010.hdf5")
data.LoadParticlesByType("gas")
data.LoadParticlesByType("star")
gas, star = data.gas, data.star

centre = np.average(star.pos, weights=star.mass, axis=0)
R = pj.orientation_from_angular_momentum(gas.pos - centre, gas.vel, gas.mass, r_max=5.0)

kw = dict(rotation=R, centre=centre, half_extent=15.0, npix=512)   # snapshot length units
fig, axes = plt.subplots(2, 2, figsize=(10, 10))
for col, view in enumerate(["face", "edge"]):
    # gas: rho-weighted density = sum(m*rho)/sum(m)
    g = pj.project(gas.pos, gas.mass, weight_by=gas.density, view=view, **kw)
    g.plot(ax=axes[0, col], cmap="magma"); axes[0, col].set_title(f"gas, {view}-on")
    # stars: surface density
    s = pj.project(star.pos, star.mass, view=view, **kw)
    s.plot(ax=axes[1, col], cmap="gray"); axes[1, col].set_title(f"stars, {view}-on")
plt.tight_layout()
```

If the disc already lies in the xy plane, drop `rotation=R`. To weight stars by a stellar density instead, use `weight_by=pj.sph_density(star.pos, star.mass)`.

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
