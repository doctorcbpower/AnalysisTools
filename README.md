# AnalysisTools

**AnalysisTools** is a Python toolkit for working with particle-based simulation data, designed for **cosmological simulations and galaxy formation studies**.

It provides tools to:

- Read and write simulation snapshots (HDF5 and GADGET binary)
- Extract and manipulate particle data by type, and build initial conditions (merge, taper, rescale)
- Select particles by geometry, and grid/smooth particle fields onto meshes
- Render publication-quality Voronoi-mesh projections of Arepo output via [vortrace](https://github.com/gusbeane/vortrace) (optional dependency)
- Drive [imbas_renderer](https://github.com/doctorcbpower/imbas_renderer), an MPI-parallel SPH/CIC renderer for GADGET4/SWIFT/Arepo snapshots, including generating SLURM batch scripts for HPC use
- Read and standardise halo catalogues (SubFind, AHF, VELOCIraptor, SWIFT FOF)
- Walk halo merger trees (SubFind-HBT, TreeFrog/VELOCIraptor, AHF MergerTree) and analyse orbits/infall
- Compute radial/vertical density and kinematic profiles
- Read FDM (fuzzy dark matter) wavefunction field snapshots
- Provide lazy, cached, multi-snapshot access to halo catalogues across redshift (`HaloModel`)
- Read and analyse SHARK semi-analytic galaxy catalogues, including star formation histories and fsps-based photometry (`shark` package)

See [DEVELOPMENT.md](DEVELOPMENT.md) for the roadmap towards a unified interface across all data sources.

---

## Repository Layout

| Path | Contents |
|------|----------|
| `analysistools/` | Core package: snapshot, halo catalogue, merger tree, profile, gridding, rendering, FDM tools; unified `api/` layer |
| `analysistools/shark/` | SHARK semi-analytic catalogue tools: `SharkModel` (lazy reader), `Analysis`, `Plotter`, CLI, fsps-based photometry |
| `shark/` | Deprecated import shim → `analysistools.shark` |
| `docs/` | Per-topic usage guides (see Documentation below) |
| `data/` | Small example data: `snap_0031.hdf5`, VELOCIraptor walkable trees, `halos/` catalogues |
| `notebooks/` | `UnifiedInterfaceDemo.ipynb` — executed end-to-end example of the unified interface (`PlotSimulationAnalysis.ipynb` is the legacy example) |

---

## Installation

```bash
git clone https://github.com/doctorcbpower/AnalysisTools.git
cd AnalysisTools
pip install -e .
```

The package is imported as `analysistools`; submodules use relative imports (e.g. `halo_tools.py` imports from `.haloio_subfind`), so it must be installed or run as a package rather than as loose scripts.

### Optional Extras

| Extra | Adds | Used by |
|-------|------|---------|
| `test` | `pytest` | Test suite |
| `notebooks` | `nbformat`, `nbclient`, `ipykernel`, `matplotlib-inline` | Executing example notebooks |
| `photometry` | `fsps` | `shark.photometry` |
| `rendering` | `vortrace` (from GitHub; not on PyPI), `pyyaml` | `render_tools.VortraceRenderer`, `render_tools.ImbasRenderer` |

```bash
pip install -e ".[rendering]"
```

`render_tools` itself imports fine without `vortrace`/`pyyaml` installed -- only instantiating `VortraceRenderer`, or calling `ImbasRenderer.write_config`/`render`, requires them. `imbas_renderer`'s `render_image` executable is a separate C build (see its own README) and is not part of this extra -- put it on `PATH` or pass `executable=<path>` to `ImbasRenderer`.

---

## Module Overview

| Module | Main class(es) / functions | Purpose |
|--------|------------------------------|---------|
| `snapshot_tools.py` | `SnapshotTools`, `select_particles`, `place_points_in_mesh` | Read/write particle snapshots, split by type, geometric selection |
| `snapio_hdf5.py` | `read_hdf5`, `write_hdf5` | Low-level HDF5 snapshot I/O (used by `SnapshotTools`) |
| `snapio_binary.py` | `read_binary` | Low-level GADGET binary format 1/2 snapshot I/O |
| `ic_tools.py` | `ICTools` (extends `SnapshotTools`) | Build initial conditions: merge components, taper halo, rescale masses, excise types |
| `gridding_tools.py` | `GriddingTools` | NGP/CIC mesh assignment, smoothing, Voronoi/nearest-neighbour slices, 3D grid plotting |
| `render_tools.py` | `VortraceRenderer`, `ImbasRenderer` | Exact Voronoi-mesh projections via `vortrace` (in-process); YAML/CLI + SLURM-script wrapper around the standalone, MPI-parallel `imbas_renderer` executable |
| `profile_tools.py` | `ProfileTools`, `MassFunctionTools` | Volume/surface density, velocity dispersion, scale height profiles; halo mass functions |
| `halo_tools.py` | `HaloTools` | Unified reader for halo catalogues (SubFind, AHF, VELOCIraptor, SWIFT FOF) |
| `haloio_subfind.py`, `haloio_ahf.py`, `haloio_velociraptor.py`, `haloio_swiftfof.py` | `read_subfind`, `read_ahf`, `read_velociraptor`, `read_swiftfof` | Format-specific halo catalogue readers (used by `HaloTools`) |
| `halo_tools_standardise_names.py` | `standardise_catalogue_names` | Maps native per-format field names onto a common schema (`mass`, `pos`, `vel`, `radius`, `halo_id`, `num_part`) |
| `halo_model.py` | `HaloModel` | Lazy, cached access to a halo catalogue across multiple snapshots/redshifts, with derived-quantity accessors (`mvir`, `rvir`, `vmax`, `position`, ...) |
| `merger_tree_tools.py` | `MergerTreeTools` | Unified interface to walk merger trees and plot/analyse tracks; dispatches to `treeio_*` readers |
| `treeio_subfind.py`, `treeio_treefrog.py`, `treeio_ahf.py` | `read_subfind_hbt`/`walk_subfind`, `read_treefrog`/`walk_treefrog`, `read_ahf_mergertree`/`walk_ahf` | Format-specific merger tree readers and main-branch walkers |
| `merger_tree_types.py` | `HaloTrack`, `OrbitAnalysis`, `MergerTreeError`, `periodic_delta` | Shared, format-agnostic containers used by `MergerTreeTools` |
| `fdm_field_tools.py` | `FDMFieldTools` | Read FDM wavefunction field snapshots (mesh data, not particles); density, radial profile, slices |
| `galaxy_tools.py` | `GalaxyTools` | Legacy SHARK catalogue reader — superseded by `shark.model.SharkModel` (see [docs/shark.md](docs/shark.md)) |

---

## Documentation

The unified interface is summarised below; everything else is a per-topic guide in `docs/`:

| Topic | Guide |
|-------|-------|
| Snapshots (`SnapshotTools`, `ICTools`) | [docs/snapshots.md](docs/snapshots.md) |
| Halo catalogues (`HaloTools`, `HaloModel`) | [docs/halo_catalogues.md](docs/halo_catalogues.md) |
| Merger trees (`MergerTreeTools`) | [docs/merger_trees.md](docs/merger_trees.md) |
| Profiles and mass functions (`ProfileTools`) | [docs/profiles.md](docs/profiles.md) |
| Gridding and rendering (`GriddingTools`, `VortraceRenderer`, `ImbasRenderer`/HPC) | [docs/rendering.md](docs/rendering.md) |
| FDM field snapshots (`FDMFieldTools`) | [docs/fdm.md](docs/fdm.md) |
| SHARK semi-analytic catalogues (`shark` package) | [docs/shark.md](docs/shark.md) |
| Unified `Dataset`/`Simulation` interface (full example) | [docs/unified_interface.md](docs/unified_interface.md) |

All of the per-class API below remains fully supported -- the unified interface is a layer on top of these classes, and `ds.backend` always exposes the underlying object. Prefer the unified interface for new work.

---

## Unified Interface (new)

`analysistools.load()` opens any supported data product with one call — kind and format sniffed from the file — and returns a `Dataset` with common syntax across sources (snapshots, halo catalogues, merger trees, and SHARK galaxy catalogues):

```python
import analysistools as at

snap  = at.load("snap_0031.hdf5")                      # SWIFT/GADGET sniffed
halos = at.load("snap_0031.VELOCIraptor.properties.0") # VELOCIraptor sniffed

halos["mass"], halos["pos"], snap["pos"]     # same field vocabulary everywhere
snap.dm, snap.gas, snap.star, snap.bh        # per-species views

big = halos.select(mass=(1e12, None))                          # chainable selections
parts = snap.dm.select(centre=big["pos"][0], size=2.0)
```

Everything is lazy (no I/O until first field access), metadata is uniform (`ds.meta["redshift"]`, `["boxsize"]`, `["h0"]`). See [docs/unified_interface.md](docs/unified_interface.md) for the full example, including dataset-aware plotting, merger-tree tracks, SHARK galaxies, and the `Simulation` cross-matching API, and [DEVELOPMENT.md](DEVELOPMENT.md) for the design and roadmap.

---

## Status

Active development. The core snapshot tools are stable; halo catalogue, merger tree, profile, gridding, FDM field, and SHARK tooling are evolving alongside ongoing research use. A unified `Dataset`/`Simulation` interface across all data sources is planned — design and roadmap in [DEVELOPMENT.md](DEVELOPMENT.md).

---

## Authors

Chris Power
University of Western Australia

Balu Sreedhar
University of Seville

---

## Notes

This code is designed for research workflows and assumes familiarity with simulation data formats.
