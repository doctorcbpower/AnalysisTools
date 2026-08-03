# AnalysisTools

**AnalysisTools** is a Python toolkit for working with particle-based simulation data, designed for **cosmological simulations and galaxy formation studies**.

It provides tools to:

- Read and write simulation snapshots (HDF5 and GADGET binary)
- Extract and manipulate particle data by type, and build initial conditions (merge, taper, rescale)
- Select particles by geometry, and grid/smooth particle fields onto meshes
- Render publication-quality Voronoi-mesh projections of Arepo output via [vortrace](https://github.com/gusbeane/vortrace) (optional dependency)
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
| `rendering` | `vortrace` (from GitHub; not on PyPI) | `render_tools.VortraceRenderer` |

```bash
pip install -e ".[rendering]"
```

`render_tools` itself imports fine without `vortrace` installed -- only instantiating `VortraceRenderer` requires it.

---

## Module Overview

| Module | Main class(es) / functions | Purpose |
|--------|------------------------------|---------|
| `snapshot_tools.py` | `SnapshotTools`, `select_particles`, `place_points_in_mesh` | Read/write particle snapshots, split by type, geometric selection |
| `snapio_hdf5.py` | `read_hdf5`, `write_hdf5` | Low-level HDF5 snapshot I/O (used by `SnapshotTools`) |
| `snapio_binary.py` | `read_binary` | Low-level GADGET binary format 1/2 snapshot I/O |
| `ic_tools.py` | `ICTools` (extends `SnapshotTools`) | Build initial conditions: merge components, taper halo, rescale masses, excise types |
| `gridding_tools.py` | `GriddingTools` | NGP/CIC mesh assignment, smoothing, Voronoi/nearest-neighbour slices, 3D grid plotting |
| `render_tools.py` | `VortraceRenderer` | Exact Voronoi-mesh projections/volume rendering via the optional `vortrace` package (publication-quality alternative to `GriddingTools`' grid deposition) |
| `profile_tools.py` | `ProfileTools`, `MassFunctionTools` | Volume/surface density, velocity dispersion, scale height profiles; halo mass functions |
| `halo_tools.py` | `HaloTools` | Unified reader for halo catalogues (SubFind, AHF, VELOCIraptor, SWIFT FOF) |
| `haloio_subfind.py`, `haloio_ahf.py`, `haloio_velociraptor.py`, `haloio_swiftfof.py` | `read_subfind`, `read_ahf`, `read_velociraptor`, `read_swiftfof` | Format-specific halo catalogue readers (used by `HaloTools`) |
| `halo_tools_standardise_names.py` | `standardise_catalogue_names` | Maps native per-format field names onto a common schema (`mass`, `pos`, `vel`, `radius`, `halo_id`, `num_part`) |
| `halo_model.py` | `HaloModel` | Lazy, cached access to a halo catalogue across multiple snapshots/redshifts, with derived-quantity accessors (`mvir`, `rvir`, `vmax`, `position`, ...) |
| `merger_tree_tools.py` | `MergerTreeTools` | Unified interface to walk merger trees and plot/analyse tracks; dispatches to `treeio_*` readers |
| `treeio_subfind.py`, `treeio_treefrog.py`, `treeio_ahf.py` | `read_subfind_hbt`/`walk_subfind`, `read_treefrog`/`walk_treefrog`, `read_ahf_mergertree`/`walk_ahf` | Format-specific merger tree readers and main-branch walkers |
| `merger_tree_types.py` | `HaloTrack`, `OrbitAnalysis`, `MergerTreeError`, `periodic_delta` | Shared, format-agnostic containers used by `MergerTreeTools` |
| `fdm_field_tools.py` | `FDMFieldTools` | Read FDM wavefunction field snapshots (mesh data, not particles); density, radial profile, slices |
| `galaxy_tools.py` | `GalaxyTools` | Legacy SHARK catalogue reader — superseded by `shark.model.SharkModel` (see below) |

---

## Unified Interface (new)

`analysistools.load()` opens any supported data product with one call — kind and format sniffed from the file — and returns a `Dataset` with common syntax across sources (snapshots, halo catalogues, merger trees, and SHARK galaxy catalogues; see DEVELOPMENT.md):

```python
import numpy as np
import analysistools as at

snap  = at.load("snap_0031.hdf5")                      # SWIFT/GADGET sniffed
halos = at.load("snap_0031.VELOCIraptor.properties.0") # VELOCIraptor sniffed

# same field vocabulary everywhere: pos, vel, mass, id, radius, ...
halos["mass"], halos["pos"], snap["pos"]

# per-species views (snapshots)
snap.dm["pos"]          # or snap.dm.pos — both work
snap.gas, snap.star, snap.bh

# chainable selections: geometric and attribute cuts
big = halos.select(mass=(1e12, None))
c   = big["pos"][np.argmax(big["mass"])]
parts = snap.dm.select(centre=c, size=2.0 * big["radius"].max())

# dataset-aware plotting (GriddingTools/ProfileTools underneath)
from analysistools.api import plotting as plt2
ax = plt2.density_map(snap.dm, centre=c, size=5.0, method="CIC")
plt2.overlay_points(ax, big, centre=c)            # halos on the map
prof = plt2.profile(snap.dm, c, kind="density", rmin=0.05, rmax=1.0)
plt2.plot_profile(prof)
plt2.mass_function([halos])                       # halo or galaxy Datasets

# merger trees: queries return TrackDataset (rows = snapshots)
tree = at.load("VELOCIraptor.walkabletree.hdf5", halos=halos)
tr = tree.from_halo(halos, index=0)
tr["mass"], tr["redshift"], tr.infall()

# SHARK galaxies: one file, or a SharkModel frozen at a redshift
gals = at.load("199/0/galaxies.hdf5")             # kind sniffed
gals["pos"], gals["mass"], gals.select(type=(0, 1))   # centrals
z0 = shark_model.at(0.0)                          # EpochModel -> Dataset
z0h = halo_model.at(0.0)                          # works for HaloModel too

# Simulation: bind matched products, cross-match between them
sim = at.Simulation(
    snapshots={0.0: "snap_0031.hdf5"},
    halos={0.0: "snap_0031.VELOCIraptor.properties.0"},
    trees="VELOCIraptor.walkabletree.hdf5",
    galaxies={0.0: "199/0/galaxies.hdf5"},
    snapnums={0.0: 31}, label="CDM",
)
epoch = sim.at(0.0)                               # matched views, lazy
parts = epoch.particles_in_halo(index=0, r_scale=2.0)
gals  = epoch.galaxies_in_halo(index=0)           # position match (default)
gals2 = epoch.galaxies_in_halo(hid, match_by="id")
track = epoch.track_of(index=0)                   # catalogue auto-linked
```

Everything is lazy (no I/O until first field access), metadata is uniform (`ds.meta["redshift"]`, `["boxsize"]`, `["h0"]`), and the underlying legacy object stays reachable via `ds.backend`. The existing APIs below are unchanged.

---

## SHARK Semi-Analytic Catalogues (`shark` package)

The `shark` package now lives at `analysistools.shark` (the old top-level `import shark` still works via a deprecation shim). It provides lazy, cached access to SHARK galaxy catalogues, mirroring the `HaloModel` API so both can feed the same analysis and plotting code:

```python
from analysistools.shark.common import _redshift_table, parse_subvolumes
from analysistools.shark.model import SharkModel

rt = _redshift_table("redshift_list.txt")
sv = parse_subvolumes("0-63")

m = SharkModel("./CDM/base_model", rt, sv, label="CDM")

mstar = m.mstars(redshift=0.0)          # mstars_disk + mstars_bulge
sfr   = m.sfr(redshift=0.0)
raw   = m.get("mvir_hosthalo", redshift=0.0)   # any native field by name
```

Available accessors include `mstars`, `mgas`, `sfr`, `ssfr`, `rstar`, `mhalo`, `msubhalo`, `mbh`, `mbulge`, `galaxy_type`, `h0`, `volume`, `age_at_z`, plus star formation history access (`sfh_disk`, `sfh_bulge`, `Z_disk_history`, ...).

Higher-level components:

| Module | Purpose |
|--------|---------|
| `shark.analysis` | `Analysis`: halo/stellar mass functions, SFR main sequence, size–mass, BH–bulge relations |
| `shark.plots` | `Plotter`: standard comparison plots over one or more `SharkModel`s |
| `shark.cli` | Command-line driver |
| `shark.photometry` | fsps-based luminosities/magnitudes from metallicities and star formation histories |

---

## Merger Tree Formats

| `treefileformat` | Tree builder | Reader/walker |
|------------------|--------------|---------------|
| `"SubFind"` | SubFind-HBT (GADGET-4) | `treeio_subfind` |
| `"TreeFrog"` | TreeFrog / VELOCIraptor | `treeio_treefrog` |
| `"MergerTree"` | AHF MergerTree | `treeio_ahf` |

All return format-agnostic `HaloTrack` objects (`merger_tree_types`), so downstream analysis (`find_infall`, `analyse_orbit`) and plotting are format-independent.

---

## Appendix: Legacy Per-Class API

Everything below documents the original per-class API (`SnapshotTools`, `HaloTools`, `MergerTreeTools`, ...). It remains fully supported -- the unified interface above is a layer on top of these classes, and `ds.backend` always exposes the underlying object. Prefer the unified interface for new work.

### Quick Examples

### 1. Load and Read a Snapshot

`SnapshotTools` is configured at construction time, then used to load and read snapshots. The two steps can be combined or kept separate:

```python
from analysistools import SnapshotTools

# Step 1: create the tools object, specifying the file format
snap = SnapshotTools(snapfileformat="HDF5")

# Step 2: load (register the file) and read (pull data into memory)
snap.load_snapshot("snap_010.hdf5")
snap.read_snapshot()

# Or combine both steps in one call:
snap.read_snapshot("snap_010.hdf5")
```

After reading, particle data is available directly as attributes:

```python
snap.pos     # positions  (N, 3)
snap.vel     # velocities (N, 3)
snap.pids    # particle IDs
snap.mass    # masses
```

---

### 2. Supported Snapshot Formats

The `snapfileformat` argument accepts either a string or an integer code:

| Value | Meaning |
|-------|---------|
| `"HDF5"` or `3` | HDF5 (GADGET, AREPO, SWIFT) |
| `"SNAP1"` or `1` | GADGET binary format 1 |
| `"SNAP2"` or `2` | GADGET binary format 2 |

Multi-file snapshots (e.g. `snap_010.0.hdf5`, `snap_010.1.hdf5`, …) are detected automatically.

---

### 3. Reading Extra Blocks

Additional data blocks can be requested via `extra_blocks` at construction time:

```python
snap = SnapshotTools(
    snapfileformat="HDF5",
    extra_blocks=["AGE", "Z", "SFR", "POT", "GROUPID"],
)
snap.read_snapshot("snap_010.hdf5")
```

---

### 4. Splitting Particles by Type

Once a snapshot has been read, `LoadParticlesByType` splits the flat arrays into per-species objects:

```python
snap.LoadParticlesByType('all')   # or 'gas', 'dm', 'star', 'bh'

snap.dm.pos      # dark matter positions
snap.gas.pos     # gas positions
snap.gas.density # gas density (if rho was in snapshot)
snap.star.pos    # stellar positions
```

Each species object exposes `pos`, `vel`, `pids`, `mass`, `potential`. Gas additionally exposes `internal_energy` and `density`; all types expose `groupid` when available.

---

### 5. Selecting Particles by Region

```python
import numpy as np
from analysistools.snapshot_tools import select_particles

centre = np.array([50.0, 50.0, 50.0])

idx = select_particles(
    snap.pos, centre,
    size=1.0,
    geometry='spherical',
    periodic=True,
    scale_length=snap.box_size,
)
```

---

### 6. Writing a Subset Snapshot

```python
ptype = np.ones(len(idx), dtype=np.int64)

snap.write_snapshot(
    filename="snap_010_subset",
    idx=idx,
    idx_type=ptype,
    convention="GADGET4",                # or "SWIFT", "AREPO", "GADGET2/3"
    blocks_to_write=["pos", "vel", "pids", "mass"],
)
```

Default blocks written when `blocks_to_write` is omitted: `pos`, `vel`, `pids`, `mass`, `groupid`.

---

### 7. Building Initial Conditions

`ICTools` extends `SnapshotTools` with operations for merging components, tapering a halo, rescaling component masses, and excising particle types:

```python
from analysistools.ic_tools import ICTools

ic = ICTools(snapfileformat="HDF5")
ic.read_snapshot("dm_halo.hdf5")

gas = ICTools(snapfileformat="HDF5")
gas.read_snapshot("gas_disc.hdf5")

combined = ic.merge([gas])
combined.taper_halo(mass_fraction=0.95)
combined.rescale_component_mass(ptype=combined.dm_type, target_mass=1e10)
combined.write_snapshot("combined.hdf5", idx=..., idx_type=...)
```

---

### Working with Halo Catalogues

```python
from analysistools.halo_tools import HaloTools

ht = HaloTools(comoving_units=True)
halos = ht.read_catalogue(
    filename="groups_010.hdf5",
    fileformat="SubFind",       # or "AHF", "VELOCIraptor", "SWIFT_FOF"
    standardise=True,
    snapnum=10,
)
ht.summary()
```

`standardise=True` runs `halo_tools_standardise_names.standardise_catalogue_names`, mapping native per-format fields onto a common schema: `mass`, `pos`, `vel`, `radius`, `halo_id`, `num_part`.

### Multi-snapshot access with `HaloModel`

`HaloModel` gives lazy, cached access to a halo catalogue across many snapshots/redshifts, with derived-quantity accessors:

```python
from analysistools.halo_model import HaloModel

catalogues = {
    0.0: "./catalogues/snap_200/fof_output.hdf5",
    1.0: "./catalogues/snap_151/fof_output.hdf5",
}

hm = HaloModel(
    catalogues=catalogues,
    fileformat="SWIFT_FOF",
    label="CDM",
    comoving=True,
    standardise=True,
)

mvir = hm.mvir(redshift=0.0)     # file is only read on first request
pos  = hm.position(redshift=0.0)
```

---

### Merger Trees

`MergerTreeTools` dispatches to format-specific readers/walkers (SubFind-HBT, TreeFrog/VELOCIraptor, AHF MergerTree) and returns format-agnostic `HaloTrack` objects for downstream analysis and plotting:

```python
from analysistools.merger_tree_tools import MergerTreeTools

mt = MergerTreeTools("trees.hdf5", treefileformat="SubFind")
track = mt.from_halo(ht, index=42, object_type="Subhalo")   # ht = a HaloTools instance
mt.plot_track(track)

host_track = mt.get_track(host_id, snapnum)
mt.plot_relative(track, host_track)

event = track.infall_snapshot()
print(event["snapshot"], event["mass"])
```

---

### Profiles and Mass Functions

```python
from analysistools.profile_tools import ProfileTools

pt = ProfileTools(numbins=40)

rho = pt.volume_density(pos, mass, centre, rmin=0.1, rmax=300)
sigma = pt.surface_density(pos, mass, centre, rmin=0.1, rmax=50)
hz = pt.scale_height(pos, mass, centre, rmin=0.1, rmax=50)

pt.plot(rho, "density", ylabel=r"$\rho(r)$")
```

---

### Gridding and Mesh Tools

```python
from analysistools.gridding_tools import GriddingTools

gt = GriddingTools()
grid = gt.smooth_to_grid(positions, values, grid_size, grid_limits)
gt.plot_3d_slice(grid, grid_limits)
```

---

### Rendering with vortrace

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

---

### FDM (Fuzzy Dark Matter) Field Snapshots

`FDMFieldTools` is a mesh-field analogue to `SnapshotTools`, for reading FDM wavefunction field snapshots (not particle data):

```python
from analysistools.fdm_field_tools import FDMFieldTools

field = FDMFieldTools("fdm_field_010.hdf5")

field.psi       # complex128 array, shape (N, N, N)
field.density   # |psi|^2, real array, shape (N, N, N)
field.N         # grid size
field.L         # box size (code length units)
field.mass      # mc^2, eV
```

---

### Supported Data Blocks (Snapshots)

| Block key | Description |
|-----------|-------------|
| `pos` | Particle positions |
| `vel` | Velocities |
| `pids` | Particle IDs |
| `mass` | Masses |
| `u` | Gas internal energy |
| `rho` | Gas density |
| `hsml` | Smoothing lengths |
| `gas_Z` | Gas metallicity |
| `stellar_Z` | Stellar metallicity |
| `sfr` | Star formation rate |
| `age` | Stellar ages |
| `initmass` | Stellar initial masses |
| `groupid` | Halo/group membership |

---

### Constructor Options (`SnapshotTools`)

```python
SnapshotTools(
    snapfileformat = "HDF5",    # "HDF5", "SNAP1", "SNAP2" (or 1, 2, 3)
    gas_type       = 0,         # GADGET particle type for gas
    dm_type        = 1,         # GADGET particle type for dark matter
    star_type      = 4,         # GADGET particle type for stars
    bh_type        = 5,         # GADGET particle type for black holes
    # keyword options:
    convention     = "GADGET2/3",
    positions_only = False,     # read only positions
    hires_only     = False,     # skip low-res particle types
    extra_blocks   = [],        # e.g. ["AGE", "Z", "SFR", "POT", "GROUPID"]
    positions_type = "float32",
    pids_type      = 32,
    not_hires_ptypes = [2, 3, 7],
)
```

### Constructor Options (`HaloTools`)

```python
HaloTools(
    comoving_units = False,
    usehalocatonly = False,
    usesubstructure_file = False,
    loglevel = logging.INFO,
)
```

Supported `fileformat` values for `read_catalogue`: `"SUBFIND"`, `"AHF"`, `"VELOCIraptor"`, `"SWIFT_FOF"` (or integer codes 1-4).

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
