# AnalysisTools

**AnalysisTools** is a Python toolkit for working with particle-based simulation data, designed for **cosmological simulations and galaxy formation studies**.

It provides tools to:

- Read and write simulation snapshots (HDF5 and GADGET binary)
- Extract and manipulate particle data by type, and build initial conditions (merge, taper, rescale)
- Select particles by geometry, and grid/smooth particle fields onto meshes
- Read and standardise halo catalogues (SubFind, AHF, VELOCIraptor, SWIFT FOF)
- Walk halo merger trees (SubFind-HBT, TreeFrog/VELOCIraptor, AHF MergerTree) and analyse orbits/infall
- Compute radial/vertical density and kinematic profiles
- Read FDM (fuzzy dark matter) wavefunction field snapshots
- Provide lazy, cached, multi-snapshot access to halo catalogues across redshift (`HaloModel`)

---

## Installation

```bash
git clone https://github.com/doctorcbpower/AnalysisTools.git
cd AnalysisTools
pip install -e .
```

The package is imported as `analysistools`; submodules use relative imports (e.g. `halo_tools.py` imports from `.haloio_subfind`), so it must be installed or run as a package rather than as loose scripts.

---

## Module Overview

| Module | Main class(es) / functions | Purpose |
|--------|------------------------------|---------|
| `snapshot_tools.py` | `SnapshotTools`, `select_particles`, `place_points_in_mesh` | Read/write particle snapshots, split by type, geometric selection |
| `snapio_hdf5.py` | `read_hdf5`, `write_hdf5` | Low-level HDF5 snapshot I/O (used by `SnapshotTools`) |
| `snapio_binary.py` | `read_binary` | Low-level GADGET binary format 1/2 snapshot I/O |
| `ic_tools.py` | `ICTools` (extends `SnapshotTools`) | Build initial conditions: merge components, taper halo, rescale masses, excise types |
| `gridding_tools.py` | `GriddingTools` | NGP/CIC mesh assignment, smoothing, Voronoi/nearest-neighbour slices, 3D grid plotting |
| `profile_tools.py` | `ProfileTools`, `MassFunctionTools` | Volume/surface density, velocity dispersion, scale height profiles; halo mass functions |
| `halo_tools.py` | `HaloTools` | Unified reader for halo catalogues (SubFind, AHF, VELOCIraptor, SWIFT FOF) |
| `haloio_subfind.py`, `haloio_ahf.py`, `haloio_velociraptor.py`, `haloio_swiftfof.py` | `read_subfind`, `read_ahf`, `read_velociraptor`, `read_swiftfof` | Format-specific halo catalogue readers (used by `HaloTools`) |
| `halo_tools_standardise_names.py` | `standardise_catalogue_names` | Maps native per-format field names onto a common schema (`mass`, `pos`, `vel`, `radius`, `halo_id`, `num_part`) |
| `halo_model.py` | `HaloModel` | Lazy, cached access to a halo catalogue across multiple snapshots/redshifts, with derived-quantity accessors (`mvir`, `rvir`, `vmax`, `position`, ...) |
| `merger_tree_tools.py` | `MergerTreeTools` | Unified interface to walk merger trees and plot/analyse tracks; dispatches to `treeio_*` readers |
| `treeio_subfind.py`, `treeio_treefrog.py`, `treeio_ahf.py` | `read_subfind_hbt`/`walk_subfind`, `read_treefrog`/`walk_treefrog`, `read_ahf_mergertree`/`walk_ahf` | Format-specific merger tree readers and main-branch walkers |
| `merger_tree_types.py` | `HaloTrack`, `OrbitAnalysis`, `MergerTreeError`, `periodic_delta` | Shared, format-agnostic containers used by `MergerTreeTools` |
| `fdm_field_tools.py` | `FDMFieldTools` | Read FDM wavefunction field snapshots (mesh data, not particles); density, radial profile, slices |

---

## Quick Examples

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

## Working with Halo Catalogues

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

## Merger Trees

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

## Profiles and Mass Functions

```python
from analysistools.profile_tools import ProfileTools

pt = ProfileTools(numbins=40)

rho = pt.volume_density(pos, mass, centre, rmin=0.1, rmax=300)
sigma = pt.surface_density(pos, mass, centre, rmin=0.1, rmax=50)
hz = pt.scale_height(pos, mass, centre, rmin=0.1, rmax=50)

pt.plot(rho, "density", ylabel=r"$\rho(r)$")
```

---

## Gridding and Mesh Tools

```python
from analysistools.gridding_tools import GriddingTools

gt = GriddingTools()
grid = gt.smooth_to_grid(positions, values, grid_size, grid_limits)
gt.plot_3d_slice(grid, grid_limits)
```

---

## FDM (Fuzzy Dark Matter) Field Snapshots

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

## Supported Data Blocks (Snapshots)

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

## Constructor Options (`SnapshotTools`)

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

## Constructor Options (`HaloTools`)

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

Active development. The core snapshot tools are stable; halo catalogue, merger tree, profile, gridding, and FDM field tooling are evolving alongside ongoing research use.

---

## Authors

Chris Power
University of Western Australia

Balu Sreedhar
University of Seville

---

## Notes

This code is designed for research workflows and assumes familiarity with simulation data formats.
