# Snapshots

Covers `SnapshotTools` (read/write/select/split particle snapshots) and `ICTools` (build initial conditions). See also [unified_interface.md](unified_interface.md) for the higher-level `Dataset` API these sit underneath, including `SnapshotDataset`'s `comoving` kwarg (little-*h* handling) and metadata access.

## 1. Load and Read a Snapshot

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

## 2. Supported Snapshot Formats

The `snapfileformat` argument accepts either a string or an integer code:

| Value | Meaning |
|-------|---------|
| `"HDF5"` or `3` | HDF5 (GADGET, AREPO, SWIFT) |
| `"SNAP1"` or `1` | GADGET binary format 1 |
| `"SNAP2"` or `2` | GADGET binary format 2 |

Multi-file snapshots (e.g. `snap_010.0.hdf5`, `snap_010.1.hdf5`, …) are detected automatically.

## 3. Reading Extra Blocks

Additional data blocks can be requested via `extra_blocks` at construction time:

```python
snap = SnapshotTools(
    snapfileformat="HDF5",
    extra_blocks=["AGE", "Z", "SFR", "POT", "GROUPID"],
)
snap.read_snapshot("snap_010.hdf5")
```

## 4. Splitting Particles by Type

Once a snapshot has been read, `LoadParticlesByType` splits the flat arrays into per-species objects:

```python
snap.LoadParticlesByType('all')   # or 'gas', 'dm', 'star', 'bh'

snap.dm.pos      # dark matter positions
snap.gas.pos     # gas positions
snap.gas.density # gas density (if rho was in snapshot)
snap.star.pos    # stellar positions
```

Each species object exposes `pos`, `vel`, `pids`, `mass`, `potential`. Gas additionally exposes `internal_energy` and `density`; all types expose `groupid` when available.

## 5. Selecting Particles by Region

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

## 6. Writing a Subset Snapshot

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

## 7. Building Initial Conditions

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

## Supported Data Blocks

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
