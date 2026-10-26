# Snapshots

Covers `SnapshotTools` (read/split/select particle snapshots) and `ICTools` (build initial conditions). For the higher-level `Dataset` API these sit underneath (`at.load("snap.hdf5")`), including `SnapshotDataset`'s `comoving`/`little_h` kwargs (two independent axes -- scale factor vs. little-*h*, don't conflate them) and metadata access, see [unified_interface.md](unified_interface.md#units-comoving-vs-physical-and-little-h).

## 1. Load and Read a Snapshot

`SnapshotTools` holds configuration (format, convention, particle-type numbers, extra blocks). `read_snapshot()` **returns** a `SnapshotData` object holding the arrays and metadata (and also mirrors them onto the tools object, see below).

```python
from analysistools import SnapshotTools

snap = SnapshotTools(snapfileformat="HDF5", convention="AREPO")
data = snap.read_snapshot("snap_010.hdf5")        # returns SnapshotData
```

The two steps can also be separated (register the file, then read it):

```python
snap.load_snapshot("snap_010.hdf5")      # detects format, single/multi-file; returns None
data = snap.read_snapshot()
```

`convention` selects which header layout the reader expects, and the match is **case-sensitive** in the reader (`snapio_hdf5.py`):

| `convention` | Cosmology read from |
|---|---|
| `"SWIFT"` | `/Header` (`Scale-factor`, `Dimension`, `BoxSize`) and `/Cosmology` |
| `"GADGET4"` or `"AREPO"` (upper case) | `/Header` (`Time`, `BoxSize`) and `/Parameters` (`Omega0`, `OmegaLambda`, `HubbleParam`) |
| anything else, e.g. the default `"GADGET2/3"` | `/Header` (`Time`, `BoxSize`, `Omega0`, `OmegaLambda`, `HubbleParam`) |

So `"Arepo"` (mixed case) silently takes the GADGET2/3 path, which is correct for Arepo files that keep `Omega0` etc. in `/Header` but not for files that only have them in `/Parameters`. Reading with the wrong layout fails with a missing-attribute `KeyError` (e.g. `Omega0` for a SWIFT file). It can also be passed to `read_snapshot(convention=...)`.

`data` then carries the raw, all-species arrays and header values:

```python
data.pos             # positions  (N, 3), all species concatenated in type order
data.vel             # velocities (N, 3)
data.pids            # particle IDs
data.mass            # masses
data.ptype           # particle type of each row
data.num_part_total  # particle count per type
data.box_size, data.scale_factor, data.hubble_param, data.omega_0, ...
```

`read_snapshot` also mirrors these arrays onto the `SnapshotTools` object itself (`snap.pos`, `snap.ptype`, ... -- the same arrays, not copies), because `write_snapshot` and `ICTools` operate on them there. `LoadParticlesByType` and the per-species objects exist only on the returned `SnapshotData`.

## 2. Supported Snapshot Formats

The `snapfileformat` argument accepts either a string or an integer code:

| Value | Meaning |
|-------|---------|
| `"HDF5"` or `3` | HDF5 (GADGET, AREPO, SWIFT) |
| `"SNAP1"` or `1` | GADGET binary format 1 |
| `"SNAP2"` or `2` | GADGET binary format 2 |

Multi-file snapshots (e.g. `snap_010.0.hdf5`, `snap_010.1.hdf5`, …) are detected automatically by `load_snapshot`.

## 3. Reading Extra Blocks

Additional data blocks can be requested via `extra_blocks` at construction time. Recognised keys: `"AGE"`, `"Z"`, `"SFR"`, `"DUST"`, `"STELLARGENS"`, `"INIT_MASS"`, `"POT"`, `"GROUPID"`.

```python
snap = SnapshotTools(
    snapfileformat="HDF5",
    convention="AREPO",
    extra_blocks=["AGE", "Z", "SFR", "POT", "GROUPID"],
)
data = snap.read_snapshot("snap_010.hdf5")
```

## 4. Splitting Particles by Type

`LoadParticlesByType` is a method of the returned `SnapshotData` (not of `SnapshotTools`). It slices the flat arrays into per-species objects, using the type numbers given to the `SnapshotTools` constructor (`gas_type=0`, `dm_type=1`, `star_type=4`, `bh_type=5` by default):

```python
data.LoadParticlesByType('all')   # or 'gas', 'dm', 'star', 'bh'

data.dm.pos      # dark matter positions
data.gas.pos     # gas positions
data.gas.density # gas density (only if the snapshot has a density block)
data.star.pos    # stellar positions
```

Each species object exposes `pos`, `vel`, `pids`, `mass`, `potential`. Gas additionally exposes `internal_energy` and `density` (when the snapshot contains them); all types expose `groupid` when available. A species with no particles in the snapshot is not created (`data.gas` will not exist for a DM-only file), and each call prints the start/end offsets it loaded.

## 5. Selecting Particles by Region

```python
import numpy as np
from analysistools.snapshot_tools import select_particles

centre = np.array([50.0, 50.0, 50.0])

idx = select_particles(
    data.pos, centre,
    size=1.0,
    geometry='spherical',
    periodic=True,
    scale_length=data.box_size,
)
```

`geometry='spherical'` treats `size` as a **radius** and returns an array of **indices**; `geometry='cubic'` treats `size` as the **full side length** and returns a **boolean mask**. Periodic wrapping is applied only when both `periodic=True` and `scale_length` are given.

## 6. Writing a Subset Snapshot

Write from the `SnapshotTools` (or `ICTools`) object that did the reading:

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

`ICTools` provides operations for merging components, tapering a halo, rescaling component masses, and excising particle types:

```python
from analysistools.ic_tools import ICTools

ic = ICTools(snapfileformat="HDF5", convention="AREPO")
ic.read_snapshot("dm_halo.hdf5")

gas = ICTools(snapfileformat="HDF5", convention="AREPO")
gas.read_snapshot("gas_disc.hdf5")

combined = ic.merge([gas])
combined.taper_halo(mass_fraction=0.95)
combined.rescale_component_mass(ptype=combined.dm_type, target_mass=1e10)
n = len(combined.pos)
combined.write_snapshot("combined.hdf5", idx=np.arange(n), idx_type=combined.ptype,
                        convention="GADGET4")
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
    convention     = "GADGET2/3",   # "SWIFT", "GADGET4", "AREPO" (case-sensitive), else GADGET2/3 layout
    positions_only = False,     # read only positions
    hires_only     = False,     # skip low-res particle types
    extra_blocks   = [],        # e.g. ["AGE", "Z", "SFR", "POT", "GROUPID"]
    positions_type = "float32",
    pids_type      = 32,
    not_hires_ptypes = [2, 3, 7],
)
```
