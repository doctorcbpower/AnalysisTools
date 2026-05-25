# AnalysisTools

**AnalysisTools** is a Python toolkit for working with particle-based simulation data, designed for **cosmological simulations and galaxy formation studies**.

It provides simple, flexible tools to:

- Read and write simulation snapshots (HDF5 and GADGET binary)
- Extract and manipulate particle data by type
- Select particles by geometry
- Interface with halo catalogues and merger trees
- Perform basic analysis workflows
- Work with outputs from the **SHARK** semi-analytic model

---

## Installation

```bash
git clone https://github.com/doctorcbpower/AnalysisTools.git
cd AnalysisTools
pip install -e .
```

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

# Now available (where present in the file):
snap.stellarage        # stellar ages
snap.gas_metallicity   # gas metallicity
snap.stellar_metallicity
snap.gas_sfr           # star formation rate
snap.potential         # gravitational potential
snap.groupid           # group/halo membership
```

---

### 4. Splitting Particles by Type

Once a snapshot has been read, `LoadParticlesByType` splits the flat arrays into per-species objects:

```python
snap.LoadParticlesByType('all')   # or 'gas', 'dm', 'star', 'bh'

# Access by species:
snap.dm.pos      # dark matter positions
snap.gas.pos     # gas positions
snap.gas.density # gas density (if rho was in snapshot)
snap.star.pos    # stellar positions
```

Each species object exposes: `pos`, `vel`, `pids`, `mass`, `potential`.  
Gas additionally exposes `internal_energy` and `density`; all types expose `groupid` when available.

---

### 5. Selecting Particles by Region

```python
import numpy as np
from analysistools import select_particles

centre = np.array([50.0, 50.0, 50.0])

# Spherical selection (returns indices)
idx = select_particles(
    snap.pos, centre,
    size=1.0,                     # radius in snapshot length units
    geometry='spherical',
    periodic=True,
    scale_length=snap.box_size,
)

# Cubic selection
idx = select_particles(
    snap.pos, centre,
    size=2.0,                     # half-width per side
    geometry='cubic',
)
```

---

### 6. Writing a Subset Snapshot

Use `write_snapshot` to write a selection of particles to a new HDF5 file.  
`idx` is an integer array of particle indices; `idx_type` assigns each particle to a GADGET particle-type (0–5).

```python
import numpy as np

# Example: write all selected particles as dark matter (type 1)
ptype = np.ones(len(idx), dtype=np.int64)

snap.write_snapshot(
    filename="snap_010_subset",          # .hdf5 is appended automatically
    idx=idx,
    idx_type=ptype,
    convention="GADGET4",                # or "SWIFT", "AREPO", "GADGET2/3"
    blocks_to_write=["pos", "vel", "pids", "mass"],
)
```

Default blocks written when `blocks_to_write` is omitted: `pos`, `vel`, `pids`, `mass`, `groupid`.

---

### 7. Adding Metadata to Output Snapshots

Any keyword argument passed to `write_snapshot` is stored as an HDF5 header attribute. Recognised metadata keywords:

```python
snap.write_snapshot(
    filename="snap_010_halo",
    idx=idx,
    idx_type=ptype,
    convention="GADGET4",
    halo_centre=np.array([50.0, 50.0, 50.0]),
    halo_systemic_velocity=np.array([0.0, 0.0, 0.0]),
    halo_extent=1.0,
    run_label="ZoomRun",
    periodic=1,
)
```

---

### 8. Output Conventions

The `convention` argument controls the HDF5 header layout and dataset names written:

| Convention | Notes |
|-----------|-------|
| `"SWIFT"` | Writes `Scale-factor`, `BoxSize` as 3-vector; adds `Cosmology` group |
| `"GADGET4"` / `"AREPO"` | Writes `Time`, `Redshift`, cosmological params; adds `Parameters` group |
| `"GADGET2/3"` | Legacy header layout (default for reading) |

---

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

---

## Constructor Options

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

---

## Working with Halo Catalogues

```python
from analysistools import HaloCatalogue

cat = HaloCatalogue("halo_catalogue_122")

masses  = cat["mass"]
centres = cat["centre"]

halo_id = 10
centre  = centres[halo_id]
```

---

## Merger Trees

```python
from analysistools import MergerTree

tree    = MergerTree("tree_file")
history = tree.get_main_branch(halo_id)
```

---

## SHARK Semi-Analytic Model

```python
from analysistools import SharkData

shark        = SharkData("shark_output.hdf5")
stellar_mass = shark["stellar_mass"]
sfr          = shark["sfr"]
```

---

## Status

Active development. The core snapshot tools are stable; additional functionality (trees, SHARK, analysis helpers) is evolving.

---

## Authors

Chris Power  
University of Western Australia

Balu Sreedhar  
University of Seville

---

## Notes

This code is designed for research workflows and assumes familiarity with simulation data formats.
