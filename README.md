# AnalysisTools

**AnalysisTools** is a Python toolkit for working with particle-based simulation data, designed for **cosmological simulations and galaxy formation studies**.

It provides simple, flexible tools to:

-   Read and write simulation snapshots (HDF5)
-   Extract and manipulate particle data
-   Interface with halo catalogues and merger trees
-   Perform basic analysis workflows
-   Work with outputs from the **SHARK** semi-analytic model

------------------------------------------------------------------------

## Installation

``` bash
git clone https://github.com/doctorcbpower/AnalysisTools.git
cd AnalysisTools
pip install -e .
```

------------------------------------------------------------------------

## Quick Examples

### 1. Load a Snapshot

``` python
from analysistools import SnapshotTools

snap = SnapshotTools("snapshot_122")
```

------------------------------------------------------------------------

### 2. Select Particles and Write a New Snapshot

``` python
import numpy as np

mask = ...  # e.g. particles within a halo
idx = np.where(mask)[0]
ptype = np.ones(len(idx), dtype=np.int64)

snap.write_snapshot(
    filename="snap_122_subset.hdf5",
    idx=idx,
    idx_type=ptype,
    convention="GADGET4",
    blocks_to_write=["pos", "vel", "pids", "mass"],
)
```

------------------------------------------------------------------------

### 3. Add Metadata to Output Snapshots

``` python
snap.write_snapshot(
    filename="snap_122_halo.hdf5",
    idx=idx,
    idx_type=ptype,
    halo_centre=halo_centre,
    halo_systemic_velocity=halo_velocity,
    halo_extent=halo_extent,
    run_label="ZoomRun",
    periodic=1,
)
```

These values are written to the **HDF5 header**.

------------------------------------------------------------------------

### 4. Working with Halo Catalogues

``` python
from analysistools import HaloCatalogue

cat = HaloCatalogue("halo_catalogue_122")

# Access halo properties
masses = cat["mass"]
centres = cat["centre"]

# Select a halo
halo_id = 10
centre = centres[halo_id]
```

------------------------------------------------------------------------

### 5. Merger Trees

``` python
from analysistools import MergerTree

tree = MergerTree("tree_file")

# Follow a halo through time
history = tree.get_main_branch(halo_id)
```

------------------------------------------------------------------------

### 6. SHARK Semi-Analytic Model

``` python
from analysistools import SharkData

shark = SharkData("shark_output.hdf5")

stellar_mass = shark["stellar_mass"]
sfr = shark["sfr"]
```

------------------------------------------------------------------------

## Snapshot Writing

The main interface is:

``` python
snap.write_snapshot(...)
```

Key features: - Select particles using `idx` - Control particle types via `idx_type` - Choose which datasets to write (`blocks_to_write`) - Attach metadata via keyword arguments

Default blocks:

``` python
["pos", "vel", "pids", "mass"]
```

------------------------------------------------------------------------

## Supported Data Blocks

-   `pos`, `vel`, `pids`, `mass`
-   `u`, `rho`, `hsml`
-   `gas_Z`, `stellar_Z`
-   `sfr`, `age`, `initmass`

------------------------------------------------------------------------

## Status

Active development. The core snapshot tools are stable; additional functionality (trees, SHARK, analysis helpers) is evolving.

------------------------------------------------------------------------

## Authors

Chris Power\
University of Western Australia

Balu Sreedhar\
University of Seville

------------------------------------------------------------------------

## Notes

This code is designed for research workflows and assumes familiarity
with simulation data formats.
