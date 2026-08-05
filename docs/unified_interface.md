# Unified Interface

`analysistools.load()` opens any supported data product with one call — kind and format sniffed from the file — and returns a `Dataset` with common syntax across sources (snapshots, halo catalogues, merger trees, and SHARK galaxy catalogues; see [DEVELOPMENT.md](../DEVELOPMENT.md) for the design and roadmap).

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

Everything is lazy (no I/O until first field access), metadata is uniform (`ds.meta["redshift"]`, `["boxsize"]`, `["h0"]`), and the underlying legacy object stays reachable via `ds.backend`. The per-class API documented in the other `docs/` pages is unchanged and fully supported — the unified interface is a layer on top of it.

## Units: little-h

Both `SnapshotDataset` and `HaloCatalogue` take a `comoving` kwarg (default `True` on both) that controls whether the little-*h* factor is stripped from length/mass fields:

```python
snap  = at.load("snap_0031.hdf5", comoving=True)                       # pos, mass, boxsize / h0
halos = at.load("halos.VELOCIraptor.properties.0", fileformat="VELOCIraptor", comoving=True)  # pos / h0
```

With the defaults, both come back with the little-*h* factor already divided out, so `snap["pos"]` and `halos["pos"]` are directly comparable. Pass `comoving=False` on either (or both) to get values exactly as stored in the file instead (typically h⁻¹ code units for GADGET/Arepo/SWIFT snapshots and h⁻¹ Mpc for halo catalogues) -- set both adapters the same way rather than mixing, or particle/halo positions will be off by a factor of `h0`. Check `ds.meta["comoving"]` (snapshots) if you need to confirm which state a given `Dataset` ended up in.

This only affects `pos`/`mass`/`boxsize`; it does not touch the separate comoving-vs-physical (scale factor) axis -- see `SnapshotTools.UnitConversion`'s `convert_to_physical`/`convert_to_comoving` in [snapshots.md](snapshots.md) for that.

## Accessing metadata / header info

Every `Dataset` (`SnapshotDataset`, `HaloCatalogue`, `MergerTree`/`TrackDataset`, `GalaxyCatalogue`, ...) exposes `ds.meta` -- a plain dict, populated on first field access (or call `ds.preload()` to force it without touching any column). Common keys: `redshift`, `boxsize`, `h0`, `scale_factor`, `units`. `ds.summary()` prints those plus row count and available fields.

```python
snap = at.load("snap_0031.hdf5").preload()
snap.meta["redshift"], snap.meta["h0"], snap.meta["boxsize"]
snap.summary()
```

Both `HaloCatalogue` and `SnapshotDataset` also keep a raw header/config dict at `ds.meta["native_meta"]`, for fields that don't have a standardised `meta` key -- e.g. `omega_dm`/`omega_bar` (SWIFT-only split of `omega_0`), `mass_table`, `dimension`, `num_files`, `ispotential`/`isgroupid`. For `SnapshotDataset` this is always the raw, un-h-stripped value regardless of the `comoving` setting -- only the standardised top-level keys (`pos`, `mass`, `boxsize`) get the little-*h* conversion.

```python
snap.meta["native_meta"]["omega_dm"], snap.meta["native_meta"]["mass_table"]
```

See also: [snapshots.md](snapshots.md), [halo_catalogues.md](halo_catalogues.md), [merger_trees.md](merger_trees.md), [shark.md](shark.md).
