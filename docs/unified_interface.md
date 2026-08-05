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

## Units: comoving vs. physical, and little-h

These are two **independent** axes -- don't conflate them, and don't assume `comoving` controls h (an earlier version of this doc, and of the code, did exactly that). Both `SnapshotDataset` and `HaloCatalogue` take both kwargs, with the same meaning on each:

- **`comoving`** (default `True` on both): scale-factor axis. Catalogues/snapshots are always stored comoving on disk, so `comoving=True` is a no-op; `comoving=False` multiplies `pos`/`boxsize` by the scale factor `a` to convert to physical coordinates. Never touches `mass`.
- **`little_h`** (default `False` on both): whether `pos`/`mass`/`boxsize` are in little-*h* units (Mpc/h-family, 1e10 Msol/h) or h-free (Mpc-family, Msol). Whether `little_h=False` (the default) actually *does* anything depends on the **native convention of the specific format/code the data came from** -- SWIFT already stores h-free values; GADGET/Arepo-family codes (SubFind, and conventionally AHF/VELOCIraptor) store h-scaled values. Mixing a SWIFT-origin snapshot against a SubFind-origin catalogue (or vice versa) without accounting for this is exactly the kind of mismatch that motivated this split.

```python
snap  = at.load("snap_0031.hdf5")                                              # comoving=True, little_h=False
halos = at.load("halos.VELOCIraptor.properties.0", fileformat="VELOCIraptor",
                 native_includes_h=False)                                      # see caveat below
```

With the defaults, `snap["pos"]` and `halos["pos"]` come back in the same (comoving, h-free) units, directly comparable -- *provided* the little-h native-convention guess is correct for both. Check `ds.meta["comoving"]`/`ds.meta["little_h"]` to confirm which state a given `Dataset` actually ended up in (it may differ from what was requested if `HubbleParam`/scale factor weren't available -- a warning is logged when that happens).

**Caveat for AHF/VELOCIraptor**: unlike SubFind (which *is* GADGET/Arepo's own group finder, so always inherits its h-scaled convention), AHF and VELOCIraptor are standalone halo finders that inherit whatever convention their *input snapshot* used -- there's no fixed answer for "AHF's native units" or "VELOCIraptor's native units" the way there is for SubFind or SWIFT_FOF. Confirmed with this repo's own bundled test data: `data/halos/snap_0031.VELOCIraptor.properties.0` was run against the SWIFT snapshot `data/snap_0031.hdf5`, and its `Period` is SWIFT's raw h-free `BoxSize` passed straight through -- not converted to a Mpc/h convention. `HaloTools`/`HaloCatalogue`/`HaloModel`'s default guess (`NATIVE_INCLUDES_LITTLE_H`, `True` for both formats) would silently double-strip h for a catalogue like this. Pass `native_includes_h=False` (or `True`) explicitly for AHF/VELOCIraptor once you've checked which convention your specific catalogue actually uses (e.g. compare its `BoxSize`/`Period` against its source snapshot's, as done above) -- don't trust the default for these two formats.

## Accessing metadata / header info

Every `Dataset` (`SnapshotDataset`, `HaloCatalogue`, `MergerTree`/`TrackDataset`, `GalaxyCatalogue`, ...) exposes `ds.meta` -- a plain dict, populated on first field access (or call `ds.preload()` to force it without touching any column). Common keys: `redshift`, `boxsize`, `h0`, `scale_factor`, `units`. `ds.summary()` prints those plus row count and available fields.

```python
snap = at.load("snap_0031.hdf5").preload()
snap.meta["redshift"], snap.meta["h0"], snap.meta["boxsize"]
snap.summary()
```

Both `HaloCatalogue` and `SnapshotDataset` also keep a raw header/config dict at `ds.meta["native_meta"]`, for fields that don't have a standardised `meta` key -- e.g. `omega_dm`/`omega_bar` (SWIFT-only split of `omega_0`), `mass_table`, `dimension`, `num_files`, `ispotential`/`isgroupid`. This is always the raw, unconverted value regardless of `comoving`/`little_h` -- only the standardised top-level keys (`pos`, `mass`, `boxsize`) get those conversions.

```python
snap.meta["native_meta"]["omega_dm"], snap.meta["native_meta"]["mass_table"]
```

See also: [snapshots.md](snapshots.md), [halo_catalogues.md](halo_catalogues.md), [merger_trees.md](merger_trees.md), [shark.md](shark.md).
