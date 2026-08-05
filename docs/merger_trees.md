# Merger Trees

`MergerTreeTools` dispatches to format-specific readers/walkers and returns format-agnostic `HaloTrack` objects (`merger_tree_types`), so downstream analysis (`find_infall`, `analyse_orbit`) and plotting are format-independent.

## Supported Formats

| `treefileformat` | Tree builder | Reader/walker |
|------------------|--------------|---------------|
| `"SubFind"` | SubFind-HBT (GADGET-4) | `treeio_subfind` |
| `"TreeFrog"` | TreeFrog / VELOCIraptor | `treeio_treefrog` |
| `"MergerTree"` | AHF MergerTree | `treeio_ahf` |

## Usage

```python
from analysistools.merger_tree_tools import MergerTreeTools

mt = MergerTreeTools("trees.hdf5", treefileformat="SubFind", comoving=True, little_h=False)
track = mt.from_halo(ht, index=42, object_type="Subhalo")   # ht = a HaloTools instance
mt.plot_track(track)

host_track = mt.get_track(host_id, snapnum)
mt.plot_relative(track, host_track)

event = track.infall_snapshot()
print(event["snapshot"], event["mass"])
```

## Units: comoving vs. physical, and little-h

Same two independent axes as [`HaloTools`](halo_catalogues.md) (see [unified_interface.md](unified_interface.md#units-comoving-vs-physical-and-little-h) for the full explanation) -- `comoving` (scale factor) and `little_h` (whether Mpc/h-family or Mpc-family), applied once when the tree is read (there's no separate `standardise=True` gate here, unlike `HaloTools`).

Two format-specific wrinkles worth knowing:

- **SubFind-HBT stores `SubhaloPos`/`SubhaloVel` *physical*, not comoving** -- the one exception to the usual "everything is comoving on disk" rule. `comoving=True` (default) correctly converts it *to* comoving (divides by the scale factor, per snapshot); `comoving=False` leaves it physical. This was already handled correctly before this split existed -- it's preserved, not new.
- **TreeFrog's little-h convention is a guess, like `HaloTools`' VELOCIraptor entry** -- and until this was split out, TreeFrog's `comoving` flag actually did `pos * HubbleParam / a` in one expression, conflating both axes (and very likely applying h in the wrong direction for h-included-native data). That's fixed now, but the underlying ambiguity about TreeFrog's *native* little-h convention remains: pass `native_includes_h=` explicitly once you've checked your tree file (or its source snapshot) rather than trusting the default.

```python
mt = MergerTreeTools("VELOCIraptor.tree.hdf5", treefileformat="TreeFrog",
                      comoving=True, little_h=False, native_includes_h=False)
```

`little_h`/`native_includes_h` no-op (with a warning) wherever `HubbleParam` isn't available -- currently TreeFrog walkable trees (topology only, no properties at all) and AHF (no properties of its own either).

See also: [halo_catalogues.md](halo_catalogues.md) for the `HaloTools` instance trees are built from, and [unified_interface.md](unified_interface.md) for the `TrackDataset` API (`tree.from_halo(...)`, `tr.infall()`).
