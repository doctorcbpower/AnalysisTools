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

mt = MergerTreeTools("trees.hdf5", treefileformat="SubFind")
track = mt.from_halo(ht, index=42, object_type="Subhalo")   # ht = a HaloTools instance
mt.plot_track(track)

host_track = mt.get_track(host_id, snapnum)
mt.plot_relative(track, host_track)

event = track.infall_snapshot()
print(event["snapshot"], event["mass"])
```

See also: [halo_catalogues.md](halo_catalogues.md) for the `HaloTools` instance trees are built from, and [unified_interface.md](unified_interface.md) for the `TrackDataset` API (`tree.from_halo(...)`, `tr.infall()`).
