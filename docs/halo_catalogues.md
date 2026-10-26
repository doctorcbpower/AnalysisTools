# Halo Catalogues

Covers `HaloTools` (single-catalogue reader) and `HaloModel` (lazy, cached, multi-snapshot access).

## Working with Halo Catalogues

```python
from analysistools.halo_tools import HaloTools

ht = HaloTools(comoving=True, little_h=False)
halos = ht.read_catalogue(
    filename="groups_010.hdf5",
    fileformat="SubFind",       # or "AHF", "VELOCIraptor", "SWIFT_FOF"
    standardise=True,
    snapnum=10,
)
ht.summary()
```

`standardise=True` runs `halo_tools_standardise_names.standardise_catalogue_names`, mapping native per-format fields onto a common schema (`mass`, `pos`, `vel`, `radius`, `halo_id`, `num_part`), then applies the `comoving`/`little_h` unit conversions -- both are no-ops with `standardise=False`, raw fields come back exactly as stored in the file. See [unified_interface.md](unified_interface.md#units-comoving-vs-physical-and-little-h) for what `comoving`/`little_h` actually mean (two independent axes -- don't conflate them) and, importantly, the `native_includes_h` override you should set explicitly for AHF/VELOCIraptor catalogues rather than trusting the per-format default guess.

## Centring SUBFIND groups on their primary subhalo

By default, a SUBFIND group's `pos`/`vel` are `GroupPos`/`GroupVel` -- the FOF group's own centre, which substructure can offset from where the halo actually is. Pass `centre_on_subhalo=True` to use the group's primary subhalo instead (`GroupFirstSub` -> `SubhaloPos`/`SubhaloVel`), which is usually the better proxy for the halo centre:

```python
ht = HaloTools(comoving=True, little_h=False, centre_on_subhalo=True)
halos = ht.read_catalogue(
    filename="groups_010.hdf5", fileformat="SUBFIND", standardise=True,
)

halos["pos"]                 # primary-subhalo position where available
halos["centred_on_subhalo"]  # bool per group: False where GroupFirstSub was
                              # invalid (no bound subhalo) and GroupPos/GroupVel
                              # were kept instead
```

Only applies to SUBFIND catalogues with a Subhalo table; other formats ignore the flag. `vel` is substituted from `SubhaloVel` too when present, but that's SubFind's own subhalo bulk velocity -- not a self-consistent recomputation from bound particles -- so treat it as an approximation rather than a definitive centre-of-mass velocity.

Works the same way through the unified interface: `HaloCatalogue(path, fileformat="SUBFIND", centre_on_subhalo=True)` or `at.load(path, fileformat="SUBFIND", centre_on_subhalo=True)` (see [unified_interface.md](unified_interface.md)).

## Multi-snapshot access with `HaloModel`

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
    little_h=False,
    standardise=True,
)

mvir = hm.mvir(redshift=0.0)     # file is only read on first request
pos  = hm.position(redshift=0.0)
vol  = hm.volume(redshift=0.0)   # respects little_h/native_includes_h too
```

## AHF catalogues

`read_ahf` (in `haloio_ahf.py`) reads the `.AHF_halos` table. Columns are selected by the names in the `#` header line (`ID(1) hostHalo(2) ...`), not by position, so a different AHF version raises an error that prints the header instead of shifting columns silently. Every row of the file is returned as the Group table; the rows with `hostHalo > 0` (inside another halo) are also returned as the Subhalo table. A `hostHalo` of 0 or -1 means the halo is not inside another halo. Values are as stored in the file, and the file does not give `HubbleParam`, so pass `native_includes_h` explicitly (see above).

Two further readers are in the same module and are not used by `HaloTools`:

```python
from analysistools.haloio_ahf import read_ahf_substructure, read_ahf_particles, AHF_PTYPE_STAR

subs = read_ahf_substructure("run.AHF_substructure")      # {host ID: array of subhalo IDs}
blocks = read_ahf_particles("run.AHF_particles", wanted=[host_id])
ids, types = blocks[host_id]["ids"], blocks[host_id]["types"]
star_ids = ids[types == AHF_PTYPE_STAR]                  # GADGET types: 0 gas, 1 DM, 4 star
```

`read_ahf_substructure` expects, for each host, a line with the host ID and the number of subhalos followed by one line listing their IDs. It raises an error naming the line if the file differs. `read_ahf_particles` skips blocks that are not wanted and refuses a file with no particle-type column. Each block carries `n_declared`, the block size in its header, which should equal `NumPart` in the halo table; comparing the two is left to the caller.

## Constructor Options (`HaloTools`)

```python
HaloTools(
    comoving = True,
    little_h = False,
    native_includes_h = None,    # override the per-format guess; see above
    centre_on_subhalo = False,   # SUBFIND only, see above
    usehalocatonly = False,
    usesubstructure_file = False,
    loglevel = logging.INFO,
)
```

Supported `fileformat` values for `read_catalogue`: `"SUBFIND"`, `"AHF"`, `"VELOCIraptor"`, `"SWIFT_FOF"` (or integer codes 1-4).

See also: [merger_trees.md](merger_trees.md) for walking trees built on top of these catalogues, and [unified_interface.md](unified_interface.md#units-comoving-vs-physical-and-little-h) for the `Dataset`-based API, the full `comoving`/`little_h` semantics (shared with `SnapshotDataset`), and metadata access.
