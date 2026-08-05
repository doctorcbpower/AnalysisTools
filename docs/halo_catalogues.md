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
