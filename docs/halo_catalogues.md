# Halo Catalogues

Covers `HaloTools` (single-catalogue reader) and `HaloModel` (lazy, cached, multi-snapshot access).

## Working with Halo Catalogues

```python
from analysistools.halo_tools import HaloTools

ht = HaloTools(comoving_units=True)
halos = ht.read_catalogue(
    filename="groups_010.hdf5",
    fileformat="SubFind",       # or "AHF", "VELOCIraptor", "SWIFT_FOF"
    standardise=True,
    snapnum=10,
)
ht.summary()
```

`standardise=True` runs `halo_tools_standardise_names.standardise_catalogue_names`, mapping native per-format fields onto a common schema: `mass`, `pos`, `vel`, `radius`, `halo_id`, `num_part`.

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
    standardise=True,
)

mvir = hm.mvir(redshift=0.0)     # file is only read on first request
pos  = hm.position(redshift=0.0)
```

## Constructor Options (`HaloTools`)

```python
HaloTools(
    comoving_units = False,
    usehalocatonly = False,
    usesubstructure_file = False,
    loglevel = logging.INFO,
)
```

Supported `fileformat` values for `read_catalogue`: `"SUBFIND"`, `"AHF"`, `"VELOCIraptor"`, `"SWIFT_FOF"` (or integer codes 1-4).

See also: [merger_trees.md](merger_trees.md) for walking trees built on top of these catalogues, and [unified_interface.md](unified_interface.md) for the `Dataset`-based API.
