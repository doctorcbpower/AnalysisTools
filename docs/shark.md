# SHARK Semi-Analytic Catalogues (`shark` package)

The `shark` package lives at `analysistools.shark` (the old top-level `import shark` still works via a deprecation shim). It provides lazy, cached access to SHARK galaxy catalogues, mirroring the [`HaloModel`](halo_catalogues.md) API so both can feed the same analysis and plotting code.

```python
from analysistools.shark.common import _redshift_table, parse_subvolumes
from analysistools.shark.model import SharkModel

rt = _redshift_table("redshift_list.txt")
sv = parse_subvolumes("0-63")

m = SharkModel("./CDM/base_model", rt, sv, label="CDM")

mstar = m.mstars(redshift=0.0)          # mstars_disk + mstars_bulge
sfr   = m.sfr(redshift=0.0)
raw   = m.get("mvir_hosthalo", redshift=0.0)   # any native field by name
```

Available accessors include `mstars`, `mgas`, `sfr`, `ssfr`, `rstar`, `mhalo`, `msubhalo`, `mbh`, `mbulge`, `galaxy_type`, `h0`, `volume`, `age_at_z`, plus star formation history access (`sfh_disk`, `sfh_bulge`, `Z_disk_history`, ...).

## Higher-level components

| Module | Purpose |
|--------|---------|
| `shark.analysis` | `Analysis`: halo/stellar mass functions, SFR main sequence, size–mass, BH–bulge relations |
| `shark.plots` | `Plotter`: standard comparison plots over one or more `SharkModel`s |
| `shark.cli` | Command-line driver |
| `shark.photometry` | fsps-based luminosities/magnitudes from metallicities and star formation histories (requires the `photometry` extra) |

`galaxy_tools.py` (`GalaxyTools`) is the legacy SHARK catalogue reader, superseded by `SharkModel` above.

See also: [unified_interface.md](unified_interface.md) for loading SHARK catalogues via `at.load()`.
