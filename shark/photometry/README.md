# shark/photometry

Standalone photometry pipeline for [Shark](https://github.com/ICRAR/shark) semi-analytic model output.

Converts Shark SFH and metallicity histories into broadband absolute magnitudes via `python-fsps` SSP convolution. Completely independent of ProSpect or any Shark internal tooling.

## Dependencies

```
python >= 3.9
numpy
h5py
astropy
fsps  (python-fsps)
```

Install python-fsps:
```bash
pip install fsps
```

Note: python-fsps requires the FSPS Fortran code to be compiled. See the [python-fsps docs](https://dfm.io/python-fsps).

## Quick start

```python
from shark.model import SharkModel
from shark.common import _redshift_table, parse_subvolumes
from shark.photometry import PhotometryPipeline
import numpy as np

rt    = _redshift_table("redshift_list.txt")
sv    = parse_subvolumes("0-63")
model = SharkModel("./CDM/base_model", rt, sv, label="CDM")

pipe = PhotometryPipeline(
    model,
    z_obs=0.1,
    bands=["v"],
    imf_type=1,      # Chabrier — match your Shark run
    add_dust=False,  # intrinsic magnitudes
)

# All galaxies
M_V = pipe.abs_mag_V()

# Subset
M_V = pipe.abs_mag_V(gal_indices=np.arange(1000))

# Multiple bands
mags = pipe.abs_mags()   # shape (n_gal, n_bands)

# Mass-to-light ratio
ML_V = pipe.mass_to_light()
```

## Module structure

```
shark/photometry/
├── __init__.py      top-level exports
├── io.py            galaxy_data(), sfh_ages_from_model() — adapters over SharkModel
├── sps.py           SPSEngine — python-fsps tabulated SFH convolution
└── photometry.py    PhotometryPipeline — wires everything together
```

There is no `cosmology.py` and no `SharkReader`. All HDF5 reading (including
`star_formation_histories.hdf5`) is owned by `shark.model.SharkModel`; this
package only consumes `SharkModel`'s cached arrays and `SharkModel.age_at_z()`.

## Normalisation

Shark tracks **surviving stellar mass** M★. FSPS normalises internally to **1 M☉ formed**. The pipeline reconciles this by:

```
mags = mags_1msun_formed - 2.5 * log10(M★_surviving / fsps.stellar_mass)
```

where `fsps.stellar_mass` is FSPS's own prediction of the surviving fraction for the same SFH, ensuring self-consistent normalisation.

## Notes

- **Disk and bulge** are convolved separately and combined in linear flux space.
- **Metallicity** is derived as Z(t) = metal mass formed per SFH bin / SFR per bin, from `star_formation_histories.hdf5`, clipped to the FSPS SSP grid range [1e-4, 0.03].
- **Dust**: disabled by default. Enable with `add_dust=True` (Calzetti 2000 law).
- **IMF**: must match the IMF assumed in the Shark run (typically Chabrier).
- Magnitudes are **AB**, **absolute** (not observer-frame).
