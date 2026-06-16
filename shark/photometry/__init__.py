"""
shark.photometry
================
Standalone photometry pipeline for Shark semi-analytic model output.
Converts Shark SFHs and metallicity histories into broadband magnitudes
via python-fsps SSP convolution.

Submodules
----------
io          : Adapter functions over SharkModel (galaxy_data, sfh_ages_from_model)
sps         : SPSEngine — python-fsps tabulated SFH convolution
photometry  : PhotometryPipeline — main entry point

Cosmology is owned entirely by SharkModel (model.age_at_z / model.lookback_at_z);
there is no separate cosmology object in this package.

Typical usage
-------------
>>> from shark.model import SharkModel
>>> from shark.common import _redshift_table, parse_subvolumes
>>> from shark.photometry import PhotometryPipeline
>>> import numpy as np

>>> rt    = _redshift_table("redshift_list.txt")
>>> sv    = parse_subvolumes("0-63")
>>> model = SharkModel("./CDM/base_model", rt, sv, label="CDM")

>>> pipe  = PhotometryPipeline(model, z_obs=0.1, bands=["v"])
>>> M_V   = pipe.abs_mag_V(gal_indices=np.arange(1000))
"""

from .photometry import PhotometryPipeline
from .sps import SPSEngine
from .io import galaxy_data, sfh_ages_from_model

__all__ = [
    "PhotometryPipeline",
    "SPSEngine",
    "galaxy_data",
    "sfh_ages_from_model",
]
