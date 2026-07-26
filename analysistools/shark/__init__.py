"""
shark — Analysis package for SHARK galaxy formation model output.

Public API
----------
SharkModel
    Lazy, cached access to a single SHARK model run's HDF5 catalogues.

Analysis
    Computes binned statistics (mass functions, scaling relations) from
    one or more SharkModel instances.

Plotter
    Renders Analysis result dicts to publication-quality matplotlib figures.

Quick-start
-----------
>>> from common import _redshift_table, parse_subvolumes
>>> from shark import SharkModel, Analysis, Plotter
>>> import numpy as np

>>> rt  = _redshift_table("redshift_list.txt")
>>> sv  = parse_subvolumes("0-63")
>>> z   = np.array([0.0, 1.0, 2.0])

>>> models = [
...     SharkModel("./CDM/base_model", rt, sv, label="CDM", colour="#e41a1c"),
...     SharkModel("./WDM/base_model", rt, sv, label="WDM", colour="#377eb8"),
... ]

>>> analysis = Analysis(models)
>>> plotter  = Plotter("./plots")

>>> plotter.plot_halo_mf(   analysis.compute_halo_mf(z))
>>> plotter.plot_stellar_mf(analysis.compute_stellar_mf(z))
>>> plotter.plot_sfr_main_sequence(analysis.compute_sfr_main_sequence(z))
>>> plotter.plot_size_mass( analysis.compute_size_mass(z))
>>> plotter.plot_bh_bulge(  analysis.compute_bh_bulge(z))

# Custom relation — cold gas fraction vs stellar mass
>>> results = analysis.compute_custom(
...     redshifts=[0.0],
...     x_func=lambda m, z: m.mstars(z),
...     y_func=lambda m, z: m.mgas(z) / (m.mstars(z) + m.mgas(z)),
...     sel_func=lambda m, z: m.mstars(z) > 1e8,
...     x_log=True,
...     y_log=False,
... )
>>> plotter.plot_custom(results, axis_cfg=dict(
...     xtit=r"$\\log_{10}(M_\\star/M_\\odot)$",
...     ytit=r"$f_{\\rm gas}$",
...     ymin=0, ymax=1,
... ))

CLI
---
python -m shark.cli \\
    -z redshift_list.txt -v '0-63' -o ./plots \\
    --model-dirs ./CDM/base_model ./WDM/base_model \\
    --redshifts 0 1 2 --labels CDM WDM \\
    --plots halo_mf stellar_mf
"""

from .model    import SharkModel, GALAXY_FIELDS
from .analysis import Analysis, make_bins, number_density
from .plots    import Plotter, DEFAULT_COLOURS

__all__ = [
    # Core classes
    "SharkModel",
    "Analysis",
    "Plotter",
    # Field catalogue
    "GALAXY_FIELDS",
    # Analysis utilities
    "make_bins",
    "number_density",
    # Plot utilities
    "DEFAULT_COLOURS",
]
