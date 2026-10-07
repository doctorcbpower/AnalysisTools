"""
analysistools.orbits
--------------------
Standalone orbit tracking: host-frame kinematics, infall/pericentre/
apocentre events, spherical host potentials, orphan continuation, and an
HDF5 sidecar (``OrbitCatalogue``) that the master-catalogue pipeline reads.

See ``docs/orbits.md`` for definitions and assumptions.
"""
from .cosmology import FlatLCDM, G_KPC_KMS2_MSUN, KMS_TO_KPC_PER_GYR
from .kinematics import (TrackUnits, HostFrameSeries, FineOrbit, Event,
                         host_frame_series, host_reference, interpolate_orbit,
                         find_crossings, find_turning_points,
                         CROSS_IN, CROSS_OUT, PERICENTRE, APOCENTRE)
from .potential import (Potential, NFWPotential, TabulatedSpherical,
                        TimeDependentPotential)
from .orphans import integrate_orphans
from .selection import (SatelliteSelection, select_satellites,
                        select_satellites_subfind)
from . import benson
from .tracker import (OrbitConfig, OrbitRecord, OrbitTracker, OrbitCatalogue,
                      tree_nfw_potential, SCALARS)

__all__ = [
    "FlatLCDM", "G_KPC_KMS2_MSUN", "KMS_TO_KPC_PER_GYR",
    "TrackUnits", "HostFrameSeries", "FineOrbit", "Event",
    "host_frame_series", "host_reference", "interpolate_orbit",
    "find_crossings", "find_turning_points",
    "CROSS_IN", "CROSS_OUT", "PERICENTRE", "APOCENTRE",
    "Potential", "NFWPotential", "TabulatedSpherical",
    "TimeDependentPotential", "integrate_orphans",
    "SatelliteSelection", "select_satellites", "select_satellites_subfind",
    "OrbitConfig", "OrbitRecord", "OrbitTracker", "OrbitCatalogue",
    "tree_nfw_potential", "SCALARS", "benson",
]
