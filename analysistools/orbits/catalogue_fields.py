"""
analysistools.orbits.catalogue_fields
-------------------------------------
Mapping from orbit-sidecar scalars (``tracker.SCALARS``) to master-catalogue
columns. Single source of truth for both ``OrbitalPropertiesStage`` and the
``schema_v1.1.yaml`` fragment (see ``schema_fragment()``).

Pre-existing schema 1.0 fields (filled from the sidecar, same names/units):

* ``OrbitalPericentre``  <- ``r_peri_min``  (smallest confirmed pericentre
  since infall; NaN until the first pericentre has been passed -- the old
  stage's snapshot-sampled ``min(distance)`` is gone)
* ``OrbitalApocentre``   <- ``r_apo_last``  (latest *confirmed* apocentre
  since infall, not the maximum over the whole history)
* ``OrbitalEccentricity``<- ``ecc_last``
* ``OrbitalPeriod``      <- ``period_peri`` (measured, needs >= 2 pericentres)
* ``OrbitalEnergy``      <- ``E_infall``    (needs a host potential)
* ``OrbitalAngularMomentum`` <- ``L_infall`` (at first R_ref crossing)

Schema 1.1 adds the rest, all prefixed ``Orbit*`` so none collides with the
FOF-membership-based ``HaloPropertiesStage`` infall fields (which keep their
meaning).
"""
from __future__ import annotations

from .tracker import SCALARS

GROUP = "Satellites/HaloProperties"

LEGACY_FIELD_MAP = {
    "OrbitalPericentre": "r_peri_min",
    "OrbitalApocentre": "r_apo_last",
    "OrbitalEccentricity": "ecc_last",
    "OrbitalPeriod": "period_peri",
    "OrbitalEnergy": "E_infall",
    "OrbitalAngularMomentum": "L_infall",
}

EXTENDED_FIELD_MAP = {
    "OrbitInfallTime": "t_infall",
    "OrbitInfallRedshift": "z_infall",
    "OrbitInfallSnapshot": "snap_infall",
    "OrbitInfallMass": "mass_infall",
    "OrbitInfallVmax": "vmax_infall",
    "OrbitInfallRadialVelocity": "vr_infall",
    "OrbitInfallTangentialVelocity": "vt_infall",
    "OrbitBornInside": "born_inside",
    "OrbitFOFInfallTime": "t_fof_infall",
    "OrbitPreprocessed": "preprocessed",
    "OrbitNumberOfInfalls": "n_infalls",
    "OrbitNumberOfOutfalls": "n_outfalls",
    "OrbitIsBacksplash": "is_backsplash",
    "OrbitInsideReference": "inside_now",
    "OrbitInCountSample": "in_count_sample",
    "OrbitInOrbitSample": "in_orbit_sample",
    "OrbitNumberOfPericentres": "n_pericentres",
    "OrbitPericentreFirst": "r_peri_first",
    "OrbitPericentreLast": "r_peri_last",
    "OrbitTimeLastPericentre": "t_peri_last",
    "OrbitApocentreFirst": "r_apo_first",
    "OrbitPericentreEL": "r_peri_EL",
    "OrbitApocentreEL": "r_apo_EL",
    "OrbitPeriodEL": "period_EL",
    "OrbitFinalRadius": "r_final",
    "OrbitPeakMass": "mass_peak",
    "OrbitFinalMass": "mass_final",
    "OrbitPeakParticles": "npart_peak",
    "OrbitFinalParticles": "npart_final",
    "OrbitIsOrphan": "is_orphan_final",
    "OrbitOrphanContinued": "orphan_continued",
}

_DTYPE = {"float64": "float64", "int32": "int32", "int16": "int16",
          "bool": "bool"}


def schema_fragment() -> str:
    """YAML field entries for ``EXTENDED_FIELD_MAP`` (schema 1.1 additions)."""
    out = []
    for col, key in EXTENDED_FIELD_MAP.items():
        dtype, _, units, desc = SCALARS[key]
        out.append(
            f"  - {{name: {col}, group: {GROUP}, dtype: {_DTYPE[dtype]}, "
            f"units: {units}, is_derived: true, "
            f"provenance: \"orbit tracker (analysistools.orbits)\",\n"
            f"     description: \"{desc}\"}}")
    return "\n".join(out) + "\n"
