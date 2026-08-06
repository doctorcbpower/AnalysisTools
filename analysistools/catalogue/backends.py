#!/usr/bin/env python3
"""
analysistools.catalogue.backends
----------------------------------
GalaxyBackend: the one piece of the catalogue pipeline that genuinely
differs between projects, made pluggable so the schema and every other
stage stay simulation-agnostic.

Dorcha gets galaxy properties from SHARK, a semi-analytic model bolted onto
a DMO merger tree. A hydrodynamic rerun of the same halo (Arepo-Solas,
SWIFT) has stars *in the snapshot* -- there is no SAM step, but the target
fields (StellarMass, MetallicityStellar, StarFormationRate, ...) are the
same physical quantities, just measured a different way. A GalaxyBackend
takes an Epoch and a halo row and returns one row of
``Satellites/GalaxyProperties`` fields; ``pipeline.CatalogueBuilder`` never
knows or cares which backend produced them.

See docs/master_catalogue.md ("The one piece that genuinely differs
between projects") and DEVELOPMENT.md Phase 6.
"""
from __future__ import annotations

from typing import Any, Dict, Optional, Protocol

import logging

import numpy as np

logger = logging.getLogger(__name__)


class GalaxyBackend(Protocol):
    """Interface every galaxy-property source implements."""

    def galaxy_properties(self, epoch: "analysistools.api.simulation.Epoch",
                          halo_row: int) -> Dict[str, Any]:
        """Return one row of Satellites/GalaxyProperties fields (native
        field names, unconverted) for the halo at ``halo_row`` in
        ``epoch.halos``. Missing quantities should be omitted (not NaN-
        filled) so ``validation.IntegrityValidator`` can distinguish
        "not applicable" from "failed to compute"."""
        ...


class SharkGalaxyBackend:
    """Galaxy properties from a SHARK semi-analytic catalogue, matched to
    the halo catalogue via ``Epoch.galaxies_in_halo`` (position matching by
    default -- see api/simulation.py).

    When more than one SHARK galaxy falls within the match aperture
    (``r_scale`` * halo radius), the most massive one (by stellar mass) is
    treated as "the" galaxy hosted by this halo -- the central, in SAM
    convention. Returns ``{}`` (per the ``GalaxyBackend`` protocol) when
    nothing matches, or when no stellar mass is resolvable at all (so no
    central can be identified).

    Fields returned, and where they come from (see
    ``analysistools.shark.model.GALAXY_FIELDS`` for the full native-name
    registry):

    - ``StellarMass``, ``GasMass_Cold``, ``StarFormationRate``,
      ``BlackHoleMass``: the ``GalaxyCatalogue``'s own standardised
      aliases (``mass``, ``mgas``, ``sfr``, ``mbh``).
    - ``GasMass_Hot``: native ``"mhot"`` (not part of the standardised
      alias set).
    - ``MetallicityStellar``/``MetallicityGas``: metal mass / total mass,
      from the native ``*_metals_disk``/``*_metals_bulge`` fields -- a
      genuinely unit- and h-independent ratio, unlike everything else
      here (see units caveat below).

    Units caveat (same as ``pipeline.HaloExtractStage``/
    ``derived.EnvironmentStage``): values are returned exactly as
    ``GalaxyCatalogue`` provides them -- SHARK's own native convention
    (``Msun/h`` masses, ``Msun/Gyr/h`` star formation rates -- note that's
    *also* a Gyr-vs-yr mismatch against schema.py's declared "Msun/yr",
    not just little-h), not schema.py's declared units. Unit
    reconciliation isn't implemented anywhere in this pipeline yet.

    Deliberately **not** returned, documented rather than guessed:
    ``LuminosityV``/``Luminosity_ugriz``/``AbsoluteMagnitude_ugriz`` (need
    ``shark.photometry``'s separate, optional fsps-based pipeline -- not
    wired into this backend), ``HalfLightRadius`` (depends on
    luminosity), ``HalfMassRadiusStellar`` (``rstar_disk``/``rstar_bulge``
    are exponential/profile scale lengths, not half-mass radii --
    converting one to the other needs an assumed disk/bulge profile
    shape), ``SersicIndex`` (SHARK produces no structural fit at all).
    """

    name = "shark"

    def __init__(self, match_by: str = "position", r_scale: float = 1.0):
        self.match_by = match_by
        self.r_scale = r_scale

    def galaxy_properties(self, epoch, halo_row: int) -> Dict[str, Any]:
        matched = epoch.galaxies_in_halo(
            index=halo_row, match_by=self.match_by, r_scale=self.r_scale)
        if len(matched) == 0 or "mass" not in matched:
            return {}

        mass = np.asarray(matched["mass"])
        central = int(np.argmax(mass))

        def scalar(field: str, default: Optional[float] = None):
            arr = matched.get(field)
            return (float(np.asarray(arr)[central]) if arr is not None
                    else default)

        def metal_mass(disk_field: str, bulge_field: str):
            """Sum of disk+bulge metal mass, or None if *neither* component
            is available -- unlike scalar()'s own default=0.0 use below,
            this must not silently produce 0.0 for "field never existed"
            (that would compute a fake MetallicityStellar/Gas=0.0 instead
            of omitting it, exactly what the GalaxyBackend protocol says
            not to do)."""
            disk, bulge = scalar(disk_field), scalar(bulge_field)
            if disk is None and bulge is None:
                return None
            return (disk or 0.0) + (bulge or 0.0)

        stellar_mass = scalar("mass")
        gas_mass_cold = scalar("mgas")
        stellar_metal_mass = metal_mass("mstars_metals_disk",
                                        "mstars_metals_bulge")
        gas_metal_mass = metal_mass("mgas_metals_disk", "mgas_metals_bulge")

        props: Dict[str, Any] = {}
        if stellar_mass is not None:
            props["StellarMass"] = stellar_mass
            if stellar_metal_mass is not None and stellar_mass > 0:
                props["MetallicityStellar"] = \
                    stellar_metal_mass / stellar_mass
        if gas_mass_cold is not None:
            props["GasMass_Cold"] = gas_mass_cold
            if gas_metal_mass is not None and gas_mass_cold > 0:
                props["MetallicityGas"] = gas_metal_mass / gas_mass_cold
        hot = scalar("mhot")
        if hot is not None:
            props["GasMass_Hot"] = hot
        sfr = scalar("sfr")
        if sfr is not None:
            props["StarFormationRate"] = sfr
        mbh = scalar("mbh")
        if mbh is not None:
            props["BlackHoleMass"] = mbh

        return props


class HydroGalaxyBackend:
    """Galaxy properties synthesised from star particles in a
    hydrodynamic snapshot, matched to the halo via
    ``Epoch.particles_in_halo(species="star", ...)``.

    Aggregation recipes (to implement in Phase 6b) should target the exact
    same field names as SharkGalaxyBackend so a hydro catalogue is
    schema-identical wherever the physics is comparable:

    - StellarMass       <- sum of star particle masses
    - MetallicityStellar <- mass-weighted mean of particle metallicity
    - MeanStellarAge     <- mass-weighted mean of (snapshot time - particle
                             formation time)
    - StarFormationRate  <- mass formed in young (< age threshold)
                             particles / threshold time window
    - HalfLightRadius / HalfMassRadiusStellar <- computed directly from
                             particle positions relative to the halo centre
    """

    name = "hydro"

    def __init__(self, r_scale: float = 1.0,
                 young_star_age_threshold: float = 0.1):
        self.r_scale = r_scale
        self.young_star_age_threshold = young_star_age_threshold

    def galaxy_properties(self, epoch, halo_row: int) -> Dict[str, Any]:
        raise NotImplementedError(
            "Phase 6b: epoch.particles_in_halo(index=halo_row, "
            "species='star', r_scale=self.r_scale) -> aggregate into the "
            "same GalaxyProperties field names as SharkGalaxyBackend."
        )


BACKENDS = {
    "shark": SharkGalaxyBackend,
    "hydro": HydroGalaxyBackend,
}


def get_backend(name: str, **kwargs) -> GalaxyBackend:
    """Resolve a project config's ``galaxy_backend: shark|hydro`` string to
    an instantiated backend."""
    try:
        cls = BACKENDS[name]
    except KeyError:
        raise ValueError(
            f"Unknown galaxy_backend '{name}'; available: "
            f"{sorted(BACKENDS)}")
    return cls(**kwargs)
