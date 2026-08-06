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

from typing import Any, Dict, Protocol

import logging

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
    """

    name = "shark"

    def __init__(self, match_by: str = "position", r_scale: float = 1.0):
        self.match_by = match_by
        self.r_scale = r_scale

    def galaxy_properties(self, epoch, halo_row: int) -> Dict[str, Any]:
        raise NotImplementedError(
            "Phase 6b: epoch.galaxies_in_halo(index=halo_row, "
            "match_by=self.match_by, r_scale=self.r_scale) -> aggregate "
            "GalaxyCatalogue fields (mass, sfr, metallicity, sfh, ...) "
            "into the GalaxyProperties field names in schema.py."
        )


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
