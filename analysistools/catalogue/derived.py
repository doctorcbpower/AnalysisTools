#!/usr/bin/env python3
"""
analysistools.catalogue.derived
----------------------------------
DerivedQuantityStage subclasses -- one per physically-related family, kept
small and single-purpose so a future paper can contribute its own stage as
a self-contained, reviewable addition rather than a patch to one giant
function. See docs/dorcha_master_catalogue_design.md section 3 for the
full recommended derived-quantity list and section 2.2 for field
definitions; each class below computes one row of that table.

A project's config lists which of these run, by name (see pipeline.py /
docs/master_catalogue.md "template mechanism"). All operate on a
``PipelineContext`` populated by the Extract/Cross-match stages plus
whichever ``GalaxyBackend`` the project selected.
"""
from __future__ import annotations

from .pipeline import PipelineStage


class HaloPropertiesStage(PipelineStage):
    """Mpeak, SnapshotAtMpeak, Vpeak, Minfall/Vinfall/RedshiftInfall,
    NumberOfInfalls, IsBacksplash -- from the merger tree's main branch."""

    name = "halo_properties"
    inputs = ("MergerTrees/main_branch",)
    outputs = ("Satellites/HaloProperties/Mpeak",
               "Satellites/HaloProperties/RedshiftInfall")

    def run(self, context):
        raise NotImplementedError("Phase 6b.")


class OrbitalPropertiesStage(PipelineStage):
    """Pericentre/apocentre/eccentricity/period/energy/angular momentum,
    TidalTrackClass -- orbit-fitting against the tree's position/velocity
    history relative to the host."""

    name = "orbital_properties"
    inputs = ("MergerTrees/main_branch",)
    outputs = ("Satellites/HaloProperties/OrbitalPericentre",
               "Satellites/HaloProperties/OrbitalApocentre")

    def run(self, context):
        raise NotImplementedError("Phase 6b.")


class StarFormationHistoryStage(PipelineStage):
    """Rebins the galaxy backend's native SFH onto the common
    ``Snapshots/time_bin_edges_sfh`` grid; derives MeanStellarAge,
    QuenchingTime, IsQuenched_z0."""

    name = "star_formation_history"
    inputs = ("Satellites/GalaxyProperties/StellarMass",)
    outputs = ("Satellites/GalaxyProperties/SFH",
               "Satellites/GalaxyProperties/QuenchingTime")

    def run(self, context):
        raise NotImplementedError("Phase 6b.")


class EnvironmentStage(PipelineStage):
    """LocalNumberDensity, DistanceToNearestMassiveNeighbour, TidalIndex,
    CosmicWebClass, host isolation/pairing copy-through."""

    name = "environment"
    inputs = ("Haloes/Position",)
    outputs = ("Satellites/Environment/LocalNumberDensity",
               "Satellites/Environment/CosmicWebClass")

    def run(self, context):
        raise NotImplementedError("Phase 6b.")


class ObservabilityStage(PipelineStage):
    """Heliocentric distance, Galactic coordinates, apparent magnitudes,
    surface brightness, Rubin detection probability, completeness weight."""

    name = "observability"
    inputs = ("Satellites/HaloProperties/",
              "Satellites/GalaxyProperties/LuminosityV")
    outputs = ("Satellites/Observability/RubinDetectionProbability",
               "Satellites/Observability/CompletenessWeight")

    def run(self, context):
        raise NotImplementedError("Phase 6b.")


class DorchaSpecificStage(PipelineStage):
    """ProgenitorParticleFraction, EarliestProgenitorRedshift,
    PeakOverdensity, FormationEnvironmentClass, FossilFraction -- requires
    particle tagging + progenitor-field analysis inputs; not applicable to
    projects without particle tagging (omit from ``derived_stages``)."""

    name = "dorcha_specific"
    inputs = ("Satellites/ParticleTags/particle_ids",)
    outputs = ("Satellites/DorchaProperties/ProgenitorParticleFraction",
               "Satellites/DorchaProperties/FossilFraction")

    def run(self, context):
        raise NotImplementedError("Phase 6b.")


# Registry used by CatalogueBuilder to resolve a config's `derived_stages`
# name list to classes.
STAGES = {
    "halo_properties": HaloPropertiesStage,
    "orbital_properties": OrbitalPropertiesStage,
    "star_formation_history": StarFormationHistoryStage,
    "environment": EnvironmentStage,
    "observability": ObservabilityStage,
    "dorcha_specific": DorchaSpecificStage,
}
