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

import logging

import numpy as np

from .pipeline import PipelineStage

logger = logging.getLogger(__name__)


class HaloPropertiesStage(PipelineStage):
    """Mpeak, SnapshotAtMpeak, Vpeak, Minfall/Vinfall/RedshiftInfall/
    SnapshotInfall, NumberOfInfalls, IsBacksplash -- from each satellite's
    own merger-tree main branch (``MergerTrees/main_branch``, one
    ``TrackDataset`` per satellite, produced by
    ``pipeline.TreeExtractStage``).

    "Infall" here is a satellite's own IsSubhalo history (the first
    central->satellite transition along its track, via
    ``HaloTrack.infall_snapshot()``) -- a simple, always-available
    definition. It is *not* the more physical "crossed inside the host's
    R200c" definition (that needs the host's own track + Rvir history via
    ``MergerTreeTools.analyse_orbit``/``OrbitAnalysis.first_crossing()``,
    which is ``OrbitalPropertiesStage``'s job, not this one).

    Satellites with no resolved tree entry (``main_branch`` entry is
    ``None``, e.g. below the tree's resolution) get NaN/-1 sentinels for
    every field here rather than raising -- flagging that as a QC issue is
    ``IntegrityValidator``'s job (Phase 6c), not this stage's.
    """

    name = "halo_properties"
    inputs = ("MergerTrees/main_branch",)
    outputs = ("Satellites/HaloProperties/Mpeak",
               "Satellites/HaloProperties/SnapshotAtMpeak",
               "Satellites/HaloProperties/Vpeak",
               "Satellites/HaloProperties/Minfall",
               "Satellites/HaloProperties/Vinfall",
               "Satellites/HaloProperties/RedshiftInfall",
               "Satellites/HaloProperties/SnapshotInfall",
               "Satellites/HaloProperties/NumberOfInfalls",
               "Satellites/HaloProperties/IsBacksplash")

    def run(self, context):
        main_branch = context.columns["MergerTrees/main_branch"]
        n = len(main_branch)

        mpeak = np.full(n, np.nan)
        snapshot_at_mpeak = np.full(n, -1, dtype=np.int32)
        vpeak = np.full(n, np.nan)
        minfall = np.full(n, np.nan)
        vinfall = np.full(n, np.nan)
        redshift_infall = np.full(n, np.nan)
        snapshot_infall = np.full(n, -1, dtype=np.int32)
        number_of_infalls = np.zeros(n, dtype=np.int16)
        is_backsplash = np.zeros(n, dtype=bool)

        n_no_track = 0
        for i, track_ds in enumerate(main_branch):
            if track_ds is None:
                n_no_track += 1
                continue
            track = track_ds.track
            if len(track) == 0:
                n_no_track += 1
                continue

            mass = track.Mass
            i_peak = int(np.argmax(mass))
            mpeak[i] = float(mass[i_peak])
            snapshot_at_mpeak[i] = int(track.SnapNum[i_peak])

            vmax_history = track.extra.get("Vmax")
            has_vmax = vmax_history is not None and len(vmax_history) == len(track)
            if has_vmax:
                vpeak[i] = float(np.max(vmax_history))

            event = track.infall_snapshot()
            if event is not None:
                minfall[i] = event["mass"]
                redshift_infall[i] = event["redshift"]
                snapshot_infall[i] = event["snapshot"]
                if has_vmax:
                    vinfall[i] = float(vmax_history[event["index"]])

            is_sub = track.IsSubhalo
            # count False->True transitions walking forward in time; a
            # track that's already a subhalo at its earliest snapshot
            # doesn't have an observed transition, so isn't counted here
            number_of_infalls[i] = int(
                np.count_nonzero(is_sub[1:] & ~is_sub[:-1]))
            is_backsplash[i] = bool(np.any(is_sub)) and not bool(is_sub[-1])

        if n_no_track:
            logger.warning(
                "%s: %d/%d satellites have no resolved tree entry "
                "(MergerTrees/main_branch is None or empty); their "
                "HaloProperties fields are NaN/-1.",
                self.name, n_no_track, n)

        context.columns["Satellites/HaloProperties/Mpeak"] = mpeak
        context.columns["Satellites/HaloProperties/SnapshotAtMpeak"] = \
            snapshot_at_mpeak
        context.columns["Satellites/HaloProperties/Vpeak"] = vpeak
        context.columns["Satellites/HaloProperties/Minfall"] = minfall
        context.columns["Satellites/HaloProperties/Vinfall"] = vinfall
        context.columns["Satellites/HaloProperties/RedshiftInfall"] = \
            redshift_infall
        context.columns["Satellites/HaloProperties/SnapshotInfall"] = \
            snapshot_infall
        context.columns["Satellites/HaloProperties/NumberOfInfalls"] = \
            number_of_infalls
        context.columns["Satellites/HaloProperties/IsBacksplash"] = \
            is_backsplash

        context.record_stage(self.name, n_satellites=n,
                             n_no_track=n_no_track)
        return context


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
