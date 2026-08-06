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

from ..merger_tree_types import MergerTreeError, periodic_delta
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
    """OrbitalPericentre/OrbitalApocentre/OrbitalEccentricity and
    OrbitalAngularMomentum (at infall) -- computed from each satellite's
    orbit relative to the host, via ``MergerTreeTools.analyse_orbit()``.

    Needs the host's own tree track, which ``TreeExtractStage`` doesn't
    produce (it only tracks satellites) -- built here directly from
    ``(epoch, host_row)``, matching ``HaloExtractStage``'s convention of
    taking host/satellite selection explicitly rather than guessing it.

    "At infall" uses the satellite's own ``IsSubhalo``-based
    ``HaloTrack.infall_snapshot()`` -- the same definition
    ``HaloPropertiesStage`` uses -- rather than
    ``OrbitAnalysis.first_crossing()``'s R200-based one. That's a deliberate
    choice: the host's virial radius isn't reliably on the tree's `extra`
    for every format (SubFind-HBT carries it as "GroupR200"; a TreeFrog
    track built from a linked halo catalogue, as ``TreeExtractStage``
    produces for the walkable-tree case, carries no radius field at all).

    Deliberately **not** computed here -- each needs an assumed physical
    model this pipeline doesn't specify anywhere, so guessing one would
    silently bake in an unstated assumption:

    - ``OrbitalEnergy``: needs a host potential model (point-mass? NFW?),
      not just kinematics.
    - ``OrbitalPeriod``: needs either a potential model, or multiple fully
      resolved pericentre passages in the snapshot cadence, which isn't
      guaranteed.
    - ``TidalTrackClass``: needs a full tidal-tracks classification model
      (e.g. Penarrubia et al. tables) -- a different kind of model
      entirely, out of scope for orbit kinematics alone.

    These three are left out of ``outputs``/uncomputed here; a future
    stage (or this one, once a modelling choice is made and documented)
    should add them rather than have this stage guess.
    """

    name = "orbital_properties"
    inputs = ("MergerTrees/main_branch",)
    outputs = ("Satellites/HaloProperties/OrbitalPericentre",
               "Satellites/HaloProperties/OrbitalApocentre",
               "Satellites/HaloProperties/OrbitalEccentricity",
               "Satellites/HaloProperties/OrbitalAngularMomentum")

    def __init__(self, epoch, host_row: int):
        self.epoch = epoch
        self.host_row = int(host_row)

    def run(self, context):
        if self.epoch.tree is None:
            raise RuntimeError(
                f"Stage '{self.name}': this Epoch has no merger tree.")
        try:
            host_track = self.epoch.track_of(index=self.host_row).track
        except MergerTreeError as exc:
            raise RuntimeError(
                f"Stage '{self.name}': host halo (row {self.host_row}) has "
                f"no resolvable tree entry, so no satellite orbit can be "
                f"computed relative to it.") from exc

        analyse_orbit = self.epoch.tree.backend.analyse_orbit
        boxsize = self.epoch.boxsize

        main_branch = context.columns["MergerTrees/main_branch"]
        n = len(main_branch)

        pericentre = np.full(n, np.nan)
        apocentre = np.full(n, np.nan)
        eccentricity = np.full(n, np.nan)
        angular_momentum = np.full((n, 3), np.nan)

        n_no_track = n_no_overlap = 0
        for i, track_ds in enumerate(main_branch):
            if track_ds is None or len(track_ds.track) == 0:
                n_no_track += 1
                continue
            track = track_ds.track

            try:
                orbit = analyse_orbit(track, host_track, boxsize=boxsize)
            except MergerTreeError:
                n_no_overlap += 1
                continue

            peri = float(np.min(orbit.Distance))
            apo = float(np.max(orbit.Distance))
            pericentre[i] = peri
            apocentre[i] = apo
            if (peri + apo) > 0:
                eccentricity[i] = (apo - peri) / (apo + peri)

            event = track.infall_snapshot()
            if event is not None:
                infall_snap = track.SnapNum[event["index"]]
                match_self = np.flatnonzero(track.SnapNum == infall_snap)
                match_host = np.flatnonzero(host_track.SnapNum == infall_snap)
                if match_self.size and match_host.size:
                    i_s, i_h = int(match_self[0]), int(match_host[0])
                    dpos = periodic_delta(
                        track.Pos[i_s] - host_track.Pos[i_h], boxsize)
                    dvel = track.Vel[i_s] - host_track.Vel[i_h]
                    angular_momentum[i] = np.cross(dpos, dvel)

        if n_no_track or n_no_overlap:
            logger.warning(
                "%s: %d/%d satellites skipped (%d no resolved tree entry, "
                "%d no snapshot overlap with the host's track).",
                self.name, n_no_track + n_no_overlap, n, n_no_track,
                n_no_overlap)

        context.columns["Satellites/HaloProperties/OrbitalPericentre"] = \
            pericentre
        context.columns["Satellites/HaloProperties/OrbitalApocentre"] = \
            apocentre
        context.columns["Satellites/HaloProperties/OrbitalEccentricity"] = \
            eccentricity
        context.columns["Satellites/HaloProperties/OrbitalAngularMomentum"] = \
            angular_momentum

        context.record_stage(self.name, n_satellites=n,
                             n_no_track=n_no_track,
                             n_no_overlap=n_no_overlap)
        return context


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
