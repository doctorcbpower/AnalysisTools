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
    """``SFH``, ``MeanStellarAge``, ``QuenchingTime``, ``IsQuenched_z0`` --
    via the configured ``GalaxyBackend``'s ``star_formation_history()``,
    rebinned onto a common ``time_bin_edges`` grid (required constructor
    argument, matching schema.py's ``Snapshots/time_bin_edges_sfh`` -- no
    default is picked here, same reasoning as ``EnvironmentStage``'s
    ``mass_threshold``/``aperture_radius``). ``time_bin_edges`` must be
    ascending *lookback* time in Gyr, index 0 = now -- both backends'
    ``star_formation_history()`` and the quenching-time walk below assume
    that ordering.

    ``MeanStellarAge`` here **supersedes** whatever
    ``galaxy_backend.galaxy_properties()`` already put in context (both
    ``SharkGalaxyBackend`` and ``HydroGalaxyBackend`` compute a quick one
    too, from their own native per-galaxy/per-particle data) with a
    version computed from the same common time grid for every satellite
    regardless of backend -- which is the point: a SHARK and a hydro
    catalogue are only directly comparable if ``MeanStellarAge`` means the
    same thing in both. ``StarFormationRate`` (schema: primary, straight
    from the backend) is *not* touched here.

    ``QuenchingTime``/``IsQuenched_z0`` deliberately deviate from the
    schema's literal "sSFR below threshold x sSFR_MS" definition: that
    needs an assumed star-forming-main-sequence parametrisation (e.g.
    Speagle et al. 2014) this pipeline doesn't specify anywhere.
    ``quenched_ssfr_threshold`` is instead a required, **absolute** sSFR
    cut (1/yr -- e.g. the classic ~1e-11/yr), not main-sequence-relative.
    Extend this stage if the MS-relative definition is actually needed
    once a specific SFMS model is chosen.

    Satellites with no computable SFH (backend's
    ``star_formation_history()`` returns ``None``) get an all-NaN ``SFH``
    row and NaN/``False`` for the derived fields, not a stage failure.
    """

    name = "star_formation_history"
    inputs = ("Satellites/_internal/halo_row",
             "Satellites/GalaxyProperties/StellarMass")
    outputs = ("Satellites/GalaxyProperties/SFH",
               "Satellites/GalaxyProperties/MeanStellarAge",
               "Satellites/GalaxyProperties/QuenchingTime",
               "Satellites/GalaxyProperties/IsQuenched_z0")

    def __init__(self, epoch, galaxy_backend, time_bin_edges,
                 quenched_ssfr_threshold: float):
        self.epoch = epoch
        self.galaxy_backend = galaxy_backend
        self.time_bin_edges = np.asarray(time_bin_edges, dtype=float)
        self.quenched_ssfr_threshold = float(quenched_ssfr_threshold)

    def run(self, context):
        halo_rows = context.columns["Satellites/_internal/halo_row"]
        n = len(halo_rows)
        n_bins = len(self.time_bin_edges) - 1
        bin_centres = 0.5 * (self.time_bin_edges[:-1]
                             + self.time_bin_edges[1:])
        bin_widths_yr = np.diff(self.time_bin_edges) * 1.0e9

        stellar_mass = context.columns.get(
            "Satellites/GalaxyProperties/StellarMass")
        sfr_z0 = context.columns.get(
            "Satellites/GalaxyProperties/StarFormationRate")

        sfh = np.full((n, n_bins), np.nan)
        mean_age = np.full(n, np.nan)
        quenching_time = np.full(n, np.nan)
        is_quenched = np.zeros(n, dtype=bool)

        n_no_sfh = 0
        for i, row in enumerate(halo_rows):
            row_sfh = self.galaxy_backend.star_formation_history(
                self.epoch, int(row), self.time_bin_edges)
            if row_sfh is None:
                n_no_sfh += 1
                continue
            sfh[i] = row_sfh

            formed_mass = row_sfh * bin_widths_yr
            total_formed = float(np.sum(formed_mass))
            if total_formed > 0:
                mean_age[i] = float(
                    np.average(bin_centres, weights=formed_mass))

            mstar_i = (stellar_mass[i] if stellar_mass is not None
                      and not np.isnan(stellar_mass[i]) else None)
            sfr_i = (sfr_z0[i] if sfr_z0 is not None
                    and not np.isnan(sfr_z0[i]) else None)
            if mstar_i is not None and mstar_i > 0 and sfr_i is not None:
                is_quenched[i] = (sfr_i / mstar_i) < self.quenched_ssfr_threshold
                if is_quenched[i]:
                    ssfr_history = row_sfh / mstar_i
                    active = ssfr_history >= self.quenched_ssfr_threshold
                    if np.any(active):
                        # first (smallest-lookback, i.e. most recent)
                        # active bin -- the last time star formation was
                        # above threshold before falling permanently below
                        quenching_time[i] = bin_centres[int(np.argmax(active))]

        if n_no_sfh:
            logger.warning(
                "%s: %d/%d satellites have no computable star formation "
                "history (galaxy_backend.star_formation_history() "
                "returned None); SFH/MeanStellarAge/QuenchingTime/"
                "IsQuenched_z0 left at their NaN/False sentinels.",
                self.name, n_no_sfh, n)

        context.columns["Satellites/GalaxyProperties/SFH"] = sfh
        context.columns["Satellites/GalaxyProperties/MeanStellarAge"] = \
            mean_age
        context.columns["Satellites/GalaxyProperties/QuenchingTime"] = \
            quenching_time
        context.columns["Satellites/GalaxyProperties/IsQuenched_z0"] = \
            is_quenched
        # so downstream stages that also need to interpret the SFH array
        # (e.g. DorchaSpecificStage's FossilFraction) use the exact same
        # bins, rather than risk a second, independently-specified and
        # possibly mismatched time_bin_edges argument.
        context.meta["time_bin_edges_sfh"] = self.time_bin_edges

        context.record_stage(self.name, n_satellites=n, n_no_sfh=n_no_sfh)
        return context


class EnvironmentStage(PipelineStage):
    """LocalNumberDensity, DistanceToNearestMassiveNeighbour -- computed
    from the *other* halos in ``Epoch.halos`` (every row except the host
    and the satellites already selected by ``HaloExtractStage``), which
    for a Dorcha-style zoom-in box are the field/companion halos
    surrounding the host system.

    ``mass_threshold``/``aperture_radius`` are required constructor
    arguments, not defaulted: the design doc specifies "above a mass
    threshold" / "within a fixed aperture", but no specific values are
    documented anywhere in this codebase -- pass the ones your project
    actually wants rather than have this stage silently pick a number.

    Units caveat (same as ``HaloExtractStage``): computed in whatever
    length/mass units the Epoch's ``HaloCatalogue`` was configured with
    (comoving/little_h), not necessarily schema.py's declared "Mpc^-3
    (comoving)"/"Mpc" -- unit reconciliation against the schema isn't
    implemented anywhere in the pipeline yet.

    Deliberately **not** computed here, documented rather than guessed:

    - ``TidalIndex``: the Karachentsev-style tidal index has more than one
      non-equivalent convention in the literature (single most tidally-
      dominant neighbour vs. a sum over neighbours; what reference mass/
      offset constant, if any, is subtracted) -- picking one silently
      would bake in an unstated convention.
    - ``CosmicWebClass``: needs an actual web classifier (T-web/V-web or
      similar) -- a substantial separate piece of machinery, not a
      neighbour-counting operation.

    ``HostIsIsolated``/``HostIsPaired`` *are* computed here, but only as a
    copy-through: if ``HostEnvironmentStage`` has already run (and
    populated ``Haloes/IsIsolated``/``IsPaired``), its single host-level
    value is broadcast onto every satellite row. If it hasn't run, these
    two columns are simply omitted (not NaN-filled -- there is nothing to
    copy).
    """

    name = "environment"
    inputs = ("Satellites/_internal/pos_z0", "Satellites/_internal/halo_row")
    outputs = ("Satellites/Environment/LocalNumberDensity",
               "Satellites/Environment/DistanceToNearestMassiveNeighbour",
               "Satellites/Environment/HostIsIsolated",
               "Satellites/Environment/HostIsPaired")

    def __init__(self, epoch, host_row: int, mass_threshold: float,
                 aperture_radius: float):
        self.epoch = epoch
        self.host_row = int(host_row)
        self.mass_threshold = float(mass_threshold)
        self.aperture_radius = float(aperture_radius)

    def run(self, context):
        halos = self.epoch.halos
        if halos is None:
            raise RuntimeError(
                f"Stage '{self.name}': this Epoch has no halo catalogue.")

        satellite_halo_rows = np.asarray(
            context.columns["Satellites/_internal/halo_row"])
        excluded = set(int(r) for r in satellite_halo_rows) | {self.host_row}

        mass = np.asarray(halos["mass"])
        pos = np.asarray(halos["pos"])
        neighbour_mask = np.ones(len(mass), dtype=bool)
        neighbour_mask[list(excluded)] = False
        neighbour_mask &= mass >= self.mass_threshold
        neighbour_pos = pos[neighbour_mask]

        boxsize = self.epoch.boxsize
        satellite_pos = np.asarray(
            context.columns["Satellites/_internal/pos_z0"])
        n = satellite_pos.shape[0]

        local_density = np.full(n, np.nan)
        distance_to_nearest = np.full(n, np.nan)
        volume = (4.0 / 3.0) * np.pi * self.aperture_radius ** 3

        if neighbour_pos.shape[0] == 0:
            logger.warning(
                "%s: no neighbours found above mass_threshold=%.3e; "
                "LocalNumberDensity=0 and "
                "DistanceToNearestMassiveNeighbour left NaN for every "
                "satellite.", self.name, self.mass_threshold)
            local_density[:] = 0.0
        else:
            for i in range(n):
                dpos = periodic_delta(neighbour_pos - satellite_pos[i],
                                      boxsize)
                dist = np.linalg.norm(dpos, axis=1)
                local_density[i] = \
                    np.count_nonzero(dist <= self.aperture_radius) / volume
                distance_to_nearest[i] = float(np.min(dist))

        context.columns["Satellites/Environment/LocalNumberDensity"] = \
            local_density
        context.columns[
            "Satellites/Environment/DistanceToNearestMassiveNeighbour"] = \
            distance_to_nearest

        host_isolated = context.columns.get("Haloes/IsIsolated")
        host_paired = context.columns.get("Haloes/IsPaired")
        if host_isolated is not None:
            context.columns["Satellites/Environment/HostIsIsolated"] = \
                np.full(n, bool(np.asarray(host_isolated)[0]))
        if host_paired is not None:
            context.columns["Satellites/Environment/HostIsPaired"] = \
                np.full(n, bool(np.asarray(host_paired)[0]))

        context.record_stage(self.name, n_satellites=n,
                             n_neighbours=int(neighbour_pos.shape[0]),
                             mass_threshold=self.mass_threshold,
                             aperture_radius=self.aperture_radius)
        return context


class HostEnvironmentStage(PipelineStage):
    """``Haloes/IsIsolated``, ``IsPaired``, ``PairedHostID``,
    ``N_satellites_total`` -- host-level fields (``Haloes/`` is a
    one-host-per-pipeline-run table, see ``pipeline.HaloExtractStage``),
    computed from the *other* halos in ``Epoch.halos``, the same neighbour
    source ``EnvironmentStage`` uses for its per-satellite fields. Distinct
    from that stage (which is per-satellite and deliberately does *not*
    compute these -- see its own docstring) and from ``Satellites/
    Environment/HostIsIsolated``/``HostIsPaired`` (a copy-through of this
    stage's output onto every satellite row, done by ``EnvironmentStage``
    when this stage has already run).

    Every definition below is explicit and driven entirely by required
    constructor arguments -- no default is picked -- per the design doc's
    "explicit ownership of ambiguous definitions" principle (isolation and
    pairing both admittedly have more than one reasonable convention in
    the literature):

    - ``IsIsolated``: True iff no *other* halo in ``Epoch.halos`` more
      massive than this host lies within
      ``isolation_radius_factor * R200c_host`` of the host's position.
    - ``IsPaired`` / ``PairedHostID``: True iff at least one *other* halo
      (any of these can also count towards ``IsIsolated``'s "more massive"
      neighbour -- the two flags are independent, not mutually exclusive)
      has mass within ``[pairing_mass_ratio_min, 1/pairing_mass_ratio_min]``
      times the host's own mass, and lies within
      ``pairing_max_separation`` -- a symmetric "roughly comparable mass,
      roughly Local-Group-analogue separation" criterion. ``PairedHostID``
      is that halo's native ``halo_id`` (from ``Epoch.halos``), or -1 if
      none. If more than one halo satisfies the criterion, the nearest is
      stored and a warning is logged (the schema field is scalar, so it
      cannot represent an ambiguous multi-way pairing).
    - ``N_satellites_total``: count of *other* halos within
      ``isolation_radius_factor * R200c_host`` with mass >=
      ``completeness_mass_threshold`` -- independent of whichever
      satellite subset the caller passed to ``HaloExtractStage`` (that is
      a user selection of *which* satellites get full per-object
      processing, not necessarily the full completeness-limited
      population this field reports).

    Units caveat (same as ``HaloExtractStage``): computed in whatever
    length/mass units the Epoch's ``HaloCatalogue`` was configured with,
    not necessarily schema.py's declared units.
    """

    name = "host_environment"
    inputs = ("Haloes/HostHaloID", "Haloes/M200c", "Haloes/R200c",
             "Haloes/Position")
    outputs = ("Haloes/IsIsolated", "Haloes/IsPaired",
              "Haloes/PairedHostID", "Haloes/N_satellites_total")

    def __init__(self, epoch, host_row: int, isolation_radius_factor: float,
                 pairing_mass_ratio_min: float, pairing_max_separation: float,
                 completeness_mass_threshold: float):
        self.epoch = epoch
        self.host_row = int(host_row)
        self.isolation_radius_factor = float(isolation_radius_factor)
        self.pairing_mass_ratio_min = float(pairing_mass_ratio_min)
        self.pairing_max_separation = float(pairing_max_separation)
        self.completeness_mass_threshold = float(completeness_mass_threshold)

    def run(self, context):
        halos = self.epoch.halos
        if halos is None:
            raise RuntimeError(
                f"Stage '{self.name}': this Epoch has no halo catalogue.")

        mass = np.asarray(halos["mass"])
        pos = np.asarray(halos["pos"])
        halo_id = np.asarray(halos["halo_id"])
        boxsize = self.epoch.boxsize

        host_mass = float(np.asarray(context.columns["Haloes/M200c"])[0])
        host_r200c = float(np.asarray(context.columns["Haloes/R200c"])[0])
        host_pos = np.asarray(context.columns["Haloes/Position"])[0]

        others = np.ones(len(mass), dtype=bool)
        others[self.host_row] = False

        dpos = periodic_delta(pos - host_pos, boxsize)
        dist = np.linalg.norm(dpos, axis=1)

        isolation_radius = self.isolation_radius_factor * host_r200c
        within_isolation = others & (dist <= isolation_radius)

        is_isolated = not bool(np.any(within_isolation & (mass > host_mass)))
        n_satellites_total = int(np.count_nonzero(
            within_isolation & (mass >= self.completeness_mass_threshold)))

        ratio_lo = self.pairing_mass_ratio_min
        ratio_hi = (1.0 / self.pairing_mass_ratio_min
                   if self.pairing_mass_ratio_min > 0 else np.inf)
        mass_ratio = mass / host_mass
        pairing_candidates = (others & (dist <= self.pairing_max_separation)
                              & (mass_ratio >= ratio_lo)
                              & (mass_ratio <= ratio_hi))
        n_candidates = int(np.count_nonzero(pairing_candidates))

        is_paired = n_candidates > 0
        paired_host_id = -1
        if is_paired:
            cand_idx = np.flatnonzero(pairing_candidates)
            nearest = cand_idx[np.argmin(dist[cand_idx])]
            paired_host_id = int(halo_id[nearest])
            if n_candidates > 1:
                logger.warning(
                    "%s: %d halos satisfy the pairing criterion for this "
                    "host; PairedHostID stores only the nearest (id=%d) "
                    "since the schema field is scalar.",
                    self.name, n_candidates, paired_host_id)

        context.columns["Haloes/IsIsolated"] = np.asarray([is_isolated])
        context.columns["Haloes/IsPaired"] = np.asarray([is_paired])
        context.columns["Haloes/PairedHostID"] = np.asarray(
            [paired_host_id], dtype=np.int64)
        context.columns["Haloes/N_satellites_total"] = np.asarray(
            [n_satellites_total], dtype=np.int32)

        context.record_stage(self.name, is_isolated=is_isolated,
                             is_paired=is_paired,
                             n_satellites_total=n_satellites_total,
                             n_pairing_candidates=n_candidates)
        return context


class ObservabilityStage(PipelineStage):
    """``HeliocentricDistance`` and ``GalactocentricRadius`` -- the only
    two ``Satellites/Observability/*`` fields computable without either an
    assumed galaxy luminosity (``SharkGalaxyBackend`` doesn't compute
    ``LuminosityV``/``Luminosity_ugriz`` -- see its docstring) or a
    Rubin/LSST survey selection-function model (footprint, depth maps --
    external data this pipeline doesn't have).

    ``HeliocentricDistance`` needs a defined "mock Sun" position -- there
    is no natural one for a simulated system, so ``observer_pos`` is a
    required constructor argument (absolute box coordinates), not guessed.
    ``GalactocentricRadius`` needs no new parameter: it's the distance from
    the host's own z=0 position, already in context via
    ``HaloExtractStage``.

    Units caveat (same as ``HaloExtractStage``/``EnvironmentStage``):
    computed in whatever length units the Epoch's ``HaloCatalogue`` was
    configured with, not necessarily schema.py's declared "kpc" -- unit
    reconciliation isn't implemented anywhere in the pipeline yet.

    Deliberately **not** computed here, documented rather than guessed:

    - ``GalacticLongitude``/``GalacticLatitude``/``RA``/``Dec``: need an
      arbitrary "mock sky orientation" convention (which direction in the
      box is Galactic north, and how that maps onto equatorial
      coordinates) -- nothing in this codebase or the design doc pins one
      down, and RA/Dec in particular conventionally implies alignment
      with the *real* sky, which doesn't apply to synthetic data without
      yet another undocumented convention.
    - ``ApparentMagnitude_ugriz``/``SurfaceBrightness_V``: need
      ``Luminosity_ugriz``/``LuminosityV``/``HalfLightRadius``, none of
      which ``SharkGalaxyBackend`` computes yet.
    - ``RubinDetectionProbability``/``RubinCoaddDepth_r``/
      ``CompletenessWeight``/``SurveyFootprintFlag``: need an actual
      Rubin/LSST selection-function model (footprint polygon, coadd depth
      maps) -- external survey data this pipeline has no access to, a
      different kind of input entirely from anything computed so far.
    """

    name = "observability"
    inputs = ("Satellites/_internal/pos_z0", "Haloes/Position")
    outputs = ("Satellites/Observability/HeliocentricDistance",
               "Satellites/Observability/GalactocentricRadius")

    def __init__(self, epoch, observer_pos):
        self.epoch = epoch
        self.observer_pos = np.asarray(observer_pos, dtype=float)

    def run(self, context):
        satellite_pos = np.asarray(
            context.columns["Satellites/_internal/pos_z0"])
        host_pos = np.asarray(context.columns["Haloes/Position"])[0]
        boxsize = self.epoch.boxsize

        d_host = periodic_delta(satellite_pos - host_pos, boxsize)
        galactocentric_radius = np.linalg.norm(d_host, axis=1)

        d_obs = periodic_delta(satellite_pos - self.observer_pos, boxsize)
        heliocentric_distance = np.linalg.norm(d_obs, axis=1)

        context.columns["Satellites/Observability/HeliocentricDistance"] = \
            heliocentric_distance
        context.columns["Satellites/Observability/GalactocentricRadius"] = \
            galactocentric_radius

        context.record_stage(self.name, n_satellites=len(satellite_pos))
        return context


class DorchaSpecificStage(PipelineStage):
    """``EarliestProgenitorRedshift`` and ``FossilFraction`` -- the only
    two ``Satellites/DorchaProperties/*`` fields computable without
    infrastructure this codebase doesn't have at all (particle tagging,
    a progenitor density-field analysis, a cosmic-web classifier) or a
    common-grid machinery this stage doesn't build (see the deferred list
    below).

    ``EarliestProgenitorRedshift`` needs no new parameter: the tree walk
    (``TreeExtractStage``) already stops at the highest redshift a
    progenitor could be resolved to, so each satellite's
    ``MergerTrees/main_branch`` track's *first* entry (earliest time) is
    already the answer.

    ``FossilFraction`` ("fraction of stellar mass formed before
    reionisation") needs ``Satellites/GalaxyProperties/SFH`` (run
    ``StarFormationHistoryStage`` first -- this stage reads its
    ``context.meta["time_bin_edges_sfh"]`` rather than taking a second,
    independently-specified grid that could mismatch it) plus a
    reionisation epoch. ``reionisation_lookback_time`` (Gyr) is a
    required constructor argument, not an assumed redshift (commonly
    z~6-10 depending on definition/model) -- picking one silently would
    be exactly the kind of unstated modelling assumption this pipeline
    has avoided everywhere else (``EnvironmentStage``'s
    ``mass_threshold``, ``StarFormationHistoryStage``'s
    ``quenched_ssfr_threshold``, ...).

    Deliberately **not** computed here, documented rather than guessed:

    - ``ProgenitorParticleFraction``: needs particle tagging
      (cross-snapshot particle-ID tracking to identify which progenitor
      halo each z=0 bound particle came from) -- genuinely new
      infrastructure, nothing in this codebase does this.
    - ``PeakOverdensity``: needs a progenitor density-field analysis (a
      density field evaluated along each progenitor's trajectory) --
      also new infrastructure.
    - ``FormationEnvironmentClass``: needs an actual cosmic-web
      classifier (T-web/V-web or similar) -- the same gap as
      ``EnvironmentStage``'s ``CosmicWebClass``.
    - ``NumberOfMergers``: needs the full progenitor list at each
      snapshot (which *other* branches merged into the main branch), not
      just the main-branch walk this pipeline currently does --
      ``treeio_subfind.py``/``treeio_treefrog.py`` don't parse a
      NextProgenitor-equivalent field at all, so this is a genuine
      missing-infrastructure gap in the *reader* layer, not just an
      unstated modelling choice at this stage.
    - ``MassAccretionRateDM``: schema declares this on a fixed
      ``[N, N_snap]`` ``Snapshots/``-indexed grid shared by every
      satellite, but nothing establishes that common grid anywhere in
      this pipeline (satellites resolve back to different earliest
      snapshots); doing this properly needs the same kind of common-grid
      machinery ``StarFormationHistoryStage`` built for SFH, applied to
      mass instead of SFR -- a reasonable follow-up, not implemented
      here.
    """

    name = "dorcha_specific"
    inputs = ("MergerTrees/main_branch",)
    outputs = ("Satellites/DorchaProperties/EarliestProgenitorRedshift",
               "Satellites/DorchaProperties/FossilFraction")

    def __init__(self, reionisation_lookback_time: float):
        self.reionisation_lookback_time = float(reionisation_lookback_time)

    def run(self, context):
        main_branch = context.columns["MergerTrees/main_branch"]
        n = len(main_branch)

        earliest_z = np.full(n, np.nan)
        n_no_track = 0
        for i, track_ds in enumerate(main_branch):
            if track_ds is None or len(track_ds.track) == 0:
                n_no_track += 1
                continue
            earliest_z[i] = float(track_ds.track.Redshift[0])

        if n_no_track:
            logger.warning(
                "%s: %d/%d satellites have no resolved tree entry; "
                "EarliestProgenitorRedshift left NaN for them.",
                self.name, n_no_track, n)

        fossil_fraction = np.full(n, np.nan)
        sfh = context.columns.get("Satellites/GalaxyProperties/SFH")
        time_bin_edges = context.meta.get("time_bin_edges_sfh")

        if sfh is None:
            pass  # StarFormationHistoryStage hasn't run -- leave NaN
        elif time_bin_edges is None:
            logger.warning(
                "%s: Satellites/GalaxyProperties/SFH is present but "
                "context.meta['time_bin_edges_sfh'] isn't (run "
                "StarFormationHistoryStage, which sets it, first) -- "
                "FossilFraction left NaN for every satellite.", self.name)
        else:
            bin_widths_yr = np.diff(time_bin_edges) * 1.0e9
            # fraction of each bin's width at lookback >= reionisation
            # time (i.e. formed *before* reionisation), by fractional
            # overlap so a bin straddling the boundary isn't just snapped
            # to its nearest edge -- same technique as backends.rebin_sfh.
            overlap_lo = np.maximum(time_bin_edges[:-1],
                                    self.reionisation_lookback_time)
            overlap_hi = time_bin_edges[1:]
            bin_widths = np.diff(time_bin_edges)
            pre_reion_overlap = np.clip(overlap_hi - overlap_lo, 0.0, None)
            pre_reion_frac = np.divide(
                pre_reion_overlap, bin_widths,
                out=np.zeros_like(bin_widths), where=bin_widths > 0)

            for i in range(n):
                if np.any(np.isnan(sfh[i])):
                    continue
                formed_mass = sfh[i] * bin_widths_yr
                total_mass = float(np.sum(formed_mass))
                if total_mass > 0:
                    pre_reion_mass = float(np.sum(formed_mass * pre_reion_frac))
                    fossil_fraction[i] = pre_reion_mass / total_mass

        context.columns[
            "Satellites/DorchaProperties/EarliestProgenitorRedshift"] = \
            earliest_z
        context.columns["Satellites/DorchaProperties/FossilFraction"] = \
            fossil_fraction

        context.record_stage(self.name, n_satellites=n,
                             n_no_track=n_no_track)
        return context


# Registry used by CatalogueBuilder to resolve a config's `derived_stages`
# name list to classes.
STAGES = {
    "halo_properties": HaloPropertiesStage,
    "orbital_properties": OrbitalPropertiesStage,
    "star_formation_history": StarFormationHistoryStage,
    "host_environment": HostEnvironmentStage,
    "environment": EnvironmentStage,
    "observability": ObservabilityStage,
    "dorcha_specific": DorchaSpecificStage,
}
