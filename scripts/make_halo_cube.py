#!/usr/bin/env python3
"""
make_halo_cube.py
------------------
Extract a spherical, halo-centred cutout of dark matter particles from
a full AREPO/SUBFIND snapshot, tag each selected particle with its
SUBFIND subhalo membership, and write the result as a standalone HDF5
snapshot.

Pipeline
--------
1. Read the full snapshot (PartType1 only) with SnapshotTools.
2. Read the FOF/SUBFIND catalogue with HaloTools to get the target
   halo's centre, radius and bulk velocity.
3. Select particles within `r_scale * halo_radius` of the halo centre
   (periodic, spherical) -- this is a radius query, NOT a match on
   existing group membership, since the raw snapshot has no per-particle
   group id block.
4. Recover each particle's SUBFIND subhalo number from
   SubhaloOffsetType/SubhaloLenType -- NOT via /IDs/ParticleIDs (this
   run's catalogues carry an empty /IDs group despite Header/Nids_Total
   being nonzero; checked across fof_subhalo_tab_000/050/121/122, all
   empty), and NOT the parent FOF GroupOffsetType/GroupLenType either:
   subhalo, not FOF group, is what has continuity across snapshots for
   this run (subhalo_desc/prog/treelink_NNN.hdf5 alongside each
   fof_subhalo_tab_NNN.hdf5 build merger trees at the subhalo level), so
   tagging by SubhaloID is what actually supports tracking a structure
   backwards/forwards through time. Subhalo{Offset,Len}Type index
   directly into the snapshot's own on-disk PartType1 order, same
   convention as Group{Offset,Len}Type (identical (Nsubhalos/Ngroups, 6)
   shape and semantics, just one level finer): AREPO writes bound
   particles first, contiguously by subhalo, then unbound/field
   particles after, per species -- no ID lookup needed at all. See
   build_dm_subhaloid_array()'s docstring for the output dataset-naming
   caveat (still called "GroupIDs" on disk).
5. Write the cutout with SnapshotTools.write_snapshot(). halo_centre /
   halo_systemic_velocity / halo_extent are passed as plain kwargs,
   which write_hdf5() stashes in self.metadata and dumps verbatim into
   Header.attrs under those exact (lowercase) names -- confirmed against
   the original snap_050.cube.hdf5 and by a synthetic round-trip test.
   (There's a second, dead code path in snapio_hdf5.py:505-512 gated on
   a `selection` kwarg that writes "Halo Centre" etc. title-case instead
   -- do NOT pass selection=True, or you'll get both sets of attrs.)
"""
from __future__ import annotations

import argparse

import numpy as np
import h5py

from analysistools import SnapshotTools
from analysistools.snapshot_tools import select_particles
from analysistools.halo_tools import HaloTools

# ---------------------------------------------------------------------------
# Defaults -- used both as the CLI flags' defaults (see parse_args() below)
# and as what you get running this with no arguments at all. Edit these for
# your own interactive/dev use; override any of them per-invocation via the
# command line for batch/array jobs (see `python make_halo_cube.py --help`).
# ---------------------------------------------------------------------------

# The REFERENCE epoch is where the halo is actually found: HaloTools reads
# the catalogue to get the target halo's centre/radius/velocity, and
# select_particles() does the one and only spatial (radius) cut, against
# the reference snapshot. Every other epoch (see OTHER_EPOCHS/--other-*
# below) is reached purely by particle-ID membership from there -- so the
# reference can be your z=0 snapshot (tracing back to z=2) or your z=2
# snapshot (tracing forward to z=0) equally well; "earlier" vs "later"
# doesn't matter to the ID-matching logic, only to how you've populated
# the other epochs.
DEFAULT_REFERENCE_SNAPNUM = 122
DEFAULT_REFERENCE_SNAPSHOT = "/Volumes/ChrisPowerHardDrive1/DorchaShark/DORCHA_01/output/snapshot_122.hdf5"
DEFAULT_REFERENCE_CATALOGUE = "/Volumes/ChrisPowerHardDrive1/DorchaShark/DORCHA_01/output/fof_subhalo_tab_122.hdf5"

# `halo_id` is the row index into the reference's standardised catalogue
# (0-based); look it up by mass/position first if you don't already know
# it (e.g. np.argsort(halos["mass"])[::-1] for the most massive halos).
DEFAULT_HALO_IDS = [0]
DEFAULT_R_SCALE = 1.0

# For each halo above, also trace its particle IDs (fixed at the
# reference's spatial selection) through these other epochs -- an
# ID-membership query, NOT a spatial one (the particles will have drifted
# from wherever they were at the reference time, so select_particles'
# radius cut doesn't apply at all here).
#
# GroupIDs written for these traced particles are the REFERENCE epoch's
# own subhalo tags (looked up by ParticleID) -- i.e. "which subhalo did
# this particle belong to at the reference time", carried onto wherever
# it sits in this other snapshot. NOT that epoch's own, independently-
# derived subhalo membership -- there is no per-epoch catalogue read at
# all, since the reference groupid array already has everything needed
# (it's a free ID lookup, not a second catalogue read).
#
# Leave --other-snapnums empty (the default) to skip this stage and only
# write the reference-halo cube.
DEFAULT_OTHER_SNAPNUMS: list[int] = [50, 100]
DEFAULT_OTHER_DIR = "/Users/00075868/Currentwork/Dorcha"
DEFAULT_OTHER_SNAPSHOT_PATTERN = "snapshot_{snapnum:03d}.hdf5"
DEFAULT_OTHER_OUTFILE_PATTERN = "snap_{snapnum:03d}.cube"

DEFAULT_RADIUS_FIELD = "radius"  # standardised name; maps to Group_R_Crit200
                                  # by default via halo_tools_standardise_names
                                  # -- use the raw catalogue field directly
                                  # (standardise=False) if you need
                                  # Group_R_Mean200/TopHat200/Crit500 instead.

DEFAULT_COMOVING = True
# little_h=True (native, keep h in the units -- i.e. Mpc/h) is required
# here, NOT False. SnapshotTools.read_snapshot() applies no unit
# conversion at all to the raw Coordinates -- they stay in whatever the
# simulation wrote (comoving Mpc/h). little_h=False divides the
# catalogue's pos/radius/BoxSize by HubbleParam (native_includes_h=True
# is SUBFIND's default -- see halo_tools.py's NATIVE_INCLUDES_LITTLE_H),
# converting THEM to h-free Mpc while the snapshot stays in Mpc/h -- a
# real unit mismatch, confirmed 2026-09: with little_h=False, halo row 0
# centred at [82.86, 63.95, 72.47] and selected 0 particles; multiplying
# by HubbleParam=0.678 gives [56.18, 43.36, 49.13], matching the original
# snap_050.cube.hdf5's halo_centre almost exactly -- i.e. little_h=True
# (no conversion) is what the original notebook's units actually were.
DEFAULT_LITTLE_H = True
DEFAULT_CENTRE_ON_SUBHALO = True


def other_epochs_from_snapnums(
    snapnums: list[int],
    snapshot_dir: str,
    snapshot_pattern: str = DEFAULT_OTHER_SNAPSHOT_PATTERN,
    outfile_pattern: str = DEFAULT_OTHER_OUTFILE_PATTERN,
) -> list[dict]:
    """Build an OTHER_EPOCHS list from a list of snapshot numbers plus a
    single directory/filename convention -- for the common case (a run
    where every other-epoch snapshot lives in the same place, numbered
    consistently). Each pattern is .format()-ed with `snapnum=`; override
    them if your naming differs (e.g. "snap_{snapnum:04d}.hdf5" for
    4-digit padding, or drop the padding: "snapshot_{snapnum}.hdf5").
    """
    return [
        dict(
            snapnum=n,
            snapshot_file=f"{snapshot_dir.rstrip('/')}/{snapshot_pattern.format(snapnum=n)}",
            outfile=outfile_pattern.format(snapnum=n),
        )
        for n in snapnums
    ]


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Extract a spherical, halo-centred DM cutout from a "
                     "reference snapshot, tag it with SUBFIND subhalo "
                     "membership, and (optionally) trace the same "
                     "particle IDs through other snapshots/epochs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ref = p.add_argument_group("reference epoch")
    ref.add_argument("--reference-snapshot", default=DEFAULT_REFERENCE_SNAPSHOT,
                      help="Full particle snapshot to select the halo from.")
    ref.add_argument("--reference-catalogue", default=DEFAULT_REFERENCE_CATALOGUE,
                      help="Matching FOF/SUBFIND catalogue for the reference snapshot.")
    ref.add_argument("--reference-snapnum", type=int, default=DEFAULT_REFERENCE_SNAPNUM,
                      help="Snapshot number of the reference epoch (metadata only).")

    halo = p.add_argument_group("halo selection")
    halo.add_argument("--halo-id", type=int, nargs="+", default=DEFAULT_HALO_IDS,
                       help="Row index/indices into the reference catalogue "
                            "(0-based). Pass multiple for several halos in "
                            "one run, e.g. --halo-id 0 3 17.")
    halo.add_argument("--r-scale", type=float, default=DEFAULT_R_SCALE,
                       help="Radius multiplier applied to each halo's "
                            "catalogue radius for the spatial cut.")
    halo.add_argument("--outfile", default=None,
                       help="Base output filename for the reference-halo "
                            "cube (no .hdf5 -- the writer appends it). "
                            "Default: 'snap_<reference-snapnum>.cube', "
                            "with a '.halo<N>' suffix per halo if "
                            "--halo-id has more than one value.")

    other = p.add_argument_group("other (traced) epochs")
    other.add_argument("--other-snapnums", type=int, nargs="*",
                        default=DEFAULT_OTHER_SNAPNUMS,
                        help="Snapshot numbers to trace the reference "
                             "halo's particle IDs through. Empty (pass "
                             "--other-snapnums with nothing after it) to "
                             "skip this stage entirely.")
    other.add_argument("--other-dir", default=DEFAULT_OTHER_DIR,
                        help="Directory containing the other-epoch "
                             "snapshots (all must follow the same naming "
                             "pattern -- see --other-snapshot-pattern).")
    other.add_argument("--other-snapshot-pattern", default=DEFAULT_OTHER_SNAPSHOT_PATTERN,
                        help="Filename pattern for other-epoch snapshots, "
                             "'.format()'-ed with snapnum=N.")
    other.add_argument("--other-outfile-pattern", default=DEFAULT_OTHER_OUTFILE_PATTERN,
                        help="Output-filename pattern (no .hdf5) for "
                             "other-epoch cubes, '.format()'-ed with "
                             "snapnum=N.")

    misc = p.add_argument_group("catalogue reading")
    misc.add_argument("--radius-field", default=DEFAULT_RADIUS_FIELD,
                       help="Standardised catalogue field to use as the "
                            "halo radius; 'radius' maps to Group_R_Crit200 "
                            "-- see halo_tools_standardise_names for the "
                            "full mapping if you need a different radius "
                            "definition.")
    misc.add_argument("--comoving", action=argparse.BooleanOptionalAction,
                       default=DEFAULT_COMOVING,
                       help="Keep catalogue positions comoving (no-op; "
                            "catalogues are always stored comoving).")
    misc.add_argument("--little-h", action=argparse.BooleanOptionalAction,
                       default=DEFAULT_LITTLE_H,
                       help="Keep h in the catalogue's length/mass units "
                            "(native Mpc/h) -- must match the raw "
                            "snapshot's own units, which are never "
                            "h-converted. Leave this True unless you know "
                            "your snapshot's Coordinates are h-free; see "
                            "the module docstring for how getting this "
                            "wrong silently selects zero particles.")
    misc.add_argument("--centre-on-subhalo", action=argparse.BooleanOptionalAction,
                       default=DEFAULT_CENTRE_ON_SUBHALO,
                       help="Centre each halo on its primary subhalo's "
                            "position/velocity rather than the raw FOF "
                            "GroupPos/GroupVel.")

    return p.parse_args(argv)


def build_dm_subhaloid_array(catalogue_file: str, n_dm: int,
                              dm_type: int = 1) -> np.ndarray:
    """Per-particle SUBFIND subhalo number for a DM-only particle array of
    length `n_dm` (already sliced down to just PartType{dm_type} -- see
    restrict_to_dm_only() below), aligned 1:1 with SubhaloOffsetType/
    SubhaloLenType[:, dm_type] -- same offset-into-snapshot-order
    convention as Group/GroupOffsetType (checked: identical shape,
    (Nsubhalos, 6), same semantics), just one level finer than FOF group.
    -1 for particles not bound to any subhalo.

    Subhalo, not FOF group, is what actually has continuity across
    snapshots here -- this run's merger trees (subhalo_desc_NNN.hdf5,
    subhalo_prog_NNN.hdf5, subhalo_treelink_NNN.hdf5, alongside each
    fof_subhalo_tab_NNN.hdf5) are built at the subhalo level, so tagging
    particles by SubhaloID rather than GroupID is what actually supports
    tracking a given (sub)structure backwards/forwards through time.

    Written to the output under the same 'GroupIDs' dataset name as
    before (SnapshotTools.write_snapshot's 'groupid' block is hardcoded
    to that name) -- the values are subhalo ids now, only the label on
    disk still says "Group". Rename write_groupids()/_write_groupids in
    snapio_hdf5.py if you want the dataset itself called "SubhaloIDs".
    """
    with h5py.File(catalogue_file, "r") as f:
        offset_type = f["Subhalo/SubhaloOffsetType"][:, dm_type]
        len_type = f["Subhalo/SubhaloLenType"][:, dm_type]

    subhaloid = np.full(n_dm, -1, dtype=np.int32)
    for s in range(len(offset_type)):
        o, l = int(offset_type[s]), int(len_type[s])
        if l:
            subhaloid[o:o + l] = s
    return subhaloid


def restrict_to_dm_only(snap: SnapshotTools, dm_type: int = 1) -> None:
    """Slice snap.pos/vel/mass/pids/ptype down to just the PartType{dm_type}
    range, in place.

    Required, not just tidy: SnapshotTools._allocate_memory() allocates
    pos/vel/mass/pids with np.ndarray(shape=...) -- uninitialized memory,
    not zero-filled -- sized for ALL particle types in the file's header,
    before hires_only/not_hires_ptypes is applied. Skipped types (e.g.
    PartType5 here) are then never written into their slice during the
    read loop, so snap.pos[<those rows>] is raw garbage, not zeros or an
    absent range. Anything operating on the full arrays (select_particles,
    the group-id builder) would silently search that garbage too. This
    must run before any of that -- confirmed against snapio_hdf5.py's
    _allocate_memory/_read_particle_data_single_file, which allocate by
    sum(num_part_total) but only fill non-skipped types' slices.
    """
    n0 = int(np.sum(snap.num_part_total[:dm_type]))
    n1 = n0 + int(snap.num_part_total[dm_type])
    for attr in ("pos", "vel", "mass", "pids", "ptype"):
        setattr(snap, attr, getattr(snap, attr)[n0:n1])


def read_dm_only_snapshot(snapshot_file: str) -> SnapshotTools:
    """Read a snapshot and return a SnapshotTools instance restricted to
    just its (valid, non-garbage) PartType1 particles. Factored out of
    main() so the earlier-snapshot ID-tracking stage can reuse the exact
    same read/merge/restrict sequence as the reference snapshot."""
    snap = SnapshotTools(
        snapfileformat="HDF5",
        convention="GADGET4",
        hires_only=True,
        not_hires_ptypes=[2, 3, 5, 7],   # add your low-res/boundary type(s)
    )
    # See main()'s comment on the same pattern: read_snapshot() returns a
    # separate SnapshotData object rather than mutating `snap`, and only
    # the actual data attributes should be copied over (not reader-side
    # config like name_of_mass_block -- see git history/earlier comments
    # for why that matters).
    _data = snap.read_snapshot(snapshot_file)
    for _attr in ("pos", "vel", "mass", "pids", "ptype", "num_part_total",
                  "box_size", "scale_factor", "omega_0", "omega_lambda",
                  "hubble_param"):
        setattr(snap, _attr, getattr(_data, _attr))
    restrict_to_dm_only(snap, dm_type=1)
    return snap


def extract_by_ids(snapshot_file: str, target_ids: np.ndarray,
                    target_groupids: np.ndarray,
                    outfile: str, provenance: dict) -> None:
    """Write a cube containing exactly the DM particles in `snapshot_file`
    whose ParticleID is in `target_ids` -- an ID-membership query, not a
    spatial one (see OTHER_EPOCHS above for why).

    `target_groupids` is REFERENCE's own per-particle subhalo array,
    already sliced down to exactly `target_ids` (i.e. target_groupids[i]
    is the subhalo target_ids[i] belonged to AT THE REFERENCE TIME) --
    it gets carried onto the matched particles here via an ID lookup, NOT
    recomputed from any catalogue local to this epoch. So the GroupIDs
    written are "which reference-time subhalo does this particle belong
    to", regardless of where/how it's currently clustered in this
    snapshot.
    """
    print(f"\nReading snapshot: {snapshot_file}")
    snap = read_dm_only_snapshot(snapshot_file)

    idx = np.flatnonzero(np.isin(snap.pids, target_ids))
    print(f"  Matched {len(idx)} of {len(target_ids)} target IDs "
          f"({len(snap.pids)} DM particles in this snapshot)")
    if len(idx) < len(target_ids):
        print(f"  NOTE: {len(target_ids) - len(idx)} target IDs not found "
              f"in this snapshot -- expected if the box changed particle "
              f"content between epochs, but double check if the shortfall "
              f"is large.")

    idx_type = snap.ptype[idx].astype(np.int64)

    # ID -> reference-time groupid lookup, applied only to the matched
    # rows. Every id in snap.pids[idx] is guaranteed (by construction of
    # idx via np.isin above) to be present in target_ids, so the
    # searchsorted lookup below can't miss.
    order = np.argsort(target_ids)
    sorted_ids = target_ids[order]
    sorted_groupids = target_groupids[order]
    lookup_pos = np.searchsorted(sorted_ids, snap.pids[idx])
    full_groupid = np.full(len(snap.pos), -1, dtype=np.int32)
    full_groupid[idx] = sorted_groupids[lookup_pos]
    snap.groupid = full_groupid

    snap.write_snapshot(
        filename=outfile,
        idx=idx,
        idx_type=idx_type,
        convention="AREPO",
        blocks_to_write=["pos", "vel", "pids", "mass", "groupid"],
        **provenance,
    )
    print(f"  Wrote {outfile}.hdf5")


def main(argv=None):
    args = parse_args(argv)

    other_epochs = other_epochs_from_snapnums(
        snapnums=args.other_snapnums,
        snapshot_dir=args.other_dir,
        snapshot_pattern=args.other_snapshot_pattern,
        outfile_pattern=args.other_outfile_pattern,
    )
    n_halos = len(args.halo_id)
    base_outfile = args.outfile or f"snap_{args.reference_snapnum}.cube"

    print(f"Reading halo catalogue: {args.reference_catalogue}")
    ht = HaloTools(
        comoving=args.comoving,
        little_h=args.little_h,
        centre_on_subhalo=args.centre_on_subhalo,
    )
    meta, halos, subhalos = ht.read_catalogue(
        filename=args.reference_catalogue,
        fileformat="SubFind",
        standardise=True,
        snapnum=args.reference_snapnum,
    )

    print(f"Reading snapshot: {args.reference_snapshot}")
    snap = read_dm_only_snapshot(args.reference_snapshot)

    print("Building per-particle subhalo-id array "
          "(one-off cost, reused for every halo below)...")
    snap.groupid = build_dm_subhaloid_array(
        args.reference_catalogue, len(snap.pos), dm_type=1)

    for row in args.halo_id:
        centre = np.asarray(halos["pos"])[row]
        velocity = np.asarray(halos["vel"])[row]
        radius = float(np.asarray(halos[args.radius_field])[row])
        extent = radius * args.r_scale

        print(f"\nHalo row {row}: centre={centre}, radius={radius:.6g}, "
              f"extent={extent:.6g}")

        idx = select_particles(
            snap.pos, centre,
            size=extent,
            geometry="spherical",
            periodic=True,
            scale_length=snap.box_size,
        )
        print(f"  Selected {len(idx)} particles within r < {extent:.6g}")

        idx_type = snap.ptype[idx].astype(np.int64)
        # snap.groupid (full-length, set once above) is sliced by idx
        # internally inside write_snapshot -- nothing more to do here.

        # Only disambiguate with a .halo{row} suffix when tracing more
        # than one halo in this run -- with a single halo (the common
        # case, including the default SLURM-array-job pattern of one
        # halo per job via --halo-id), base_outfile is used exactly as
        # given, e.g. "snap_122.cube" -> "snap_122.cube.hdf5".
        halo_suffix = f".halo{row}" if n_halos > 1 else ""
        outfile = f"{base_outfile}{halo_suffix}"  # no ".hdf5" -- the writer appends it
        snap.write_snapshot(
            filename=outfile,
            idx=idx,
            idx_type=idx_type,
            convention="AREPO",
            blocks_to_write=["pos", "vel", "pids", "mass", "groupid"],
            halo_centre=centre,
            halo_systemic_velocity=velocity,
            halo_extent=extent,
        )

        print(f"  Wrote {outfile}.hdf5")

        if other_epochs:
            target_ids = snap.pids[idx]
            # REFERENCE's own subhalo tags for exactly these particles --
            # this is what gets carried onto every other epoch below, not
            # anything recomputed locally there. snap.groupid is the
            # full-length array built once via build_dm_subhaloid_array()
            # above, sliced by idx the same way write_snapshot does.
            target_groupids = snap.groupid[idx]
            provenance = dict(
                halo_centre=centre,
                halo_systemic_velocity=velocity,
                halo_extent=extent,
            )
            for epoch in other_epochs:
                extract_by_ids(
                    epoch["snapshot_file"], target_ids, target_groupids,
                    outfile=f"{epoch['outfile']}{halo_suffix}",
                    provenance=provenance,
                )


if __name__ == "__main__":
    main()
