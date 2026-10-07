#!/usr/bin/env python3
"""
make_halo_cube_backward.py
---------------------------
Inverse of make_halo_cube.py: select particles in a final snapshot within
some radius, then tag them with their SUBFIND subhalo membership from an
ARBITRARY earlier epoch via ParticleID matching.

Use case: find which historical subhalo members have drifted/been stripped
away by the present day. Supports selecting particles from ANY snapshot
(not just spatially adjacent ones) to track halo evolution across multiple
epochs.

Pipeline
--------
1. Read the REFERENCE catalogue at an ARBITRARY earlier snapshot (source of
   truth for subhalo membership).
2. Read the REFERENCE snapshot at that same earlier time.
3. Select particles within `r_scale * halo_radius` of the halo centre
   (spatial cut at the reference epoch).
4. Extract their ParticleIDs and reference-epoch subhalo membership.
5. Read the FINAL snapshot (at any later time).
6. Match all particles in the final snapshot against the reference-epoch IDs.
7. Tag matched particles with their reference-epoch subhalo index (-1 if not
   in the reference set).
8. Write the final snapshot with:
   - `groupid`: reference-epoch subhalo ID.
   - Metadata: halo_centre_ref, halo_systemic_velocity_ref, halo_extent_ref
     (the spatial selection at the reference time).
   - Metadata: reference_snapnum, final_snapnum (for tracking provenance).
"""
from __future__ import annotations

import argparse

import numpy as np
import h5py

from analysistools import SnapshotTools
from analysistools.snapshot_tools import select_particles
from analysistools.halo_tools import HaloTools

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------

# REFERENCE epoch: where we select spatially and get the subhalo tags.
# This can be ANY snapshot in the simulation (not necessarily the earliest).
DEFAULT_REFERENCE_SNAPNUM = 50
DEFAULT_REFERENCE_SNAPSHOT = "/Volumes/ChrisPowerHardDrive1/DorchaShark/DORCHA_01/output/snapshot_050.hdf5"
DEFAULT_REFERENCE_CATALOGUE = "/Volumes/ChrisPowerHardDrive1/DorchaShark/DORCHA_01/output/fof_subhalo_tab_050.hdf5"

# FINAL epoch: where we apply the tags to all particles.
# Can be earlier or later than the reference epoch.
DEFAULT_FINAL_SNAPNUM = 122
DEFAULT_FINAL_SNAPSHOT = "/Volumes/ChrisPowerHardDrive1/DorchaShark/DORCHA_01/output/snapshot_122.hdf5"

# Halo selection (at the REFERENCE time).
DEFAULT_HALO_IDS = [0]
DEFAULT_R_SCALE = 1.0

# Output
DEFAULT_OUTFILE = None  # Default: 'snap_<final-snapnum>.tagged_from_<ref-snapnum>.cube'

DEFAULT_RADIUS_FIELD = "radius"
DEFAULT_COMOVING = True
DEFAULT_LITTLE_H = True
DEFAULT_CENTRE_ON_SUBHALO = True


def build_dm_subhaloid_array(catalogue_file: str, n_dm: int,
                              dm_type: int = 1) -> np.ndarray:
    """Per-particle SUBFIND subhalo number for a DM-only particle array of
    length `n_dm`, aligned with SubhaloOffsetType/SubhaloLenType[:, dm_type].
    -1 for particles not bound to any subhalo.
    
    (Copied from make_halo_cube.py)
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
    """Slice snap.pos/vel/mass/pids/ptype down to just PartType{dm_type}.
    (Copied from make_halo_cube.py)
    """
    n0 = int(np.sum(snap.num_part_total[:dm_type]))
    n1 = n0 + int(snap.num_part_total[dm_type])
    for attr in ("pos", "vel", "mass", "pids", "ptype"):
        setattr(snap, attr, getattr(snap, attr)[n0:n1])


def read_dm_only_snapshot(snapshot_file: str) -> SnapshotTools:
    """Read a snapshot and return a SnapshotTools instance restricted to
    just its PartType1 particles.
    (Copied from make_halo_cube.py)
    """
    snap = SnapshotTools(
        snapfileformat="HDF5",
        convention="GADGET4",
        hires_only=True,
        not_hires_ptypes=[2, 3, 5, 7],
    )
    _data = snap.read_snapshot(snapshot_file)
    for _attr in ("pos", "vel", "mass", "pids", "ptype", "num_part_total",
                  "box_size", "scale_factor", "omega_0", "omega_lambda",
                  "hubble_param"):
        setattr(snap, _attr, getattr(_data, _attr))
    restrict_to_dm_only(snap, dm_type=1)
    return snap


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Select particles in a FINAL snapshot, then tag them "
                     "with their SUBFIND subhalo membership from an ARBITRARY "
                     "REFERENCE epoch via ParticleID matching. The reference "
                     "can be earlier or later than the final epoch. "
                     "Use to track halo evolution and see where historical "
                     "subhalo members end up, or where final-epoch particles "
                     "came from.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    ref = p.add_argument_group("reference epoch (spatial selection + subhalo tags)")
    ref.add_argument("--reference-snapshot", default=DEFAULT_REFERENCE_SNAPSHOT,
                        help="Particle snapshot at the reference epoch.")
    ref.add_argument("--reference-catalogue", default=DEFAULT_REFERENCE_CATALOGUE,
                        help="Matching FOF/SUBFIND catalogue for the reference snapshot.")
    ref.add_argument("--reference-snapnum", type=int, default=DEFAULT_REFERENCE_SNAPNUM,
                        help="Snapshot number of the reference epoch (metadata only).")

    final = p.add_argument_group("final epoch (particles to tag)")
    final.add_argument("--final-snapshot", default=DEFAULT_FINAL_SNAPSHOT,
                        help="Snapshot at the final time to tag with reference-epoch "
                             "subhalo IDs.")
    final.add_argument("--final-snapnum", type=int, default=DEFAULT_FINAL_SNAPNUM,
                        help="Snapshot number of the final epoch (metadata only).")

    halo = p.add_argument_group("halo selection (at reference time)")
    halo.add_argument("--halo-id", type=int, nargs="+", default=DEFAULT_HALO_IDS,
                       help="Row index/indices in the reference catalogue (0-based). "
                            "Pass multiple for several halos in one run.")
    halo.add_argument("--r-scale", type=float, default=DEFAULT_R_SCALE,
                       help="Radius multiplier applied to the halo's "
                            "catalogue radius for the spatial cut.")
    halo.add_argument("--outfile", default=None,
                       help="Base output filename (no .hdf5). "
                        "Default: 'snap_<final-snapnum>.tagged_from_<ref-snapnum>.cube'")

    misc = p.add_argument_group("catalogue reading")
    misc.add_argument("--radius-field", default=DEFAULT_RADIUS_FIELD,
                       help="Standardised catalogue field for halo radius.")
    misc.add_argument("--comoving", action=argparse.BooleanOptionalAction,
                       default=DEFAULT_COMOVING,
                       help="Keep catalogue positions comoving.")
    misc.add_argument("--little-h", action=argparse.BooleanOptionalAction,
                       default=DEFAULT_LITTLE_H,
                       help="Keep h in catalogue units (Mpc/h).")
    misc.add_argument("--centre-on-subhalo", action=argparse.BooleanOptionalAction,
                       default=DEFAULT_CENTRE_ON_SUBHALO,
                       help="Centre on primary subhalo rather than raw FOF group.")

    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    n_halos = len(args.halo_id)
    base_outfile = args.outfile or f"snap_{args.final_snapnum}.tagged_from_{args.reference_snapnum}.cube"

    # =====================================================================
    # REFERENCE EPOCH: read catalogue, snapshot, select spatially, get subhalo tags
    # =====================================================================
    print(f"\n{'='*70}")
    print(f"REFERENCE EPOCH (snapnum={args.reference_snapnum})")
    print(f"{'='*70}")

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
    snap_ref = read_dm_only_snapshot(args.reference_snapshot)

    print("Building per-particle subhalo-id array...")
    snap_ref.groupid = build_dm_subhaloid_array(
        args.reference_catalogue, len(snap_ref.pos), dm_type=1)

    # Collect particle IDs and tags from all selected halos at the reference time.
    all_target_ids = []
    all_target_groupids = []
    halo_metadata = []

    for row in args.halo_id:
        centre = np.asarray(halos["pos"])[row]
        velocity = np.asarray(halos["vel"])[row]
        radius = float(np.asarray(halos[args.radius_field])[row])
        extent = radius * args.r_scale

        print(f"\nHalo row {row}: centre={centre}, radius={radius:.6g}, "
              f"extent={extent:.6g}")

        idx = select_particles(
            snap_ref.pos, centre,
            size=extent,
            geometry="spherical",
            periodic=True,
            scale_length=snap_ref.box_size,
        )
        print(f"  Selected {len(idx)} particles at reference time within r < {extent:.6g}")

        # Extract IDs and subhalo tags.
        target_ids = snap_ref.pids[idx]
        target_groupids = snap_ref.groupid[idx]
        all_target_ids.append(target_ids)
        all_target_groupids.append(target_groupids)

        halo_metadata.append({
            "row": row,
            "centre": centre,
            "velocity": velocity,
            "radius": radius,
            "extent": extent,
        })

    # Concatenate all halos' IDs and tags.
    all_target_ids = np.concatenate(all_target_ids)
    all_target_groupids = np.concatenate(all_target_groupids)
    print(f"\nTotal reference-epoch particles from {n_halos} halo(s): {len(all_target_ids)}")

    # =====================================================================
    # FINAL EPOCH: read snapshot, match IDs, tag with reference-time subhalo IDs
    # =====================================================================
    print(f"\n{'='*70}")
    print(f"FINAL EPOCH (snapnum={args.final_snapnum})")
    print(f"{'='*70}")

    print(f"Reading snapshot: {args.final_snapshot}")
    snap_final = read_dm_only_snapshot(args.final_snapshot)

    # ID -> reference-epoch groupid lookup.
    print("Building ID lookup table for reference-epoch subhalo tags...")
    order = np.argsort(all_target_ids)
    sorted_ids = all_target_ids[order]
    sorted_groupids = all_target_groupids[order]

    # For every particle in the final snapshot, check if it was in the
    # reference set. If yes, assign its reference-time subhalo tag; else -1.
    full_groupid = np.full(len(snap_final.pos), -1, dtype=np.int32)
    idx_final = np.flatnonzero(np.isin(snap_final.pids, all_target_ids))
    lookup_pos = np.searchsorted(sorted_ids, snap_final.pids[idx_final])
    full_groupid[idx_final] = sorted_groupids[lookup_pos]

    n_tagged = np.sum(full_groupid != -1)
    print(f"Tagged {n_tagged} of {len(snap_final.pos)} particles as "
          f"members of reference-epoch halo(s)")
    print(f"  Lost/stripped: {len(all_target_ids) - n_tagged} of "
          f"{len(all_target_ids)} reference particles not found in final snapshot")
    if len(all_target_ids) > 0:
        pct_retained = 100.0 * n_tagged / len(all_target_ids)
        print(f"  Retention rate: {pct_retained:.1f}%")

    snap_final.groupid = full_groupid

    # =====================================================================
    # Write output
    # =====================================================================
    print(f"\nWriting output...")
    
    # Write one output per halo (or a single output if only one halo).
    for halo_idx, meta in enumerate(halo_metadata):
        halo_suffix = f".halo{meta['row']}" if n_halos > 1 else ""
        outfile = f"{base_outfile}{halo_suffix}"

        snap_final.write_snapshot(
            filename=outfile,
            convention="AREPO",
            blocks_to_write=["pos", "vel", "pids", "mass", "groupid"],
            # Reference-epoch halo info (where particles were selected)
            halo_centre_ref=meta["centre"],
            halo_systemic_velocity_ref=meta["velocity"],
            halo_extent_ref=meta["extent"],
            # Metadata about the tagging
            reference_snapnum=args.reference_snapnum,
            final_snapnum=args.final_snapnum,
        )
        print(f"  Wrote {outfile}.hdf5")


if __name__ == "__main__":
    main()
