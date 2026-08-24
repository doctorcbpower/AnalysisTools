#!/usr/bin/env python3
"""
analysistools.catalogue.diagnostics
-------------------------------------
Sanity-check helpers for the raw merger-tree tracks and SHARK SFH that
feed the Dorcha catalogue pipeline -- report concrete numbers for a human
to judge (mass histories, tree-vs-Group-table agreement, SFH-vs-stellar-
mass ratios), rather than reducing everything to a validator's pass/fail.

Not part of ``CatalogueBuilder``'s pipeline -- meant to be run standalone
against real data, e.g. from a notebook, before trusting a build. See
``notebooks/DorchaCatalogueEndToEnd.ipynb``'s diagnostics section.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional

import numpy as np

from ..merger_tree_types import MergerTreeError


@dataclass
class TreeSanityReport:
    """One satellite's main-branch tree walk, checked for continuity and
    cross-checked against a direct ``epoch.halos`` lookup at the same row.
    """
    halo_row: int
    n_snapshots: int
    snapnum_min: int
    snapnum_max: int
    snapnum_gaps: List[int]           # snapshots the main branch skips over
    mass_min: float
    mass_max: float
    n_nonpositive_mass: int
    n_nan_mass: int
    infall_snapnum: Optional[int]     # first snapshot IsSubhalo goes True
    group_m200_z0_from_tree: Optional[float]
    group_r200_z0_from_tree: Optional[float]
    group_m200_z0_from_catalogue: float
    group_r200_z0_from_catalogue: float
    mass_agreement_ratio: Optional[float]     # tree / catalogue, want ~1
    radius_agreement_ratio: Optional[float]

    def print_summary(self) -> None:
        print(f"halo_row={self.halo_row}: {self.n_snapshots} snapshots "
             f"(SnapNum {self.snapnum_min}-{self.snapnum_max}), "
             f"{len(self.snapnum_gaps)} gap(s) {self.snapnum_gaps or ''}")
        print(f"  Mass: [{self.mass_min:.3e}, {self.mass_max:.3e}], "
             f"{self.n_nonpositive_mass} non-positive, "
             f"{self.n_nan_mass} NaN")
        print(f"  infall snapshot: {self.infall_snapnum}")
        if self.mass_agreement_ratio is not None:
            print(f"  M200c_z0 tree/catalogue ratio: "
                 f"{self.mass_agreement_ratio:.6f} "
                 f"(tree={self.group_m200_z0_from_tree:.4e}, "
                 f"catalogue={self.group_m200_z0_from_catalogue:.4e})")
        else:
            print(f"  M200c_z0: no tree GroupM200 -- catalogue value "
                 f"{self.group_m200_z0_from_catalogue:.4e} only")
        if self.radius_agreement_ratio is not None:
            print(f"  R200c_z0 tree/catalogue ratio: "
                 f"{self.radius_agreement_ratio:.6f} "
                 f"(tree={self.group_r200_z0_from_tree:.4e}, "
                 f"catalogue={self.group_r200_z0_from_catalogue:.4e})")
        else:
            print(f"  R200c_z0: no tree GroupR200 -- catalogue value "
                 f"{self.group_r200_z0_from_catalogue:.4e} only")


def tree_sanity_check(epoch, halo_row: int) -> Optional[TreeSanityReport]:
    """Walks ``epoch.track_of(index=halo_row)``'s main branch and reports
    concrete numbers to eyeball: SnapNum continuity, mass sign/finiteness,
    the infall snapshot, and whether the tree's own resolved
    ``Group_M_Crit200``/``Group_R_Crit200`` (the same values
    ``TreeExtractStage`` refines ``M200c_z0``/``R200c_z0`` from) agree
    with a direct lookup against ``epoch.halos`` at this row -- close to
    1.0 is the expected/healthy case; a large or systematic deviation
    would mean the tree resolved a *different* group than the row this
    diagnostic was asked to check. `None` if this halo has no resolvable
    tree entry."""
    try:
        track_ds = epoch.track_of(index=int(halo_row))
    except MergerTreeError:
        return None
    track = track_ds.track

    snapnum = np.asarray(track.SnapNum)
    mass = np.asarray(track.Mass, dtype=float)
    is_subhalo = np.asarray(track.IsSubhalo, dtype=bool)

    diffs = np.diff(snapnum)
    gaps = [int(snapnum[i]) for i in range(len(diffs)) if diffs[i] != 1]

    infall_snapnum = None
    if is_subhalo.size:
        transitions = np.where(~is_subhalo[:-1] & is_subhalo[1:])[0]
        if transitions.size:
            infall_snapnum = int(snapnum[transitions[0] + 1])
        elif is_subhalo[0]:
            infall_snapnum = int(snapnum[0])

    extra = track.extra
    group_m200_tree = float(extra["GroupM200"][-1]) \
        if "GroupM200" in extra else None
    group_r200_tree = float(extra["GroupR200"][-1]) \
        if "GroupR200" in extra else None

    cat_mass = float(np.asarray(epoch.halos["mass"])[halo_row])
    cat_radius = float(np.asarray(epoch.halos["radius"])[halo_row])

    mass_ratio = (group_m200_tree / cat_mass) \
        if group_m200_tree is not None and cat_mass else None
    radius_ratio = (group_r200_tree / cat_radius) \
        if group_r200_tree is not None and cat_radius else None

    return TreeSanityReport(
        halo_row=int(halo_row),
        n_snapshots=int(snapnum.size),
        snapnum_min=int(snapnum.min()) if snapnum.size else -1,
        snapnum_max=int(snapnum.max()) if snapnum.size else -1,
        snapnum_gaps=gaps,
        mass_min=float(np.nanmin(mass)) if mass.size else float("nan"),
        mass_max=float(np.nanmax(mass)) if mass.size else float("nan"),
        n_nonpositive_mass=int(np.count_nonzero(mass <= 0)),
        n_nan_mass=int(np.count_nonzero(np.isnan(mass))),
        infall_snapnum=infall_snapnum,
        group_m200_z0_from_tree=group_m200_tree,
        group_r200_z0_from_tree=group_r200_tree,
        group_m200_z0_from_catalogue=cat_mass,
        group_r200_z0_from_catalogue=cat_radius,
        mass_agreement_ratio=mass_ratio,
        radius_agreement_ratio=radius_ratio,
    )


@dataclass
class ShaSfhSanityReport:
    """One satellite's matched SHARK galaxy: does its raw (native-grid,
    un-rebinned) SFH integrate to something in the right ballpark of its
    own current stellar mass, and is ``delta_t``/``lbt_mean`` internally
    consistent. Deliberately bypasses ``time_bin_edges`` rebinning/
    truncation (unlike the pipeline's own ``rebin_sfh`` path) to isolate
    "is SHARK's own raw output sane" from "did this project's config pick
    a wide enough time grid"."""
    halo_row: int
    galaxy_row: Optional[int]
    stellar_mass: Optional[float]
    formed_mass: Optional[float]
    ratio: Optional[float]           # stellar_mass / formed_mass
    delta_t_sum_gyr: Optional[float]
    lbt_span_gyr: Optional[float]
    n_sfh_bins: Optional[int]
    sfr_all_zero: Optional[bool]
    sfr_has_negative: Optional[bool]

    def print_summary(self) -> None:
        if self.galaxy_row is None:
            print(f"halo_row={self.halo_row}: no matched SHARK galaxy "
                 f"(stellar_mass={self.stellar_mass})")
            return
        print(f"halo_row={self.halo_row} -> galaxy_row={self.galaxy_row}: "
             f"StellarMass={self.stellar_mass:.4e}, "
             f"SFH-formed mass={self.formed_mass:.4e}, "
             f"ratio(StellarMass/formed)={self.ratio}")
        print(f"  {self.n_sfh_bins} SFH bins, delta_t sum="
             f"{self.delta_t_sum_gyr:.4f} Gyr, lbt span="
             f"{self.lbt_span_gyr} Gyr")
        if self.sfr_all_zero:
            print("  WARNING: entire SFH is zero")
        if self.sfr_has_negative:
            print("  WARNING: SFH has negative entries")


def shark_sfh_sanity_check(epoch, halo_row: int,
                           match_by: str = "hostHaloID",
                           r_scale: Optional[float] = None
                           ) -> Optional[ShaSfhSanityReport]:
    """Direct check of a satellite's matched SHARK galaxy against its own
    raw ``SharkModel.sfh_disk``/``sfh_bulge``/``get_sfh_meta`` output --
    the same matching (``Epoch.galaxies_in_halo``, most-massive-in-
    aperture as "the" galaxy) ``SharkGalaxyBackend`` uses internally, but
    exposed here so the raw numbers for one satellite can be inspected
    directly. `None` if nothing matches this halo row at all."""
    kwargs = {"match_by": match_by}
    if r_scale is not None:
        kwargs["r_scale"] = r_scale
    matched = epoch.galaxies_in_halo(index=int(halo_row), **kwargs)
    if matched is None or len(matched) == 0 or "mass" not in matched:
        return None
    mass = np.asarray(matched["mass"])
    central = int(np.argmax(mass))
    stellar_mass = float(mass[central])

    model = getattr(matched, "model", None)
    if model is None:
        return ShaSfhSanityReport(
            halo_row=int(halo_row), galaxy_row=None,
            stellar_mass=stellar_mass, formed_mass=None, ratio=None,
            delta_t_sum_gyr=None, lbt_span_gyr=None, n_sfh_bins=None,
            sfr_all_zero=None, sfr_has_negative=None)

    row = int(matched.index[central])
    redshift = getattr(epoch, "redshift", 0.0) or 0.0
    sfh_disk = np.asarray(model.sfh_disk(redshift))[row]
    sfh_bulge = np.asarray(model.sfh_bulge(redshift))[row]
    meta = model.get_sfh_meta(redshift)
    delta_t = np.asarray(meta["delta_t"], dtype=float)          # Gyr
    lbt_mean = np.asarray(meta["lbt_mean"], dtype=float)        # Gyr

    sfr_total = sfh_disk + sfh_bulge                            # Msun/yr
    formed_mass = float(np.sum(sfr_total * delta_t * 1.0e9))    # -> Msun
    ratio = (stellar_mass / formed_mass) if formed_mass > 0 else None
    lbt_span = float(lbt_mean.max() - lbt_mean.min()) \
        if lbt_mean.size > 1 else None

    return ShaSfhSanityReport(
        halo_row=int(halo_row), galaxy_row=row,
        stellar_mass=stellar_mass, formed_mass=formed_mass, ratio=ratio,
        delta_t_sum_gyr=float(delta_t.sum()), lbt_span_gyr=lbt_span,
        n_sfh_bins=int(sfr_total.size),
        sfr_all_zero=bool(np.all(sfr_total == 0)),
        sfr_has_negative=bool(np.any(sfr_total < 0)))
