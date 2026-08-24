"""
Tests for analysistools.catalogue.diagnostics -- sanity-check helpers for
merger-tree tracks and SHARK SFH, run standalone against real data (see
notebooks/DorchaCatalogueEndToEnd.ipynb's diagnostics section).

Uses hand-built fakes, same pattern as
tests/test_catalogue_pipeline_halo_extract_subhalo.py.
"""
from types import SimpleNamespace

import numpy as np

from analysistools.catalogue.diagnostics import (
    shark_sfh_sanity_check, tree_sanity_check,
)
from analysistools.merger_tree_types import MergerTreeError


# ---------------------------------------------------------------------------
# tree_sanity_check
# ---------------------------------------------------------------------------

def _fake_track(snapnum, mass, is_subhalo, extra=None):
    ht = SimpleNamespace(
        SnapNum=np.asarray(snapnum), Mass=np.asarray(mass, dtype=float),
        IsSubhalo=np.asarray(is_subhalo, dtype=bool), extra=extra or {})
    return SimpleNamespace(track=ht)


class _FakeHalos:
    def __init__(self, mass, radius):
        self.columns = {"mass": np.asarray(mass, dtype=float),
                        "radius": np.asarray(radius, dtype=float)}

    def __getitem__(self, field):
        return self.columns[field]


def test_tree_sanity_reports_continuous_healthy_track():
    track = _fake_track(
        snapnum=[10, 11, 12], mass=[5.0, 8.0, 10.0],
        is_subhalo=[False, False, True],
        extra={"GroupM200": np.array([50.0, 80.0, 100.0]),
              "GroupR200": np.array([1.0, 1.2, 1.5])})
    epoch = SimpleNamespace(
        track_of=lambda index: track,
        halos=_FakeHalos(mass=[0, 100.0], radius=[0, 1.5]))

    report = tree_sanity_check(epoch, halo_row=1)

    assert report.n_snapshots == 3
    assert report.snapnum_gaps == []
    assert report.n_nonpositive_mass == 0
    assert report.n_nan_mass == 0
    assert report.infall_snapnum == 12
    assert report.group_m200_z0_from_tree == 100.0
    assert report.mass_agreement_ratio == 1.0
    assert report.group_r200_z0_from_tree == 1.5
    assert report.radius_agreement_ratio == 1.0


def test_tree_sanity_detects_snapnum_gap():
    track = _fake_track(snapnum=[10, 11, 14], mass=[5.0, 8.0, 10.0],
                        is_subhalo=[False, False, False])
    epoch = SimpleNamespace(
        track_of=lambda index: track,
        halos=_FakeHalos(mass=[10.0], radius=[1.0]))

    report = tree_sanity_check(epoch, halo_row=0)
    assert report.snapnum_gaps == [11]  # jump happens right after snap 11


def test_tree_sanity_detects_nonpositive_and_nan_mass():
    track = _fake_track(snapnum=[10, 11, 12], mass=[5.0, -1.0, np.nan],
                        is_subhalo=[False, False, False])
    epoch = SimpleNamespace(
        track_of=lambda index: track,
        halos=_FakeHalos(mass=[10.0], radius=[1.0]))

    report = tree_sanity_check(epoch, halo_row=0)
    assert report.n_nonpositive_mass == 1
    assert report.n_nan_mass == 1


def test_tree_sanity_no_infall_for_always_central_track():
    track = _fake_track(snapnum=[10, 11, 12], mass=[5.0, 8.0, 10.0],
                        is_subhalo=[False, False, False])
    epoch = SimpleNamespace(
        track_of=lambda index: track,
        halos=_FakeHalos(mass=[10.0], radius=[1.0]))

    report = tree_sanity_check(epoch, halo_row=0)
    assert report.infall_snapnum is None


def test_tree_sanity_already_subhalo_at_tree_root():
    track = _fake_track(snapnum=[10, 11, 12], mass=[5.0, 8.0, 10.0],
                        is_subhalo=[True, True, True])
    epoch = SimpleNamespace(
        track_of=lambda index: track,
        halos=_FakeHalos(mass=[10.0], radius=[1.0]))

    report = tree_sanity_check(epoch, halo_row=0)
    assert report.infall_snapnum == 10


def test_tree_sanity_handles_missing_group_extra_fields():
    track = _fake_track(snapnum=[10], mass=[5.0], is_subhalo=[False],
                        extra={})
    epoch = SimpleNamespace(
        track_of=lambda index: track,
        halos=_FakeHalos(mass=[10.0], radius=[1.0]))

    report = tree_sanity_check(epoch, halo_row=0)
    assert report.group_m200_z0_from_tree is None
    assert report.mass_agreement_ratio is None
    assert report.group_m200_z0_from_catalogue == 10.0


def test_tree_sanity_returns_none_when_unresolvable():
    def raiser(index):
        raise MergerTreeError("not found")

    epoch = SimpleNamespace(track_of=raiser, halos=_FakeHalos([1.0], [1.0]))
    assert tree_sanity_check(epoch, halo_row=0) is None


def test_tree_sanity_flags_large_tree_catalogue_mass_disagreement():
    # a big deviation from ratio=1 would mean the tree resolved a
    # different group than the row being checked -- exactly what this
    # diagnostic exists to surface.
    track = _fake_track(snapnum=[10], mass=[5.0], is_subhalo=[False],
                        extra={"GroupM200": np.array([500.0])})
    epoch = SimpleNamespace(
        track_of=lambda index: track,
        halos=_FakeHalos(mass=[10.0], radius=[1.0]))

    report = tree_sanity_check(epoch, halo_row=0)
    assert report.mass_agreement_ratio == 50.0


# ---------------------------------------------------------------------------
# shark_sfh_sanity_check
# ---------------------------------------------------------------------------

class _FakeMatched:
    def __init__(self, mass, index, model=None):
        self.columns = {"mass": np.asarray(mass, dtype=float)}
        self.index = np.asarray(index, dtype=np.int64)
        self.model = model

    def __len__(self):
        return len(self.columns["mass"])

    def __contains__(self, field):
        return field in self.columns

    def __getitem__(self, field):
        return self.columns[field]


class _FakeModel:
    def __init__(self, sfh_disk, sfh_bulge, delta_t, lbt_mean):
        self._sfh_disk = np.asarray(sfh_disk, dtype=float)
        self._sfh_bulge = np.asarray(sfh_bulge, dtype=float)
        self._delta_t = np.asarray(delta_t, dtype=float)
        self._lbt_mean = np.asarray(lbt_mean, dtype=float)

    def sfh_disk(self, redshift):
        return self._sfh_disk

    def sfh_bulge(self, redshift):
        return self._sfh_bulge

    def get_sfh_meta(self, redshift):
        return {"delta_t": self._delta_t, "lbt_mean": self._lbt_mean}


def test_shark_sfh_sanity_reports_ratio_near_one():
    # SFR=1 Msun/yr for 1 Gyr -> formed mass = 1e9 Msun; StellarMass=1e9.
    model = _FakeModel(sfh_disk=[[1.0]], sfh_bulge=[[0.0]],
                       delta_t=[1.0], lbt_mean=[0.5])
    matched = _FakeMatched(mass=[1e9], index=[0], model=model)
    epoch = SimpleNamespace(
        galaxies_in_halo=lambda index, **kw: matched, redshift=0.0)

    report = shark_sfh_sanity_check(epoch, halo_row=0)

    assert report.galaxy_row == 0
    assert report.stellar_mass == 1e9
    assert report.formed_mass == 1e9
    assert report.ratio == 1.0
    assert report.sfr_all_zero is False
    assert report.sfr_has_negative is False


def test_shark_sfh_sanity_flags_all_zero_sfh():
    model = _FakeModel(sfh_disk=[[0.0, 0.0]], sfh_bulge=[[0.0, 0.0]],
                       delta_t=[0.5, 0.5], lbt_mean=[0.25, 0.75])
    matched = _FakeMatched(mass=[1e9], index=[0], model=model)
    epoch = SimpleNamespace(
        galaxies_in_halo=lambda index, **kw: matched, redshift=0.0)

    report = shark_sfh_sanity_check(epoch, halo_row=0)
    assert report.sfr_all_zero is True
    assert report.formed_mass == 0.0
    assert report.ratio is None   # division by zero avoided


def test_shark_sfh_sanity_flags_negative_sfh():
    model = _FakeModel(sfh_disk=[[-1.0]], sfh_bulge=[[0.0]],
                       delta_t=[1.0], lbt_mean=[0.5])
    matched = _FakeMatched(mass=[1e9], index=[0], model=model)
    epoch = SimpleNamespace(
        galaxies_in_halo=lambda index, **kw: matched, redshift=0.0)

    report = shark_sfh_sanity_check(epoch, halo_row=0)
    assert report.sfr_has_negative is True


def test_shark_sfh_sanity_returns_stellar_mass_only_for_file_backed_catalogue():
    matched = _FakeMatched(mass=[1e9], index=[0], model=None)
    epoch = SimpleNamespace(
        galaxies_in_halo=lambda index, **kw: matched, redshift=0.0)

    report = shark_sfh_sanity_check(epoch, halo_row=0)
    assert report.galaxy_row is None
    assert report.stellar_mass == 1e9
    assert report.formed_mass is None


def test_shark_sfh_sanity_returns_none_when_no_match():
    empty = _FakeMatched(mass=[], index=[])
    epoch = SimpleNamespace(
        galaxies_in_halo=lambda index, **kw: empty, redshift=0.0)

    assert shark_sfh_sanity_check(epoch, halo_row=0) is None


def test_shark_sfh_sanity_picks_most_massive_as_central():
    # sfh_disk/sfh_bulge are indexed by SHARK's own global galaxy row
    # (matched.index's values), not by position within `matched` -- rows
    # 0-2 exist so index=[0, 2] (row 2 = the more massive match) is valid.
    model = _FakeModel(sfh_disk=[[1.0], [0.0], [2.0]],
                       sfh_bulge=[[0.0], [0.0], [0.0]],
                       delta_t=[1.0], lbt_mean=[0.5])
    matched = _FakeMatched(mass=[1e8, 5e8], index=[0, 2], model=model)
    epoch = SimpleNamespace(
        galaxies_in_halo=lambda index, **kw: matched, redshift=0.0)

    report = shark_sfh_sanity_check(epoch, halo_row=0)
    assert report.galaxy_row == 2            # most massive -> index[1]
    assert report.stellar_mass == 5e8
    assert report.formed_mass == 2e9         # sfh_disk[2]=2.0 * 1 Gyr
