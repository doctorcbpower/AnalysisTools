"""
Tests for HaloExtractStage's Group-vs-Subhalo handling (GADGET/Arepo-
family catalogues, e.g. SubFind): satellite M200c_z0/R200c_z0/Vmax_z0/
SubhaloID_z0 (and the host's Vmax_host) should be sourced from each
object's *primary subhalo* (via GroupFirstSub) when available, not the
FOF Group table -- see HaloExtractStage's own docstring for the full
rationale (this is what makes PhysicalValidator's mpeak_below_m200c_z0
check apples-to-apples again: Mpeak, from the merger tree, is already a
bound-subhalo mass).

Uses a hand-built fake Group/Subhalo catalogue pair (no bundled SubFind
test fixture with a real GroupFirstSub column exists yet).
"""
import numpy as np
import pytest

from analysistools.catalogue.pipeline import HaloExtractStage, PipelineContext


class _FakeSubhaloTable:
    def __init__(self, mass, radius, halo_id, vmax=None):
        self.columns = {
            "mass": np.asarray(mass, dtype=float),
            "radius": np.asarray(radius, dtype=float),
            "halo_id": np.asarray(halo_id, dtype=np.int64),
        }
        if vmax is not None:
            self.columns["Vmax"] = np.asarray(vmax, dtype=float)

    def __len__(self):
        return len(self.columns["mass"])

    def __contains__(self, field):
        return field in self.columns

    def __getitem__(self, field):
        return self.columns[field]


class _FakeGroupTable:
    def __init__(self, mass, radius, pos, halo_id, group_first_sub=None,
                subhalos=None, vel=None):
        self.columns = {
            "mass": np.asarray(mass, dtype=float),
            "radius": np.asarray(radius, dtype=float),
            "pos": np.asarray(pos, dtype=float),
            "halo_id": np.asarray(halo_id, dtype=np.int64),
        }
        if vel is not None:
            self.columns["vel"] = np.asarray(vel, dtype=float)
        if group_first_sub is not None:
            self.columns["GroupFirstSub"] = \
                np.asarray(group_first_sub, dtype=np.int64)
        self._subhalos = subhalos
        self.fileformat = "SUBFIND"
        self.meta = {}

    def __len__(self):
        return len(self.columns["mass"])

    def __contains__(self, field):
        return field in self.columns

    def __getitem__(self, field):
        return self.columns[field]

    @property
    def has_subhalos(self):
        return self._subhalos is not None

    @property
    def subhalos(self):
        if self._subhalos is None:
            raise AttributeError("No subhalo table in this catalogue.")
        return self._subhalos


class _FakeEpoch:
    def __init__(self, halos, boxsize=None, snapnum=0):
        self.halos = halos
        self.boxsize = boxsize
        self.snapnum = snapnum


# 4 groups (host=0, satellites=1,2,3), each group's Group_M_Crit200
# deliberately LARGER than its primary subhalo's own bound SubhaloMass
# (mirroring the real relationship: spherical-overdensity mass >=
# bound-substructure mass).
GROUP_MASS = [100.0, 20.0, 30.0, 40.0]
GROUP_RADIUS = [10.0, 4.0, 5.0, 6.0]
GROUP_POS = [[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0]]
GROUP_HALO_ID = [0, 1, 2, 3]           # ROW_INDEX_SENTINEL convention

SUB_MASS = [80.0, 15.0, 22.0, 33.0]    # < corresponding GROUP_MASS
SUB_RADIUS = [8.0, 3.0, 4.0, 5.0]
SUB_VMAX = [200.0, 90.0, 110.0, 130.0]
SUB_HALO_ID = [0, 1, 2, 3]             # SubhaloNr, also row-index convention


def _catalogue_with_subhalos(group_first_sub):
    subs = _FakeSubhaloTable(SUB_MASS, SUB_RADIUS, SUB_HALO_ID, SUB_VMAX)
    return _FakeGroupTable(GROUP_MASS, GROUP_RADIUS, GROUP_POS, GROUP_HALO_ID,
                          group_first_sub=group_first_sub, subhalos=subs)


def test_satellite_mass_radius_vmax_sourced_from_primary_subhalo():
    halos = _catalogue_with_subhalos(group_first_sub=[0, 1, 2, 3])
    epoch = _FakeEpoch(halos)
    context = HaloExtractStage(epoch, host_row=0, satellite_rows=[1, 2, 3]).run(
        PipelineContext())

    np.testing.assert_allclose(
        context.columns["Satellites/HaloProperties/M200c_z0"],
        [SUB_MASS[1], SUB_MASS[2], SUB_MASS[3]])
    np.testing.assert_allclose(
        context.columns["Satellites/HaloProperties/R200c_z0"],
        [SUB_RADIUS[1], SUB_RADIUS[2], SUB_RADIUS[3]])
    np.testing.assert_allclose(
        context.columns["Satellites/HaloProperties/Vmax_z0"],
        [SUB_VMAX[1], SUB_VMAX[2], SUB_VMAX[3]])
    np.testing.assert_array_equal(
        context.columns["Satellites/Identification/SubhaloID_z0"],
        [SUB_HALO_ID[1], SUB_HALO_ID[2], SUB_HALO_ID[3]])


def test_host_m200c_r200c_stay_group_level_not_subhalo():
    # Haloes/M200c/R200c must remain the Group's own values -- that's
    # what the field is supposed to mean (host virial mass/radius), not
    # a per-object bound-subhalo quantity.
    halos = _catalogue_with_subhalos(group_first_sub=[0, 1, 2, 3])
    epoch = _FakeEpoch(halos)
    context = HaloExtractStage(epoch, host_row=0, satellite_rows=[1, 2, 3]).run(
        PipelineContext())

    assert context.columns["Haloes/M200c"][0] == pytest.approx(GROUP_MASS[0])
    assert context.columns["Haloes/R200c"][0] == pytest.approx(GROUP_RADIUS[0])


def test_host_vmax_sourced_from_primary_subhalo():
    # Group table has no Vmax at all -- Haloes/Vmax_host must come from
    # the host's own primary subhalo, not be omitted.
    halos = _catalogue_with_subhalos(group_first_sub=[0, 1, 2, 3])
    epoch = _FakeEpoch(halos)
    context = HaloExtractStage(epoch, host_row=0, satellite_rows=[1, 2, 3]).run(
        PipelineContext())

    assert context.columns["Haloes/Vmax_host"][0] == pytest.approx(SUB_VMAX[0])


def test_mpeak_consistent_relationship_mass_smaller_than_group():
    # the whole point: with the fix, Satellites/HaloProperties/M200c_z0
    # for a satellite is now *smaller* than that satellite's own
    # Group_M_Crit200 would have been -- consistent with a tree-tracked
    # bound-subhalo Mpeak, rather than systematically larger than it.
    halos = _catalogue_with_subhalos(group_first_sub=[0, 1, 2, 3])
    epoch = _FakeEpoch(halos)
    context = HaloExtractStage(epoch, host_row=0, satellite_rows=[1, 2, 3]).run(
        PipelineContext())

    m200c_z0 = context.columns["Satellites/HaloProperties/M200c_z0"]
    for i, sat_row in enumerate([1, 2, 3]):
        assert m200c_z0[i] < GROUP_MASS[sat_row]


def test_satellite_with_no_valid_primary_subhalo_gets_nan_not_group_fallback(
        caplog):
    # row 2's GroupFirstSub is -1 (sentinel: no identified subhalo) --
    # must be NaN/-1, not silently fall back to the Group's own value
    # (which would reintroduce the exact inconsistency being fixed).
    halos = _catalogue_with_subhalos(group_first_sub=[0, 1, -1, 3])
    epoch = _FakeEpoch(halos)
    with caplog.at_level("WARNING"):
        context = HaloExtractStage(
            epoch, host_row=0, satellite_rows=[1, 2, 3]).run(PipelineContext())

    m200c_z0 = context.columns["Satellites/HaloProperties/M200c_z0"]
    r200c_z0 = context.columns["Satellites/HaloProperties/R200c_z0"]
    vmax_z0 = context.columns["Satellites/HaloProperties/Vmax_z0"]
    subhalo_id = context.columns["Satellites/Identification/SubhaloID_z0"]

    assert np.isnan(m200c_z0[1])    # row 2 -> local index 1
    assert np.isnan(r200c_z0[1])
    assert np.isnan(vmax_z0[1])
    assert subhalo_id[1] == -1
    # the other two satellites are unaffected
    assert not np.isnan(m200c_z0[0])
    assert not np.isnan(m200c_z0[2])
    assert any("no valid primary subhalo" in r.message for r in caplog.records)


def test_out_of_range_group_first_sub_treated_as_invalid():
    # a GroupFirstSub value >= the subhalo table length is also invalid
    # (defensive against a corrupt/out-of-sync sentinel), not an IndexError.
    halos = _catalogue_with_subhalos(group_first_sub=[0, 1, 999, 3])
    epoch = _FakeEpoch(halos)
    context = HaloExtractStage(epoch, host_row=0, satellite_rows=[1, 2, 3]).run(
        PipelineContext())

    m200c_z0 = context.columns["Satellites/HaloProperties/M200c_z0"]
    assert np.isnan(m200c_z0[1])   # row 2 -> local index 1


def test_falls_back_to_group_table_when_no_subhalo_table():
    halos = _FakeGroupTable(GROUP_MASS, GROUP_RADIUS, GROUP_POS, GROUP_HALO_ID,
                            group_first_sub=None, subhalos=None)
    epoch = _FakeEpoch(halos)
    context = HaloExtractStage(epoch, host_row=0, satellite_rows=[1, 2, 3]).run(
        PipelineContext())

    np.testing.assert_allclose(
        context.columns["Satellites/HaloProperties/M200c_z0"],
        [GROUP_MASS[1], GROUP_MASS[2], GROUP_MASS[3]])
    np.testing.assert_array_equal(
        context.columns["Satellites/Identification/SubhaloID_z0"],
        [GROUP_HALO_ID[1], GROUP_HALO_ID[2], GROUP_HALO_ID[3]])
    assert "Haloes/Vmax_host" not in context.columns


def test_falls_back_to_group_table_when_subhalos_present_but_no_group_first_sub():
    subs = _FakeSubhaloTable(SUB_MASS, SUB_RADIUS, SUB_HALO_ID, SUB_VMAX)
    halos = _FakeGroupTable(GROUP_MASS, GROUP_RADIUS, GROUP_POS, GROUP_HALO_ID,
                            group_first_sub=None, subhalos=subs)
    epoch = _FakeEpoch(halos)
    context = HaloExtractStage(epoch, host_row=0, satellite_rows=[1, 2, 3]).run(
        PipelineContext())

    np.testing.assert_allclose(
        context.columns["Satellites/HaloProperties/M200c_z0"],
        [GROUP_MASS[1], GROUP_MASS[2], GROUP_MASS[3]])
