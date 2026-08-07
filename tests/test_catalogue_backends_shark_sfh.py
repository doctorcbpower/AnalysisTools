"""
Tests for analysistools.catalogue.backends.rebin_sfh and
SharkGalaxyBackend.star_formation_history() (Phase 6b).
"""
import numpy as np
import pytest

from analysistools.catalogue.backends import SharkGalaxyBackend, rebin_sfh


# ---------------------------------------------------------------------------
# rebin_sfh
# ---------------------------------------------------------------------------

def test_rebin_identity():
    edges = np.array([0.0, 1.0, 2.0])
    sfr = np.array([10.0, 20.0])
    np.testing.assert_allclose(rebin_sfh(edges, sfr, edges), sfr)


def test_rebin_splitting_one_bin_into_two():
    native_edges = np.array([0.0, 2.0])
    native_sfr = np.array([10.0])
    output_edges = np.array([0.0, 1.0, 2.0])
    np.testing.assert_allclose(rebin_sfh(native_edges, native_sfr,
                                        output_edges), [10.0, 10.0])


def test_rebin_merging_two_bins_into_one_conserves_mass():
    native_edges = np.array([0.0, 1.0, 2.0])
    native_sfr = np.array([10.0, 20.0])  # mass = 10 + 20 = 30 over width 2
    output_edges = np.array([0.0, 2.0])
    np.testing.assert_allclose(rebin_sfh(native_edges, native_sfr,
                                        output_edges), [15.0])


def test_rebin_partial_overlap():
    native_edges = np.array([0.0, 1.0, 2.0, 3.0])
    native_sfr = np.array([10.0, 20.0, 30.0])
    output_edges = np.array([0.5, 2.5])
    # mass = 10*0.5 + 20*1 + 30*0.5 = 5 + 20 + 15 = 40, width = 2 -> 20
    np.testing.assert_allclose(rebin_sfh(native_edges, native_sfr,
                                        output_edges), [20.0])


def test_rebin_output_beyond_native_range_gets_zero_contribution():
    native_edges = np.array([0.0, 1.0])
    native_sfr = np.array([10.0])
    output_edges = np.array([0.0, 2.0])  # half outside native coverage
    # mass = 10*1 (only [0,1] overlaps) = 10, width = 2 -> 5
    np.testing.assert_allclose(rebin_sfh(native_edges, native_sfr,
                                        output_edges), [5.0])


# ---------------------------------------------------------------------------
# SharkGalaxyBackend.star_formation_history
# ---------------------------------------------------------------------------

class _FakeModel:
    def __init__(self, sfh_disk, sfh_bulge, delta_t_myr, lbt_mean_gyr):
        self._sfh_disk = np.asarray(sfh_disk, dtype=float)
        self._sfh_bulge = np.asarray(sfh_bulge, dtype=float)
        self._meta = {"delta_t": np.asarray(delta_t_myr, dtype=float),
                     "lbt_mean": np.asarray(lbt_mean_gyr, dtype=float)}

    def sfh_disk(self, z):
        return self._sfh_disk

    def sfh_bulge(self, z):
        return self._sfh_bulge

    def get_sfh_meta(self, z):
        return self._meta


class _FakeGalaxyView:
    def __init__(self, columns, model=None, index=None):
        self.columns = {k: np.asarray(v) for k, v in columns.items()}
        n = len(next(iter(self.columns.values()))) if self.columns else 0
        self.model = model
        self.index = (np.asarray(index) if index is not None
                      else np.arange(n))

    def __len__(self):
        return len(next(iter(self.columns.values()))) if self.columns else 0

    def __contains__(self, field):
        return field in self.columns

    def __getitem__(self, field):
        return self.columns[field]

    def get(self, field, default=None):
        return self.columns.get(field, default)


class _FakeEpoch:
    def __init__(self, matched, redshift=0.0):
        self._matched = matched
        self.redshift = redshift

    def galaxies_in_halo(self, index, match_by, r_scale):
        return self._matched


def test_returns_none_when_no_match():
    epoch = _FakeEpoch(_FakeGalaxyView({}))
    result = SharkGalaxyBackend().star_formation_history(
        epoch, halo_row=0, time_bin_edges=[0.0, 1.0])
    assert result is None


def test_returns_none_when_file_backed_catalogue():
    matched = _FakeGalaxyView({"mass": [1e9]}, model=None)
    epoch = _FakeEpoch(matched)
    result = SharkGalaxyBackend().star_formation_history(
        epoch, halo_row=0, time_bin_edges=[0.0, 1.0])
    assert result is None


def test_returns_none_when_epoch_has_no_redshift():
    model = _FakeModel([[1.0]], [[0.0]], [1000.0], [0.5])
    matched = _FakeGalaxyView({"mass": [1e9]}, model=model, index=[0])
    epoch = _FakeEpoch(matched)
    epoch.redshift = None
    result = SharkGalaxyBackend().star_formation_history(
        epoch, halo_row=0, time_bin_edges=[0.0, 1.0])
    assert result is None


def test_disk_plus_bulge_summed_and_identity_rebin():
    # one native bin covering the full grid exactly [0, 1] Gyr
    model = _FakeModel(sfh_disk=[[5.0]], sfh_bulge=[[1.0]],
                       delta_t_myr=[1000.0], lbt_mean_gyr=[0.5])
    matched = _FakeGalaxyView({"mass": [1e9]}, model=model, index=[0])
    epoch = _FakeEpoch(matched, redshift=0.0)

    result = SharkGalaxyBackend().star_formation_history(
        epoch, halo_row=0, time_bin_edges=[0.0, 1.0])
    np.testing.assert_allclose(result, [6.0])  # 5 + 1


def test_two_native_bins_matching_output_grid_exactly():
    model = _FakeModel(sfh_disk=[[5.0, 10.0]], sfh_bulge=[[1.0, 2.0]],
                       delta_t_myr=[500.0, 500.0], lbt_mean_gyr=[0.25, 0.75])
    matched = _FakeGalaxyView({"mass": [1e9]}, model=model, index=[0])
    epoch = _FakeEpoch(matched, redshift=0.0)

    result = SharkGalaxyBackend().star_formation_history(
        epoch, halo_row=0, time_bin_edges=[0.0, 0.5, 1.0])
    np.testing.assert_allclose(result, [6.0, 12.0])


def test_uses_correct_row_via_index_into_full_model_table():
    # 4 galaxies in the model's own tables; the matched view represents
    # only row 3 -- star_formation_history must use model row 3, not row 0.
    disk = np.zeros((4, 1))
    disk[3] = 7.0
    bulge = np.zeros((4, 1))
    model = _FakeModel(sfh_disk=disk, sfh_bulge=bulge,
                       delta_t_myr=[1000.0], lbt_mean_gyr=[0.5])
    matched = _FakeGalaxyView({"mass": [1e9]}, model=model, index=[3])
    epoch = _FakeEpoch(matched, redshift=0.0)

    result = SharkGalaxyBackend().star_formation_history(
        epoch, halo_row=0, time_bin_edges=[0.0, 1.0])
    np.testing.assert_allclose(result, [7.0])


def test_picks_most_massive_matched_galaxy_as_central():
    disk = np.array([[1.0], [9.0]])
    bulge = np.zeros((2, 1))
    model = _FakeModel(sfh_disk=disk, sfh_bulge=bulge,
                       delta_t_myr=[1000.0], lbt_mean_gyr=[0.5])
    matched = _FakeGalaxyView({"mass": [1e9, 5e9]}, model=model,
                              index=[0, 1])
    epoch = _FakeEpoch(matched, redshift=0.0)

    result = SharkGalaxyBackend().star_formation_history(
        epoch, halo_row=0, time_bin_edges=[0.0, 1.0])
    np.testing.assert_allclose(result, [9.0])  # row 1 is more massive
