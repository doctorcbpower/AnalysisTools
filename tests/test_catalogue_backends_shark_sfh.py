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


def test_rebin_warns_when_native_mass_falls_outside_output_range(caplog):
    # native covers [0,3], output only [0,2] -- the [2,3] native bin's
    # mass (10*1=10, of a 10+10+10=30 total) is silently unplaceable.
    native_edges = np.array([0.0, 1.0, 2.0, 3.0])
    native_sfr = np.array([10.0, 10.0, 10.0])
    output_edges = np.array([0.0, 2.0])
    with caplog.at_level("WARNING"):
        rebin_sfh(native_edges, native_sfr, output_edges)
    assert any("dropped" in r.message for r in caplog.records)
    assert any("33.33%" in r.message for r in caplog.records)


def test_rebin_does_not_warn_when_native_fully_covered(caplog):
    native_edges = np.array([0.0, 1.0, 2.0])
    native_sfr = np.array([10.0, 20.0])
    output_edges = np.array([0.0, 2.0])
    with caplog.at_level("WARNING"):
        rebin_sfh(native_edges, native_sfr, output_edges)
    assert not any("dropped" in r.message for r in caplog.records)


def test_rebin_does_not_warn_for_negligible_loss(caplog):
    # loss well below the 1e-6 relative threshold
    native_edges = np.array([0.0, 1.0, 1.0 + 1e-10])
    native_sfr = np.array([10.0, 10.0])
    output_edges = np.array([0.0, 1.0])
    with caplog.at_level("WARNING"):
        rebin_sfh(native_edges, native_sfr, output_edges)
    assert not any("dropped" in r.message for r in caplog.records)


# ---------------------------------------------------------------------------
# SharkGalaxyBackend.star_formation_history
# ---------------------------------------------------------------------------

class _FakeModel:
    def __init__(self, sfh_disk, sfh_bulge, delta_t_gyr, lbt_mean_gyr):
        # delta_t is in the *same* unit as lbt_mean (Gyr) -- confirmed
        # against a real SHARK run, see backends.py's own comment on
        # this. Was wrongly assumed to be Myr (needing /1e3) previously.
        self._sfh_disk = np.asarray(sfh_disk, dtype=float)
        self._sfh_bulge = np.asarray(sfh_bulge, dtype=float)
        self._meta = {"delta_t": np.asarray(delta_t_gyr, dtype=float),
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
    model = _FakeModel([[1.0]], [[0.0]], [1.0], [0.5])
    matched = _FakeGalaxyView({"mass": [1e9]}, model=model, index=[0])
    epoch = _FakeEpoch(matched)
    epoch.redshift = None
    result = SharkGalaxyBackend().star_formation_history(
        epoch, halo_row=0, time_bin_edges=[0.0, 1.0])
    assert result is None


def test_disk_plus_bulge_summed_and_identity_rebin():
    # one native bin covering the full grid exactly [0, 1] Gyr
    model = _FakeModel(sfh_disk=[[5.0]], sfh_bulge=[[1.0]],
                       delta_t_gyr=[1.0], lbt_mean_gyr=[0.5])
    matched = _FakeGalaxyView({"mass": [1e9]}, model=model, index=[0])
    epoch = _FakeEpoch(matched, redshift=0.0)

    result = SharkGalaxyBackend().star_formation_history(
        epoch, halo_row=0, time_bin_edges=[0.0, 1.0])
    np.testing.assert_allclose(result, [6.0])  # 5 + 1


def test_two_native_bins_matching_output_grid_exactly():
    model = _FakeModel(sfh_disk=[[5.0, 10.0]], sfh_bulge=[[1.0, 2.0]],
                       delta_t_gyr=[0.5, 0.5], lbt_mean_gyr=[0.25, 0.75])
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
                       delta_t_gyr=[1.0], lbt_mean_gyr=[0.5])
    matched = _FakeGalaxyView({"mass": [1e9]}, model=model, index=[3])
    epoch = _FakeEpoch(matched, redshift=0.0)

    result = SharkGalaxyBackend().star_formation_history(
        epoch, halo_row=0, time_bin_edges=[0.0, 1.0])
    np.testing.assert_allclose(result, [7.0])


def test_delta_t_myr_instead_of_gyr_triggers_units_warning(caplog):
    # regression test for the real bug: if delta_t genuinely were Myr
    # (1000x smaller than the Gyr this code now assumes), the mismatch
    # against lbt_mean's span must be flagged, not silently miscomputed
    # again the way the old /1e3 assumption was.
    model = _FakeModel(sfh_disk=[[5.0, 10.0]], sfh_bulge=[[1.0, 2.0]],
                       delta_t_gyr=[0.0005, 0.0005],  # i.e. 0.5 Myr, mislabelled
                       lbt_mean_gyr=[0.25, 0.75])
    matched = _FakeGalaxyView({"mass": [1e9]}, model=model, index=[0])
    epoch = _FakeEpoch(matched, redshift=0.0)

    with caplog.at_level("WARNING"):
        SharkGalaxyBackend().star_formation_history(
            epoch, halo_row=0, time_bin_edges=[0.0, 0.5, 1.0])

    assert any("delta_t" in r.message and "may not be in Gyr" in r.message
              for r in caplog.records)


def test_consistent_delta_t_does_not_warn(caplog):
    model = _FakeModel(sfh_disk=[[5.0, 10.0]], sfh_bulge=[[1.0, 2.0]],
                       delta_t_gyr=[0.5, 0.5], lbt_mean_gyr=[0.25, 0.75])
    matched = _FakeGalaxyView({"mass": [1e9]}, model=model, index=[0])
    epoch = _FakeEpoch(matched, redshift=0.0)

    with caplog.at_level("WARNING"):
        SharkGalaxyBackend().star_formation_history(
            epoch, halo_row=0, time_bin_edges=[0.0, 0.5, 1.0])

    assert not any("may not be in Gyr" in r.message for r in caplog.records)


def test_single_native_bin_skips_units_sanity_check():
    # len(lbt_mean) == 1 -> no span to compare against, must not crash
    model = _FakeModel(sfh_disk=[[5.0]], sfh_bulge=[[1.0]],
                       delta_t_gyr=[1.0], lbt_mean_gyr=[0.5])
    matched = _FakeGalaxyView({"mass": [1e9]}, model=model, index=[0])
    epoch = _FakeEpoch(matched, redshift=0.0)

    result = SharkGalaxyBackend().star_formation_history(
        epoch, halo_row=0, time_bin_edges=[0.0, 1.0])
    np.testing.assert_allclose(result, [6.0])


def test_picks_most_massive_matched_galaxy_as_central():
    disk = np.array([[1.0], [9.0]])
    bulge = np.zeros((2, 1))
    model = _FakeModel(sfh_disk=disk, sfh_bulge=bulge,
                       delta_t_gyr=[1.0], lbt_mean_gyr=[0.5])
    matched = _FakeGalaxyView({"mass": [1e9, 5e9]}, model=model,
                              index=[0, 1])
    epoch = _FakeEpoch(matched, redshift=0.0)

    result = SharkGalaxyBackend().star_formation_history(
        epoch, halo_row=0, time_bin_edges=[0.0, 1.0])
    np.testing.assert_allclose(result, [9.0])  # row 1 is more massive
