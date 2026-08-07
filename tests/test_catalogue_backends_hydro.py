"""
Tests for analysistools.catalogue.backends.HydroGalaxyBackend.galaxy_properties()
(Phase 6b).

Uses a lightweight fake Epoch/SnapshotDataset-view so matched particle
values, cosmology, and geometry are all hand-checkable. MeanStellarAge/
StarFormationRate assertions derive their expected values by calling
_cosmic_time_gyr() directly (checking the *wiring* -- mass-weighting,
young-star threshold, initmass-preference, Gyr->yr conversion -- rather
than re-deriving the cosmology formula by hand); a separate test checks
that helper against a known real-world benchmark age of the universe.
"""
import numpy as np
import pytest

from analysistools.catalogue.backends import HydroGalaxyBackend, _cosmic_time_gyr

COSMO = {"h0": 0.7, "omega_0": 0.3, "omega_lambda": 0.7, "scale_factor": 1.0}


class _FakeStars:
    def __init__(self, columns, meta=None):
        self.columns = {k: np.asarray(v) for k, v in columns.items()}
        self.meta = dict(meta) if meta else {}

    def __len__(self):
        return len(next(iter(self.columns.values()))) if self.columns else 0

    def __contains__(self, field):
        return field in self.columns

    def __getitem__(self, field):
        return self.columns[field]

    def get(self, field, default=None):
        return self.columns.get(field, default)


class _FakeHalos:
    def __init__(self, pos):
        self.pos = np.asarray(pos, dtype=float)

    def __getitem__(self, key):
        return getattr(self, key)


class _FakeEpoch:
    def __init__(self, stars, halo_pos, boxsize=None):
        self._stars = stars
        self.halos = _FakeHalos(halo_pos)
        self.boxsize = boxsize
        self.calls = []

    def particles_in_halo(self, index, r_scale, species):
        self.calls.append({"index": index, "r_scale": r_scale,
                           "species": species})
        return self._stars


# ---------------------------------------------------------------------------
# _cosmic_time_gyr sanity check against a known real-world benchmark
# ---------------------------------------------------------------------------

def test_cosmic_time_gyr_matches_known_age_of_universe_today():
    # Planck-like flat LCDM (h=0.7, Om=0.3, OL=0.7): age today ~13.5 Gyr.
    assert _cosmic_time_gyr(1.0, 0.7, 0.3, 0.7) == pytest.approx(13.47, abs=0.05)


# ---------------------------------------------------------------------------
# HydroGalaxyBackend
# ---------------------------------------------------------------------------

def test_no_stars_returns_empty_dict():
    epoch = _FakeEpoch(_FakeStars({}), halo_pos=[[0, 0, 0]])
    props = HydroGalaxyBackend().galaxy_properties(epoch, halo_row=0)
    assert props == {}


def test_no_mass_field_returns_empty_dict():
    epoch = _FakeEpoch(_FakeStars({"pos": [[0, 0, 0]]}), halo_pos=[[0, 0, 0]])
    props = HydroGalaxyBackend().galaxy_properties(epoch, halo_row=0)
    assert props == {}


def test_zero_total_mass_returns_empty_dict():
    stars = _FakeStars({"pos": [[0, 0, 0]], "mass": [0.0]})
    epoch = _FakeEpoch(stars, halo_pos=[[0, 0, 0]])
    props = HydroGalaxyBackend().galaxy_properties(epoch, halo_row=0)
    assert props == {}


def test_stellar_mass_and_metallicity():
    stars = _FakeStars({
        "pos": [[1, 0, 0], [0, 1, 0]],
        "mass": [1e5, 3e5],
        "stellar_Z": [0.01, 0.02],
    })
    epoch = _FakeEpoch(stars, halo_pos=[[0, 0, 0]])
    props = HydroGalaxyBackend().galaxy_properties(epoch, halo_row=0)

    assert props["StellarMass"] == pytest.approx(4e5)
    assert props["MetallicityStellar"] == pytest.approx(
        np.average([0.01, 0.02], weights=[1e5, 3e5]))


def test_mean_stellar_age_and_sfr_with_cosmology():
    a_form = np.array([0.5, 0.999])
    mass = np.array([1e5, 3e5])
    initmass = np.array([1.2e5, 3.5e5])
    stars = _FakeStars({
        "pos": [[0, 0, 0]] * 2, "mass": mass, "initmass": initmass,
        "age": a_form,
    }, meta=COSMO)
    epoch = _FakeEpoch(stars, halo_pos=[[0, 0, 0]])
    props = HydroGalaxyBackend(young_star_age_threshold=0.1).galaxy_properties(
        epoch, halo_row=0)

    t_form = _cosmic_time_gyr(a_form, 0.7, 0.3, 0.7)
    t_now = _cosmic_time_gyr(1.0, 0.7, 0.3, 0.7)
    age_gyr = t_now - t_form
    expected_mean_age = np.average(age_gyr, weights=mass)
    assert props["MeanStellarAge"] == pytest.approx(expected_mean_age)

    young = age_gyr < 0.1
    assert bool(young[1]) and not bool(young[0])  # only a_form=0.999 is young
    expected_sfr = float(initmass[young].sum()) / (0.1 * 1.0e9)
    assert props["StarFormationRate"] == pytest.approx(expected_sfr)


def test_sfr_zero_when_no_young_stars():
    stars = _FakeStars({
        "pos": [[0, 0, 0]] * 2, "mass": [1e5, 2e5],
        "age": [0.5, 0.6],  # both well before z=0
    }, meta=COSMO)
    epoch = _FakeEpoch(stars, halo_pos=[[0, 0, 0]])
    props = HydroGalaxyBackend(young_star_age_threshold=0.1).galaxy_properties(
        epoch, halo_row=0)

    assert props["StarFormationRate"] == 0.0
    assert "MeanStellarAge" in props  # age itself is still computable


def test_sfr_falls_back_to_mass_when_no_initmass():
    stars = _FakeStars({
        "pos": [[0, 0, 0]], "mass": [1e5], "age": [0.999],
    }, meta=COSMO)
    epoch = _FakeEpoch(stars, halo_pos=[[0, 0, 0]])
    props = HydroGalaxyBackend(young_star_age_threshold=0.1).galaxy_properties(
        epoch, halo_row=0)

    assert props["StarFormationRate"] == pytest.approx(1e5 / (0.1 * 1.0e9))


def test_missing_cosmology_skips_age_and_sfr_but_keeps_mass():
    stars = _FakeStars({
        "pos": [[0, 0, 0]], "mass": [1e5], "age": [0.999],
    }, meta={"h0": None, "omega_0": 0.3, "omega_lambda": 0.7,
            "scale_factor": 1.0})
    epoch = _FakeEpoch(stars, halo_pos=[[0, 0, 0]])
    props = HydroGalaxyBackend().galaxy_properties(epoch, halo_row=0)

    assert props["StellarMass"] == pytest.approx(1e5)
    assert "MeanStellarAge" not in props
    assert "StarFormationRate" not in props


def test_missing_age_field_skips_age_and_sfr():
    stars = _FakeStars({"pos": [[0, 0, 0]], "mass": [1e5]}, meta=COSMO)
    epoch = _FakeEpoch(stars, halo_pos=[[0, 0, 0]])
    props = HydroGalaxyBackend().galaxy_properties(epoch, halo_row=0)

    assert "MeanStellarAge" not in props
    assert "StarFormationRate" not in props


def test_half_mass_radius_stellar_exact_geometry():
    # masses [1,1,2] at distances [1,2,3] from centre -> cumulative
    # [1,2,4], half of total (4) is 2, first reached at distance=2.
    stars = _FakeStars({
        "pos": [[1, 0, 0], [2, 0, 0], [3, 0, 0]],
        "mass": [1.0, 1.0, 2.0],
    })
    epoch = _FakeEpoch(stars, halo_pos=[[0, 0, 0]])
    props = HydroGalaxyBackend().galaxy_properties(epoch, halo_row=0)

    assert props["HalfMassRadiusStellar"] == pytest.approx(2.0)


def test_forwards_r_scale_and_species_to_epoch():
    stars = _FakeStars({"pos": [[0, 0, 0]], "mass": [1e5]})
    # halo_row indexes epoch.halos too (for the HalfMassRadiusStellar
    # centre lookup), so give it enough rows to be a valid index.
    epoch = _FakeEpoch(stars, halo_pos=[[0, 0, 0]] * 4)
    HydroGalaxyBackend(r_scale=2.5).galaxy_properties(epoch, halo_row=3)

    assert epoch.calls == [{"index": 3, "r_scale": 2.5, "species": "star"}]


# ---------------------------------------------------------------------------
# star_formation_history
# ---------------------------------------------------------------------------

def test_sfh_returns_none_without_stars():
    epoch = _FakeEpoch(_FakeStars({}), halo_pos=[[0, 0, 0]])
    result = HydroGalaxyBackend().star_formation_history(
        epoch, halo_row=0, time_bin_edges=[0.0, 1.0])
    assert result is None


def test_sfh_returns_none_without_age_field():
    stars = _FakeStars({"pos": [[0, 0, 0]], "mass": [1e5]})  # no "age"
    epoch = _FakeEpoch(stars, halo_pos=[[0, 0, 0]])
    result = HydroGalaxyBackend().star_formation_history(
        epoch, halo_row=0, time_bin_edges=[0.0, 1.0])
    assert result is None


def test_sfh_returns_none_without_cosmology():
    stars = _FakeStars({"pos": [[0, 0, 0]], "mass": [1e5], "age": [0.9]},
                       meta={"h0": None, "omega_0": 0.3,
                            "omega_lambda": 0.7, "scale_factor": 1.0})
    epoch = _FakeEpoch(stars, halo_pos=[[0, 0, 0]])
    result = HydroGalaxyBackend().star_formation_history(
        epoch, halo_row=0, time_bin_edges=[0.0, 1.0])
    assert result is None


def test_sfh_bins_formation_mass_by_lookback_age():
    a_form = np.array([0.999, 0.5])  # one very recent, one ancient
    initmass = np.array([1.0e5, 2.0e5])
    stars = _FakeStars({
        "pos": [[0, 0, 0]] * 2, "mass": [9.0e4, 1.9e5],
        "initmass": initmass, "age": a_form,
    }, meta=COSMO)
    epoch = _FakeEpoch(stars, halo_pos=[[0, 0, 0]])

    t_now = _cosmic_time_gyr(1.0, 0.7, 0.3, 0.7)
    lookback = t_now - _cosmic_time_gyr(a_form, 0.7, 0.3, 0.7)
    # bin edges chosen so the two ages fall in different, known bins
    edges = np.array([0.0, 1.0, 15.0])
    assert lookback[0] < 1.0 < lookback[1]  # sanity: recent vs ancient

    result = HydroGalaxyBackend().star_formation_history(
        epoch, halo_row=0, time_bin_edges=edges)

    expected = np.zeros(2)
    expected[0] = initmass[0] / ((edges[1] - edges[0]) * 1.0e9)
    expected[1] = initmass[1] / ((edges[2] - edges[1]) * 1.0e9)
    np.testing.assert_allclose(result, expected)


def test_sfh_falls_back_to_mass_when_no_initmass():
    stars = _FakeStars({
        "pos": [[0, 0, 0]], "mass": [1.0e5], "age": [0.999],
    }, meta=COSMO)
    epoch = _FakeEpoch(stars, halo_pos=[[0, 0, 0]])

    result = HydroGalaxyBackend().star_formation_history(
        epoch, halo_row=0, time_bin_edges=[0.0, 1.0])
    assert result[0] == pytest.approx(1.0e5 / (1.0 * 1.0e9))
