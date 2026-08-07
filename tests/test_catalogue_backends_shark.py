"""
Tests for analysistools.catalogue.backends.SharkGalaxyBackend.galaxy_properties()
(Phase 6b).

Uses a lightweight fake Epoch/GalaxyCatalogue-view so the exact matched
values, and the most-massive-is-central selection, are hand-checkable.
"""
import numpy as np
import pytest

from analysistools.catalogue.backends import SharkGalaxyBackend


class _FakeMatched:
    """Stand-in for the GalaxyCatalogue view epoch.galaxies_in_halo()
    returns: supports len(), 'in', [], and .get(field, default)."""

    def __init__(self, columns):
        self.columns = {k: np.asarray(v) for k, v in columns.items()}

    def __len__(self):
        return len(next(iter(self.columns.values()))) if self.columns else 0

    def __contains__(self, field):
        return field in self.columns

    def __getitem__(self, field):
        return self.columns[field]

    def get(self, field, default=None):
        return self.columns.get(field, default)


class _FakeEpoch:
    def __init__(self, matched):
        self._matched = matched
        self.calls = []

    def galaxies_in_halo(self, index, match_by, r_scale):
        self.calls.append({"index": index, "match_by": match_by,
                           "r_scale": r_scale})
        return self._matched


_TWO_GALAXIES = {
    "mass": [1e9, 5e9],
    "mgas": [2e8, 1e9],
    "mstars_metals_disk": [1e7, 5e7],
    "mstars_metals_bulge": [0.0, 0.0],
    "mgas_metals_disk": [1e6, 5e6],
    "mgas_metals_bulge": [0.0, 0.0],
    "mhot": [3e9, 4e9],
    "sfr": [0.5, 1.2],
    "mbh": [1e6, 2e6],
    "id_galaxy": [42, 99],
}


def test_no_match_returns_empty_dict():
    epoch = _FakeEpoch(_FakeMatched({}))
    props = SharkGalaxyBackend().galaxy_properties(epoch, halo_row=0)
    assert props == {}


def test_no_mass_field_returns_empty_dict():
    epoch = _FakeEpoch(_FakeMatched({"sfr": [0.5]}))
    props = SharkGalaxyBackend().galaxy_properties(epoch, halo_row=0)
    assert props == {}


def test_most_massive_galaxy_treated_as_central():
    epoch = _FakeEpoch(_FakeMatched(_TWO_GALAXIES))
    props = SharkGalaxyBackend().galaxy_properties(epoch, halo_row=0)

    # index 1 has mass=5e9, the larger of the two
    assert props["StellarMass"] == pytest.approx(5e9)
    assert props["GasMass_Cold"] == pytest.approx(1e9)
    assert props["GasMass_Hot"] == pytest.approx(4e9)
    assert props["StarFormationRate"] == pytest.approx(1.2)
    assert props["BlackHoleMass"] == pytest.approx(2e6)
    assert props["SharkGalaxyID"] == 99  # id_galaxy of the central (index 1)


def test_metallicities_computed_as_metal_mass_ratio():
    epoch = _FakeEpoch(_FakeMatched(_TWO_GALAXIES))
    props = SharkGalaxyBackend().galaxy_properties(epoch, halo_row=0)

    assert props["MetallicityStellar"] == pytest.approx(5e7 / 5e9)
    assert props["MetallicityGas"] == pytest.approx(5e6 / 1e9)


def test_missing_optional_fields_are_omitted_not_nan():
    # only mass -- no mgas/mhot/sfr/mbh/metals at all
    epoch = _FakeEpoch(_FakeMatched({"mass": [3e9]}))
    props = SharkGalaxyBackend().galaxy_properties(epoch, halo_row=0)

    assert props == {"StellarMass": 3e9}
    assert "GasMass_Cold" not in props
    assert "GasMass_Hot" not in props
    assert "StarFormationRate" not in props
    assert "BlackHoleMass" not in props
    assert "MetallicityStellar" not in props  # no metal-mass fields either


def test_zero_mass_central_skips_metallicity_no_divide_by_zero():
    epoch = _FakeEpoch(_FakeMatched({
        "mass": [0.0], "mgas": [0.0],
        "mstars_metals_disk": [0.0], "mstars_metals_bulge": [0.0],
        "mgas_metals_disk": [0.0], "mgas_metals_bulge": [0.0],
    }))
    props = SharkGalaxyBackend().galaxy_properties(epoch, halo_row=0)

    assert props["StellarMass"] == 0.0
    assert props["GasMass_Cold"] == 0.0
    assert "MetallicityStellar" not in props
    assert "MetallicityGas" not in props


def test_forwards_match_by_and_r_scale_to_epoch():
    epoch = _FakeEpoch(_FakeMatched(_TWO_GALAXIES))
    backend = SharkGalaxyBackend(match_by="id", r_scale=2.5)
    backend.galaxy_properties(epoch, halo_row=7)

    assert epoch.calls == [{"index": 7, "match_by": "id", "r_scale": 2.5}]


def test_default_match_by_and_r_scale():
    epoch = _FakeEpoch(_FakeMatched(_TWO_GALAXIES))
    SharkGalaxyBackend().galaxy_properties(epoch, halo_row=0)

    assert epoch.calls == [{"index": 0, "match_by": "position", "r_scale": 1.0}]


def test_shark_galaxy_id_omitted_when_absent():
    epoch = _FakeEpoch(_FakeMatched({"mass": [3e9]}))  # no id_galaxy
    props = SharkGalaxyBackend().galaxy_properties(epoch, halo_row=0)
    assert "SharkGalaxyID" not in props


def test_native_comoving_little_h_is_always_true_true():
    assert SharkGalaxyBackend().native_comoving_little_h() == (True, True)
    # epoch argument is accepted but ignored -- SHARK's convention is fixed
    assert SharkGalaxyBackend().native_comoving_little_h(epoch=object()) == \
        (True, True)
