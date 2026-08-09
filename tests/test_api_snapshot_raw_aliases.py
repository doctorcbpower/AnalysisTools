"""
Tests for analysistools.api.snapshot.SnapshotDataset's _RAW_ATTR_ALIASES
resolution -- the documented public field names (gas_Z/stellar_Z/
initmass/sfr/age, per docs/snapshots.md) don't match the attribute names
analysistools.snapio_hdf5's reader actually sets (gas_metallicity/
stellar_metallicity/stellarinitmass/gas_sfr/stellarage). Uses a fake
backend (monkeypatching SnapshotTools.read_snapshot) so this is exactly
hand-checkable without needing real hydro snapshot data with these extra
blocks (the bundled test snapshot is DMO-only, so this codepath was never
exercised against real data until this fix).
"""
from types import SimpleNamespace

import numpy as np

from analysistools.api.snapshot import SnapshotDataset
from analysistools.snapshot_tools import SnapshotTools


def _fake_data(**extra_attrs):
    n_gas, n_star = 3, 2
    ptype = np.array([0, 0, 0, 4, 4], dtype=np.int32)
    return SimpleNamespace(
        pos=np.zeros((n_gas + n_star, 3)), vel=np.zeros((n_gas + n_star, 3)),
        mass=np.ones(n_gas + n_star), ptype=ptype,
        scale_factor=1.0, hubble_param=0.7, box_size=100.0,
        **extra_attrs,
    )


def _load_snapshot(monkeypatch, fake_data):
    monkeypatch.setattr(SnapshotTools, "read_snapshot",
                        lambda self, filename=None, convention=None: fake_data)
    ds = SnapshotDataset("unused.hdf5")
    ds.preload()
    return ds


def test_sfr_resolves_via_gas_sfr_alias(monkeypatch):
    fake = _fake_data(gas_sfr=np.array([1.0, 2.0, 3.0]))
    ds = _load_snapshot(monkeypatch, fake)
    np.testing.assert_array_equal(ds["sfr"], [1.0, 2.0, 3.0])


def test_age_resolves_via_stellarage_alias(monkeypatch):
    fake = _fake_data(stellarage=np.array([0.5, 0.9]))
    ds = _load_snapshot(monkeypatch, fake)
    np.testing.assert_array_equal(ds["age"], [0.5, 0.9])


def test_gas_z_resolves_via_gas_metallicity_alias(monkeypatch):
    fake = _fake_data(gas_metallicity=np.array([0.01, 0.02, 0.03]))
    ds = _load_snapshot(monkeypatch, fake)
    np.testing.assert_array_equal(ds["gas_Z"], [0.01, 0.02, 0.03])


def test_stellar_z_resolves_via_stellar_metallicity_alias(monkeypatch):
    fake = _fake_data(stellar_metallicity=np.array([0.01, 0.02]))
    ds = _load_snapshot(monkeypatch, fake)
    np.testing.assert_array_equal(ds["stellar_Z"], [0.01, 0.02])


def test_initmass_resolves_via_stellarinitmass_alias(monkeypatch):
    fake = _fake_data(stellarinitmass=np.array([1e5, 2e5], dtype=np.float32))
    ds = _load_snapshot(monkeypatch, fake)
    np.testing.assert_array_equal(ds["initmass"], [1e5, 2e5])


def test_direct_attribute_still_takes_priority_over_alias(monkeypatch):
    # if a future reader version sets "sfr" directly, that should win
    # over the gas_sfr alias fallback -- getattr(d, name) is tried first.
    fake = _fake_data(sfr=np.array([9.0, 9.0, 9.0]),
                      gas_sfr=np.array([1.0, 2.0, 3.0]))
    ds = _load_snapshot(monkeypatch, fake)
    np.testing.assert_array_equal(ds["sfr"], [9.0, 9.0, 9.0])


def test_field_absent_when_neither_name_present(monkeypatch):
    fake = _fake_data()  # no gas_sfr/sfr at all
    ds = _load_snapshot(monkeypatch, fake)
    assert "sfr" not in ds
