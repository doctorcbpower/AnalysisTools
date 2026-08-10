"""
Tests for analysistools.shark.model.SharkModel's SFH accessors
(sfh_bulge/sfh_metals_bulge), specifically that they sum the
disk-instability-driven and merger-driven bulge growth channels
(bulges_diskins/bulges_mergers) -- see SFH_FIELDS' docstring.

Regression coverage for a real bug: mstars_bulge (used for StellarMass
elsewhere) is already the sum of both channels, but sfh_bulge() used to
read only bulges_diskins, so any galaxy with real merger-driven bulge
growth would show StellarMass exceeding the SFH-integrated formed mass
-- caught via PhysicalValidator's stellar_mass_exceeds_formed_mass check
against a real Dorcha SHARK run.

Builds a small synthetic on-disk SHARK-style directory tree (real group/
dataset names: disks/bulges_diskins/bulges_mergers under
star_formation_histories.hdf5) -- no actual SHARK run needed.
"""
import os

import h5py
import numpy as np
import pytest

from analysistools.shark.common import _redshift_table
from analysistools.shark.model import SharkModel

N_GAL = 4
N_SFH = 5
SNAPSHOT = 1        # the snapshot actually written with data (z=0.0)
SUBVOL = 0


@pytest.fixture
def model_dir_with_sfh(tmp_path):
    rng = np.random.default_rng(3)

    # _redshift_table needs >= 2 rows (np.loadtxt collapses a single-row
    # file to a 0-d array) and strictly increasing snapshot / decreasing
    # redshift -- snapshot 0 is an unused earlier (higher-z) dummy row.
    rt_path = tmp_path / "redshift_list.txt"
    rt_path.write_text(f"0 1.0\n{SNAPSHOT} 0.0\n")

    model_dir = tmp_path / "model"
    subdir = model_dir / str(SNAPSHOT) / str(SUBVOL)
    os.makedirs(subdir, exist_ok=True)

    delta_t = np.full(N_SFH, 100.0)   # Myr
    lbt_mean = np.linspace(0.05, 4.0, N_SFH)  # Gyr

    sfh_disk = rng.lognormal(0.0, 0.3, size=(N_GAL, N_SFH))
    sfh_bulge_diskins = rng.lognormal(-1.0, 0.3, size=(N_GAL, N_SFH))
    sfh_bulge_mergers = rng.lognormal(-1.5, 0.3, size=(N_GAL, N_SFH))
    mz_disk = 0.02 * sfh_disk
    mz_bulge_diskins = 0.02 * sfh_bulge_diskins
    mz_bulge_mergers = 0.02 * sfh_bulge_mergers

    with h5py.File(subdir / "star_formation_histories.hdf5", "w") as f:
        f.create_dataset("delta_t", data=delta_t)
        f.create_dataset("lbt_mean", data=lbt_mean)

        disks = f.create_group("disks")
        disks.create_dataset("star_formation_rate_histories", data=sfh_disk)
        disks.create_dataset("metallicity_histories", data=mz_disk)

        diskins = f.create_group("bulges_diskins")
        diskins.create_dataset("star_formation_rate_histories",
                               data=sfh_bulge_diskins)
        diskins.create_dataset("metallicity_histories", data=mz_bulge_diskins)

        mergers = f.create_group("bulges_mergers")
        mergers.create_dataset("star_formation_rate_histories",
                               data=sfh_bulge_mergers)
        mergers.create_dataset("metallicity_histories", data=mz_bulge_mergers)

    return dict(rt_path=str(rt_path), model_dir=str(model_dir),
               sfh_disk=sfh_disk, sfh_bulge_diskins=sfh_bulge_diskins,
               sfh_bulge_mergers=sfh_bulge_mergers, mz_disk=mz_disk,
               mz_bulge_diskins=mz_bulge_diskins,
               mz_bulge_mergers=mz_bulge_mergers)


def _model(fixture):
    rt = _redshift_table(fixture["rt_path"])
    return SharkModel(fixture["model_dir"], rt, subvols={SUBVOL})


def test_sfh_bulge_sums_diskins_and_mergers_channels(model_dir_with_sfh):
    model = _model(model_dir_with_sfh)
    result = model.sfh_bulge(0.0)
    expected = (model_dir_with_sfh["sfh_bulge_diskins"]
               + model_dir_with_sfh["sfh_bulge_mergers"])
    np.testing.assert_allclose(result, expected)


def test_sfh_metals_bulge_sums_diskins_and_mergers_channels(
        model_dir_with_sfh):
    model = _model(model_dir_with_sfh)
    result = model.sfh_metals_bulge(0.0)
    expected = (model_dir_with_sfh["mz_bulge_diskins"]
               + model_dir_with_sfh["mz_bulge_mergers"])
    np.testing.assert_allclose(result, expected)


def test_sfh_disk_unaffected_single_channel(model_dir_with_sfh):
    model = _model(model_dir_with_sfh)
    result = model.sfh_disk(0.0)
    np.testing.assert_allclose(result, model_dir_with_sfh["sfh_disk"])


def test_sfh_bulge_channels_individually_accessible(model_dir_with_sfh):
    # the raw per-channel keys stay available for callers that want only
    # one component (see SFH_FIELDS' docstring)
    model = _model(model_dir_with_sfh)
    np.testing.assert_allclose(
        model.get("sfh_bulge_diskins", 0.0),
        model_dir_with_sfh["sfh_bulge_diskins"])
    np.testing.assert_allclose(
        model.get("sfh_bulge_mergers", 0.0),
        model_dir_with_sfh["sfh_bulge_mergers"])


def test_sfh_bulge_total_mass_matches_stellar_mass_style_sum(
        model_dir_with_sfh):
    # the whole point of the fix: integrating sfh_bulge() over time
    # should be consistent with mstars_bulge being diskins+mergers, not
    # systematically smaller than it (see this file's module docstring).
    model = _model(model_dir_with_sfh)
    bulge = model.sfh_bulge(0.0)
    diskins_only = model.get("sfh_bulge_diskins", 0.0)
    assert np.all(bulge >= diskins_only)
    assert np.any(bulge > diskins_only)  # mergers channel actually adds mass
