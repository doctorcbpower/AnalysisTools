"""
Tests for analysistools.shark.model.SharkModel's SFH accessors.

Two real bugs, both caught via PhysicalValidator's
stellar_mass_exceeds_formed_mass check against a real Dorcha SHARK run:

1. sfh_bulge()/sfh_metals_bulge() must sum the disk-instability-driven
   and merger-driven bulge growth channels (bulges_diskins/
   bulges_mergers) -- mstars_bulge (used for StellarMass elsewhere) is
   already their sum.
2. star_formation_histories.hdf5's arrays are *not* row-aligned with
   galaxies.hdf5's -- the SFH file commonly covers only a subset of
   galaxies (confirmed on a real run: 124 SFH-file rows vs. >28,000
   galaxies.hdf5 rows in the same subvolume). Both files carry their own
   galaxies/id_galaxy; SharkModel must join on that ID, not assume row
   `i` means the same galaxy in both files -- naively indexing raised
   IndexError for out-of-range rows and silently returned a *different*
   galaxy's SFH for in-range ones.

Builds small synthetic on-disk SHARK-style directory trees (real group/
dataset names) -- no actual SHARK run needed.
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


def _write_sfh_file(path, delta_t, lbt_mean, sfh_disk, sfh_bulge_diskins,
                    sfh_bulge_mergers, mz_disk, mz_bulge_diskins,
                    mz_bulge_mergers, sfh_id_galaxy):
    with h5py.File(path, "w") as f:
        f.create_dataset("delta_t", data=delta_t)
        f.create_dataset("lbt_mean", data=lbt_mean)

        galaxies = f.create_group("galaxies")
        galaxies.create_dataset("id_galaxy",
                                data=np.asarray(sfh_id_galaxy, dtype=np.int32))

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


def _write_galaxies_file(path, id_galaxy):
    with h5py.File(path, "w") as f:
        galaxies = f.create_group("galaxies")
        galaxies.create_dataset("id_galaxy",
                                data=np.asarray(id_galaxy, dtype=np.int32))


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

    # every galaxy present in both files, same order -- row_map is the
    # identity, exercising the "happy path" (no misalignment)
    ids = [10, 11, 12, 13]
    _write_galaxies_file(subdir / "galaxies.hdf5", ids)
    _write_sfh_file(subdir / "star_formation_histories.hdf5",
                    delta_t, lbt_mean, sfh_disk, sfh_bulge_diskins,
                    sfh_bulge_mergers, mz_disk, mz_bulge_diskins,
                    mz_bulge_mergers, sfh_id_galaxy=ids)

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


# ---------------------------------------------------------------------------
# Row misalignment between galaxies.hdf5 and star_formation_histories.hdf5
# ---------------------------------------------------------------------------

@pytest.fixture
def model_dir_with_misaligned_sfh(tmp_path):
    # galaxies.hdf5: 6 galaxies, ids 10..15 (row i -> id 10+i)
    # star_formation_histories.hdf5: only 4 of them, in a shuffled order
    # and with two (12, 14) entirely absent -- mirrors the real bug
    # (SFH file covering a strict, differently-ordered subset).
    rng = np.random.default_rng(5)

    rt_path = tmp_path / "redshift_list.txt"
    rt_path.write_text(f"0 1.0\n{SNAPSHOT} 0.0\n")

    model_dir = tmp_path / "model"
    subdir = model_dir / str(SNAPSHOT) / str(SUBVOL)
    os.makedirs(subdir, exist_ok=True)

    galaxy_ids = [10, 11, 12, 13, 14, 15]
    sfh_ids = [13, 10, 15, 11]   # subset, shuffled -- rows 12 and 14 absent

    n_sfh_gal = len(sfh_ids)
    delta_t = np.full(N_SFH, 100.0)
    lbt_mean = np.linspace(0.05, 4.0, N_SFH)
    sfh_disk = rng.lognormal(0.0, 0.3, size=(n_sfh_gal, N_SFH))
    sfh_bulge_diskins = np.zeros((n_sfh_gal, N_SFH))
    sfh_bulge_mergers = np.zeros((n_sfh_gal, N_SFH))
    mz_disk = 0.02 * sfh_disk
    mz_bulge_diskins = np.zeros_like(sfh_bulge_diskins)
    mz_bulge_mergers = np.zeros_like(sfh_bulge_mergers)

    _write_galaxies_file(subdir / "galaxies.hdf5", galaxy_ids)
    _write_sfh_file(subdir / "star_formation_histories.hdf5",
                    delta_t, lbt_mean, sfh_disk, sfh_bulge_diskins,
                    sfh_bulge_mergers, mz_disk, mz_bulge_diskins,
                    mz_bulge_mergers, sfh_id_galaxy=sfh_ids)

    return dict(rt_path=str(rt_path), model_dir=str(model_dir),
               galaxy_ids=galaxy_ids, sfh_ids=sfh_ids, sfh_disk=sfh_disk)


def test_sfh_reindexed_onto_galaxies_hdf5_row_order(
        model_dir_with_misaligned_sfh):
    fixture = model_dir_with_misaligned_sfh
    model = _model(fixture)
    result = model.sfh_disk(0.0)

    assert result.shape == (len(fixture["galaxy_ids"]), N_SFH)
    for gal_row, gid in enumerate(fixture["galaxy_ids"]):
        if gid in fixture["sfh_ids"]:
            sfh_row = fixture["sfh_ids"].index(gid)
            np.testing.assert_allclose(result[gal_row],
                                       fixture["sfh_disk"][sfh_row])
        else:
            assert np.all(np.isnan(result[gal_row]))


def test_sfh_no_longer_raises_indexerror_for_out_of_range_galaxy_row(
        model_dir_with_misaligned_sfh):
    # regression test for the exact real-world crash: indexing the raw
    # (shorter) SFH array directly with a galaxies.hdf5 row used to
    # raise IndexError once that row exceeded the SFH array's length.
    model = _model(model_dir_with_misaligned_sfh)
    result = model.sfh_disk(0.0)
    assert result.shape[0] == 6   # does not raise, full galaxies.hdf5 length


def test_sfh_missing_galaxies_are_nan_not_another_galaxys_data(
        model_dir_with_misaligned_sfh):
    # rows for ids 12 and 14 (absent from the SFH file) must be NaN, not
    # silently reuse whatever raw row happened to share that index.
    fixture = model_dir_with_misaligned_sfh
    model = _model(fixture)
    result = model.sfh_disk(0.0)

    row_of_12 = fixture["galaxy_ids"].index(12)
    row_of_14 = fixture["galaxy_ids"].index(14)
    assert np.all(np.isnan(result[row_of_12]))
    assert np.all(np.isnan(result[row_of_14]))


def test_sfh_all_unavailable_when_sfh_file_lacks_id_galaxy(tmp_path):
    # older/different SHARK output with no id_galaxy in the SFH file at
    # all -- no safe row correspondence can be established, so every SFH
    # field must be explicitly unavailable rather than guessing row
    # alignment (the exact assumption that caused the real bug).
    rt_path = tmp_path / "redshift_list.txt"
    rt_path.write_text(f"0 1.0\n{SNAPSHOT} 0.0\n")
    model_dir = tmp_path / "model"
    subdir = model_dir / str(SNAPSHOT) / str(SUBVOL)
    os.makedirs(subdir, exist_ok=True)

    _write_galaxies_file(subdir / "galaxies.hdf5", [10, 11])
    with h5py.File(subdir / "star_formation_histories.hdf5", "w") as f:
        f.create_dataset("delta_t", data=np.full(N_SFH, 100.0))
        f.create_dataset("lbt_mean", data=np.linspace(0.05, 4.0, N_SFH))
        disks = f.create_group("disks")
        disks.create_dataset("star_formation_rate_histories",
                             data=np.ones((2, N_SFH)))
        disks.create_dataset("metallicity_histories",
                             data=np.ones((2, N_SFH)))
        # deliberately no "galaxies" group at all

    rt = _redshift_table(str(rt_path))
    model = SharkModel(str(model_dir), rt, subvols={SUBVOL})
    with pytest.raises(AttributeError, match="not found"):
        model.get("sfh_disk", 0.0)
