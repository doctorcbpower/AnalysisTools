"""
test_photometry.py
==================
Smoke tests for shark.photometry, driven through a real SharkModel built
from a synthetic on-disk Shark-style directory tree (galaxies.hdf5 +
star_formation_histories.hdf5 + redshift_list.txt). No actual Shark run
needed; no mocking of SharkModel internals — this exercises the real
read path end to end.

Run with:
    python test_photometry.py
"""

import sys
import os
import shutil
import tempfile

import numpy as np
import h5py

# Allow running from shark/photometry/ directly
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from shark.common import _redshift_table, parse_subvolumes
from shark.model import SharkModel
from shark.photometry.io import galaxy_data, sfh_ages_from_model
from shark.photometry.sps import SPSEngine, _combine_mags


# ------------------------------------------------------------------
# Mock Shark directory tree builder
# ------------------------------------------------------------------

N_GAL  = 10
N_SNAP = 3          # galaxies.hdf5 snapshots (only the last one matters here)
N_SFH  = 15         # SFH bins
SNAPSHOT_OBS = N_SNAP - 1   # "current" snapshot used for z_obs


def make_mock_shark_tree(root: str):
    """
    Build a directory tree of the form:

      root/redshift_list.txt
      root/model/<snapshot>/<subvol>/galaxies.hdf5
      root/model/<snapshot>/<subvol>/star_formation_histories.hdf5

    with one subvolume (0) and N_SNAP snapshots.
    """
    rng = np.random.default_rng(7)

    model_dir = os.path.join(root, "model")
    os.makedirs(model_dir, exist_ok=True)

    # --- redshift_list.txt: snapshot, z (z descending with snapshot) ---
    snaps = np.arange(N_SNAP)
    zvals = np.linspace(2.0, 0.0, N_SNAP)   # descending
    rt_path = os.path.join(root, "redshift_list.txt")
    with open(rt_path, "w") as f:
        for s, z in zip(snaps, zvals):
            f.write(f"{s} {z:.6f}\n")

    # --- per-snapshot galaxies.hdf5 + star_formation_histories.hdf5 ---
    for snap in snaps:
        subdir = os.path.join(model_dir, str(snap), "0")
        os.makedirs(subdir, exist_ok=True)

        mstars_disk  = rng.lognormal(9.0, 0.5, N_GAL)
        mstars_bulge = rng.lognormal(8.5, 0.5, N_GAL)
        sfr_disk     = rng.lognormal(0.0, 0.5, N_GAL)
        sfr_bulge    = rng.lognormal(-0.5, 0.5, N_GAL)

        with h5py.File(os.path.join(subdir, "galaxies.hdf5"), "w") as f:
            cosmo = f.create_group("cosmology")
            cosmo.create_dataset("h",      data=0.678)
            cosmo.create_dataset("OmegaM", data=0.308)
            cosmo.create_dataset("OmegaB", data=0.0486)

            ri = f.create_group("run_info")
            ri.create_dataset("effective_volume", data=100.0)

            gal = f.create_group("galaxies")
            gal.create_dataset("mstars_disk",  data=mstars_disk)
            gal.create_dataset("mstars_bulge", data=mstars_bulge)
            gal.create_dataset("sfr_disk",     data=sfr_disk)
            gal.create_dataset("sfr_burst",    data=sfr_bulge)
            gal.create_dataset("type",         data=np.zeros(N_GAL, dtype=int))

        # SFH bins: ascending lookback time in Shark's own storage order
        # (oldest bin first), spanning 0 -> 13 Gyr
        lbt_mean = np.linspace(13.0, 0.05, N_SFH)
        delta_t  = np.full(N_SFH, 13.0 / N_SFH * 1e3)  # Myr, arbitrary

        sfh_disk  = rng.lognormal(0.0, 0.4, size=(N_GAL, N_SFH))
        sfh_bulge = rng.lognormal(-0.5, 0.4, size=(N_GAL, N_SFH))
        # Metal mass formed per bin: solar-ish metallicity * mass formed proxy
        mz_disk   = 0.02 * sfh_disk
        mz_bulge  = 0.02 * sfh_bulge

        with h5py.File(
            os.path.join(subdir, "star_formation_histories.hdf5"), "w"
        ) as f:
            f.create_dataset("delta_t",  data=delta_t)
            f.create_dataset("lbt_mean", data=lbt_mean)

            sfh_grp = f.create_group("star_formation_history")
            sfh_grp.create_dataset("disk",  data=sfh_disk)
            sfh_grp.create_dataset("burst", data=sfh_bulge)

            z_grp = f.create_group("metallicity_history")
            z_grp.create_dataset("disk",  data=mz_disk)
            z_grp.create_dataset("burst", data=mz_bulge)

    return rt_path, model_dir


# ------------------------------------------------------------------
# Tests
# ------------------------------------------------------------------

def test_model_sfh_access(model):
    print("--- test_model_sfh_access ---")
    z_obs = 0.0

    sfh_d = model.sfh_disk(z_obs)
    sfh_b = model.sfh_bulge(z_obs)
    print(f"  sfh_disk shape  = {sfh_d.shape}")
    assert sfh_d.shape == (N_GAL, N_SFH)
    assert sfh_b.shape == (N_GAL, N_SFH)

    Z_d = model.Z_disk_history(z_obs)
    Z_b = model.Z_bulge_history(z_obs)
    assert Z_d.shape == (N_GAL, N_SFH)
    assert np.all(Z_d >= 1e-4) and np.all(Z_d <= 0.03)
    print(f"  Z_disk range = [{Z_d.min():.4f}, {Z_d.max():.4f}]")

    meta = model.get_sfh_meta(z_obs)
    assert "delta_t" in meta and "lbt_mean" in meta
    assert len(meta["lbt_mean"]) == N_SFH
    print("  PASS")


def test_age_at_z(model):
    print("--- test_age_at_z ---")
    age0 = model.age_at_z(0.0)
    age1 = model.age_at_z(1.0)
    print(f"  age(z=0)={age0:.3f} Gyr, age(z=1)={age1:.3f} Gyr")
    assert age0 > age1 > 0, "age should decrease with increasing z"
    print("  PASS")


def test_sfh_ages_from_model(model):
    print("--- test_sfh_ages_from_model ---")
    ages = sfh_ages_from_model(model, 0.0)
    print(f"  ages[:3]  = {ages[:3]}")
    print(f"  ages[-3:] = {ages[-3:]}")
    assert len(ages) == N_SFH
    assert np.all(np.diff(ages) > 0), "ages must be strictly ascending"
    assert ages[0] > 0
    print("  PASS")


def test_galaxy_data(model):
    print("--- test_galaxy_data ---")
    idx = np.arange(5)
    data = galaxy_data(model, 0.0, idx)
    for key in ("sfr_disk", "sfr_bulge", "Z_disk", "Z_bulge"):
        assert data[key].shape == (5, N_SFH), f"{key} shape mismatch"
    for key in ("mstar_disk", "mstar_bulge"):
        assert data[key].shape == (5,), f"{key} shape mismatch"
    assert np.all(data["Z_disk"] >= 1e-4)
    assert np.all(data["mstar_disk"] >= 0)
    print("  PASS")


def test_combine_mags():
    print("--- test_combine_mags ---")
    a = np.array([4.83])
    b = np.array([4.83])
    combined = _combine_mags(a, b)
    expected = 4.83 - 2.5 * np.log10(2.0)
    assert abs(combined[0] - expected) < 1e-6
    nan = np.array([np.nan])
    assert np.isclose(_combine_mags(nan, a)[0], a[0])
    assert np.isnan(_combine_mags(nan, nan)[0])
    print("  PASS")


def test_sps_engine(model):
    print("--- test_sps_engine (requires python-fsps) ---")
    try:
        import fsps  # noqa: F401
    except ImportError:
        print("  SKIP (python-fsps not installed)")
        return

    z_obs = 0.0
    ages  = sfh_ages_from_model(model, z_obs)
    tage  = model.age_at_z(z_obs)
    data  = galaxy_data(model, z_obs, np.array([0]))

    engine = SPSEngine(imf_type=1, bands=["v"])
    mag = engine.mag_galaxy(
        sfr_disk    = data["sfr_disk"][0],
        sfr_bulge   = data["sfr_bulge"][0],
        Z_disk      = data["Z_disk"][0],
        Z_bulge     = data["Z_bulge"][0],
        snap_ages   = ages,
        mstar_disk  = float(data["mstar_disk"][0]),
        mstar_bulge = float(data["mstar_bulge"][0]),
        tage        = tage,
        z_obs       = z_obs,
    )
    print(f"  M_V (galaxy 0) = {mag[0]:.3f}")
    assert np.isfinite(mag[0])
    assert -30 < mag[0] < 0
    print("  PASS")


def test_pipeline(model):
    print("--- test_pipeline (requires python-fsps) ---")
    try:
        import fsps  # noqa: F401
    except ImportError:
        print("  SKIP (python-fsps not installed)")
        return

    from shark.photometry import PhotometryPipeline

    pipe = PhotometryPipeline(model, z_obs=0.0, bands=["v"], progress=False)
    print(f"  n_galaxies = {pipe.n_galaxies}")

    M_V = pipe.abs_mag_V(gal_indices=np.arange(5))
    print(f"  M_V[:5] = {M_V}")
    assert M_V.shape == (5,)
    assert np.sum(np.isfinite(M_V)) > 0

    ML = pipe.mass_to_light(gal_indices=np.arange(5))
    print(f"  (M/L)_V[:5] = {ML}")
    print("  PASS")


# ------------------------------------------------------------------
# Main
# ------------------------------------------------------------------

if __name__ == "__main__":
    root = tempfile.mkdtemp()
    try:
        rt_path, model_dir = make_mock_shark_tree(root)

        rt = _redshift_table(rt_path)
        sv = parse_subvolumes("0")
        model = SharkModel(model_dir, rt, sv, label="mock")

        test_model_sfh_access(model)
        test_age_at_z(model)
        test_sfh_ages_from_model(model)
        test_galaxy_data(model)
        test_combine_mags()
        test_sps_engine(model)
        test_pipeline(model)

        print("\nAll tests passed.")
    finally:
        shutil.rmtree(root, ignore_errors=True)
