#!/usr/bin/env python3
"""
Phase-2 tests: MergerTree/TrackDataset, GalaxyCatalogue, EpochModel.at(z),
and the shark package move + compatibility shim.
"""
import os
import warnings

import h5py
import numpy as np
import pytest

import analysistools as at

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TREE = os.path.join(HERE, "data", "VELOCIraptor.walkabletree.hdf5")
VRCAT = os.path.join(HERE, "data", "halos",
                     "snap_0031.VELOCIraptor.properties.0")

needs_tree = pytest.mark.skipif(not os.path.exists(TREE),
                                reason="example tree missing")
needs_vr = pytest.mark.skipif(not os.path.exists(VRCAT),
                              reason="example VR catalogue missing")


# ---------------------------------------------------------------------------
# Merger trees
# ---------------------------------------------------------------------------

@needs_tree
def test_load_sniffs_walkable_tree():
    tree = at.load(TREE)
    assert isinstance(tree, at.MergerTree)
    assert tree.fileformat == "TreeFrog"


@needs_tree
@needs_vr
def test_track_dataset_with_linked_catalogue():
    cat = at.load(VRCAT, snapnum=31)
    tree = at.load(TREE, halos=cat)
    tr = tree.from_halo(cat, index=0)          # catalogue row -> track
    assert isinstance(tr, at.TrackDataset)
    assert len(tr) > 1
    # Dataset contract: standardised names work
    assert tr["mass"].shape == (len(tr),)
    assert tr["pos"].shape == (len(tr), 3)
    assert tr["redshift"][0] > tr["redshift"][-1]      # earliest -> latest
    # properties filled at snap 31 (the only linked catalogue epoch)
    assert np.isfinite(tr["mass"][tr["snapnum"] == 31]).all()
    # HaloTrack analysis remains reachable
    assert tr.infall() is None or "snapshot" in tr.infall()


@needs_tree
def test_track_by_temporal_id():
    tree = at.load(TREE)
    with h5py.File(TREE) as f:
        tid = int(f["Snapshots/Snap_031/ID"][0])
    tr = tree.track(tid, snapnum=31)
    assert len(tr) > 1
    assert (tr["snapnum"][-1]) == 31
    # no catalogue linked -> topology walks, properties are NaN
    assert np.isnan(tr["mass"]).all()


# ---------------------------------------------------------------------------
# SHARK galaxies (synthetic fixture: no SHARK output ships in data/)
# ---------------------------------------------------------------------------

@pytest.fixture
def shark_file(tmp_path):
    n = 100
    rng = np.random.default_rng(42)
    path = tmp_path / "galaxies.hdf5"
    with h5py.File(path, "w") as f:
        g = f.create_group("galaxies")
        for name in ("position_x", "position_y", "position_z"):
            g.create_dataset(name, data=rng.uniform(0, 100, n))
        for name in ("velocity_x", "velocity_y", "velocity_z"):
            g.create_dataset(name, data=rng.normal(0, 200, n))
        g.create_dataset("mstars_disk", data=rng.lognormal(22, 1, n))
        g.create_dataset("mstars_bulge", data=rng.lognormal(21, 1, n))
        g.create_dataset("mgas_disk", data=rng.lognormal(21, 1, n))
        g.create_dataset("mgas_bulge", data=rng.lognormal(20, 1, n))
        g.create_dataset("sfr_disk", data=rng.lognormal(0, 1, n))
        g.create_dataset("sfr_burst", data=rng.lognormal(-1, 1, n))
        g.create_dataset("id_galaxy", data=np.arange(n, dtype=np.int64))
        g.create_dataset("type", data=rng.integers(0, 3, n))
        g.create_dataset("mvir_hosthalo", data=rng.lognormal(26, 1, n))
        f.create_dataset("cosmology/h", data=0.6751)
        f.create_dataset("run_info/effective_volume", data=(100.0 / 0.6751)**3)
        f.create_dataset("run_info/redshift", data=0.0)
    return str(path)


def test_galaxy_catalogue_from_file(shark_file):
    gals = at.load(shark_file)
    assert isinstance(gals, at.GalaxyCatalogue)
    assert gals.kind == "galaxies"
    assert len(gals) == 100
    assert gals["pos"].shape == (100, 3)
    assert gals["vel"].shape == (100, 3)
    np.testing.assert_allclose(
        gals["mass"], gals["mstars_disk"] + gals["mstars_bulge"])
    np.testing.assert_allclose(
        gals["sfr"], gals["sfr_disk"] + gals["sfr_burst"])
    assert gals.meta["redshift"] == 0.0
    assert gals.meta["h0"] == pytest.approx(0.6751)


def test_galaxy_catalogue_select(shark_file):
    gals = at.load(shark_file)
    lo = float(np.median(gals["mass"]))
    big = gals.select(mass=(lo, None))
    assert 0 < len(big) < len(gals)
    assert big["mass"].min() >= lo
    assert big["pos"].shape == (len(big), 3)
    centrals = gals.select(type=(0, 1))
    assert np.all(centrals["type"] == 0)


# ---------------------------------------------------------------------------
# EpochModel + at(z)
# ---------------------------------------------------------------------------

@needs_vr
def test_halo_model_at_z():
    hm = at.HaloModel({0.0: VRCAT}, fileformat="VELOCIraptor",
                      label="test", standardise=True)
    assert isinstance(hm, at.EpochModel)        # virtual subclass
    view = hm.at(0.0)
    assert view.kind == "halos"
    assert len(view) == 522
    assert view["mass"].shape == (522,)
    assert view["pos"].shape == (522, 3)
    sub = view.select(mass=(float(np.median(view["mass"])), None))
    assert 0 < len(sub) < len(view)
    assert view.meta["redshift"] == 0.0


def test_shark_model_is_epoch_model():
    from analysistools.shark.model import SharkModel
    assert issubclass(SharkModel, at.EpochModel)
    assert hasattr(SharkModel, "at")


# ---------------------------------------------------------------------------
# shark package move + shim
# ---------------------------------------------------------------------------

def test_shark_moved_into_analysistools():
    from analysistools.shark.model import SharkModel
    from analysistools.shark.analysis import Analysis
    assert SharkModel and Analysis


def test_old_shark_import_warns_and_aliases():
    import sys
    for mod in [m for m in sys.modules if m == "shark"
                or m.startswith("shark.")]:
        del sys.modules[mod]
    from analysistools.shark.model import SharkModel
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        import shark
        from shark.model import SharkModel as S2
    assert any(issubclass(x.category, DeprecationWarning) for x in w)
    assert S2 is SharkModel                    # same module object


def test_galaxy_tools_deprecated():
    from analysistools.galaxy_tools import GalaxyTools
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        GalaxyTools("nonexistent.hdf5", "SHARK")   # warns at construction
    assert any(issubclass(x.category, DeprecationWarning) for x in w)
