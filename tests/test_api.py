#!/usr/bin/env python3
"""
Phase-1 tests for the unified analysistools.api layer, run against the
small example data shipped in data/.

Run from the repository root:  pytest tests/ -q
"""
import os

import numpy as np
import pytest

import analysistools as at

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SNAP = os.path.join(HERE, "data", "snap_0031.hdf5")
VRCAT = os.path.join(HERE, "data", "halos",
                     "snap_0031.VELOCIraptor.properties.0")

needs_snap = pytest.mark.skipif(not os.path.exists(SNAP),
                                reason="example snapshot missing")
needs_vr = pytest.mark.skipif(not os.path.exists(VRCAT),
                              reason="example VR catalogue missing")


# ---------------------------------------------------------------------------
# load() sniffing
# ---------------------------------------------------------------------------

@needs_snap
def test_load_sniffs_swift_snapshot():
    ds = at.load(SNAP)
    assert isinstance(ds, at.SnapshotDataset)
    assert ds.kind == "snapshot"
    assert ds.meta == {"units": {}} or not ds._loaded  # lazy: nothing read yet


@needs_vr
def test_load_sniffs_velociraptor_catalogue():
    ds = at.load(VRCAT)
    assert isinstance(ds, at.HaloCatalogue)
    assert ds.kind == "halos"


def test_load_missing_file_raises():
    with pytest.raises(FileNotFoundError):
        at.load("no_such_file.hdf5")


def test_load_future_kind_raises_helpfully():
    with pytest.raises(NotImplementedError):
        at.load(SNAP, kind="tree")


# ---------------------------------------------------------------------------
# SnapshotDataset
# ---------------------------------------------------------------------------

@needs_snap
def test_snapshot_fields_and_meta():
    ds = at.load(SNAP)
    n = 262144
    assert len(ds) == n
    assert ds["pos"].shape == (n, 3)
    assert ds["vel"].shape == (n, 3)
    assert ds["mass"].shape == (n,)
    assert ds["id"].shape == (n,)                # alias for pids
    np.testing.assert_array_equal(ds["id"], ds["pids"])
    assert ds.meta["redshift"] == pytest.approx(0.0)
    assert ds.meta["boxsize"] == pytest.approx(142.2475, rel=1e-4)
    assert "pos" in ds and "no_such_field" not in ds


@needs_snap
def test_snapshot_species_views():
    ds = at.load(SNAP)
    dm = ds.dm
    assert len(dm) == 262144                     # DM-only snapshot
    assert dm["pos"].shape == (262144, 3)
    np.testing.assert_array_equal(dm.pos, dm["pos"])   # attribute sugar
    assert len(ds.gas) == 0                      # empty species -> empty view
    assert ds.select(species="dm", mass=(None, None)).is_view


@needs_snap
def test_snapshot_geometric_selection():
    ds = at.load(SNAP)
    centre = ds["pos"].mean(axis=0)
    sub = ds.select(centre=centre, size=20.0)
    assert 0 < len(sub) < len(ds)
    r = np.linalg.norm(sub["pos"] - centre, axis=1)
    # periodic wrap can place selected particles across the box edge;
    # distances in the wrapped frame must all be < size
    box = ds.meta["boxsize"]
    d = (sub["pos"] - centre + 0.5 * box) % box - 0.5 * box
    assert np.all(np.linalg.norm(d, axis=1) < 20.0)


@needs_snap
def test_snapshot_selection_chains():
    ds = at.load(SNAP)
    centre = ds["pos"].mean(axis=0)
    a = ds.select(centre=centre, size=30.0)
    b = a.select(centre=centre, size=10.0)
    assert len(b) <= len(a)
    assert b.is_view


@needs_snap
def test_snapshot_backend_escape_hatch():
    ds = at.load(SNAP)
    assert isinstance(ds.backend, at.SnapshotTools)


# ---------------------------------------------------------------------------
# HaloCatalogue
# ---------------------------------------------------------------------------

@needs_vr
def test_halo_catalogue_standardised_fields():
    cat = at.load(VRCAT)
    assert len(cat) == 522
    assert cat["mass"].shape == (522,)
    assert cat["pos"].shape == (522, 3)
    assert cat["vel"].shape == (522, 3)
    assert cat["radius"].shape == (522,)
    assert cat["num_part"].shape == (522,)
    # native fields pass through
    assert "Vmax" in cat.fields or "npart" in cat.fields


@needs_vr
def test_halo_catalogue_subhalos_and_cuts():
    cat = at.load(VRCAT)
    assert cat.has_subhalos
    subs = cat.subhalos
    assert len(subs) == 52
    lo = float(np.median(cat["mass"]))
    massive = cat.select(mass=(lo, None))
    assert 0 < len(massive) < len(cat)
    assert massive["mass"].min() >= lo


@needs_vr
def test_halo_catalogue_string_format_case_insensitive():
    # regression: "VELOCIraptor" vs "VELOCIRAPTOR" vs integer 3
    for fmt in ("VELOCIraptor", "velociraptor", 3):
        cat = at.HaloCatalogue(VRCAT, fileformat=fmt)
        assert len(cat) == 522


# ---------------------------------------------------------------------------
# Cross-source: same syntax everywhere
# ---------------------------------------------------------------------------

@needs_snap
@needs_vr
def test_common_syntax_across_sources():
    snap = at.load(SNAP)
    cat = at.load(VRCAT)
    for ds in (snap, cat):
        assert ds["pos"].shape[1] == 3
        assert len(ds["mass"]) == len(ds)
        assert "pos" in ds.fields
        sub = ds.select(mass=(np.median(ds["mass"]), None))
        assert isinstance(sub, type(ds))
