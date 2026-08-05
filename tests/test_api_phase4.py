#!/usr/bin/env python3
"""
Phase-4 tests: the dataset-aware plotting layer (api.plotting), run
headless (Agg) against the example data.
"""
import os

import matplotlib
matplotlib.use("Agg")

import numpy as np
import pytest

import analysistools as at
from analysistools.api import plotting as plt2

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SNAP = os.path.join(HERE, "data", "snap_0031.hdf5")
VRCAT = os.path.join(HERE, "data", "halos",
                     "snap_0031.VELOCIraptor.properties.0")

have_data = all(os.path.exists(p) for p in (SNAP, VRCAT))
needs_data = pytest.mark.skipif(not have_data, reason="example data missing")


@pytest.fixture(scope="module")
def snap():
    return at.load(SNAP)


@pytest.fixture(scope="module")
def halos():
    # native_includes_h=False: this bundled VELOCIraptor catalogue was run
    # against the SWIFT snapshot above and inherits its h-free convention,
    # not the h-included default HaloCatalogue otherwise guesses for
    # VELOCIraptor -- see docs/unified_interface.md#units-comoving-vs-physical-and-little-h.
    return at.load(VRCAT, snapnum=31, native_includes_h=False)


# ---------------------------------------------------------------------------
# density_map + overlay_points
# ---------------------------------------------------------------------------

@needs_data
def test_density_map_full_box(snap):
    ax = plt2.density_map(snap.dm, grid_size=64, method="NGP")
    grid = ax._at_grid
    assert grid.shape == (64, 64)
    box = snap.meta["boxsize"]
    cell = (box / 64) ** 2
    # projected mass is conserved: sum(sigma * cell area) == total mass
    np.testing.assert_allclose(grid.sum() * cell, snap["mass"].sum(),
                               rtol=1e-4)


@needs_data
def test_density_map_region_with_overlay(snap, halos):
    row = int(np.argmax(halos["mass"]))
    c = halos["pos"][row]
    ax = plt2.density_map(snap.dm, centre=c, size=3.0, grid_size=128,
                          method="CIC")
    near = halos.select(centre=c, size=3.0, geometry="cubic")
    sc = plt2.overlay_points(ax, near, centre=c)
    assert sc.get_offsets().shape[0] == len(near)
    # wrapped frame: overlay coordinates within +/- size
    assert np.all(np.abs(sc.get_offsets()) <= 3.0 + 1e-9)


@needs_data
def test_density_map_bad_projection(snap):
    with pytest.raises(ValueError):
        plt2.density_map(snap, projection="qq")


# ---------------------------------------------------------------------------
# profile + plot_profile
# ---------------------------------------------------------------------------

@needs_data
def test_profile_density_and_plot(snap, halos):
    row = int(np.argmax(halos["mass"]))
    c = halos["pos"][row]
    prof = plt2.profile(snap.dm, c, kind="density", rmin=0.05, rmax=2.0,
                        nbins=15)
    assert prof["kind"] == "density"
    assert len(prof["r"]) == 15
    assert np.nanmax(prof["density"]) > 0
    # density should broadly decline outwards for a halo centre
    valid = prof["density"] > 0
    assert prof["density"][valid][0] > prof["density"][valid][-1]
    ax = plt2.plot_profile(prof)
    assert ax.get_yscale() == "log"


@needs_data
def test_profile_velocity_dispersion(snap, halos):
    row = int(np.argmax(halos["mass"]))
    prof = plt2.profile(snap.dm, halos["pos"][row],
                        kind="velocity_dispersion",
                        rmin=0.05, rmax=1.0, nbins=10)
    assert "sigma" in prof
    assert np.nanmax(prof["sigma"]) > 0


@needs_data
def test_profile_unknown_kind(snap):
    with pytest.raises(ValueError):
        plt2.profile(snap, [0, 0, 0], kind="nonsense", rmin=0.1, rmax=1.0)


# ---------------------------------------------------------------------------
# mass_function
# ---------------------------------------------------------------------------

@needs_data
def test_mass_function_counts_and_volume(halos):
    res = plt2.mass_function(halos)
    key = [k for k in res if k != "_ax"][0]
    r = res[key]
    assert r["counts"].sum() == len(halos)          # every halo binned
    assert r["volume"] == pytest.approx(142.2475**3, rel=1e-3)
    assert np.all(r["phi"][r["counts"] > 0] > 0)


@needs_data
def test_mass_function_multiple_datasets_shared_axes(halos):
    big = halos.select(mass=(float(np.median(halos["mass"])), None))
    big.label = "massive half"
    res = plt2.mass_function([halos, big], cumulative=True)
    labels = [k for k in res if k != "_ax"]
    assert len(labels) == 2
    # cumulative counts at the low-mass end equal the sample sizes
    r_all = res[halos.label]
    assert r_all["counts"].max() == len(halos)
