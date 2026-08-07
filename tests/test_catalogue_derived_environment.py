"""
Tests for analysistools.catalogue.derived.EnvironmentStage (Phase 6b).

Uses a lightweight fake halo catalogue/Epoch for exact, hand-checkable
neighbour-counting/nearest-neighbour values, plus one integration test
chaining it after HaloExtractStage/CrossMatchStage against the bundled
real halo-catalogue data.
"""
import os

import numpy as np
import pytest

import analysistools as at
from analysistools.catalogue.derived import EnvironmentStage
from analysistools.catalogue.pipeline import (
    CrossMatchStage, HaloExtractStage, PipelineContext,
)

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
VRCAT = os.path.join(HERE, "data", "halos",
                     "snap_0031.VELOCIraptor.properties.0")
have_data = os.path.exists(VRCAT)
needs_data = pytest.mark.skipif(not have_data, reason="example data missing")


class _FakeHalos:
    def __init__(self, mass, pos):
        self.mass = np.asarray(mass, dtype=float)
        self.pos = np.asarray(pos, dtype=float)

    def __getitem__(self, key):
        return getattr(self, key)

    def __len__(self):
        return len(self.mass)


class _FakeEpoch:
    def __init__(self, halos, boxsize=None):
        self.halos = halos
        self.boxsize = boxsize


def _context(halo_row, pos_z0):
    context = PipelineContext()
    context.columns["Satellites/_internal/halo_row"] = np.asarray(halo_row)
    context.columns["Satellites/_internal/pos_z0"] = np.asarray(pos_z0)
    return context


# ---------------------------------------------------------------------------
# Deterministic, hand-built halo catalogue
# ---------------------------------------------------------------------------

def _five_halo_catalogue():
    # row0 = host, row1 = the one satellite being evaluated, rows 2-4 are
    # candidate neighbours at known distances/masses from the satellite.
    mass = [1e12, 1e9, 1e11, 1e11, 1e8]
    pos = [[0, 0, 0], [0, 0, 0], [1, 0, 0], [3, 0, 0], [0.5, 0, 0]]
    return _FakeHalos(mass, pos)


def test_local_density_and_nearest_neighbour():
    halos = _five_halo_catalogue()
    epoch = _FakeEpoch(halos)
    context = _context(halo_row=[1], pos_z0=[[0, 0, 0]])

    stage = EnvironmentStage(epoch, host_row=0, mass_threshold=1e10,
                             aperture_radius=2.0)
    result = stage.run(context)

    # neighbours passing mass_threshold=1e10 (excluding host row0, sat
    # row1): row2 (mass=1e11, dist=1), row3 (mass=1e11, dist=3). row4
    # (mass=1e8) is below threshold despite being closest in raw distance.
    density = result.columns["Satellites/Environment/LocalNumberDensity"][0]
    nearest = result.columns[
        "Satellites/Environment/DistanceToNearestMassiveNeighbour"][0]

    expected_volume = (4.0 / 3.0) * np.pi * 2.0 ** 3
    assert density == pytest.approx(1 / expected_volume)  # only row2 is <= 2.0
    assert nearest == pytest.approx(1.0)  # row2, the closer of the two


def test_host_and_self_excluded_from_neighbours():
    # host is far more massive than any threshold and coincides in
    # position with the satellite -- if it weren't excluded it would
    # dominate both outputs.
    halos = _FakeHalos(mass=[1e14, 1e9], pos=[[0, 0, 0], [0, 0, 0]])
    epoch = _FakeEpoch(halos)
    context = _context(halo_row=[1], pos_z0=[[0, 0, 0]])

    stage = EnvironmentStage(epoch, host_row=0, mass_threshold=1e10,
                             aperture_radius=5.0)
    result = stage.run(context)

    assert result.columns["Satellites/Environment/LocalNumberDensity"][0] == 0.0
    assert np.isnan(result.columns[
        "Satellites/Environment/DistanceToNearestMassiveNeighbour"][0])


def test_no_neighbours_above_threshold_warns_and_zeros_density(caplog):
    halos = _FakeHalos(mass=[1e12, 1e9, 1e8], pos=[[0, 0, 0]] * 3)
    epoch = _FakeEpoch(halos)
    context = _context(halo_row=[1], pos_z0=[[0, 0, 0]])

    stage = EnvironmentStage(epoch, host_row=0, mass_threshold=1e10,
                             aperture_radius=5.0)
    with caplog.at_level("WARNING"):
        result = stage.run(context)

    assert result.columns["Satellites/Environment/LocalNumberDensity"][0] == 0.0
    assert np.isnan(result.columns[
        "Satellites/Environment/DistanceToNearestMassiveNeighbour"][0])
    assert any("no neighbours" in r.message for r in caplog.records)


def test_periodic_wrap_finds_true_nearest_distance():
    # satellite near one edge, neighbour near the opposite edge -- the true
    # (wrapped) separation is 0.7, not the naive 9.3.
    halos = _FakeHalos(mass=[1e12, 1e9, 1e11],
                       pos=[[5, 0, 0], [9.5, 0, 0], [0.2, 0, 0]])
    epoch = _FakeEpoch(halos, boxsize=10.0)
    context = _context(halo_row=[1], pos_z0=[[9.5, 0, 0]])

    stage = EnvironmentStage(epoch, host_row=0, mass_threshold=1e10,
                             aperture_radius=1.0)
    result = stage.run(context)

    assert result.columns[
        "Satellites/Environment/DistanceToNearestMassiveNeighbour"][0] == \
        pytest.approx(0.7)


def test_multiple_satellites_independent_results():
    mass = [1e12, 1e9, 1e9, 1e11]
    pos = [[0, 0, 0], [0, 0, 0], [10, 0, 0], [1, 0, 0]]
    halos = _FakeHalos(mass, pos)
    epoch = _FakeEpoch(halos)
    context = _context(halo_row=[1, 2], pos_z0=[[0, 0, 0], [10, 0, 0]])

    stage = EnvironmentStage(epoch, host_row=0, mass_threshold=1e10,
                             aperture_radius=20.0)
    result = stage.run(context)

    nearest = result.columns[
        "Satellites/Environment/DistanceToNearestMassiveNeighbour"]
    assert nearest[0] == pytest.approx(1.0)   # satellite at origin -> row3
    assert nearest[1] == pytest.approx(9.0)   # satellite at [10,0,0] -> row3


def test_host_flags_copied_through_when_present():
    halos = _five_halo_catalogue()
    epoch = _FakeEpoch(halos)
    context = _context(halo_row=[1, 2], pos_z0=[[0, 0, 0], [1, 0, 0]])
    context.columns["Haloes/IsIsolated"] = np.array([True])
    context.columns["Haloes/IsPaired"] = np.array([False])

    stage = EnvironmentStage(epoch, host_row=0, mass_threshold=1e10,
                             aperture_radius=2.0)
    result = stage.run(context)

    np.testing.assert_array_equal(
        result.columns["Satellites/Environment/HostIsIsolated"],
        [True, True])
    np.testing.assert_array_equal(
        result.columns["Satellites/Environment/HostIsPaired"],
        [False, False])


def test_host_flags_omitted_when_host_environment_stage_did_not_run():
    halos = _five_halo_catalogue()
    epoch = _FakeEpoch(halos)
    context = _context(halo_row=[1], pos_z0=[[0, 0, 0]])

    stage = EnvironmentStage(epoch, host_row=0, mass_threshold=1e10,
                             aperture_radius=2.0)
    result = stage.run(context)

    assert "Satellites/Environment/HostIsIsolated" not in result.columns
    assert "Satellites/Environment/HostIsPaired" not in result.columns


def test_raises_when_no_halo_catalogue():
    epoch = _FakeEpoch(None)
    stage = EnvironmentStage(epoch, host_row=0, mass_threshold=1e10,
                             aperture_radius=1.0)
    with pytest.raises(RuntimeError, match="no halo catalogue"):
        stage.run(_context(halo_row=[], pos_z0=np.empty((0, 3))))


def test_check_inputs_requires_halo_extract_columns():
    epoch = _FakeEpoch(_five_halo_catalogue())
    stage = EnvironmentStage(epoch, host_row=0, mass_threshold=1e10,
                             aperture_radius=1.0)
    with pytest.raises(RuntimeError, match="environment"):
        stage.check_inputs(PipelineContext())


# ---------------------------------------------------------------------------
# Integration: real halo catalogue
# ---------------------------------------------------------------------------

@needs_data
def test_full_pipeline_produces_sane_shapes():
    cat = at.load(VRCAT, kind="halos", fileformat="VELOCIraptor",
                  native_includes_h=False, snapnum=31)
    cat.preload()  # meta isn't populated until first load

    class _RealHaloEpoch:
        halos = cat
        boxsize = cat.meta["boxsize"]
        snapnum = cat.meta.get("snapnum")

    epoch = _RealHaloEpoch()

    mass = np.asarray(cat["mass"])
    order = np.argsort(mass)[::-1]
    host_row = int(order[0])
    satellite_rows = [int(r) for r in order[1:8]]

    context = HaloExtractStage(epoch, host_row, satellite_rows).run(
        PipelineContext())
    context = CrossMatchStage().run(context)

    threshold = float(np.median(mass))
    stage = EnvironmentStage(epoch, host_row, mass_threshold=threshold,
                             aperture_radius=5.0)
    result = stage.run(context)

    n = len(satellite_rows)
    density = result.columns["Satellites/Environment/LocalNumberDensity"]
    nearest = result.columns[
        "Satellites/Environment/DistanceToNearestMassiveNeighbour"]
    assert density.shape == (n,)
    assert nearest.shape == (n,)
    assert np.all(density >= 0.0)
