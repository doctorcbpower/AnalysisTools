"""
Tests for analysistools.catalogue.derived.ObservabilityStage (Phase 6b).

Uses a lightweight fake Epoch (just needs .boxsize) and hand-built context
columns for exact, hand-checkable distances, plus one integration test
chaining it after HaloExtractStage/CrossMatchStage against the bundled
real halo-catalogue data.
"""
import os

import numpy as np
import pytest

import analysistools as at
from analysistools.catalogue.derived import ObservabilityStage
from analysistools.catalogue.pipeline import (
    CrossMatchStage, HaloExtractStage, PipelineContext,
)

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
VRCAT = os.path.join(HERE, "data", "halos",
                     "snap_0031.VELOCIraptor.properties.0")
have_data = os.path.exists(VRCAT)
needs_data = pytest.mark.skipif(not have_data, reason="example data missing")


class _FakeEpoch:
    def __init__(self, boxsize=None):
        self.boxsize = boxsize


def _context(host_pos, satellite_pos):
    context = PipelineContext()
    context.columns["Haloes/Position"] = np.asarray([host_pos], dtype=float)
    context.columns["Satellites/_internal/pos_z0"] = \
        np.asarray(satellite_pos, dtype=float)
    return context


def _col(result, name):
    return result.columns[f"Satellites/Observability/{name}"]


# ---------------------------------------------------------------------------
# Deterministic geometry
# ---------------------------------------------------------------------------

def test_galactocentric_and_heliocentric_distance():
    # host at origin, satellite forms a 3-4-5 triangle with it
    context = _context(host_pos=[0, 0, 0], satellite_pos=[[3, 4, 0]])
    epoch = _FakeEpoch()
    stage = ObservabilityStage(epoch, observer_pos=[8, 0, 0])
    result = stage.run(context)

    assert _col(result, "GalactocentricRadius")[0] == pytest.approx(5.0)
    # distance from [8,0,0] to [3,4,0] = sqrt(5^2 + 4^2) = sqrt(41)
    assert _col(result, "HeliocentricDistance")[0] == pytest.approx(
        np.sqrt(41.0))


def test_observer_at_host_gives_equal_distances():
    context = _context(host_pos=[1, 1, 1], satellite_pos=[[4, 5, 1]])
    epoch = _FakeEpoch()
    stage = ObservabilityStage(epoch, observer_pos=[1, 1, 1])
    result = stage.run(context)

    assert _col(result, "GalactocentricRadius")[0] == pytest.approx(
        _col(result, "HeliocentricDistance")[0])


def test_periodic_wrap_applies_to_both_distances():
    # box=10: host near one edge, satellite near the opposite edge, and the
    # observer also near the opposite edge -- true wrapped separations are
    # small, naive ones would be ~9.
    context = _context(host_pos=[0.5, 0, 0], satellite_pos=[[9.5, 0, 0]])
    epoch = _FakeEpoch(boxsize=10.0)
    stage = ObservabilityStage(epoch, observer_pos=[9.8, 0, 0])
    result = stage.run(context)

    assert _col(result, "GalactocentricRadius")[0] == pytest.approx(1.0)
    assert _col(result, "HeliocentricDistance")[0] == pytest.approx(0.3)


def test_multiple_satellites_independent_results():
    context = _context(host_pos=[0, 0, 0],
                       satellite_pos=[[3, 4, 0], [0, 0, 6]])
    epoch = _FakeEpoch()
    stage = ObservabilityStage(epoch, observer_pos=[0, 0, 0])
    result = stage.run(context)

    np.testing.assert_allclose(_col(result, "GalactocentricRadius"),
                               [5.0, 6.0])


def test_check_inputs_requires_halo_extract_columns():
    epoch = _FakeEpoch()
    stage = ObservabilityStage(epoch, observer_pos=[0, 0, 0])
    with pytest.raises(RuntimeError, match="observability"):
        stage.check_inputs(PipelineContext())


# ---------------------------------------------------------------------------
# Integration: real halo catalogue
# ---------------------------------------------------------------------------

@needs_data
def test_full_pipeline_produces_sane_shapes():
    cat = at.load(VRCAT, kind="halos", fileformat="VELOCIraptor",
                  native_includes_h=False, snapnum=31)
    cat.preload()

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

    host_pos = context.columns["Haloes/Position"][0]
    # arbitrary mock observer, offset from the host
    observer_pos = host_pos + np.array([8.0, 0.0, 0.0])
    stage = ObservabilityStage(epoch, observer_pos=observer_pos)
    result = stage.run(context)

    n = len(satellite_rows)
    assert _col(result, "GalactocentricRadius").shape == (n,)
    assert _col(result, "HeliocentricDistance").shape == (n,)
    assert np.all(_col(result, "GalactocentricRadius") >= 0.0)
    assert np.all(_col(result, "HeliocentricDistance") >= 0.0)
