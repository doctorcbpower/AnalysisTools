"""
Tests for analysistools.catalogue.derived.HostEnvironmentStage (Phase 6b).

Uses a lightweight fake halo catalogue/Epoch (same style as
test_catalogue_derived_environment.py) so isolation/pairing/completeness
counts are exactly hand-checkable, plus one integration test chaining it
after HaloExtractStage against the bundled real halo-catalogue data.
"""
import os

import numpy as np
import pytest

import analysistools as at
from analysistools.catalogue.derived import HostEnvironmentStage
from analysistools.catalogue.pipeline import HaloExtractStage, PipelineContext

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
VRCAT = os.path.join(HERE, "data", "halos",
                     "snap_0031.VELOCIraptor.properties.0")
have_data = os.path.exists(VRCAT)
needs_data = pytest.mark.skipif(not have_data, reason="example data missing")


class _FakeHalos:
    def __init__(self, mass, pos, halo_id=None, radius=None):
        self.mass = np.asarray(mass, dtype=float)
        self.pos = np.asarray(pos, dtype=float)
        n = len(self.mass)
        self.halo_id = (np.asarray(halo_id, dtype=np.int64) if halo_id
                        is not None else np.arange(n, dtype=np.int64))
        self.radius = (np.asarray(radius, dtype=float) if radius is not None
                       else np.ones(n))

    def __getitem__(self, key):
        return getattr(self, key)

    def __len__(self):
        return len(self.mass)


class _FakeEpoch:
    def __init__(self, halos, boxsize=None):
        self.halos = halos
        self.boxsize = boxsize


def _host_context(host_row, mass, radius, pos):
    context = PipelineContext()
    context.columns["Haloes/HostHaloID"] = np.asarray([host_row])
    context.columns["Haloes/M200c"] = np.asarray([mass])
    context.columns["Haloes/R200c"] = np.asarray([radius])
    context.columns["Haloes/Position"] = np.asarray([pos])
    return context


# ---------------------------------------------------------------------------
# IsIsolated / N_satellites_total
# ---------------------------------------------------------------------------

def test_isolated_when_no_more_massive_neighbour_within_radius():
    # host row0, one less-massive neighbour at dist=1 well within R200c*3
    halos = _FakeHalos(mass=[1e12, 1e10], pos=[[0, 0, 0], [1, 0, 0]],
                       radius=[1.0, 1.0])
    epoch = _FakeEpoch(halos)
    context = _host_context(0, 1e12, 1.0, [0, 0, 0])

    stage = HostEnvironmentStage(
        epoch, host_row=0, isolation_radius_factor=3.0,
        pairing_mass_ratio_min=0.3, pairing_max_separation=0.5,
        completeness_mass_threshold=1e9)
    result = stage.run(context)

    assert result.columns["Haloes/IsIsolated"][0] == True  # noqa: E712


def test_not_isolated_when_more_massive_neighbour_within_radius():
    halos = _FakeHalos(mass=[1e12, 2e12], pos=[[0, 0, 0], [1, 0, 0]],
                       radius=[1.0, 1.0])
    epoch = _FakeEpoch(halos)
    context = _host_context(0, 1e12, 1.0, [0, 0, 0])

    stage = HostEnvironmentStage(
        epoch, host_row=0, isolation_radius_factor=3.0,
        pairing_mass_ratio_min=0.3, pairing_max_separation=0.5,
        completeness_mass_threshold=1e9)
    result = stage.run(context)

    assert result.columns["Haloes/IsIsolated"][0] == False  # noqa: E712


def test_more_massive_neighbour_outside_isolation_radius_still_isolated():
    halos = _FakeHalos(mass=[1e12, 2e12], pos=[[0, 0, 0], [10, 0, 0]],
                       radius=[1.0, 1.0])
    epoch = _FakeEpoch(halos)
    context = _host_context(0, 1e12, 1.0, [0, 0, 0])

    stage = HostEnvironmentStage(
        epoch, host_row=0, isolation_radius_factor=3.0,  # radius=3.0
        pairing_mass_ratio_min=0.3, pairing_max_separation=0.5,
        completeness_mass_threshold=1e9)
    result = stage.run(context)

    assert result.columns["Haloes/IsIsolated"][0] == True  # noqa: E712


def test_n_satellites_total_counts_above_completeness_within_isolation_radius():
    mass = [1e12, 1e10, 1e9, 1e7]         # host, above, above, below thresh
    pos = [[0, 0, 0], [1, 0, 0], [2, 0, 0], [1.5, 0, 0]]
    halos = _FakeHalos(mass=mass, pos=pos, radius=[1.0, 1.0, 1.0, 1.0])
    epoch = _FakeEpoch(halos)
    context = _host_context(0, 1e12, 1.0, [0, 0, 0])

    stage = HostEnvironmentStage(
        epoch, host_row=0, isolation_radius_factor=3.0,
        pairing_mass_ratio_min=0.3, pairing_max_separation=0.5,
        completeness_mass_threshold=1e8)
    result = stage.run(context)

    assert result.columns["Haloes/N_satellites_total"][0] == 2


def test_self_row_excluded_from_all_counts():
    halos = _FakeHalos(mass=[1e12], pos=[[0, 0, 0]])
    epoch = _FakeEpoch(halos)
    context = _host_context(0, 1e12, 1.0, [0, 0, 0])

    stage = HostEnvironmentStage(
        epoch, host_row=0, isolation_radius_factor=3.0,
        pairing_mass_ratio_min=0.3, pairing_max_separation=0.5,
        completeness_mass_threshold=1e9)
    result = stage.run(context)

    assert result.columns["Haloes/IsIsolated"][0] == True  # noqa: E712
    assert result.columns["Haloes/N_satellites_total"][0] == 0
    assert result.columns["Haloes/IsPaired"][0] == False  # noqa: E712


# ---------------------------------------------------------------------------
# IsPaired / PairedHostID
# ---------------------------------------------------------------------------

def test_paired_with_comparable_mass_companion_within_separation():
    halos = _FakeHalos(mass=[1e12, 0.8e12], pos=[[0, 0, 0], [0.5, 0, 0]],
                       halo_id=[100, 200])
    epoch = _FakeEpoch(halos)
    context = _host_context(0, 1e12, 1.0, [0, 0, 0])

    stage = HostEnvironmentStage(
        epoch, host_row=0, isolation_radius_factor=0.5,
        pairing_mass_ratio_min=0.3, pairing_max_separation=1.0,
        completeness_mass_threshold=1e9)
    result = stage.run(context)

    assert result.columns["Haloes/IsPaired"][0] == True  # noqa: E712
    assert result.columns["Haloes/PairedHostID"][0] == 200


def test_not_paired_when_mass_ratio_outside_range():
    # companion is 100x less massive -> ratio 0.01 < pairing_mass_ratio_min
    halos = _FakeHalos(mass=[1e12, 1e10], pos=[[0, 0, 0], [0.5, 0, 0]])
    epoch = _FakeEpoch(halos)
    context = _host_context(0, 1e12, 1.0, [0, 0, 0])

    stage = HostEnvironmentStage(
        epoch, host_row=0, isolation_radius_factor=0.5,
        pairing_mass_ratio_min=0.3, pairing_max_separation=1.0,
        completeness_mass_threshold=1e9)
    result = stage.run(context)

    assert result.columns["Haloes/IsPaired"][0] == False  # noqa: E712
    assert result.columns["Haloes/PairedHostID"][0] == -1


def test_not_paired_when_companion_outside_max_separation():
    halos = _FakeHalos(mass=[1e12, 1e12], pos=[[0, 0, 0], [5.0, 0, 0]])
    epoch = _FakeEpoch(halos)
    context = _host_context(0, 1e12, 1.0, [0, 0, 0])

    stage = HostEnvironmentStage(
        epoch, host_row=0, isolation_radius_factor=0.5,
        pairing_mass_ratio_min=0.3, pairing_max_separation=1.0,
        completeness_mass_threshold=1e9)
    result = stage.run(context)

    assert result.columns["Haloes/IsPaired"][0] == False  # noqa: E712


def test_multiple_candidates_picks_nearest_and_warns(caplog):
    halos = _FakeHalos(mass=[1e12, 1e12, 1e12],
                       pos=[[0, 0, 0], [0.8, 0, 0], [0.3, 0, 0]],
                       halo_id=[10, 20, 30])
    epoch = _FakeEpoch(halos)
    context = _host_context(0, 1e12, 1.0, [0, 0, 0])

    stage = HostEnvironmentStage(
        epoch, host_row=0, isolation_radius_factor=0.5,
        pairing_mass_ratio_min=0.3, pairing_max_separation=1.0,
        completeness_mass_threshold=1e9)
    with caplog.at_level("WARNING"):
        result = stage.run(context)

    assert result.columns["Haloes/PairedHostID"][0] == 30  # nearer of the two
    assert any("nearest" in r.message for r in caplog.records)


def test_raises_when_no_halo_catalogue():
    epoch = _FakeEpoch(None)
    stage = HostEnvironmentStage(
        epoch, host_row=0, isolation_radius_factor=3.0,
        pairing_mass_ratio_min=0.3, pairing_max_separation=1.0,
        completeness_mass_threshold=1e9)
    with pytest.raises(RuntimeError, match="no halo catalogue"):
        stage.run(_host_context(0, 1e12, 1.0, [0, 0, 0]))


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
    satellite_rows = [int(r) for r in order[1:6]]

    context = HaloExtractStage(epoch, host_row, satellite_rows).run(
        PipelineContext())

    stage = HostEnvironmentStage(
        epoch, host_row, isolation_radius_factor=3.0,
        pairing_mass_ratio_min=0.3, pairing_max_separation=1.0,
        completeness_mass_threshold=float(np.median(mass)))
    result = stage.run(context)

    assert result.columns["Haloes/IsIsolated"].shape == (1,)
    assert result.columns["Haloes/IsPaired"].shape == (1,)
    assert result.columns["Haloes/PairedHostID"].shape == (1,)
    assert result.columns["Haloes/N_satellites_total"].shape == (1,)
    assert result.columns["Haloes/N_satellites_total"][0] >= 0
