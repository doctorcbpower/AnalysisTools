"""
Tests for analysistools.catalogue.derived.HaloPropertiesStage (Phase 6b).

Uses hand-built HaloTrack/TrackDataset fixtures for deterministic control
over infall/backsplash/peak-timing edge cases (real tree data doesn't let
you engineer "exactly 2 infalls" on demand), plus one integration test
chaining it after HaloExtractStage/TreeExtractStage/CrossMatchStage against
the bundled real snapshot/halo-catalogue/tree data.
"""
import os

import numpy as np
import pytest

import analysistools as at
from analysistools.api.trees import TrackDataset
from analysistools.catalogue.derived import HaloPropertiesStage
from analysistools.catalogue.pipeline import (
    CrossMatchStage, HaloExtractStage, PipelineContext, TreeExtractStage,
)
from analysistools.merger_tree_types import HaloTrack

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SNAP = os.path.join(HERE, "data", "snap_0031.hdf5")
VRCAT = os.path.join(HERE, "data", "halos",
                     "snap_0031.VELOCIraptor.properties.0")
TREE = os.path.join(HERE, "data", "VELOCIraptor.walkabletree.hdf5")
have_data = all(os.path.exists(p) for p in (SNAP, VRCAT, TREE))
needs_data = pytest.mark.skipif(not have_data, reason="example data missing")


def _track(snapnums, mass, is_subhalo, vmax=None, halo_id=1):
    """Build a TrackDataset from hand-specified per-snapshot arrays,
    ascending time (matching HaloTrack's ordering convention)."""
    snapnums = np.asarray(snapnums, dtype=np.int32)
    n = len(snapnums)
    redshift = np.linspace(1.0, 0.0, n)  # ascending time -> descending z
    extra = {}
    if vmax is not None:
        extra["Vmax"] = np.asarray(vmax, dtype=np.float64)

    track = HaloTrack(
        halo_id=halo_id, query_snapnum=int(snapnums[-1]),
        treefileformat="test",
        SnapNum=snapnums,
        Redshift=redshift,
        Time=1.0 / (1.0 + redshift),
        Mass=np.asarray(mass, dtype=np.float64),
        Pos=np.zeros((n, 3)),
        Vel=np.zeros((n, 3)),
        IsSubhalo=np.asarray(is_subhalo, dtype=bool),
        HostID=np.full(n, -1, dtype=np.int64),
        extra=extra,
    )
    return TrackDataset(track)


def _run(main_branch):
    context = PipelineContext()
    context.columns["MergerTrees/main_branch"] = main_branch
    return HaloPropertiesStage().run(context)


def _col(result, name):
    return result.columns[f"Satellites/HaloProperties/{name}"]


# ---------------------------------------------------------------------------
# Deterministic, hand-built tracks
# ---------------------------------------------------------------------------

def test_mpeak_and_snapshot_at_mpeak():
    track = _track(snapnums=[10, 11, 12, 13],
                   mass=[1e10, 5e10, 3e10, 2e10],
                   is_subhalo=[False, False, False, False])
    result = _run([track])
    assert _col(result, "Mpeak")[0] == pytest.approx(5e10)
    assert _col(result, "SnapshotAtMpeak")[0] == 11


def test_never_a_subhalo_has_no_infall_and_is_not_backsplash():
    track = _track(snapnums=[10, 11, 12], mass=[1e10, 2e10, 3e10],
                   is_subhalo=[False, False, False])
    result = _run([track])
    assert np.isnan(_col(result, "Minfall")[0])
    assert np.isnan(_col(result, "RedshiftInfall")[0])
    assert _col(result, "SnapshotInfall")[0] == -1
    assert _col(result, "NumberOfInfalls")[0] == 0
    assert _col(result, "IsBacksplash")[0] == False  # noqa: E712


def test_single_infall_currently_a_subhalo():
    track = _track(snapnums=[10, 11, 12, 13],
                   mass=[1e10, 2e10, 1.8e10, 1.5e10],
                   is_subhalo=[False, False, True, True],
                   vmax=[20.0, 30.0, 28.0, 25.0])
    result = _run([track])
    assert _col(result, "SnapshotInfall")[0] == 12
    assert _col(result, "RedshiftInfall")[0] == pytest.approx(
        track.track.Redshift[2])
    assert _col(result, "Minfall")[0] == pytest.approx(1.8e10)
    assert _col(result, "Vinfall")[0] == pytest.approx(28.0)
    assert _col(result, "NumberOfInfalls")[0] == 1
    assert _col(result, "IsBacksplash")[0] == False  # noqa: E712
    assert _col(result, "Vpeak")[0] == pytest.approx(30.0)


def test_backsplash_currently_not_a_subhalo_but_was():
    track = _track(snapnums=[10, 11, 12, 13],
                   mass=[1e10, 2e10, 1.5e10, 1.2e10],
                   is_subhalo=[False, True, True, False])
    result = _run([track])
    assert _col(result, "IsBacksplash")[0] == True  # noqa: E712
    assert _col(result, "SnapshotInfall")[0] == 11


def test_multiple_infalls_counted():
    track = _track(snapnums=[10, 11, 12, 13, 14, 15], mass=[1e10] * 6,
                   is_subhalo=[False, True, False, True, False, True])
    result = _run([track])
    assert _col(result, "NumberOfInfalls")[0] == 3


def test_already_subhalo_at_earliest_snapshot_not_counted_as_infall_event():
    # infall_snapshot() still reports the earliest index as "infall", but
    # NumberOfInfalls only counts *observed* False->True transitions.
    track = _track(snapnums=[10, 11, 12], mass=[1e10, 1e10, 1e10],
                   is_subhalo=[True, True, True])
    result = _run([track])
    assert _col(result, "SnapshotInfall")[0] == 10
    assert _col(result, "NumberOfInfalls")[0] == 0


def test_missing_vmax_history_leaves_vpeak_and_vinfall_nan():
    track = _track(snapnums=[10, 11], mass=[1e10, 2e10],
                   is_subhalo=[False, True])  # no vmax= given
    result = _run([track])
    assert np.isnan(_col(result, "Vpeak")[0])
    assert np.isnan(_col(result, "Vinfall")[0])


def test_unresolved_track_gets_sentinels_not_crash():
    good = _track(snapnums=[10, 11], mass=[1e10, 2e10],
                  is_subhalo=[False, True])
    result = _run([good, None])

    assert _col(result, "Mpeak")[0] == pytest.approx(2e10)
    assert np.isnan(_col(result, "Mpeak")[1])
    assert _col(result, "SnapshotAtMpeak")[1] == -1
    assert _col(result, "SnapshotInfall")[1] == -1
    assert _col(result, "IsBacksplash")[1] == False  # noqa: E712
    assert result.provenance[0]["n_no_track"] == 1


def test_check_inputs_requires_main_branch():
    with pytest.raises(RuntimeError, match="halo_properties"):
        HaloPropertiesStage().check_inputs(PipelineContext())


# ---------------------------------------------------------------------------
# Integration: full real-data pipeline through HaloPropertiesStage
# ---------------------------------------------------------------------------

@needs_data
def test_full_pipeline_produces_sane_shapes_and_values():
    sim = at.Simulation(
        snapshots={0.0: SNAP}, halos={0.0: VRCAT}, trees=TREE,
        snapnums={0.0: 31}, label="test",
        load_kwargs={"halos": {"native_includes_h": False}},
    )
    epoch = sim.at(0.0)
    mass = np.asarray(epoch.halos["mass"])
    order = np.argsort(mass)[::-1]
    host_row = int(order[0])
    satellite_rows = [int(r) for r in order[1:8]]

    context = HaloExtractStage(epoch, host_row, satellite_rows).run(
        PipelineContext())
    context = TreeExtractStage(epoch).run(context)
    context = CrossMatchStage().run(context)
    context = HaloPropertiesStage().run(context)

    n = len(satellite_rows)
    mpeak = _col(context, "Mpeak")
    m200c_z0 = _col(context, "M200c_z0")
    assert mpeak.shape == (n,)
    for name in ("SnapshotAtMpeak", "Vpeak", "Minfall", "Vinfall",
                 "RedshiftInfall", "SnapshotInfall", "NumberOfInfalls",
                 "IsBacksplash"):
        assert _col(context, name).shape == (n,)

    # Mpeak (max mass over history) should never be below the z=0 mass,
    # for satellites with a resolved tree.
    resolved = ~np.isnan(mpeak)
    assert np.all(mpeak[resolved] >= m200c_z0[resolved] - 1e-6)
