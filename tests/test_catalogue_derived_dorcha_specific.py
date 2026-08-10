"""
Tests for analysistools.catalogue.derived.DorchaSpecificStage (Phase 6b):
EarliestProgenitorRedshift, FossilFraction, and NumberOfMergers --
everything else in DorchaProperties is documented as deferred (see the
class docstring).
"""
import os

import numpy as np
import pytest

import analysistools as at
from analysistools.api.trees import TrackDataset
from analysistools.catalogue.backends import HydroGalaxyBackend
from analysistools.catalogue.derived import (
    DorchaSpecificStage, HaloPropertiesStage, StarFormationHistoryStage,
)
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

TIME_BIN_EDGES = np.array([0.0, 1.0, 2.0, 3.0, 4.0])  # Gyr lookback


def _track(redshifts, halo_id=1):
    n = len(redshifts)
    track = HaloTrack(
        halo_id=halo_id, query_snapnum=0, treefileformat="test",
        SnapNum=np.arange(n, dtype=np.int32),
        Redshift=np.asarray(redshifts, dtype=float),
        Time=1.0 / (1.0 + np.asarray(redshifts, dtype=float)),
        Mass=np.ones(n), Pos=np.zeros((n, 3)), Vel=np.zeros((n, 3)),
        IsSubhalo=np.zeros(n, dtype=bool), HostID=np.full(n, -1, dtype=np.int64),
        extra={},
    )
    return TrackDataset(track)


def _context(main_branch, sfh=None, time_bin_edges=None):
    context = PipelineContext()
    context.columns["MergerTrees/main_branch"] = main_branch
    if sfh is not None:
        context.columns["Satellites/GalaxyProperties/SFH"] = np.asarray(sfh)
    if time_bin_edges is not None:
        context.meta["time_bin_edges_sfh"] = np.asarray(time_bin_edges)
    return context


def _col(result, name):
    return result.columns[f"Satellites/DorchaProperties/{name}"]


class _FakeTreeBackend:
    """Stand-in for MergerTreeTools -- only count_mergers() is used by
    DorchaSpecificStage."""

    def __init__(self, counts_by_id=None):
        self.counts_by_id = counts_by_id or {}
        self.calls = []

    def count_mergers(self, track):
        self.calls.append(track.halo_id)
        return self.counts_by_id.get(track.halo_id)


class _FakeTree:
    def __init__(self, backend=None):
        self.backend = backend or _FakeTreeBackend()


class _FakeEpoch:
    """No-tree by default (the plain reionisation-only test's usual
    case); pass tree=_FakeTree(...) for the NumberOfMergers tests."""

    def __init__(self, tree=None):
        self.tree = tree


# ---------------------------------------------------------------------------
# EarliestProgenitorRedshift
# ---------------------------------------------------------------------------

def test_earliest_progenitor_redshift_is_tracks_first_entry():
    # ascending time -> descending redshift; first entry = highest z
    track = _track(redshifts=[8.0, 4.0, 1.0, 0.0])
    context = _context([track])
    stage = DorchaSpecificStage(_FakeEpoch(), reionisation_lookback_time=13.0)
    result = stage.run(context)

    assert _col(result, "EarliestProgenitorRedshift")[0] == pytest.approx(8.0)


def test_earliest_progenitor_redshift_nan_for_unresolved_track():
    context = _context([None])
    stage = DorchaSpecificStage(_FakeEpoch(), reionisation_lookback_time=13.0)
    result = stage.run(context)

    assert np.isnan(_col(result, "EarliestProgenitorRedshift")[0])
    assert result.provenance[0]["n_no_track"] == 1


def test_check_inputs_requires_main_branch():
    stage = DorchaSpecificStage(_FakeEpoch(), reionisation_lookback_time=13.0)
    with pytest.raises(RuntimeError, match="dorcha_specific"):
        stage.check_inputs(PipelineContext())


# ---------------------------------------------------------------------------
# FossilFraction
# ---------------------------------------------------------------------------

def test_fossil_fraction_nan_without_sfh():
    track = _track(redshifts=[8.0, 0.0])
    context = _context([track])  # no SFH at all
    stage = DorchaSpecificStage(_FakeEpoch(), reionisation_lookback_time=2.0)
    result = stage.run(context)

    assert np.isnan(_col(result, "FossilFraction")[0])


def test_fossil_fraction_nan_and_warns_without_time_bin_edges_meta(caplog):
    track = _track(redshifts=[8.0, 0.0])
    context = _context([track], sfh=[[10.0, 10.0, 10.0, 10.0]])
    # SFH present but meta["time_bin_edges_sfh"] deliberately not set
    stage = DorchaSpecificStage(_FakeEpoch(), reionisation_lookback_time=2.0)
    with caplog.at_level("WARNING"):
        result = stage.run(context)

    assert np.isnan(_col(result, "FossilFraction")[0])
    assert any("time_bin_edges_sfh" in r.message for r in caplog.records)


def test_fossil_fraction_exact_fractional_overlap():
    # constant SFR=10 Msun/yr in all 4 bins (widths=1 Gyr each);
    # reionisation at lookback=2.5 straddles bin index 2 ([2,3]).
    track = _track(redshifts=[8.0, 0.0])
    context = _context([track], sfh=[[10.0, 10.0, 10.0, 10.0]],
                       time_bin_edges=TIME_BIN_EDGES)
    stage = DorchaSpecificStage(_FakeEpoch(), reionisation_lookback_time=2.5)
    result = stage.run(context)

    # bin2: 0.5 Gyr of its 1 Gyr width is pre-reionisation; bin3: all of it.
    # pre_mass = 1e10*0.5 + 1e10*1.0 = 1.5e10; total = 4e10 -> 0.375
    assert _col(result, "FossilFraction")[0] == pytest.approx(0.375)


def test_fossil_fraction_all_before_reionisation():
    track = _track(redshifts=[8.0, 0.0])
    context = _context([track], sfh=[[10.0, 10.0, 10.0, 10.0]],
                       time_bin_edges=TIME_BIN_EDGES)
    stage = DorchaSpecificStage(_FakeEpoch(), reionisation_lookback_time=0.0)
    result = stage.run(context)
    assert _col(result, "FossilFraction")[0] == pytest.approx(1.0)


def test_fossil_fraction_all_after_reionisation():
    track = _track(redshifts=[8.0, 0.0])
    context = _context([track], sfh=[[10.0, 10.0, 10.0, 10.0]],
                       time_bin_edges=TIME_BIN_EDGES)
    stage = DorchaSpecificStage(_FakeEpoch(), reionisation_lookback_time=10.0)
    result = stage.run(context)
    assert _col(result, "FossilFraction")[0] == pytest.approx(0.0)


def test_fossil_fraction_zero_total_mass_stays_nan():
    track = _track(redshifts=[8.0, 0.0])
    context = _context([track], sfh=[[0.0, 0.0, 0.0, 0.0]],
                       time_bin_edges=TIME_BIN_EDGES)
    stage = DorchaSpecificStage(_FakeEpoch(), reionisation_lookback_time=2.0)
    result = stage.run(context)
    assert np.isnan(_col(result, "FossilFraction")[0])


def test_fossil_fraction_skips_satellite_with_nan_sfh_row():
    track_a = _track(redshifts=[8.0, 0.0])
    track_b = _track(redshifts=[6.0, 0.0])
    sfh = np.array([[10.0, 10.0, 10.0, 10.0], [np.nan] * 4])
    context = _context([track_a, track_b], sfh=sfh,
                       time_bin_edges=TIME_BIN_EDGES)
    stage = DorchaSpecificStage(_FakeEpoch(), reionisation_lookback_time=2.5)
    result = stage.run(context)

    assert not np.isnan(_col(result, "FossilFraction")[0])
    assert np.isnan(_col(result, "FossilFraction")[1])


# ---------------------------------------------------------------------------
# NumberOfMergers
# ---------------------------------------------------------------------------

def test_number_of_mergers_uses_epoch_tree_backend_count_mergers():
    track_a = _track(redshifts=[8.0, 0.0], halo_id=1)
    track_b = _track(redshifts=[6.0, 0.0], halo_id=2)
    backend = _FakeTreeBackend({1: 3, 2: 0})
    context = _context([track_a, track_b])
    stage = DorchaSpecificStage(_FakeEpoch(tree=_FakeTree(backend)),
                                reionisation_lookback_time=13.0)
    result = stage.run(context)

    np.testing.assert_array_equal(_col(result, "NumberOfMergers"), [3, 0])
    assert sorted(backend.calls) == [1, 2]


def test_number_of_mergers_minus_one_sentinel_when_backend_returns_none():
    # e.g. a SubFind-HBT tree, where count_mergers() always returns None
    track = _track(redshifts=[8.0, 0.0], halo_id=1)
    backend = _FakeTreeBackend({})  # no entry -> count_mergers returns None
    context = _context([track])
    stage = DorchaSpecificStage(_FakeEpoch(tree=_FakeTree(backend)),
                                reionisation_lookback_time=13.0)
    result = stage.run(context)

    assert _col(result, "NumberOfMergers")[0] == -1
    assert result.provenance[0]["n_no_merger_info"] == 1


def test_number_of_mergers_minus_one_sentinel_when_epoch_has_no_tree(caplog):
    track = _track(redshifts=[8.0, 0.0], halo_id=1)
    context = _context([track])
    stage = DorchaSpecificStage(_FakeEpoch(tree=None),
                                reionisation_lookback_time=13.0)
    with caplog.at_level("WARNING"):
        result = stage.run(context)

    assert _col(result, "NumberOfMergers")[0] == -1
    assert result.provenance[0]["n_no_merger_info"] == 1
    assert any("no merger tree" in r.message for r in caplog.records)


def test_number_of_mergers_minus_one_sentinel_for_unresolved_track():
    context = _context([None])
    stage = DorchaSpecificStage(_FakeEpoch(tree=_FakeTree()),
                                reionisation_lookback_time=13.0)
    result = stage.run(context)

    assert _col(result, "NumberOfMergers")[0] == -1


def test_number_of_mergers_dtype_is_int16():
    track = _track(redshifts=[8.0, 0.0], halo_id=1)
    backend = _FakeTreeBackend({1: 2})
    context = _context([track])
    stage = DorchaSpecificStage(_FakeEpoch(tree=_FakeTree(backend)),
                                reionisation_lookback_time=13.0)
    result = stage.run(context)

    assert _col(result, "NumberOfMergers").dtype == np.int16


# ---------------------------------------------------------------------------
# Integration: full real-data pipeline
# ---------------------------------------------------------------------------

@needs_data
def test_full_pipeline_does_not_crash():
    sim = at.Simulation(
        snapshots={0.0: SNAP}, halos={0.0: VRCAT}, trees=TREE,
        snapnums={0.0: 31}, label="test",
        load_kwargs={"halos": {"native_includes_h": False}},
    )
    epoch = sim.at(0.0)
    mass = np.asarray(epoch.halos["mass"])
    order = np.argsort(mass)[::-1]
    host_row = int(order[0])
    satellite_rows = [int(r) for r in order[1:6]]

    backend = HydroGalaxyBackend()
    context = HaloExtractStage(epoch, host_row, satellite_rows).run(
        PipelineContext())
    context = TreeExtractStage(epoch).run(context)
    context = CrossMatchStage(galaxy_backend=backend, epoch=epoch).run(context)
    context = HaloPropertiesStage().run(context)
    context = StarFormationHistoryStage(
        epoch, backend, time_bin_edges=TIME_BIN_EDGES,
        quenched_ssfr_threshold=1e-11).run(context)
    context = DorchaSpecificStage(epoch, reionisation_lookback_time=13.0).run(context)

    n = len(satellite_rows)
    assert _col(context, "EarliestProgenitorRedshift").shape == (n,)
    assert _col(context, "FossilFraction").shape == (n,)
    # DMO snapshot -> no stars -> SFH all NaN -> FossilFraction all NaN too
    assert np.all(np.isnan(_col(context, "FossilFraction")))

    # bundled tree is a TreeFrog walkable tree -> NumberOfMergers should
    # be genuinely computed (not left at the -1 sentinel), and match an
    # independent call to the same backend method per satellite.
    number_of_mergers = _col(context, "NumberOfMergers")
    assert number_of_mergers.shape == (n,)
    assert number_of_mergers.dtype == np.int16
    main_branch = context.columns["MergerTrees/main_branch"]
    for i, track_ds in enumerate(main_branch):
        expected = epoch.tree.backend.count_mergers(track_ds.track)
        assert number_of_mergers[i] == expected
