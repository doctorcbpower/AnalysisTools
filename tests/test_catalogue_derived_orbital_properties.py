"""
Tests for analysistools.catalogue.derived.OrbitalPropertiesStage (Phase 6b).

Uses hand-built HaloTrack/TrackDataset fixtures (host + satellite) for
deterministic control over pericentre/apocentre/angular-momentum values --
real tree data doesn't let you engineer an exact orbit on demand -- plus
one integration test chaining it after HaloExtractStage/TreeExtractStage/
CrossMatchStage/HaloPropertiesStage against the bundled real data.
"""
import os

import numpy as np
import pytest

import analysistools as at
from analysistools.api.trees import TrackDataset
from analysistools.catalogue.derived import (
    HaloPropertiesStage, OrbitalPropertiesStage,
)
from analysistools.catalogue.pipeline import (
    CrossMatchStage, HaloExtractStage, PipelineContext, TreeExtractStage,
)
from analysistools.merger_tree_tools import MergerTreeTools
from analysistools.merger_tree_types import HaloTrack, MergerTreeError

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SNAP = os.path.join(HERE, "data", "snap_0031.hdf5")
VRCAT = os.path.join(HERE, "data", "halos",
                     "snap_0031.VELOCIraptor.properties.0")
TREE = os.path.join(HERE, "data", "VELOCIraptor.walkabletree.hdf5")
have_data = all(os.path.exists(p) for p in (SNAP, VRCAT, TREE))
needs_data = pytest.mark.skipif(not have_data, reason="example data missing")


def _track(snapnums, pos, vel=None, is_subhalo=None, halo_id=1):
    snapnums = np.asarray(snapnums, dtype=np.int32)
    n = len(snapnums)
    redshift = np.linspace(1.0, 0.0, n)
    pos = np.asarray(pos, dtype=float)
    vel = np.zeros((n, 3)) if vel is None else np.asarray(vel, dtype=float)
    is_subhalo = (np.zeros(n, dtype=bool) if is_subhalo is None
                 else np.asarray(is_subhalo, dtype=bool))
    track = HaloTrack(
        halo_id=halo_id, query_snapnum=int(snapnums[-1]),
        treefileformat="test",
        SnapNum=snapnums, Redshift=redshift, Time=1.0 / (1.0 + redshift),
        Mass=np.full(n, 1e10), Pos=pos, Vel=vel,
        IsSubhalo=is_subhalo, HostID=np.full(n, -1, dtype=np.int64),
        extra={},
    )
    return TrackDataset(track)


class _MockTree:
    def __init__(self):
        self.backend = object.__new__(MergerTreeTools)


class _MockEpoch:
    """analyse_orbit() only touches self.BoxSize on MergerTreeTools when
    boxsize=None is passed through -- OrbitalPropertiesStage always passes
    epoch.boxsize explicitly, so a bare MergerTreeTools (no __init__) is a
    safe stand-in here."""

    def __init__(self, host_track_ds, boxsize=None):
        self._host_track_ds = host_track_ds
        self.tree = _MockTree()
        self.boxsize = boxsize

    def track_of(self, index):
        if self._host_track_ds is None:
            raise MergerTreeError("host not found")
        return self._host_track_ds


def _run(main_branch, host_track_ds, boxsize=None, host_row=0):
    context = PipelineContext()
    context.columns["MergerTrees/main_branch"] = main_branch
    epoch = _MockEpoch(host_track_ds, boxsize=boxsize)
    return OrbitalPropertiesStage(epoch, host_row).run(context)


def _col(result, name):
    return result.columns[f"Satellites/HaloProperties/{name}"]


# ---------------------------------------------------------------------------
# Deterministic, hand-built host + satellite tracks
# ---------------------------------------------------------------------------

def test_pericentre_apocentre_eccentricity():
    snaps = [10, 11, 12, 13]
    host = _track(snaps, pos=[[0, 0, 0]] * 4)
    sat = _track(snaps, pos=[[10, 0, 0], [5, 0, 0], [8, 0, 0], [12, 0, 0]])

    result = _run([sat], host)

    assert _col(result, "OrbitalPericentre")[0] == pytest.approx(5.0)
    assert _col(result, "OrbitalApocentre")[0] == pytest.approx(12.0)
    assert _col(result, "OrbitalEccentricity")[0] == pytest.approx(7.0 / 17.0)


def test_angular_momentum_at_infall():
    snaps = [10, 11, 12]
    host = _track(snaps, pos=[[0, 0, 0]] * 3)
    sat = _track(snaps, pos=[[10, 0, 0], [5, 0, 0], [5, 0, 0]],
                vel=[[0, 0, 0], [0, 3, 0], [0, 3, 0]],
                is_subhalo=[False, True, True])

    result = _run([sat], host)

    # infall at snap 11: r=[5,0,0], v=[0,3,0] -> L = r x v = [0,0,15]
    np.testing.assert_allclose(_col(result, "OrbitalAngularMomentum")[0],
                               [0.0, 0.0, 15.0])


def test_never_a_subhalo_has_nan_angular_momentum_but_real_peri_apo():
    snaps = [10, 11, 12]
    host = _track(snaps, pos=[[0, 0, 0]] * 3)
    sat = _track(snaps, pos=[[10, 0, 0], [6, 0, 0], [8, 0, 0]])

    result = _run([sat], host)

    assert not np.isnan(_col(result, "OrbitalPericentre")[0])
    assert np.all(np.isnan(_col(result, "OrbitalAngularMomentum")[0]))


def test_no_common_snapshots_leaves_nan_and_counts_no_overlap():
    host = _track([10, 11, 12], pos=[[0, 0, 0]] * 3)
    sat = _track([20, 21, 22], pos=[[1, 0, 0]] * 3)  # disjoint SnapNum

    result = _run([sat], host)

    assert np.isnan(_col(result, "OrbitalPericentre")[0])
    assert result.provenance[0]["n_no_overlap"] == 1


def test_unresolved_satellite_track_skipped():
    host = _track([10, 11], pos=[[0, 0, 0]] * 2)
    result = _run([None], host)

    assert np.isnan(_col(result, "OrbitalPericentre")[0])
    assert result.provenance[0]["n_no_track"] == 1


def test_no_tree_raises():
    epoch = _MockEpoch(None)
    epoch.tree = None
    context = PipelineContext()
    context.columns["MergerTrees/main_branch"] = []
    with pytest.raises(RuntimeError, match="no merger tree"):
        OrbitalPropertiesStage(epoch, 0).run(context)


def test_unresolvable_host_raises():
    epoch = _MockEpoch(None)  # track_of raises MergerTreeError
    context = PipelineContext()
    context.columns["MergerTrees/main_branch"] = []
    with pytest.raises(RuntimeError, match="host halo"):
        OrbitalPropertiesStage(epoch, 0).run(context)


def test_check_inputs_requires_main_branch():
    epoch = _MockEpoch(None)
    with pytest.raises(RuntimeError, match="orbital_properties"):
        OrbitalPropertiesStage(epoch, 0).check_inputs(PipelineContext())


def test_boxsize_wraps_periodic_separation():
    # host at x=1, satellite at x=99, box=100 -- true separation is 2
    # (wrapping around the periodic boundary), not 98.
    host = _track([10], pos=[[1.0, 0.0, 0.0]])
    sat = _track([10], pos=[[99.0, 0.0, 0.0]])

    result = _run([sat], host, boxsize=100.0)
    assert _col(result, "OrbitalPericentre")[0] == pytest.approx(2.0)


# ---------------------------------------------------------------------------
# Integration: full real-data pipeline through OrbitalPropertiesStage
# ---------------------------------------------------------------------------

@needs_data
def test_full_pipeline_produces_sane_shapes_and_ordering():
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
    context = OrbitalPropertiesStage(epoch, host_row).run(context)

    n = len(satellite_rows)
    for name in ("OrbitalPericentre", "OrbitalApocentre",
                 "OrbitalEccentricity", "OrbitalAngularMomentum"):
        col = _col(context, name)
        assert col.shape[0] == n

    peri = _col(context, "OrbitalPericentre")
    apo = _col(context, "OrbitalApocentre")
    ecc = _col(context, "OrbitalEccentricity")
    resolved = ~np.isnan(peri)
    assert np.all(peri[resolved] <= apo[resolved])
    assert np.all(ecc[resolved] >= 0.0) and np.all(ecc[resolved] <= 1.0)
