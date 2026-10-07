"""
Tests for analysistools.catalogue.derived.OrbitalPropertiesStage.

The stage no longer measures orbits itself: it maps an
``analysistools.orbits`` sidecar (or a tracker run) into catalogue columns.
The orbit physics is tested in tests/test_orbits.py against analytic and
brute-force answers; these tests cover the wiring -- id matching, fill
values, schema consistency, the explicit-configuration requirement -- plus
one end-to-end run against the bundled example data.

(The previous version of this file pinned the old behaviour, e.g.
``OrbitalApocentre == max(distance over the whole tree overlap)``, which
reported a satellite's pre-infall distance as its apocentre; that
behaviour was wrong and is deliberately not preserved.)
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
from analysistools.catalogue.schema import default_schema
from analysistools.merger_tree_types import MergerTreeError
from analysistools.orbits import OrbitConfig, OrbitTracker, SCALARS
from analysistools.orbits.catalogue_fields import (
    EXTENDED_FIELD_MAP, GROUP, LEGACY_FIELD_MAP)

import test_orbits as TO

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SNAP = os.path.join(HERE, "data", "snap_0031.hdf5")
VRCAT = os.path.join(HERE, "data", "halos",
                     "snap_0031.VELOCIraptor.properties.0")
TREE = os.path.join(HERE, "data", "VELOCIraptor.walkabletree.hdf5")
have_data = all(os.path.exists(p) for p in (SNAP, VRCAT, TREE))
needs_data = pytest.mark.skipif(not have_data, reason="example data missing")


class _MockEpoch:
    def __init__(self, host_track, boxsize=None):
        self._host = host_track
        self.tree = object() if host_track is not None else None
        self.boxsize = boxsize

    def track_of(self, index):
        if self._host is None:
            raise MergerTreeError("host not found")
        return TrackDataset(self._host)


def _sidecar_and_stage_inputs(extended=False):
    """Two resolved satellites (one backsplash flyby, one bound orbit) and
    one unresolved, plus the tracker that made them."""
    t, xs, vs = TO._truth_orbit(np.array([300.0, 0, 0]),
                                np.array([-20.0, 70.0, 0]), 9.0, 0.17)
    host = TO._host(t)
    bound = TO._sat(host, xs, vs, halo_id=101)
    n = len(t)
    fly = np.zeros((n, 3)); fly[:, 0] = -700.0 + 200.0 * TO.KMS * t; fly[:, 1] = 50
    flyby = TO._sat(host, fly, np.tile([200.0, 0, 0], (n, 1)), halo_id=102)
    pot = TO.TimeDependentPotential([TO.T0 + t[0], TO.T0 + t[-1]],
                                    [TO.NFW, TO.NFW])
    tracker = OrbitTracker(TO.COSMO, TO.PHYS, OrbitConfig(), potential=pot)
    return host, tracker, bound, flyby


def _context(*tracks):
    c = PipelineContext()
    c.columns["MergerTrees/main_branch"] = [
        None if tr is None else TrackDataset(tr) for tr in tracks]
    return c


def _col(result, name):
    return result.columns[f"{GROUP}/{name}"]


def test_requires_exactly_one_source():
    host, tracker, bound, _ = _sidecar_and_stage_inputs()
    with pytest.raises(ValueError, match="exactly one"):
        OrbitalPropertiesStage(_MockEpoch(host), 0)
    with pytest.raises(ValueError, match="exactly one"):
        OrbitalPropertiesStage(_MockEpoch(host), 0, tracker=tracker,
                               tracker_options={"length_kpc": 1.0})


def test_tracker_mode_fills_columns_and_matches_by_halo_id():
    host, tracker, bound, flyby = _sidecar_and_stage_inputs()
    ctx = _context(flyby, None, bound)
    result = OrbitalPropertiesStage(_MockEpoch(host), 0, tracker=tracker) \
        .run(ctx)

    peri = _col(result, "OrbitalPericentre")
    # flyby: straight line, impact parameter 50 kpc
    assert peri[0] == pytest.approx(50.0, rel=0.01)
    assert np.isnan(peri[1])                         # unresolved -> fill
    # bound orbit: one pericentre ~84 kpc (see test_orbits)
    assert peri[2] == pytest.approx(84.0, rel=0.03)
    assert _col(result, "OrbitalAngularMomentum").shape == (3, 3)
    assert np.all(np.isnan(_col(result, "OrbitalAngularMomentum")[1]))
    assert result.provenance[0]["n_no_track"] == 1


def test_apocentre_is_confirmed_turning_point_not_history_maximum():
    """Regression for the old stage: a satellite that was 5 Mpc away early
    on must not report that as its apocentre."""
    host, tracker, bound, _ = _sidecar_and_stage_inputs()
    bound.Pos[:3] = [[5000.0, 0, 0], [3000.0, 0, 0], [1200.0, 0, 0]]
    result = OrbitalPropertiesStage(_MockEpoch(host), 0, tracker=tracker) \
        .run(_context(bound))
    apo = _col(result, "OrbitalApocentre")[0]
    assert np.isnan(apo) or apo < 400.0


def test_orbit_catalogue_mode_and_unknown_ids(tmp_path):
    host, tracker, bound, flyby = _sidecar_and_stage_inputs()
    cat = tracker.run([bound], host, ids=[bound.halo_id])
    path = str(tmp_path / "orbits.h5")
    cat.write_hdf5(path)
    result = OrbitalPropertiesStage(None, 0, orbit_catalogue=path) \
        .run(_context(bound, flyby))
    assert not np.isnan(_col(result, "OrbitalPericentre")[0])
    assert np.isnan(_col(result, "OrbitalPericentre")[1])   # id 102 not in sidecar
    assert result.provenance[0]["n_not_in_sidecar"] == 1


def test_extended_fields_are_declared_in_schema_1_1_and_legacy_in_1_0():
    s10, s11 = default_schema("1.0"), default_schema("1.1")
    for col in LEGACY_FIELD_MAP:
        assert f"{GROUP}/{col}" in s10 and f"{GROUP}/{col}" in s11
    for col in EXTENDED_FIELD_MAP:
        assert f"{GROUP}/{col}" not in s10
        assert f"{GROUP}/{col}" in s11
    assert set(s10.fields) < set(s11.fields)          # additive only
    # declared dtype matches what the stage produces
    for col, key in EXTENDED_FIELD_MAP.items():
        assert s11[f"{GROUP}/{col}"].dtype == SCALARS[key][0]


def test_extended_stage_writes_every_declared_column():
    host, tracker, bound, flyby = _sidecar_and_stage_inputs()
    stage = OrbitalPropertiesStage(_MockEpoch(host), 0, tracker=tracker,
                                   extended_fields=True)
    result = stage.run(_context(bound, flyby))
    for path in stage.outputs:
        assert path in result.columns
    assert bool(_col(result, "OrbitIsBacksplash")[1])
    assert not bool(_col(result, "OrbitIsBacksplash")[0])
    assert _col(result, "OrbitNumberOfInfalls").dtype == np.int16


def test_no_tree_raises():
    host, tracker, *_ = _sidecar_and_stage_inputs()
    epoch = _MockEpoch(None)
    ctx = PipelineContext(); ctx.columns["MergerTrees/main_branch"] = []
    with pytest.raises(RuntimeError, match="no merger tree"):
        OrbitalPropertiesStage(epoch, 0, tracker=tracker).run(ctx)


def test_unresolvable_host_raises():
    host, tracker, *_ = _sidecar_and_stage_inputs()
    epoch = _MockEpoch(host); epoch._host = None; epoch.tree = object()
    ctx = PipelineContext(); ctx.columns["MergerTrees/main_branch"] = []
    with pytest.raises(RuntimeError, match="host halo"):
        OrbitalPropertiesStage(epoch, 0, tracker=tracker).run(ctx)


def test_check_inputs_requires_main_branch():
    host, tracker, *_ = _sidecar_and_stage_inputs()
    with pytest.raises(RuntimeError, match="orbital_properties"):
        OrbitalPropertiesStage(_MockEpoch(host), 0, tracker=tracker) \
            .check_inputs(PipelineContext())


# ---------------------------------------------------------------------------
# Integration: bundled real data (TreeFrog walkable tree, no host radius in
# the tree -> an explicit fixed reference radius and cosmology are supplied,
# which is exactly the explicitness the stage insists on)
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
    stage = OrbitalPropertiesStage(
        epoch, host_row,
        tracker_options={"length_kpc": 1000.0, "mass_msun": 1.0e10,
                         "potential": None,
                         "cosmology": {"h": 0.7, "omega_m": 0.3},
                         "fixed_radius_kpc": 300.0, "orphans": False},
        extended_fields=True)
    context = stage.run(context)

    n = len(satellite_rows)
    for path in stage.outputs:
        assert context.columns[path].shape[0] == n

    # wiring check only: the bundled "satellites" are just the next-most-
    # massive unrelated halos (see HaloExtractStage), not real subhaloes of
    # the host, so no physical orbit is expected -- but every track must at
    # least have produced a finite host distance
    assert np.all(np.isfinite(_col(context, "OrbitFinalRadius")))

    peri = _col(context, "OrbitalPericentre")
    apo = _col(context, "OrbitalApocentre")
    ecc = _col(context, "OrbitalEccentricity")
    both = ~np.isnan(peri) & ~np.isnan(apo)
    assert np.all(peri[both] <= apo[both])
    e_ok = ~np.isnan(ecc)
    assert np.all(ecc[e_ok] >= 0.0) and np.all(ecc[e_ok] <= 1.0)
