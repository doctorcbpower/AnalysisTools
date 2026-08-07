"""
Tests for the implemented parts of Phase 6b: HaloExtractStage,
TreeExtractStage, and CrossMatchStage in analysistools.catalogue.pipeline.

Uses the bundled real snapshot/halo-catalogue/tree test data (same fixture
pattern as test_api_phase3.py's make_sim(), including the
native_includes_h=False override -- this VELOCIraptor catalogue was built
against the bundled SWIFT snapshot and is h-free natively, not h-included
like HaloCatalogue's per-format default guess assumes). Host/satellite
selection is arbitrary row picks for structural testing -- there's no
physically meaningful "host" in this generic test data, only enough
structure to exercise the pipeline plumbing.
"""
import os
from types import SimpleNamespace

import numpy as np
import pytest

import analysistools as at
from analysistools.catalogue.pipeline import (
    CrossMatchStage, HaloExtractStage, PipelineContext, TreeExtractStage,
)
from analysistools.merger_tree_types import MergerTreeError

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SNAP = os.path.join(HERE, "data", "snap_0031.hdf5")
VRCAT = os.path.join(HERE, "data", "halos",
                     "snap_0031.VELOCIraptor.properties.0")
TREE = os.path.join(HERE, "data", "VELOCIraptor.walkabletree.hdf5")

have_data = all(os.path.exists(p) for p in (SNAP, VRCAT, TREE))
needs_data = pytest.mark.skipif(not have_data, reason="example data missing")


@pytest.fixture(scope="module")
def epoch():
    sim = at.Simulation(
        snapshots={0.0: SNAP},
        halos={0.0: VRCAT},
        trees=TREE,
        snapnums={0.0: 31},
        label="test",
        load_kwargs={"halos": {"native_includes_h": False}},
    )
    return sim.at(0.0)


@pytest.fixture
def host_and_satellite_rows(epoch):
    mass = np.asarray(epoch.halos["mass"])
    order = np.argsort(mass)[::-1]  # most massive first
    host_row = int(order[0])
    satellite_rows = [int(r) for r in order[1:6]]  # arbitrary next 5
    return host_row, satellite_rows


# ---------------------------------------------------------------------------
# HaloExtractStage
# ---------------------------------------------------------------------------

@needs_data
def test_halo_extract_populates_host_and_satellite_columns(
        epoch, host_and_satellite_rows):
    host_row, satellite_rows = host_and_satellite_rows
    context = HaloExtractStage(epoch, host_row, satellite_rows).run(
        PipelineContext())

    halo_id = np.asarray(epoch.halos["halo_id"])
    mass = np.asarray(epoch.halos["mass"])
    radius = np.asarray(epoch.halos["radius"])

    assert context.columns["Haloes/HostHaloID"] == [halo_id[host_row]]
    assert context.columns["Haloes/M200c"] == pytest.approx([mass[host_row]])
    assert context.columns["Haloes/R200c"] == pytest.approx([radius[host_row]])
    assert context.columns["Haloes/Position"].shape == (1, 3)

    n = len(satellite_rows)
    np.testing.assert_array_equal(
        context.columns["Satellites/Identification/SubhaloID_z0"],
        halo_id[satellite_rows])
    np.testing.assert_allclose(
        context.columns["Satellites/HaloProperties/M200c_z0"],
        mass[satellite_rows])
    np.testing.assert_allclose(
        context.columns["Satellites/HaloProperties/R200c_z0"],
        radius[satellite_rows])
    assert context.columns["Satellites/Identification/Snapshot"].shape == (n,)
    assert (context.columns["Satellites/Identification/Snapshot"] == 31).all()

    np.testing.assert_array_equal(
        context.columns["Satellites/_internal/halo_row"], satellite_rows)


@needs_data
def test_halo_extract_finds_vmax_and_velocity(epoch, host_and_satellite_rows):
    # this bundled catalogue does have Vmax/vel fields -- confirm the
    # candidate-list lookup finds them rather than silently omitting.
    host_row, satellite_rows = host_and_satellite_rows
    context = HaloExtractStage(epoch, host_row, satellite_rows).run(
        PipelineContext())

    assert "Satellites/HaloProperties/Vmax_z0" in context.columns
    assert context.columns["Satellites/HaloProperties/Vmax_z0"].shape == \
        (len(satellite_rows),)
    assert "Haloes/Velocity" in context.columns
    assert "Satellites/_internal/vel_z0" in context.columns


@needs_data
def test_halo_extract_records_provenance(epoch, host_and_satellite_rows):
    host_row, satellite_rows = host_and_satellite_rows
    context = HaloExtractStage(epoch, host_row, satellite_rows).run(
        PipelineContext())
    assert context.provenance == [
        {"stage": "halo_extract", "host_row": host_row,
         "n_satellites": len(satellite_rows)}
    ]


@needs_data
def test_halo_extract_records_comoving_little_h_meta(
        epoch, host_and_satellite_rows):
    # read by validation.SchemaValidator's little-h/comoving cross-check
    host_row, satellite_rows = host_and_satellite_rows
    context = HaloExtractStage(epoch, host_row, satellite_rows).run(
        PipelineContext())

    state = context.meta["comoving_little_h"]["Haloes"]
    assert state["comoving"] == epoch.halos.meta.get("comoving")
    assert state["little_h"] == epoch.halos.meta.get("little_h")


def test_halo_extract_raises_on_empty_satellite_rows():
    stage = HaloExtractStage(SimpleNamespace(halos=object()), 0, [])
    with pytest.raises(RuntimeError, match="satellite_rows is empty"):
        stage.run(PipelineContext())


def test_halo_extract_raises_when_no_halo_catalogue():
    stage = HaloExtractStage(SimpleNamespace(halos=None), 0, [1, 2])
    with pytest.raises(RuntimeError, match="no halo catalogue"):
        stage.run(PipelineContext())


# ---------------------------------------------------------------------------
# TreeExtractStage
# ---------------------------------------------------------------------------

@needs_data
def test_tree_extract_after_halo_extract(epoch, host_and_satellite_rows):
    host_row, satellite_rows = host_and_satellite_rows
    context = HaloExtractStage(epoch, host_row, satellite_rows).run(
        PipelineContext())
    context = TreeExtractStage(epoch).run(context)

    n = len(satellite_rows)
    tree_ids = context.columns["Satellites/Identification/MergerTreeID"]
    main_branch = context.columns["MergerTrees/main_branch"]

    assert tree_ids.shape == (n,)
    assert tree_ids.dtype == np.int64
    assert len(main_branch) == n
    # every entry is either a resolved track or the documented sentinel
    for tid, track in zip(tree_ids, main_branch):
        if track is None:
            assert tid == -1
        else:
            assert tid == int(track.track.halo_id)


@needs_data
def test_tree_extract_requires_halo_extract_first(epoch):
    with pytest.raises(RuntimeError, match="tree_extract"):
        TreeExtractStage(epoch).check_inputs(PipelineContext())


def test_tree_extract_raises_when_no_tree():
    stage = TreeExtractStage(SimpleNamespace(tree=None))
    with pytest.raises(RuntimeError, match="no merger tree"):
        stage.run(PipelineContext())


def test_tree_extract_handles_unresolvable_halo_gracefully():
    class _FailingTree:
        def track_of(self, index):
            raise MergerTreeError("not found")

    epoch = SimpleNamespace(tree=object(), track_of=_FailingTree().track_of)
    context = PipelineContext()
    context.columns["Satellites/_internal/halo_row"] = np.array([0, 1, 2])

    result = TreeExtractStage(epoch).run(context)
    tree_ids = result.columns["Satellites/Identification/MergerTreeID"]
    main_branch = result.columns["MergerTrees/main_branch"]

    assert (tree_ids == -1).all()
    assert main_branch == [None, None, None]
    assert result.provenance[0]["n_missing"] == 3


# ---------------------------------------------------------------------------
# CrossMatchStage
# ---------------------------------------------------------------------------

@needs_data
def test_cross_match_sorts_by_subhalo_id_consistently(
        epoch, host_and_satellite_rows):
    host_row, satellite_rows = host_and_satellite_rows
    context = HaloExtractStage(epoch, host_row, satellite_rows).run(
        PipelineContext())
    context = TreeExtractStage(epoch).run(context)

    # snapshot the pre-sort (subhalo_id -> mass, -> tree_id) mapping so we
    # can check the permutation kept every column in sync, not just each
    # one independently sorted
    pre_subhalo_id = np.array(
        context.columns["Satellites/Identification/SubhaloID_z0"])
    pre_mass = np.array(context.columns["Satellites/HaloProperties/M200c_z0"])
    pre_tree_id = np.array(
        context.columns["Satellites/Identification/MergerTreeID"])
    mass_by_id = dict(zip(pre_subhalo_id, pre_mass))
    tree_id_by_id = dict(zip(pre_subhalo_id, pre_tree_id))

    result = CrossMatchStage().run(context)

    subhalo_id = result.columns["Satellites/Identification/SubhaloID_z0"]
    mass = result.columns["Satellites/HaloProperties/M200c_z0"]
    tree_id = result.columns["Satellites/Identification/MergerTreeID"]

    assert (np.diff(subhalo_id) >= 0).all()  # sorted
    for sid, m, tid in zip(subhalo_id, mass, tree_id):
        assert m == pytest.approx(mass_by_id[sid])
        assert tid == tree_id_by_id[sid]

    n = len(satellite_rows)
    np.testing.assert_array_equal(
        result.columns["Satellites/Identification/SatelliteID"], np.arange(n))

    host_id = int(np.asarray(epoch.halos["halo_id"])[host_row])
    assert (result.columns["Satellites/Identification/HostHaloID"]
           == host_id).all()

    # host-level columns untouched
    assert len(result.columns["Haloes/HostHaloID"]) == 1


@needs_data
def test_cross_match_permutes_main_branch_list_too(
        epoch, host_and_satellite_rows):
    host_row, satellite_rows = host_and_satellite_rows
    context = HaloExtractStage(epoch, host_row, satellite_rows).run(
        PipelineContext())
    context = TreeExtractStage(epoch).run(context)

    pre_subhalo_id = list(
        context.columns["Satellites/Identification/SubhaloID_z0"])
    pre_main_branch = list(context.columns["MergerTrees/main_branch"])
    track_by_id = dict(zip(pre_subhalo_id, pre_main_branch))

    result = CrossMatchStage().run(context)
    subhalo_id = result.columns["Satellites/Identification/SubhaloID_z0"]
    main_branch = result.columns["MergerTrees/main_branch"]

    assert len(main_branch) == len(subhalo_id)
    for sid, track in zip(subhalo_id, main_branch):
        assert track is track_by_id[sid]


@needs_data
def test_cross_match_requires_halo_extract_first(epoch):
    with pytest.raises(RuntimeError, match="cross_match"):
        CrossMatchStage().check_inputs(PipelineContext())


@needs_data
def test_cross_match_warns_but_does_not_crash_with_galaxy_backend(
        epoch, host_and_satellite_rows, caplog):
    host_row, satellite_rows = host_and_satellite_rows
    context = HaloExtractStage(epoch, host_row, satellite_rows).run(
        PipelineContext())

    with caplog.at_level("WARNING"):
        CrossMatchStage(galaxy_backend=object()).run(context)

    assert any("galaxy_backend" in r.message for r in caplog.records)


# ---------------------------------------------------------------------------
# CrossMatchStage: galaxy_backend.galaxy_properties() wiring
# ---------------------------------------------------------------------------

class _FakeGalaxyBackend:
    def __init__(self, props_by_row):
        self.props_by_row = props_by_row
        self.calls = []

    def galaxy_properties(self, epoch, halo_row):
        self.calls.append(halo_row)
        return self.props_by_row.get(halo_row, {})


@needs_data
def test_galaxy_properties_populate_columns_in_final_row_order(
        epoch, host_and_satellite_rows):
    host_row, satellite_rows = host_and_satellite_rows
    props_by_row = {row: {"StellarMass": float(1000 + row)}
                    for row in satellite_rows}
    backend = _FakeGalaxyBackend(props_by_row)

    context = HaloExtractStage(epoch, host_row, satellite_rows).run(
        PipelineContext())
    result = CrossMatchStage(galaxy_backend=backend, epoch=epoch).run(context)

    halo_rows = result.columns["Satellites/_internal/halo_row"]
    stellar_mass = result.columns["Satellites/GalaxyProperties/StellarMass"]
    for row, mass in zip(halo_rows, stellar_mass):
        assert mass == pytest.approx(props_by_row[int(row)]["StellarMass"])

    # called once per satellite, with the (post-sort) halo rows
    assert sorted(backend.calls) == sorted(int(r) for r in halo_rows)


@needs_data
def test_galaxy_properties_missing_fields_become_nan(
        epoch, host_and_satellite_rows):
    host_row, satellite_rows = host_and_satellite_rows
    # only the first satellite has GasMass_Cold
    props_by_row = {row: {"StellarMass": 1e9} for row in satellite_rows}
    props_by_row[satellite_rows[0]]["GasMass_Cold"] = 5e7
    backend = _FakeGalaxyBackend(props_by_row)

    context = HaloExtractStage(epoch, host_row, satellite_rows).run(
        PipelineContext())
    result = CrossMatchStage(galaxy_backend=backend, epoch=epoch).run(context)

    halo_rows = result.columns["Satellites/_internal/halo_row"]
    gas_mass = result.columns["Satellites/GalaxyProperties/GasMass_Cold"]
    n_with_gas = sum(1 for row in halo_rows if int(row) == satellite_rows[0])
    assert n_with_gas == 1
    assert np.sum(~np.isnan(gas_mass)) == 1
    assert np.sum(np.isnan(gas_mass)) == len(satellite_rows) - 1


@needs_data
def test_galaxy_backend_without_epoch_warns_and_skips(
        epoch, host_and_satellite_rows, caplog):
    host_row, satellite_rows = host_and_satellite_rows
    backend = _FakeGalaxyBackend({})
    context = HaloExtractStage(epoch, host_row, satellite_rows).run(
        PipelineContext())

    with caplog.at_level("WARNING"):
        result = CrossMatchStage(galaxy_backend=backend).run(context)  # no epoch=

    assert backend.calls == []
    assert not any(k.startswith("Satellites/GalaxyProperties/")
                   for k in result.columns)


@needs_data
def test_cross_match_records_comoving_little_h_meta_when_backend_supports_it(
        epoch, host_and_satellite_rows):
    class _FakeBackendWithNativeState(_FakeGalaxyBackend):
        def native_comoving_little_h(self, epoch):
            return True, False

    host_row, satellite_rows = host_and_satellite_rows
    backend = _FakeBackendWithNativeState(
        {row: {"StellarMass": 1e9} for row in satellite_rows})
    context = HaloExtractStage(epoch, host_row, satellite_rows).run(
        PipelineContext())
    result = CrossMatchStage(galaxy_backend=backend, epoch=epoch).run(context)

    state = result.meta["comoving_little_h"]["Satellites/GalaxyProperties"]
    assert state == {"comoving": True, "little_h": False}


@needs_data
def test_cross_match_skips_comoving_little_h_meta_when_backend_lacks_method(
        epoch, host_and_satellite_rows):
    host_row, satellite_rows = host_and_satellite_rows
    backend = _FakeGalaxyBackend(  # no native_comoving_little_h method
        {row: {"StellarMass": 1e9} for row in satellite_rows})
    context = HaloExtractStage(epoch, host_row, satellite_rows).run(
        PipelineContext())
    result = CrossMatchStage(galaxy_backend=backend, epoch=epoch).run(context)

    assert "Satellites/GalaxyProperties" not in \
        result.meta.get("comoving_little_h", {})


def test_galaxy_backend_without_halo_extract_first_warns_and_skips(caplog):
    backend = _FakeGalaxyBackend({})
    context = PipelineContext()
    context.columns["Satellites/Identification/SubhaloID_z0"] = \
        np.array([2, 1])
    context.columns["Haloes/HostHaloID"] = np.array([99])

    with caplog.at_level("WARNING"):
        result = CrossMatchStage(
            galaxy_backend=backend, epoch=object()).run(context)

    assert backend.calls == []
    assert any("halo_row" in r.message for r in caplog.records)
    assert not any(k.startswith("Satellites/GalaxyProperties/")
                   for k in result.columns)
