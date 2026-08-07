"""
Tests for analysistools.catalogue.derived.StarFormationHistoryStage
(Phase 6b).

Uses a fake GalaxyBackend exposing only star_formation_history() (real
per-backend SFH computation is tested separately in
test_catalogue_backends_shark_sfh.py / test_catalogue_backends_hydro.py)
so the stage's own logic -- SFH stacking, MeanStellarAge weighting, the
IsQuenched_z0/QuenchingTime walk -- is exactly hand-checkable.
"""
import os

import numpy as np
import pytest

import analysistools as at
from analysistools.catalogue.backends import HydroGalaxyBackend
from analysistools.catalogue.derived import (
    HaloPropertiesStage, StarFormationHistoryStage,
)
from analysistools.catalogue.pipeline import (
    CrossMatchStage, HaloExtractStage, PipelineContext, TreeExtractStage,
)

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SNAP = os.path.join(HERE, "data", "snap_0031.hdf5")
VRCAT = os.path.join(HERE, "data", "halos",
                     "snap_0031.VELOCIraptor.properties.0")
TREE = os.path.join(HERE, "data", "VELOCIraptor.walkabletree.hdf5")
have_data = all(os.path.exists(p) for p in (SNAP, VRCAT, TREE))
needs_data = pytest.mark.skipif(not have_data, reason="example data missing")

TIME_BIN_EDGES = np.array([0.0, 1.0, 2.0, 3.0])  # Gyr lookback
BIN_CENTRES = np.array([0.5, 1.5, 2.5])


class _FakeSFHBackend:
    def __init__(self, sfh_by_row):
        self.sfh_by_row = sfh_by_row
        self.calls = []

    def star_formation_history(self, epoch, halo_row, time_bin_edges):
        self.calls.append((halo_row, tuple(time_bin_edges)))
        return self.sfh_by_row.get(halo_row)


def _context(halo_rows, stellar_mass=None, sfr=None):
    context = PipelineContext()
    context.columns["Satellites/_internal/halo_row"] = np.asarray(halo_rows)
    if stellar_mass is not None:
        context.columns["Satellites/GalaxyProperties/StellarMass"] = \
            np.asarray(stellar_mass, dtype=float)
    if sfr is not None:
        context.columns["Satellites/GalaxyProperties/StarFormationRate"] = \
            np.asarray(sfr, dtype=float)
    return context


def _col(result, name):
    return result.columns[f"Satellites/GalaxyProperties/{name}"]


def test_sfh_stacked_with_correct_shape_and_values():
    backend = _FakeSFHBackend({0: [1.0, 2.0, 3.0], 1: [4.0, 5.0, 6.0]})
    context = _context([0, 1])
    stage = StarFormationHistoryStage(epoch=object(), galaxy_backend=backend,
                                      time_bin_edges=TIME_BIN_EDGES,
                                      quenched_ssfr_threshold=1e-11)
    result = stage.run(context)

    sfh = _col(result, "SFH")
    assert sfh.shape == (2, 3)
    np.testing.assert_allclose(sfh[0], [1.0, 2.0, 3.0])
    np.testing.assert_allclose(sfh[1], [4.0, 5.0, 6.0])

    # forwarded the exact time_bin_edges given to the stage
    for _, edges in backend.calls:
        np.testing.assert_allclose(edges, TIME_BIN_EDGES)


def test_time_bin_edges_recorded_in_context_meta():
    # DorchaSpecificStage's FossilFraction reads this rather than taking
    # its own, independently-specified (and possibly mismatched) grid.
    backend = _FakeSFHBackend({0: [1.0, 2.0, 3.0]})
    context = _context([0])
    stage = StarFormationHistoryStage(epoch=object(), galaxy_backend=backend,
                                      time_bin_edges=TIME_BIN_EDGES,
                                      quenched_ssfr_threshold=1e-11)
    result = stage.run(context)

    np.testing.assert_array_equal(result.meta["time_bin_edges_sfh"],
                                  TIME_BIN_EDGES)


def test_mean_stellar_age_is_mass_weighted_over_bins():
    # all formed mass in the middle bin -> mean age = that bin's centre
    backend = _FakeSFHBackend({0: [0.0, 5.0, 0.0]})
    context = _context([0])
    stage = StarFormationHistoryStage(epoch=object(), galaxy_backend=backend,
                                      time_bin_edges=TIME_BIN_EDGES,
                                      quenched_ssfr_threshold=1e-11)
    result = stage.run(context)

    assert _col(result, "MeanStellarAge")[0] == pytest.approx(1.5)


def test_no_sfh_leaves_nan_sentinels_and_is_counted():
    backend = _FakeSFHBackend({0: None})
    context = _context([0])
    stage = StarFormationHistoryStage(epoch=object(), galaxy_backend=backend,
                                      time_bin_edges=TIME_BIN_EDGES,
                                      quenched_ssfr_threshold=1e-11)
    result = stage.run(context)

    assert np.all(np.isnan(_col(result, "SFH")[0]))
    assert np.isnan(_col(result, "MeanStellarAge")[0])
    assert np.isnan(_col(result, "QuenchingTime")[0])
    assert _col(result, "IsQuenched_z0")[0] == False  # noqa: E712
    assert result.provenance[0]["n_no_sfh"] == 1


def test_is_quenched_true_when_current_sfr_below_absolute_threshold():
    backend = _FakeSFHBackend({0: [0.01, 0.02, 5.0]})
    context = _context([0], stellar_mass=[1e10], sfr=[1e-3])
    stage = StarFormationHistoryStage(epoch=object(), galaxy_backend=backend,
                                      time_bin_edges=TIME_BIN_EDGES,
                                      quenched_ssfr_threshold=1e-11)
    result = stage.run(context)

    # sSFR_z0 = 1e-3 / 1e10 = 1e-13 < 1e-11 -> quenched
    assert _col(result, "IsQuenched_z0")[0] == True  # noqa: E712
    # ssfr per bin = [0.01,0.02,5.0]/1e10 -> threshold*mstar = 0.1 Msun/yr
    # only bin 2 (5.0) is >= 0.1 -> most recent active bin is index 2
    assert _col(result, "QuenchingTime")[0] == pytest.approx(BIN_CENTRES[2])


def test_is_quenched_false_when_above_threshold_no_quenching_time():
    backend = _FakeSFHBackend({0: [1.0, 1.0, 1.0]})
    context = _context([0], stellar_mass=[1e10], sfr=[1.0])
    stage = StarFormationHistoryStage(epoch=object(), galaxy_backend=backend,
                                      time_bin_edges=TIME_BIN_EDGES,
                                      quenched_ssfr_threshold=1e-11)
    result = stage.run(context)

    assert _col(result, "IsQuenched_z0")[0] == False  # noqa: E712
    assert np.isnan(_col(result, "QuenchingTime")[0])


def test_quenching_skipped_without_stellar_mass_or_sfr_in_context():
    backend = _FakeSFHBackend({0: [0.01, 0.01, 0.01]})
    context = _context([0])  # no StellarMass/StarFormationRate columns
    stage = StarFormationHistoryStage(epoch=object(), galaxy_backend=backend,
                                      time_bin_edges=TIME_BIN_EDGES,
                                      quenched_ssfr_threshold=1e-11)
    result = stage.run(context)

    assert _col(result, "IsQuenched_z0")[0] == False  # noqa: E712
    assert np.isnan(_col(result, "QuenchingTime")[0])
    # SFH/MeanStellarAge are still computed independently of quenching
    assert not np.all(np.isnan(_col(result, "SFH")[0]))


def test_check_inputs_requires_halo_row_and_stellar_mass():
    backend = _FakeSFHBackend({})
    stage = StarFormationHistoryStage(epoch=object(), galaxy_backend=backend,
                                      time_bin_edges=TIME_BIN_EDGES,
                                      quenched_ssfr_threshold=1e-11)
    with pytest.raises(RuntimeError, match="star_formation_history"):
        stage.check_inputs(PipelineContext())


# ---------------------------------------------------------------------------
# Integration: full real-data pipeline with a real backend (HydroGalaxyBackend)
# ---------------------------------------------------------------------------

@needs_data
def test_full_pipeline_with_real_hydro_backend_does_not_crash():
    # The bundled snapshot is DMO (no star particles), so HydroGalaxyBackend
    # will find nothing for every satellite -- this checks the full chain
    # wires together and degrades to all-NaN/False rather than crashing,
    # not that it produces non-trivial SFH values (there's no star data in
    # this repo's test fixtures to produce those from).
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

    n = len(satellite_rows)
    sfh = _col(context, "SFH")
    assert sfh.shape == (n, len(TIME_BIN_EDGES) - 1)
    assert np.all(np.isnan(sfh))  # no star particles anywhere in this snapshot
    assert np.all(_col(context, "IsQuenched_z0") == False)  # noqa: E712
