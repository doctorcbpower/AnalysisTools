"""
Tests for analysistools.catalogue.pipeline.CatalogueBuilder,
HostJob, and _concatenate_contexts (Phase 6c).

_load_config/_concatenate_contexts are tested directly against small
hand-built fixtures; run()/run_from_stage() are exercised end to end
against the bundled real SWIFT/VELOCIraptor/tree data with
HydroGalaxyBackend (no bundled SHARK catalogue -- see
docs/phase6_remaining_work.md), using two disjoint host systems carved out
of the same halo catalogue to exercise multi-host concatenation.
"""
import os

import h5py
import numpy as np
import pytest
import yaml

import analysistools as at
from analysistools.catalogue.pipeline import (
    CatalogueBuilder, HostJob, PipelineContext, QualityControlStage,
    _concatenate_contexts,
)

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SNAP = os.path.join(HERE, "data", "snap_0031.hdf5")
VRCAT = os.path.join(HERE, "data", "halos",
                     "snap_0031.VELOCIraptor.properties.0")
TREE = os.path.join(HERE, "data", "VELOCIraptor.walkabletree.hdf5")
have_data = all(os.path.exists(p) for p in (SNAP, VRCAT, TREE))
needs_data = pytest.mark.skipif(not have_data, reason="example data missing")


# ---------------------------------------------------------------------------
# _load_config
# ---------------------------------------------------------------------------

def _write_config(tmp_path, **overrides):
    config = {
        "schema_version": "1.0",
        "galaxy_backend": "hydro",
        "galaxy_backend_options": {},
        "derived_stages": ["halo_properties"],
    }
    config.update(overrides)
    path = tmp_path / "project.yaml"
    path.write_text(yaml.safe_dump(config))
    return str(path)


def test_load_config_resolves_schema_backend_and_stages(tmp_path):
    path = _write_config(tmp_path, derived_stages=["halo_properties",
                                                    "orbital_properties"])
    builder = CatalogueBuilder(path)
    builder._load_config()

    assert builder.schema.version == "1.0"
    assert builder._galaxy_backend.name == "hydro"
    assert builder.stages == ["halo_properties", "orbital_properties"]


def test_load_config_raises_on_missing_required_key(tmp_path):
    path = _write_config(tmp_path)
    config = yaml.safe_load(open(path))
    del config["galaxy_backend"]
    open(path, "w").write(yaml.safe_dump(config))

    with pytest.raises(ValueError, match="galaxy_backend"):
        CatalogueBuilder(path)._load_config()


def test_load_config_raises_on_unknown_derived_stage(tmp_path):
    path = _write_config(tmp_path, derived_stages=["not_a_real_stage"])
    with pytest.raises(ValueError, match="not_a_real_stage"):
        CatalogueBuilder(path)._load_config()


def test_load_config_reads_stage_options(tmp_path):
    path = _write_config(
        tmp_path, derived_stages=["dorcha_specific"],
        stage_options={"dorcha_specific": {"reionisation_lookback_time": 13.0}})
    builder = CatalogueBuilder(path)
    builder._load_config()

    assert builder._stage_options == {
        "dorcha_specific": {"reionisation_lookback_time": 13.0}}


def test_load_config_defaults_stage_options_to_empty():
    builder = CatalogueBuilder("configs/dorcha.yaml")
    # the repo's real example config -- just check it loads without error
    # and every listed stage resolves.
    builder._load_config()
    assert builder.schema.version == "1.0"
    assert "environment" in builder.stages


# ---------------------------------------------------------------------------
# _concatenate_contexts
# ---------------------------------------------------------------------------

def _job_context(host_id, satellite_ids, mass, extra_galaxy_field=None):
    context = PipelineContext()
    context.columns["Haloes/HostHaloID"] = np.array([host_id])
    context.columns["Haloes/M200c"] = np.array([mass])
    n = len(satellite_ids)
    context.columns["Satellites/Identification/SatelliteID"] = \
        np.arange(n, dtype=np.int64)
    context.columns["Satellites/Identification/SubhaloID_z0"] = \
        np.asarray(satellite_ids)
    context.columns["MergerTrees/main_branch"] = [None] * n
    if extra_galaxy_field is not None:
        context.columns["Satellites/GalaxyProperties/BlackHoleMass"] = \
            np.asarray(extra_galaxy_field)
    context.record_stage("halo_extract", host_row=0, n_satellites=n)
    return context


def test_concatenate_stacks_haloes_and_satellites():
    c1 = _job_context(host_id=1, satellite_ids=[10, 11], mass=1e12)
    c2 = _job_context(host_id=2, satellite_ids=[20, 21, 22], mass=2e12)

    combined = _concatenate_contexts([c1, c2])

    np.testing.assert_array_equal(combined.columns["Haloes/HostHaloID"],
                                  [1, 2])
    np.testing.assert_allclose(combined.columns["Haloes/M200c"],
                               [1e12, 2e12])
    np.testing.assert_array_equal(
        combined.columns["Satellites/Identification/SubhaloID_z0"],
        [10, 11, 20, 21, 22])


def test_concatenate_renumbers_satellite_id_globally():
    c1 = _job_context(host_id=1, satellite_ids=[10, 11], mass=1e12)
    c2 = _job_context(host_id=2, satellite_ids=[20, 21, 22], mass=2e12)

    combined = _concatenate_contexts([c1, c2])
    np.testing.assert_array_equal(
        combined.columns["Satellites/Identification/SatelliteID"],
        [0, 1, 2, 3, 4])


def test_concatenate_merges_list_columns():
    c1 = _job_context(host_id=1, satellite_ids=[10, 11], mass=1e12)
    c2 = _job_context(host_id=2, satellite_ids=[20, 21, 22], mass=2e12)

    combined = _concatenate_contexts([c1, c2])
    assert combined.columns["MergerTrees/main_branch"] == [None] * 5


def test_concatenate_merges_provenance_and_meta():
    c1 = _job_context(host_id=1, satellite_ids=[10], mass=1e12)
    c1.meta["comoving_little_h"] = {"Haloes": {"comoving": True}}
    c2 = _job_context(host_id=2, satellite_ids=[20], mass=2e12)

    combined = _concatenate_contexts([c1, c2])
    assert len(combined.provenance) == 2
    assert combined.meta["comoving_little_h"] == {"Haloes": {"comoving": True}}


def test_concatenate_backfills_missing_satellite_field_with_nan():
    c1 = _job_context(host_id=1, satellite_ids=[10, 11], mass=1e12,
                      extra_galaxy_field=[5e6, 6e6])
    c2 = _job_context(host_id=2, satellite_ids=[20], mass=2e12)  # no BH field

    combined = _concatenate_contexts([c1, c2])
    bh = combined.columns["Satellites/GalaxyProperties/BlackHoleMass"]
    np.testing.assert_allclose(bh[:2], [5e6, 6e6])
    assert np.isnan(bh[2])


def test_concatenate_raises_on_inconsistent_haloes_columns():
    c1 = _job_context(host_id=1, satellite_ids=[10], mass=1e12)
    c2 = _job_context(host_id=2, satellite_ids=[20], mass=2e12)
    c2.columns["Haloes/R200c"] = np.array([100.0])  # c1 doesn't have this

    with pytest.raises(RuntimeError, match="Haloes/R200c"):
        _concatenate_contexts([c1, c2])


def test_concatenate_raises_on_empty_list():
    with pytest.raises(ValueError, match="empty"):
        _concatenate_contexts([])


def test_concatenate_raises_without_satellite_id():
    context = PipelineContext()
    context.columns["Haloes/HostHaloID"] = np.array([1])
    with pytest.raises(RuntimeError, match="SatelliteID"):
        _concatenate_contexts([context])


def test_concatenate_single_job_still_renumbers():
    c1 = _job_context(host_id=1, satellite_ids=[10, 11], mass=1e12)
    combined = _concatenate_contexts([c1])
    np.testing.assert_array_equal(
        combined.columns["Satellites/Identification/SatelliteID"], [0, 1])


# ---------------------------------------------------------------------------
# run() / run_from_stage() -- against real bundled data
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def epoch():
    sim = at.Simulation(
        snapshots={0.0: SNAP}, halos={0.0: VRCAT}, trees=TREE,
        snapnums={0.0: 31}, label="test",
        load_kwargs={"halos": {"native_includes_h": False}},
    )
    return sim.at(0.0)


@pytest.fixture
def two_host_jobs(epoch):
    mass = np.asarray(epoch.halos["mass"])
    order = np.argsort(mass)[::-1]
    job_a = HostJob(epoch=epoch, host_row=int(order[0]),
                    satellite_rows=[int(r) for r in order[1:4]])
    job_b = HostJob(epoch=epoch, host_row=int(order[4]),
                    satellite_rows=[int(r) for r in order[5:8]])
    return [job_a, job_b]


def _minimal_config(tmp_path, extra_stages=()):
    config = {
        "schema_version": "1.0",
        "galaxy_backend": "hydro",
        "galaxy_backend_options": {},
        "derived_stages": ["halo_properties", "orbital_properties",
                          *extra_stages],
        "stage_options": {
            "environment": {"mass_threshold": 1e8, "aperture_radius": 5.0},
            "host_environment": {
                "isolation_radius_factor": 3.0,
                "pairing_mass_ratio_min": 0.3,
                "pairing_max_separation": 1.0,
                "completeness_mass_threshold": 1e8},
            "observability": {"observer_pos": [0.0, 0.0, 0.0]},
            "star_formation_history": {
                "time_bin_edges": [0.0, 1.0, 2.0],
                "quenched_ssfr_threshold": 1e-11},
            "dorcha_specific": {"reionisation_lookback_time": 13.0},
        },
    }
    path = tmp_path / "project.yaml"
    path.write_text(yaml.safe_dump(config))
    return str(path)


@needs_data
def test_run_writes_catalogue_with_both_hosts(tmp_path, two_host_jobs):
    config_path = _minimal_config(tmp_path)
    builder = CatalogueBuilder(config_path)
    out_path = str(tmp_path / "out.h5")

    report = builder.run(two_host_jobs, out_path=out_path)

    assert report.passed
    with h5py.File(out_path, "r") as f:
        assert f["Haloes/HostHaloID"].shape == (2,)
        n_sat = f["Satellites/Identification/SatelliteID"].shape[0]
        assert n_sat == 6  # 3 satellites per host x 2 hosts
        np.testing.assert_array_equal(
            f["Satellites/Identification/SatelliteID"][()], np.arange(6))


@needs_data
def test_run_writes_config_metadata(tmp_path, two_host_jobs):
    config_path = _minimal_config(tmp_path)
    config = yaml.safe_load(open(config_path))
    config["metadata"] = {"Metadata/catalogue_version": "1.2.3"}
    open(config_path, "w").write(yaml.safe_dump(config))

    out_path = str(tmp_path / "out.h5")
    CatalogueBuilder(config_path).run(two_host_jobs, out_path=out_path)

    with h5py.File(out_path, "r") as f:
        assert f["Metadata"].attrs["catalogue_version"] == "1.2.3"


@needs_data
def test_run_refuses_to_overwrite(tmp_path, two_host_jobs):
    config_path = _minimal_config(tmp_path)
    out_path = str(tmp_path / "out.h5")
    CatalogueBuilder(config_path).run(two_host_jobs, out_path=out_path)

    with pytest.raises(FileExistsError):
        CatalogueBuilder(config_path).run(two_host_jobs, out_path=out_path)


@needs_data
def test_run_uses_output_path_template_when_out_path_omitted(
        tmp_path, two_host_jobs, monkeypatch):
    monkeypatch.chdir(tmp_path)
    config = {
        "schema_version": "1.0", "galaxy_backend": "hydro",
        "derived_stages": ["halo_properties"],
        "output": {"path_template": "my_catalogue_v{version}.h5"},
    }
    (tmp_path / "project.yaml").write_text(yaml.safe_dump(config))

    CatalogueBuilder(str(tmp_path / "project.yaml")).run(two_host_jobs)
    assert (tmp_path / "my_catalogue_v1.0.h5").exists()


@needs_data
def test_run_raises_before_writing_on_validation_error(
        tmp_path, two_host_jobs, monkeypatch):
    def _boom(self, context):
        raise RuntimeError("Stage 'quality_control': 1 validation error(s)")
    monkeypatch.setattr(QualityControlStage, "run", _boom)

    config_path = _minimal_config(tmp_path)
    out_path = str(tmp_path / "out.h5")
    with pytest.raises(RuntimeError, match="validation error"):
        CatalogueBuilder(config_path).run(two_host_jobs, out_path=out_path)

    assert not os.path.exists(out_path)


@needs_data
def test_run_from_stage_resumes_single_job(tmp_path, epoch, two_host_jobs):
    config_path = _minimal_config(tmp_path)
    builder = CatalogueBuilder(config_path)
    builder._load_config()

    job = two_host_jobs[0]
    # cached context as if halo_extract/tree_extract/cross_match already
    # ran (the "edited orbital_properties.py, don't want to redo
    # extraction" dev-iteration scenario run_from_stage exists for)
    from analysistools.catalogue.pipeline import (
        CrossMatchStage, HaloExtractStage, TreeExtractStage,
    )
    context = HaloExtractStage(
        job.epoch, job.host_row, job.satellite_rows).run(PipelineContext())
    context = TreeExtractStage(job.epoch).run(context)
    cached = CrossMatchStage(
        galaxy_backend=builder._galaxy_backend, epoch=job.epoch
    ).run(context)

    result = builder.run_from_stage("halo_properties", cached, [job])
    assert "Satellites/HaloProperties/Mpeak" in result.columns
    assert "Satellites/HaloProperties/OrbitalPericentre" in result.columns


@needs_data
def test_run_from_stage_rejects_multiple_jobs(tmp_path, two_host_jobs):
    config_path = _minimal_config(tmp_path)
    builder = CatalogueBuilder(config_path)
    with pytest.raises(ValueError, match="single-host"):
        builder.run_from_stage("halo_properties", PipelineContext(),
                               two_host_jobs)


@needs_data
def test_run_from_stage_rejects_unknown_stage_name(tmp_path, two_host_jobs):
    config_path = _minimal_config(tmp_path)
    builder = CatalogueBuilder(config_path)
    with pytest.raises(ValueError, match="not_a_stage"):
        builder.run_from_stage("not_a_stage", PipelineContext(),
                               [two_host_jobs[0]])
