"""
Tests for analysistools.catalogue.pipeline.WriteStage (Phase 6c).

Uses small hand-built CatalogueSchema/PipelineContext fixtures and reads
the written file back with h5py directly, so every check is against the
actual on-disk bytes rather than a mocked writer.
"""
import json

import h5py
import numpy as np
import pytest

from analysistools.catalogue.pipeline import PipelineContext, WriteStage
from analysistools.catalogue.schema import CatalogueSchema, FieldSpec
from analysistools.catalogue.validation import ValidationReport


def _schema():
    schema = CatalogueSchema(version="1.0")
    schema.add(FieldSpec(name="M200c", group="Haloes", dtype="float64",
                         units="Msun", description="Host virial mass.",
                         provenance="halo catalogue"))
    schema.add(FieldSpec(name="SatelliteID", group="Satellites/Identification",
                         dtype="int64", units="dimensionless"))
    schema.add(FieldSpec(name="Mpeak", group="Satellites/HaloProperties",
                         dtype="float64", units="Msun", is_derived=True,
                         provenance="derived from merger tree"))
    schema.add(FieldSpec(name="SatelliteName", group="Satellites/Identification",
                         dtype="str", units="dimensionless"))
    schema.add(FieldSpec(name="IsBacksplash", group="Satellites/HaloProperties",
                         dtype="bool", units="dimensionless"))
    schema.add(FieldSpec(name="Position", group="Haloes", dtype="float64",
                         units="Mpc", shape=(3,)))
    schema.add(FieldSpec(name="catalogue_version", group="Metadata",
                         dtype="str", units="dimensionless",
                         is_attribute=True))
    return schema


def _context():
    context = PipelineContext()
    context.columns["Haloes/HostHaloID"] = np.array([42])
    context.columns["Haloes/M200c"] = np.array([1e12])
    context.columns["Haloes/Position"] = np.array([[1.0, 2.0, 3.0]])
    context.columns["Satellites/Identification/SatelliteID"] = \
        np.arange(3, dtype=np.int64)
    context.columns["Satellites/HaloProperties/Mpeak"] = \
        np.array([1e9, 2e9, 3e9])
    context.columns["Satellites/HaloProperties/IsBacksplash"] = \
        np.array([True, False, True])
    context.columns["Satellites/Identification/SatelliteName"] = \
        np.array(["", "Fornax", ""])
    context.record_stage("halo_extract", host_row=0, n_satellites=3)
    return context


# ---------------------------------------------------------------------------
# Basic write + round trip
# ---------------------------------------------------------------------------

def test_write_creates_file_with_declared_datasets(tmp_path):
    out_path = str(tmp_path / "cat.h5")
    WriteStage(out_path, _schema()).run(_context())

    with h5py.File(out_path, "r") as f:
        assert f["Haloes/M200c"][()] == pytest.approx([1e12])
        np.testing.assert_allclose(f["Haloes/Position"][()], [[1.0, 2.0, 3.0]])
        np.testing.assert_array_equal(
            f["Satellites/Identification/SatelliteID"][()], [0, 1, 2])
        np.testing.assert_allclose(
            f["Satellites/HaloProperties/Mpeak"][()], [1e9, 2e9, 3e9])
        np.testing.assert_array_equal(
            f["Satellites/HaloProperties/IsBacksplash"][()],
            [True, False, True])


def test_string_dataset_round_trips(tmp_path):
    out_path = str(tmp_path / "cat.h5")
    WriteStage(out_path, _schema()).run(_context())

    with h5py.File(out_path, "r") as f:
        names = f["Satellites/Identification/SatelliteName"].asstr()[()]
        assert list(names) == ["", "Fornax", ""]


def test_attribute_field_written_as_group_attribute(tmp_path):
    out_path = str(tmp_path / "cat.h5")
    context = _context()
    context.columns["Metadata/catalogue_version"] = "1.0.0"
    WriteStage(out_path, _schema()).run(context)

    with h5py.File(out_path, "r") as f:
        assert f["Metadata"].attrs["catalogue_version"] == "1.0.0"
        assert "catalogue_version" not in f  # not a dataset


def test_dataset_attrs_carry_units_description_provenance_is_derived(tmp_path):
    out_path = str(tmp_path / "cat.h5")
    WriteStage(out_path, _schema()).run(_context())

    with h5py.File(out_path, "r") as f:
        dset = f["Satellites/HaloProperties/Mpeak"]
        assert dset.attrs["units"] == "Msun"
        assert dset.attrs["is_derived"] == True  # noqa: E712
        assert dset.attrs["provenance"] == "derived from merger tree"
        assert f["Haloes/M200c"].attrs["is_derived"] == False  # noqa: E712


# ---------------------------------------------------------------------------
# Omission of undeclared / absent fields
# ---------------------------------------------------------------------------

def test_undeclared_columns_are_not_written(tmp_path):
    out_path = str(tmp_path / "cat.h5")
    context = _context()
    context.columns["Satellites/_internal/halo_row"] = np.array([0, 1, 2])
    context.columns["MergerTrees/main_branch"] = [None, None, None]
    WriteStage(out_path, _schema()).run(context)

    with h5py.File(out_path, "r") as f:
        assert "Satellites" in f
        assert "_internal" not in f["Satellites"]
        assert "MergerTrees" not in f


def test_absent_declared_field_is_simply_omitted(tmp_path):
    out_path = str(tmp_path / "cat.h5")
    context = _context()
    del context.columns["Haloes/Position"]  # declared in schema, not given
    WriteStage(out_path, _schema()).run(context)

    with h5py.File(out_path, "r") as f:
        assert "Position" not in f["Haloes"]


# ---------------------------------------------------------------------------
# dtype casting
# ---------------------------------------------------------------------------

def test_dtype_is_cast_to_schema_declared_dtype(tmp_path):
    out_path = str(tmp_path / "cat.h5")
    context = _context()
    context.columns["Satellites/HaloProperties/Mpeak"] = \
        np.array([1, 2, 3], dtype=np.int32)  # wrong dtype, castable
    WriteStage(out_path, _schema()).run(context)

    with h5py.File(out_path, "r") as f:
        assert f["Satellites/HaloProperties/Mpeak"].dtype == np.float64


# ---------------------------------------------------------------------------
# Documentation / Provenance
# ---------------------------------------------------------------------------

def test_documentation_schema_json_matches_schema(tmp_path):
    out_path = str(tmp_path / "cat.h5")
    schema = _schema()
    WriteStage(out_path, schema).run(_context())

    with h5py.File(out_path, "r") as f:
        written = json.loads(f["Documentation/schema.json"][()])
        assert written == schema.to_json_schema_dict()


def test_provenance_pipeline_stages_recorded(tmp_path):
    out_path = str(tmp_path / "cat.h5")
    WriteStage(out_path, _schema()).run(_context())

    with h5py.File(out_path, "r") as f:
        stages = json.loads(f["Provenance/pipeline_stages"][()])
        # WriteStage's own record_stage() call happens after the file is
        # closed (it needs n_fields_written, only known once writing is
        # done), so the file's own provenance log ends at whatever ran
        # before WriteStage -- it does not record itself finishing.
        assert stages == [{"stage": "halo_extract", "host_row": 0,
                          "n_satellites": 3}]


def test_provenance_validation_report_written_when_attached(tmp_path):
    out_path = str(tmp_path / "cat.h5")
    context = _context()
    report = ValidationReport()
    report.add("warning", "physical", "borderline value")
    context.meta["validation_report"] = report
    WriteStage(out_path, _schema()).run(context)

    with h5py.File(out_path, "r") as f:
        written = json.loads(f["Provenance/validation_report"][()])
        assert written["n_warnings"] == 1


def test_provenance_validation_report_absent_when_not_attached(tmp_path):
    out_path = str(tmp_path / "cat.h5")
    WriteStage(out_path, _schema()).run(_context())

    with h5py.File(out_path, "r") as f:
        assert "validation_report" not in f["Provenance"]


# ---------------------------------------------------------------------------
# Automatic Metadata
# ---------------------------------------------------------------------------

def test_n_satellites_n_hosts_schema_version_populated_automatically(tmp_path):
    schema = _schema()
    schema.add(FieldSpec(name="n_satellites", group="Metadata",
                         dtype="int64", units="dimensionless",
                         is_attribute=True))
    schema.add(FieldSpec(name="n_hosts", group="Metadata", dtype="int64",
                         units="dimensionless", is_attribute=True))
    schema.add(FieldSpec(name="schema_version", group="Metadata",
                         dtype="str", units="dimensionless",
                         is_attribute=True))
    out_path = str(tmp_path / "cat.h5")
    WriteStage(out_path, schema).run(_context())

    with h5py.File(out_path, "r") as f:
        assert f["Metadata"].attrs["n_satellites"] == 3
        assert f["Metadata"].attrs["n_hosts"] == 1
        assert f["Metadata"].attrs["schema_version"] == "1.0"


def test_automatic_metadata_does_not_override_explicit_value(tmp_path):
    schema = _schema()
    schema.add(FieldSpec(name="n_satellites", group="Metadata",
                         dtype="int64", units="dimensionless",
                         is_attribute=True))
    context = _context()
    context.columns["Metadata/n_satellites"] = 999  # caller override
    out_path = str(tmp_path / "cat.h5")
    WriteStage(out_path, schema).run(context)

    with h5py.File(out_path, "r") as f:
        assert f["Metadata"].attrs["n_satellites"] == 999


# ---------------------------------------------------------------------------
# time_bin_edges_sfh meta bridge
# ---------------------------------------------------------------------------

def test_time_bin_edges_sfh_bridged_from_meta(tmp_path):
    schema = _schema()
    schema.add(FieldSpec(name="time_bin_edges_sfh", group="Snapshots",
                         dtype="float64", units="Gyr", shape=("n_bins_plus_1",)))
    context = _context()
    context.meta["time_bin_edges_sfh"] = np.array([0.0, 1.0, 2.0])
    out_path = str(tmp_path / "cat.h5")
    WriteStage(out_path, schema).run(context)

    with h5py.File(out_path, "r") as f:
        np.testing.assert_allclose(
            f["Snapshots/time_bin_edges_sfh"][()], [0.0, 1.0, 2.0])


# ---------------------------------------------------------------------------
# Overwrite protection / atomic write
# ---------------------------------------------------------------------------

def test_refuses_to_overwrite_existing_file(tmp_path):
    out_path = str(tmp_path / "cat.h5")
    WriteStage(out_path, _schema()).run(_context())

    with pytest.raises(FileExistsError, match="refusing to overwrite"):
        WriteStage(out_path, _schema()).run(_context())


def test_no_tmp_file_left_behind_after_success(tmp_path):
    out_path = str(tmp_path / "cat.h5")
    WriteStage(out_path, _schema()).run(_context())
    assert not (tmp_path / "cat.h5.tmp").exists()


def test_failed_write_does_not_leave_out_path_or_tmp_file(tmp_path, monkeypatch):
    out_path = str(tmp_path / "cat.h5")

    def _boom(self, f, context):
        raise RuntimeError("boom")

    monkeypatch.setattr(WriteStage, "_write_provenance", _boom)
    with pytest.raises(RuntimeError, match="boom"):
        WriteStage(out_path, _schema()).run(_context())

    assert not (tmp_path / "cat.h5").exists()
    assert not (tmp_path / "cat.h5.tmp").exists()


def test_record_stage_reports_out_path_and_field_count(tmp_path):
    out_path = str(tmp_path / "cat.h5")
    context = WriteStage(out_path, _schema()).run(_context())

    write_entry = context.provenance[-1]
    assert write_entry["stage"] == "write"
    assert write_entry["out_path"] == out_path
    # M200c, Position, SatelliteID, Mpeak, IsBacksplash, SatelliteName
    assert write_entry["n_fields_written"] == 6
