"""
Tests for analysistools.catalogue.validation.SchemaValidator (Phase 6b:
pulled forward from Phase 6c specifically to catch the little-h/comoving
"unit reconciliation isn't implemented anywhere in the pipeline yet"
caveat repeated across derived.py/backends.py docstrings).

Uses small hand-built CatalogueSchema/PipelineContext fixtures rather than
the full bundled schema_v1.0.yaml, so each check is exactly hand-checkable
against a minimal, known field set.
"""
import numpy as np
import pytest

from analysistools.catalogue.pipeline import PipelineContext
from analysistools.catalogue.schema import CatalogueSchema, FieldSpec
from analysistools.catalogue.validation import SchemaValidator


def _schema(*specs):
    schema = CatalogueSchema(version="test")
    for spec in specs:
        schema.add(spec)
    return schema


def _issues(report, check=None):
    return [i for i in report.issues if check is None or i.check == check]


# ---------------------------------------------------------------------------
# Presence
# ---------------------------------------------------------------------------

def test_no_issues_when_columns_exactly_match_schema():
    schema = _schema(
        FieldSpec(name="M200c", group="Haloes", dtype="float64",
                  units="Msun"))
    context = PipelineContext()
    context.columns["Haloes/M200c"] = np.array([1e12])

    report = SchemaValidator().check(context, schema)
    assert report.issues == []


def test_missing_declared_field_is_a_warning():
    schema = _schema(
        FieldSpec(name="M200c", group="Haloes", dtype="float64",
                  units="Msun"),
        FieldSpec(name="R200c", group="Haloes", dtype="float64",
                  units="kpc"))
    context = PipelineContext()
    context.columns["Haloes/M200c"] = np.array([1e12])

    report = SchemaValidator().check(context, schema)
    presence = _issues(report, "schema_presence")
    assert len(presence) == 1
    assert presence[0].severity == "warning"
    assert "missing declared field 'Haloes/R200c'" in presence[0].message


def test_undeclared_field_in_known_group_is_an_error():
    schema = _schema(
        FieldSpec(name="M200c", group="Haloes", dtype="float64",
                  units="Msun"))
    context = PipelineContext()
    context.columns["Haloes/M200c"] = np.array([1e12])
    context.columns["Haloes/Typo"] = np.array([1.0])

    report = SchemaValidator().check(context, schema)
    presence = _issues(report, "schema_presence")
    assert len(presence) == 1
    assert presence[0].severity == "error"
    assert "undeclared field 'Haloes/Typo'" in presence[0].message


def test_undeclared_group_is_an_error():
    schema = _schema(
        FieldSpec(name="M200c", group="Haloes", dtype="float64",
                  units="Msun"))
    context = PipelineContext()
    context.columns["Haloes/M200c"] = np.array([1e12])
    context.columns["NotAGroup/Field"] = np.array([1.0])

    report = SchemaValidator().check(context, schema)
    groups = _issues(report, "schema_group")
    assert len(groups) == 1
    assert groups[0].severity == "error"
    assert "'NotAGroup'" in groups[0].message


def test_internal_and_mergertrees_columns_are_ignored():
    schema = _schema(
        FieldSpec(name="M200c", group="Haloes", dtype="float64",
                  units="Msun"))
    context = PipelineContext()
    context.columns["Haloes/M200c"] = np.array([1e12])
    context.columns["Satellites/_internal/halo_row"] = np.array([0, 1])
    context.columns["MergerTrees/main_branch"] = [None, None]

    report = SchemaValidator().check(context, schema)
    assert report.issues == []


# ---------------------------------------------------------------------------
# dtype
# ---------------------------------------------------------------------------

def test_dtype_mismatch_is_a_warning():
    schema = _schema(
        FieldSpec(name="SnapshotAtMpeak", group="Satellites/HaloProperties",
                  dtype="int32", units="dimensionless"))
    context = PipelineContext()
    context.columns["Satellites/HaloProperties/SnapshotAtMpeak"] = \
        np.array([1.0, 2.0])  # float, not int

    report = SchemaValidator().check(context, schema)
    dtype_issues = _issues(report, "schema_dtype")
    assert len(dtype_issues) == 1
    assert dtype_issues[0].severity == "warning"
    assert "does not match schema-declared 'int32'" in dtype_issues[0].message


def test_dtype_match_is_silent():
    schema = _schema(
        FieldSpec(name="IsBacksplash", group="Satellites/HaloProperties",
                  dtype="bool", units="dimensionless"))
    context = PipelineContext()
    context.columns["Satellites/HaloProperties/IsBacksplash"] = \
        np.array([True, False])

    report = SchemaValidator().check(context, schema)
    assert _issues(report, "schema_dtype") == []


def test_str_dtype_is_never_checked():
    schema = _schema(
        FieldSpec(name="SatelliteName", group="Satellites/Identification",
                  dtype="str", units="dimensionless"))
    context = PipelineContext()
    context.columns["Satellites/Identification/SatelliteName"] = \
        np.array(["a", "b"])

    report = SchemaValidator().check(context, schema)
    assert _issues(report, "schema_dtype") == []


def test_list_valued_columns_skip_dtype_check_without_crashing():
    # e.g. MergerTrees/main_branch-style list columns, if ever declared
    schema = _schema(
        FieldSpec(name="Thing", group="Weird", dtype="float64",
                  units="dimensionless"))
    context = PipelineContext()
    context.columns["Weird/Thing"] = [None, None]

    report = SchemaValidator().check(context, schema)
    assert _issues(report, "schema_dtype") == []


# ---------------------------------------------------------------------------
# little-h / comoving vs. declared units
# ---------------------------------------------------------------------------

def test_little_h_mismatch_is_flagged():
    schema = _schema(
        FieldSpec(name="M200c", group="Haloes", dtype="float64",
                  units="Msun"))
    context = PipelineContext()
    context.columns["Haloes/M200c"] = np.array([1e12])
    context.meta["comoving_little_h"] = {
        "Haloes": {"comoving": False, "little_h": True}}

    report = SchemaValidator().check(context, schema)
    units_issues = _issues(report, "schema_units")
    assert len(units_issues) == 1
    assert units_issues[0].severity == "warning"
    assert "little-h-scaled" in units_issues[0].message


def test_little_h_false_does_not_flag_plain_units():
    schema = _schema(
        FieldSpec(name="M200c", group="Haloes", dtype="float64",
                  units="Msun"))
    context = PipelineContext()
    context.columns["Haloes/M200c"] = np.array([1e12])
    context.meta["comoving_little_h"] = {
        "Haloes": {"comoving": False, "little_h": False}}

    report = SchemaValidator().check(context, schema)
    assert _issues(report, "schema_units") == []


def test_little_h_true_with_h_suffix_already_declared_is_silent():
    schema = _schema(
        FieldSpec(name="M200c", group="Haloes", dtype="float64",
                  units="Msun/h"))
    context = PipelineContext()
    context.columns["Haloes/M200c"] = np.array([1e12])
    context.meta["comoving_little_h"] = {
        "Haloes": {"comoving": False, "little_h": True}}

    report = SchemaValidator().check(context, schema)
    assert _issues(report, "schema_units") == []


def test_comoving_mismatch_is_flagged_for_length_field():
    schema = _schema(
        FieldSpec(name="Position", group="Haloes", dtype="float64",
                  units="Mpc", shape=(3,)))
    context = PipelineContext()
    context.columns["Haloes/Position"] = np.array([[1.0, 2.0, 3.0]])
    context.meta["comoving_little_h"] = {
        "Haloes": {"comoving": True, "little_h": False}}

    report = SchemaValidator().check(context, schema)
    units_issues = _issues(report, "schema_units")
    assert len(units_issues) == 1
    assert "comoving" in units_issues[0].message


def test_comoving_does_not_flag_non_length_fields():
    schema = _schema(
        FieldSpec(name="M200c", group="Haloes", dtype="float64",
                  units="Msun"))
    context = PipelineContext()
    context.columns["Haloes/M200c"] = np.array([1e12])
    context.meta["comoving_little_h"] = {
        "Haloes": {"comoving": True, "little_h": False}}

    report = SchemaValidator().check(context, schema)
    assert _issues(report, "schema_units") == []


def test_no_recorded_native_state_skips_unit_check_entirely():
    schema = _schema(
        FieldSpec(name="M200c", group="Haloes", dtype="float64",
                  units="Msun"))
    context = PipelineContext()
    context.columns["Haloes/M200c"] = np.array([1e12])
    # no context.meta["comoving_little_h"] at all

    report = SchemaValidator().check(context, schema)
    assert _issues(report, "schema_units") == []


def test_sfr_field_h_suffix_is_checked_like_mass():
    schema = _schema(
        FieldSpec(name="StarFormationRate", group="Satellites/GalaxyProperties",
                  dtype="float64", units="Msun/yr"))
    context = PipelineContext()
    context.columns["Satellites/GalaxyProperties/StarFormationRate"] = \
        np.array([1.5])
    context.meta["comoving_little_h"] = {
        "Satellites/GalaxyProperties": {"comoving": True, "little_h": True}}

    report = SchemaValidator().check(context, schema)
    units_issues = _issues(report, "schema_units")
    assert len(units_issues) == 1
    assert "little-h-scaled" in units_issues[0].message


# ---------------------------------------------------------------------------
# Integration with the real bundled schema
# ---------------------------------------------------------------------------

def test_runs_against_default_schema_without_crashing():
    from analysistools.catalogue.schema import default_schema
    schema = default_schema("1.0")
    context = PipelineContext()
    context.columns["Haloes/M200c"] = np.array([1e12])
    context.meta["comoving_little_h"] = {
        "Haloes": {"comoving": False, "little_h": True}}

    report = SchemaValidator().check(context, schema)
    # every other Haloes/* field is legitimately "missing" here -> warnings
    assert all(i.severity in ("warning", "error") for i in report.issues)
    assert any(i.check == "schema_units" for i in report.issues)
