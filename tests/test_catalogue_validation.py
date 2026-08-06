"""
Tests for analysistools.catalogue.validation (Phase 6a).

ValidationIssue/ValidationReport are fully implemented and covered here.
Validator is abstract; SchemaValidator/IntegrityValidator/PhysicalValidator
(and ValidationReport.to_json()) are still Phase 6c stubs -- pinned the
same way test_catalogue_pipeline.py pins pipeline.py's stubs.
"""
import pytest

from analysistools.catalogue.validation import (
    DEFAULT_VALIDATORS, IntegrityValidator, PhysicalValidator,
    SchemaValidator, Validator, ValidationIssue, ValidationReport,
)


# ---------------------------------------------------------------------------
# ValidationIssue
# ---------------------------------------------------------------------------

def test_validation_issue_default_field_is_empty():
    issue = ValidationIssue(severity="error", check="schema", message="boom")
    assert issue.field == ""


def test_validation_issue_stores_all_fields():
    issue = ValidationIssue(severity="warning", check="physical",
                            message="Mpeak < M200c_z0", field="Mpeak")
    assert issue.severity == "warning"
    assert issue.check == "physical"
    assert issue.message == "Mpeak < M200c_z0"
    assert issue.field == "Mpeak"


# ---------------------------------------------------------------------------
# ValidationReport
# ---------------------------------------------------------------------------

def test_report_starts_empty_and_passes():
    report = ValidationReport()
    assert report.issues == []
    assert report.errors == []
    assert report.warnings == []
    assert report.passed is True


def test_add_appends_issue_and_defaults_field():
    report = ValidationReport()
    report.add("error", "schema", "missing field")
    assert len(report.issues) == 1
    issue = report.issues[0]
    assert issue.severity == "error"
    assert issue.check == "schema"
    assert issue.message == "missing field"
    assert issue.field == ""


def test_errors_and_warnings_filter_by_severity():
    report = ValidationReport()
    report.add("error", "schema", "missing Mpeak field", field="Mpeak")
    report.add("warning", "physical", "Mpeak slightly below M200c_z0",
              field="Mpeak")
    report.add("error", "integrity", "duplicate SatelliteID")

    assert len(report.errors) == 2
    assert len(report.warnings) == 1
    assert {i.check for i in report.errors} == {"schema", "integrity"}
    assert report.warnings[0].check == "physical"


def test_passed_is_false_when_any_error_present():
    report = ValidationReport()
    report.add("warning", "physical", "borderline value")
    assert report.passed is True  # warnings alone don't fail
    report.add("error", "schema", "missing field")
    assert report.passed is False


def test_reports_are_independent_instances():
    # dataclass default_factory=list -- guards against an accidental
    # mutable-default bug across instances.
    a = ValidationReport()
    b = ValidationReport()
    a.add("error", "schema", "boom")
    assert b.issues == []


def test_summary_prints_error_and_warning_counts(capsys):
    report = ValidationReport()
    report.add("error", "schema", "missing Mpeak field", field="Mpeak")
    report.add("warning", "physical", "borderline value")
    report.summary()

    out = capsys.readouterr().out
    assert "1 error(s)" in out
    assert "1 warning(s)" in out
    assert "ERROR" in out
    assert "WARNING" in out
    assert "Mpeak" in out  # field is appended in parens when present


def test_summary_omits_field_parens_when_absent(capsys):
    report = ValidationReport()
    report.add("error", "schema", "top-level failure")
    report.summary()
    out = capsys.readouterr().out
    assert "top-level failure" in out
    assert "()" not in out


def test_to_json_not_yet_implemented():
    with pytest.raises(NotImplementedError, match="Phase 6c"):
        ValidationReport().to_json()


# ---------------------------------------------------------------------------
# Validator / concrete stub validators
# ---------------------------------------------------------------------------

def test_validator_is_abstract():
    with pytest.raises(TypeError):
        Validator()


def test_default_validators_are_the_three_concrete_classes():
    assert [type(v) for v in DEFAULT_VALIDATORS] == [
        SchemaValidator, IntegrityValidator, PhysicalValidator,
    ]


@pytest.mark.parametrize("validator_cls,expected_name", [
    (SchemaValidator, "schema"),
    (IntegrityValidator, "integrity"),
    (PhysicalValidator, "physical"),
])
def test_validator_name(validator_cls, expected_name):
    assert validator_cls().name == expected_name


@pytest.mark.parametrize("validator_cls", [
    SchemaValidator, IntegrityValidator, PhysicalValidator,
])
def test_validator_check_not_yet_implemented(validator_cls):
    with pytest.raises(NotImplementedError, match="Phase 6c"):
        validator_cls().check(context=None, schema=None)
