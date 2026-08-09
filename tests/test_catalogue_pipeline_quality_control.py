"""
Tests for analysistools.catalogue.pipeline.QualityControlStage (Phase 6c).

Uses small hand-built schema/context fixtures plus fake Validators so the
halt-vs-warn behaviour and validation_report attachment are exactly
hand-checkable, independent of the concrete Validator subclasses' own
logic (covered separately in test_catalogue_validation*.py).
"""
import pytest

from analysistools.catalogue.pipeline import PipelineContext, QualityControlStage
from analysistools.catalogue.schema import CatalogueSchema, FieldSpec
from analysistools.catalogue.validation import (
    IntegrityValidator, PhysicalValidator, SchemaValidator, ValidationReport,
)


class _FakeValidator:
    def __init__(self, name, issues):
        self.name = name
        self._issues = issues

    def check(self, context, schema):
        report = ValidationReport()
        for severity, check, message in self._issues:
            report.add(severity, check, message)
        return report


def _schema():
    return CatalogueSchema(version="test")


def test_no_issues_passes_and_attaches_empty_report():
    context = PipelineContext()
    stage = QualityControlStage(_schema(), validators=[
        _FakeValidator("a", []), _FakeValidator("b", []),
    ])
    result = stage.run(context)

    report = result.meta["validation_report"]
    assert report.issues == []
    assert report.passed is True


def test_warnings_only_do_not_raise():
    context = PipelineContext()
    stage = QualityControlStage(_schema(), validators=[
        _FakeValidator("a", [("warning", "physical", "borderline")]),
    ])
    result = stage.run(context)  # must not raise

    report = result.meta["validation_report"]
    assert len(report.warnings) == 1
    assert report.passed is True


def test_errors_raise_and_still_attach_report():
    context = PipelineContext()
    stage = QualityControlStage(_schema(), validators=[
        _FakeValidator("a", [("error", "integrity", "duplicate SatelliteID")]),
    ])
    with pytest.raises(RuntimeError, match="duplicate SatelliteID"):
        stage.run(context)

    # attached before the raise, so a caller catching the error can inspect it
    report = context.meta["validation_report"]
    assert len(report.errors) == 1


def test_issues_combined_across_all_validators():
    context = PipelineContext()
    stage = QualityControlStage(_schema(), validators=[
        _FakeValidator("a", [("warning", "schema", "missing field")]),
        _FakeValidator("b", [("warning", "physical", "borderline")]),
    ])
    result = stage.run(context)

    report = result.meta["validation_report"]
    assert len(report.issues) == 2
    assert {i.check for i in report.issues} == {"schema", "physical"}


def test_error_message_previews_up_to_five_errors():
    context = PipelineContext()
    issues = [("error", "integrity", f"problem {i}") for i in range(7)]
    stage = QualityControlStage(_schema(), validators=[
        _FakeValidator("a", issues),
    ])
    with pytest.raises(RuntimeError) as exc_info:
        stage.run(context)

    msg = str(exc_info.value)
    assert "7 validation error(s)" in msg
    assert "problem 0" in msg
    assert "problem 4" in msg
    assert "problem 5" not in msg
    assert "+2 more" in msg


def test_record_stage_reports_error_and_warning_counts():
    context = PipelineContext()
    stage = QualityControlStage(_schema(), validators=[
        _FakeValidator("a", [("warning", "physical", "borderline")]),
    ])
    result = stage.run(context)

    assert result.provenance == [
        {"stage": "quality_control", "n_errors": 0, "n_warnings": 1}]


def test_default_validators_are_the_three_concrete_classes():
    stage = QualityControlStage(_schema())
    assert [type(v) for v in stage.validators] == [
        SchemaValidator, IntegrityValidator, PhysicalValidator]


# ---------------------------------------------------------------------------
# Integration: real SchemaValidator via the default validator set
# ---------------------------------------------------------------------------

def test_default_validators_flag_missing_schema_fields_as_warnings():
    schema = CatalogueSchema(version="test")
    schema.add(FieldSpec(name="M200c", group="Haloes", dtype="float64",
                         units="Msun"))
    context = PipelineContext()  # Haloes/M200c never populated

    stage = QualityControlStage(schema)
    result = stage.run(context)  # missing declared field -> warning, not error

    report = result.meta["validation_report"]
    assert report.passed is True
    assert any(i.check == "schema_presence" for i in report.warnings)
