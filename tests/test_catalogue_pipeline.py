"""
Tests for analysistools.catalogue.pipeline's PipelineContext/PipelineStage
scaffolding, plus the ExtractStage base class and the still-unimplemented
QualityControlStage/WriteStage/CatalogueBuilder (Phase 6c).
HaloExtractStage/TreeExtractStage/CrossMatchStage are implemented (Phase
6b) and tested separately in test_catalogue_pipeline_extract_crossmatch.py.
"""
import pytest

from analysistools.catalogue.pipeline import (
    CatalogueBuilder, ExtractStage, PipelineContext,
    PipelineStage, QualityControlStage, WriteStage,
)


# ---------------------------------------------------------------------------
# PipelineContext
# ---------------------------------------------------------------------------

def test_context_starts_empty():
    context = PipelineContext()
    assert context.columns == {}
    assert context.meta == {}
    assert context.provenance == []


def test_record_stage_appends_provenance_entry():
    context = PipelineContext()
    context.record_stage("extract")
    assert context.provenance == [{"stage": "extract"}]


def test_record_stage_carries_extra_info():
    context = PipelineContext()
    context.record_stage("extract", n_rows=42, source="halo_catalogue")
    assert context.provenance == [
        {"stage": "extract", "n_rows": 42, "source": "halo_catalogue"}
    ]


def test_record_stage_multiple_calls_accumulate_in_order():
    context = PipelineContext()
    context.record_stage("extract")
    context.record_stage("cross_match")
    assert [p["stage"] for p in context.provenance] == ["extract", "cross_match"]


def test_context_state_is_independent_per_instance():
    # append-only, per-instance state -- mutating one context must not leak
    # into a fresh one (guards against an accidental mutable-default bug).
    a = PipelineContext()
    b = PipelineContext()
    a.columns["Satellites/HaloProperties/Mpeak"] = [1, 2, 3]
    a.meta["schema_version"] = "1.0"
    a.record_stage("extract")
    assert b.columns == {}
    assert b.meta == {}
    assert b.provenance == []


# ---------------------------------------------------------------------------
# PipelineStage
# ---------------------------------------------------------------------------

def test_pipeline_stage_is_abstract():
    with pytest.raises(TypeError):
        PipelineStage()


class _BareStage(PipelineStage):
    """A stage that overrides nothing but the required abstract run()."""

    def run(self, context):
        return context


class _DummyStage(PipelineStage):
    name = "dummy"
    inputs = ("Satellites/Identification/SatelliteID",
              "Satellites/HaloProperties/Mpeak")

    def run(self, context):
        context.record_stage(self.name)
        return context


def test_stage_default_name_inputs_outputs():
    stage = _BareStage()
    assert stage.name == "stage"
    assert stage.inputs == ()
    assert stage.outputs == ()


def test_check_inputs_passes_when_all_present():
    context = PipelineContext()
    context.columns["Satellites/Identification/SatelliteID"] = [0, 1, 2]
    context.columns["Satellites/HaloProperties/Mpeak"] = [1e10, 1e11, 1e12]

    _DummyStage().check_inputs(context)  # must not raise


def test_check_inputs_raises_on_missing_with_stage_name_and_fields():
    context = PipelineContext()
    context.columns["Satellites/Identification/SatelliteID"] = [0, 1, 2]
    # Mpeak deliberately not populated

    with pytest.raises(RuntimeError, match="dummy") as exc_info:
        _DummyStage().check_inputs(context)
    assert "Satellites/HaloProperties/Mpeak" in str(exc_info.value)


def test_check_inputs_lists_all_missing_not_just_first():
    context = PipelineContext()  # neither input populated
    with pytest.raises(RuntimeError) as exc_info:
        _DummyStage().check_inputs(context)
    msg = str(exc_info.value)
    assert "Satellites/Identification/SatelliteID" in msg
    assert "Satellites/HaloProperties/Mpeak" in msg


def test_run_executes_and_records_provenance():
    context = PipelineContext()
    context.columns["Satellites/Identification/SatelliteID"] = [0, 1, 2]
    context.columns["Satellites/HaloProperties/Mpeak"] = [1e10, 1e11, 1e12]

    stage = _DummyStage()
    stage.check_inputs(context)
    result = stage.run(context)

    assert result is context  # passed through, not copied
    assert result.provenance == [{"stage": "dummy"}]


# ---------------------------------------------------------------------------
# Concrete stub stages / CatalogueBuilder -- still Phase 6b/6c, pin that
# they fail loudly rather than silently. Update when implemented.
#
# HaloExtractStage/TreeExtractStage/CrossMatchStage are now implemented --
# see test_catalogue_pipeline_extract_crossmatch.py for their tests. The
# ExtractStage *base class* stays an unimplemented stub (particle tagging/
# selection functions/Rubin detectability have no Epoch-side machinery yet).
# ---------------------------------------------------------------------------

def test_extract_stage_base_not_yet_implemented():
    with pytest.raises(NotImplementedError, match="Phase 6b"):
        ExtractStage().run(PipelineContext())


def test_quality_control_stage_not_yet_implemented():
    with pytest.raises(NotImplementedError, match="Phase 6c"):
        QualityControlStage().run(PipelineContext())


def test_write_stage_not_yet_implemented(tmp_path):
    stage = WriteStage(str(tmp_path / "out.h5"), schema=None)
    with pytest.raises(NotImplementedError, match="Phase 6c"):
        stage.run(PipelineContext())


def test_catalogue_builder_load_config_not_yet_implemented():
    with pytest.raises(NotImplementedError, match="Phase 6a"):
        CatalogueBuilder("configs/dorcha.yaml")._load_config()


def test_catalogue_builder_run_not_yet_implemented():
    with pytest.raises(NotImplementedError, match=r"Phase 6b/6c"):
        CatalogueBuilder("configs/dorcha.yaml").run(simulations=[])


def test_catalogue_builder_run_from_stage_not_yet_implemented():
    builder = CatalogueBuilder("configs/dorcha.yaml")
    with pytest.raises(NotImplementedError, match="Phase 6b"):
        builder.run_from_stage("extract", PipelineContext(), simulations=[])
