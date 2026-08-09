"""
Tests for analysistools.catalogue.validation.IntegrityValidator (Phase 6c).

Uses hand-built PipelineContext fixtures -- each check is exactly
hand-checkable against a small, known column set. `schema` is accepted
by the interface but unused by this validator, so tests pass `schema=None`.
"""
import numpy as np

from analysistools.catalogue.pipeline import PipelineContext
from analysistools.catalogue.validation import IntegrityValidator


def _issues(report, check=None):
    return [i for i in report.issues if check is None or i.check == check]


# ---------------------------------------------------------------------------
# Row-count consistency
# ---------------------------------------------------------------------------

def test_no_issues_when_all_satellite_columns_agree_in_length():
    context = PipelineContext()
    context.columns["Satellites/Identification/SatelliteID"] = \
        np.arange(3, dtype=np.int64)
    context.columns["Satellites/HaloProperties/Mpeak"] = np.array(
        [1.0, 2.0, 3.0])

    report = IntegrityValidator().check(context, schema=None)
    assert report.issues == []


def test_mismatched_row_count_is_an_error():
    context = PipelineContext()
    context.columns["Satellites/Identification/SatelliteID"] = \
        np.arange(3, dtype=np.int64)
    context.columns["Satellites/HaloProperties/Mpeak"] = np.array([1.0, 2.0])

    report = IntegrityValidator().check(context, schema=None)
    row_count = _issues(report, "row_count")
    assert len(row_count) == 1
    assert row_count[0].severity == "error"
    assert "length 2" in row_count[0].message


def test_row_count_falls_back_to_mode_without_satellite_id():
    context = PipelineContext()
    context.columns["Satellites/HaloProperties/Mpeak"] = np.array(
        [1.0, 2.0, 3.0])
    context.columns["Satellites/HaloProperties/Vpeak"] = np.array(
        [1.0, 2.0, 3.0])
    context.columns["Satellites/GalaxyProperties/StellarMass"] = \
        np.array([1.0])  # the odd one out

    report = IntegrityValidator().check(context, schema=None)
    row_count = _issues(report, "row_count")
    assert len(row_count) == 1
    assert row_count[0].field == "Satellites/GalaxyProperties/StellarMass"


def test_haloes_and_internal_columns_are_not_row_count_checked():
    context = PipelineContext()
    context.columns["Satellites/Identification/SatelliteID"] = \
        np.arange(3, dtype=np.int64)
    context.columns["Haloes/HostHaloID"] = np.array([99])  # length 1, fine

    report = IntegrityValidator().check(context, schema=None)
    assert _issues(report, "row_count") == []


def test_list_valued_mergertree_column_length_is_checked():
    context = PipelineContext()
    context.columns["Satellites/Identification/SatelliteID"] = \
        np.arange(3, dtype=np.int64)
    context.columns["MergerTrees/main_branch"] = [None, None]  # length 2

    report = IntegrityValidator().check(context, schema=None)
    row_count = _issues(report, "row_count")
    assert len(row_count) == 1
    assert row_count[0].field == "MergerTrees/main_branch"


# ---------------------------------------------------------------------------
# SatelliteID uniqueness
# ---------------------------------------------------------------------------

def test_duplicate_satellite_id_is_an_error():
    context = PipelineContext()
    context.columns["Satellites/Identification/SatelliteID"] = \
        np.array([0, 1, 1, 2], dtype=np.int64)

    report = IntegrityValidator().check(context, schema=None)
    dup = _issues(report, "satellite_id_uniqueness")
    assert len(dup) == 1
    assert dup[0].severity == "error"
    assert "[1]" in dup[0].message


def test_unique_satellite_id_is_silent():
    context = PipelineContext()
    context.columns["Satellites/Identification/SatelliteID"] = \
        np.array([0, 1, 2], dtype=np.int64)

    report = IntegrityValidator().check(context, schema=None)
    assert _issues(report, "satellite_id_uniqueness") == []


# ---------------------------------------------------------------------------
# Orphaned foreign keys
# ---------------------------------------------------------------------------

def test_orphaned_host_halo_id_is_an_error():
    context = PipelineContext()
    context.columns["Satellites/Identification/HostHaloID"] = \
        np.array([1, 1, 99])
    context.columns["Haloes/HostHaloID"] = np.array([1])

    report = IntegrityValidator().check(context, schema=None)
    orphaned = _issues(report, "orphaned_foreign_key")
    assert len(orphaned) == 1
    assert "[99]" in orphaned[0].message


def test_no_orphaned_foreign_key_when_all_resolve():
    context = PipelineContext()
    context.columns["Satellites/Identification/HostHaloID"] = \
        np.array([1, 1, 1])
    context.columns["Haloes/HostHaloID"] = np.array([1])

    report = IntegrityValidator().check(context, schema=None)
    assert _issues(report, "orphaned_foreign_key") == []


# ---------------------------------------------------------------------------
# inf detection
# ---------------------------------------------------------------------------

def test_inf_value_is_an_error():
    context = PipelineContext()
    context.columns["Satellites/HaloProperties/Vpeak"] = \
        np.array([1.0, np.inf, 3.0])

    report = IntegrityValidator().check(context, schema=None)
    inf_issues = _issues(report, "unexpected_inf")
    assert len(inf_issues) == 1
    assert "1 inf value" in inf_issues[0].message


def test_nan_alone_is_not_flagged_as_inf():
    context = PipelineContext()
    context.columns["Satellites/HaloProperties/Vpeak"] = \
        np.array([1.0, np.nan, 3.0])

    report = IntegrityValidator().check(context, schema=None)
    assert _issues(report, "unexpected_inf") == []


def test_int_and_bool_columns_are_not_inf_checked():
    context = PipelineContext()
    context.columns["Satellites/Identification/SatelliteID"] = \
        np.array([0, 1, 2], dtype=np.int64)
    context.columns["Satellites/HaloProperties/IsBacksplash"] = \
        np.array([True, False])

    report = IntegrityValidator().check(context, schema=None)
    assert _issues(report, "unexpected_inf") == []


# ---------------------------------------------------------------------------
# CSR offsets
# ---------------------------------------------------------------------------

def test_non_monotonic_csr_offsets_is_an_error():
    context = PipelineContext()
    context.columns["Satellites/ParticleTags/particle_id_offsets"] = \
        np.array([0, 10, 5, 20])

    report = IntegrityValidator().check(context, schema=None)
    csr = _issues(report, "csr_offsets_not_monotonic")
    assert len(csr) == 1


def test_monotonic_csr_offsets_is_silent():
    context = PipelineContext()
    context.columns["Satellites/ParticleTags/particle_id_offsets"] = \
        np.array([0, 10, 10, 20])  # non-decreasing (flat run allowed)

    report = IntegrityValidator().check(context, schema=None)
    assert _issues(report, "csr_offsets_not_monotonic") == []


def test_absent_csr_offsets_is_not_checked():
    context = PipelineContext()
    context.columns["Satellites/Identification/SatelliteID"] = \
        np.array([0, 1])

    report = IntegrityValidator().check(context, schema=None)
    assert _issues(report, "csr_offsets_not_monotonic") == []


def test_empty_context_produces_no_issues():
    report = IntegrityValidator().check(PipelineContext(), schema=None)
    assert report.issues == []
