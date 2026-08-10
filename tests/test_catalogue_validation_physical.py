"""
Tests for analysistools.catalogue.validation.PhysicalValidator (Phase 6c).

Uses hand-built PipelineContext fixtures. Two checks deliberately deviate
from the validator's original one-line description -- see the class
docstring -- and are tested against their actual (substituted)
definitions here, not the literal original wording.
"""
import numpy as np

from analysistools.catalogue.pipeline import PipelineContext
from analysistools.catalogue.validation import PhysicalValidator


def _issues(report, check=None):
    return [i for i in report.issues if check is None or i.check == check]


# ---------------------------------------------------------------------------
# Mpeak >= M200c_z0
# ---------------------------------------------------------------------------

def test_mpeak_below_m200c_z0_is_an_error():
    context = PipelineContext()
    context.columns["Satellites/HaloProperties/Mpeak"] = np.array([1e9])
    context.columns["Satellites/HaloProperties/M200c_z0"] = np.array([1e10])

    report = PhysicalValidator().check(context, schema=None)
    issues = _issues(report, "mpeak_below_m200c_z0")
    assert len(issues) == 1
    assert issues[0].severity == "error"


def test_mpeak_at_or_above_m200c_z0_is_silent():
    context = PipelineContext()
    context.columns["Satellites/HaloProperties/Mpeak"] = np.array([1e10, 1e10])
    context.columns["Satellites/HaloProperties/M200c_z0"] = np.array([1e10, 5e9])

    report = PhysicalValidator().check(context, schema=None)
    assert _issues(report, "mpeak_below_m200c_z0") == []


def test_nan_entries_skipped_for_mpeak_check():
    context = PipelineContext()
    context.columns["Satellites/HaloProperties/Mpeak"] = np.array([np.nan])
    context.columns["Satellites/HaloProperties/M200c_z0"] = np.array([1e10])

    report = PhysicalValidator().check(context, schema=None)
    assert _issues(report, "mpeak_below_m200c_z0") == []


def test_mpeak_check_skipped_when_either_column_absent():
    context = PipelineContext()
    context.columns["Satellites/HaloProperties/Mpeak"] = np.array([1e9])

    report = PhysicalValidator().check(context, schema=None)
    assert _issues(report, "mpeak_below_m200c_z0") == []


# ---------------------------------------------------------------------------
# OrbitalApocentre >= OrbitalPericentre
# ---------------------------------------------------------------------------

def test_apocentre_below_pericentre_is_an_error():
    context = PipelineContext()
    context.columns["Satellites/HaloProperties/OrbitalApocentre"] = \
        np.array([5.0])
    context.columns["Satellites/HaloProperties/OrbitalPericentre"] = \
        np.array([10.0])

    report = PhysicalValidator().check(context, schema=None)
    assert len(_issues(report, "apocentre_below_pericentre")) == 1


def test_apocentre_at_or_above_pericentre_is_silent():
    context = PipelineContext()
    context.columns["Satellites/HaloProperties/OrbitalApocentre"] = \
        np.array([10.0, 5.0])
    context.columns["Satellites/HaloProperties/OrbitalPericentre"] = \
        np.array([5.0, 5.0])

    report = PhysicalValidator().check(context, schema=None)
    assert _issues(report, "apocentre_below_pericentre") == []


# ---------------------------------------------------------------------------
# HeliocentricDistance > 0
# ---------------------------------------------------------------------------

def test_non_positive_heliocentric_distance_is_an_error():
    context = PipelineContext()
    context.columns["Satellites/Observability/HeliocentricDistance"] = \
        np.array([10.0, 0.0, -1.0])

    report = PhysicalValidator().check(context, schema=None)
    issues = _issues(report, "non_positive_heliocentric_distance")
    assert len(issues) == 1
    assert "2 satellite" in issues[0].message


def test_positive_heliocentric_distance_is_silent():
    context = PipelineContext()
    context.columns["Satellites/Observability/HeliocentricDistance"] = \
        np.array([10.0, 20.0])

    report = PhysicalValidator().check(context, schema=None)
    assert _issues(report, "non_positive_heliocentric_distance") == []


# ---------------------------------------------------------------------------
# CompletenessWeight >= 1
# ---------------------------------------------------------------------------

def test_completeness_weight_below_one_is_an_error():
    context = PipelineContext()
    context.columns["Satellites/Observability/CompletenessWeight"] = \
        np.array([1.0, 0.5])

    report = PhysicalValidator().check(context, schema=None)
    assert len(_issues(report, "completeness_weight_below_one")) == 1


def test_completeness_weight_absent_is_not_checked():
    report = PhysicalValidator().check(PipelineContext(), schema=None)
    assert _issues(report, "completeness_weight_below_one") == []


# ---------------------------------------------------------------------------
# Backsplash implies at least one infall (substituted check)
# ---------------------------------------------------------------------------

def test_backsplash_without_infall_is_an_error():
    context = PipelineContext()
    context.columns["Satellites/HaloProperties/IsBacksplash"] = \
        np.array([True, False])
    context.columns["Satellites/HaloProperties/NumberOfInfalls"] = \
        np.array([0, 0])

    report = PhysicalValidator().check(context, schema=None)
    issues = _issues(report, "backsplash_without_infall")
    assert len(issues) == 1


def test_backsplash_with_infall_is_silent():
    context = PipelineContext()
    context.columns["Satellites/HaloProperties/IsBacksplash"] = \
        np.array([True, False])
    context.columns["Satellites/HaloProperties/NumberOfInfalls"] = \
        np.array([1, 0])

    report = PhysicalValidator().check(context, schema=None)
    assert _issues(report, "backsplash_without_infall") == []


# ---------------------------------------------------------------------------
# Total SFH-formed mass >= StellarMass (substituted check)
# ---------------------------------------------------------------------------

def _sfh_context(sfh, stellar_mass, edges):
    context = PipelineContext()
    context.columns["Satellites/GalaxyProperties/SFH"] = np.asarray(sfh)
    context.columns["Satellites/GalaxyProperties/StellarMass"] = \
        np.asarray(stellar_mass)
    context.meta["time_bin_edges_sfh"] = np.asarray(edges)
    return context


def test_stellar_mass_exceeding_formed_mass_is_a_warning_not_an_error():
    # SFR=1 Msun/yr for 1 Gyr -> formed mass = 1e9 Msun; claim 2e9 today.
    # Warning, not error: StellarMass can legitimately exceed a galaxy's
    # own in-situ SFH integral due to ex-situ/merger-accreted mass a SAM
    # like SHARK doesn't record in the per-galaxy SFH -- see this
    # validator's docstring for the real case that downgraded this from
    # a hard error.
    context = _sfh_context(sfh=[[1.0]], stellar_mass=[2e9], edges=[0.0, 1.0])

    report = PhysicalValidator().check(context, schema=None)
    issues = _issues(report, "stellar_mass_exceeds_formed_mass")
    assert len(issues) == 1
    assert issues[0].severity == "warning"
    assert report.passed is True   # warnings alone don't fail a build


def test_stellar_mass_exceeding_formed_mass_message_reports_worst_case_ratio():
    # two offenders: row 0 ratio=2x, row 1 ratio=4x -- message should
    # single out the worse of the two (row 1) with its actual values.
    context = _sfh_context(sfh=[[1.0], [1.0]], stellar_mass=[2e9, 4e9],
                           edges=[0.0, 1.0])

    report = PhysicalValidator().check(context, schema=None)
    issues = _issues(report, "stellar_mass_exceeds_formed_mass")
    assert len(issues) == 1
    msg = issues[0].message
    assert "row 1" in msg
    assert "4.000e+09" in msg          # StellarMass
    assert "1.000e+09" in msg          # formed_mass
    assert "4.000e+00" in msg          # ratio


def test_stellar_mass_below_formed_mass_is_silent_mass_loss_allowed():
    context = _sfh_context(sfh=[[1.0]], stellar_mass=[5e8], edges=[0.0, 1.0])

    report = PhysicalValidator().check(context, schema=None)
    assert _issues(report, "stellar_mass_exceeds_formed_mass") == []


def test_stellar_mass_equal_to_formed_mass_is_silent():
    context = _sfh_context(sfh=[[1.0]], stellar_mass=[1e9], edges=[0.0, 1.0])

    report = PhysicalValidator().check(context, schema=None)
    assert _issues(report, "stellar_mass_exceeds_formed_mass") == []


def test_all_nan_sfh_row_is_skipped_not_flagged():
    context = _sfh_context(sfh=[[np.nan]], stellar_mass=[1e9], edges=[0.0, 1.0])

    report = PhysicalValidator().check(context, schema=None)
    assert _issues(report, "stellar_mass_exceeds_formed_mass") == []


def test_sfh_check_skipped_without_time_bin_edges_meta():
    context = PipelineContext()
    context.columns["Satellites/GalaxyProperties/SFH"] = np.array([[1.0]])
    context.columns["Satellites/GalaxyProperties/StellarMass"] = \
        np.array([2e9])
    # no context.meta["time_bin_edges_sfh"]

    report = PhysicalValidator().check(context, schema=None)
    assert _issues(report, "stellar_mass_exceeds_formed_mass") == []


def test_empty_context_produces_no_issues():
    report = PhysicalValidator().check(PipelineContext(), schema=None)
    assert report.issues == []
