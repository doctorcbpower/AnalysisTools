"""
Tests for analysistools.catalogue.webreport (Phase 6a).

FigureSpec, WebReport.add_figure(), and WebReport.from_catalogue() are
fully implemented and covered here. WebReport.build() is still a Phase 6c
stub -- pinned the same way other catalogue submodule tests pin stubs.
"""
from types import SimpleNamespace

import pytest

from analysistools.catalogue.webreport import FigureSpec, WebReport


# ---------------------------------------------------------------------------
# FigureSpec
# ---------------------------------------------------------------------------

def test_figure_spec_defaults():
    spec = FigureSpec(name="mass_function", title="Mass Function", figure="fig")
    assert spec.description == ""
    assert spec.group == "default"


def test_figure_spec_stores_all_fields():
    spec = FigureSpec(name="orbits", title="Orbital Tracks", figure="fig",
                      description="Pericentre/apocentre distribution",
                      group="orbits")
    assert spec.name == "orbits"
    assert spec.title == "Orbital Tracks"
    assert spec.figure == "fig"
    assert spec.description == "Pericentre/apocentre distribution"
    assert spec.group == "orbits"


# ---------------------------------------------------------------------------
# WebReport
# ---------------------------------------------------------------------------

def test_webreport_defaults():
    report = WebReport(catalogue_version="v1.0.0")
    assert report.title == "Catalogue Report"
    assert report.figures == []


def test_add_figure_appends_spec_and_defaults_title_to_name():
    report = WebReport(catalogue_version="v1.0.0")
    report.add_figure("mass_function", "fig1")

    assert len(report.figures) == 1
    spec = report.figures[0]
    assert spec.name == "mass_function"
    assert spec.title == "mass_function"  # title defaults to name
    assert spec.figure == "fig1"
    assert spec.description == ""
    assert spec.group == "default"


def test_add_figure_accepts_explicit_title_description_group():
    report = WebReport(catalogue_version="v1.0.0")
    report.add_figure("mass_function", "fig1", title="Mass Function",
                      description="dn/dlogM", group="halo_properties")

    spec = report.figures[0]
    assert spec.title == "Mass Function"
    assert spec.description == "dn/dlogM"
    assert spec.group == "halo_properties"


def test_add_figure_is_chainable_and_accumulates():
    report = WebReport(catalogue_version="v1.0.0")
    result = report.add_figure("a", "fig_a").add_figure("b", "fig_b")

    assert result is report
    assert [s.name for s in report.figures] == ["a", "b"]


def test_reports_figures_are_independent_instances():
    # dataclass default_factory=list -- guards against an accidental
    # mutable-default bug across instances.
    a = WebReport(catalogue_version="v1")
    b = WebReport(catalogue_version="v2")
    a.add_figure("only_in_a", "fig")
    assert b.figures == []


def test_from_catalogue_reads_version_from_meta():
    catalogue = SimpleNamespace(meta={"catalogue_version": "v2.3.4"})
    report = WebReport.from_catalogue(catalogue, title="My Report")

    assert report.catalogue_version == "v2.3.4"
    assert report.title == "My Report"
    assert report.figures == []


def test_from_catalogue_defaults_version_when_absent():
    catalogue = SimpleNamespace(meta={})
    report = WebReport.from_catalogue(catalogue)

    assert report.catalogue_version == "unknown"
    assert report.title == "Catalogue Report"


def test_build_not_yet_implemented(tmp_path):
    report = WebReport(catalogue_version="v1.0.0")
    with pytest.raises(NotImplementedError, match="Phase 6c"):
        report.build(str(tmp_path))
