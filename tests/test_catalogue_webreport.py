"""
Tests for analysistools.catalogue.webreport (Phase 6a/6c).

FigureSpec, WebReport.add_figure(), WebReport.from_catalogue(), and
WebReport.build() are all implemented and covered here.
"""
from types import SimpleNamespace

import plotly.graph_objects as go
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


def _fig(y):
    return go.Figure(data=[go.Scatter(y=y)])


def test_build_writes_index_and_group_pages(tmp_path):
    report = WebReport(catalogue_version="v1.0.0", title="My Report")
    report.add_figure("mf", _fig([1, 2, 3]), title="Mass Function",
                      description="dn/dlogM", group="halo_properties")
    report.add_figure("orbits", _fig([4, 5, 6]), title="Orbits",
                      group="orbits")

    out = report.build(str(tmp_path))

    assert out == tmp_path
    assert (tmp_path / "index.html").exists()
    assert (tmp_path / "halo_properties.html").exists()
    assert (tmp_path / "orbits.html").exists()


def test_build_index_links_every_group_page(tmp_path):
    report = WebReport(catalogue_version="v1.0.0")
    report.add_figure("a", _fig([1]), group="group_a")
    report.add_figure("b", _fig([2]), group="group_b")

    report.build(str(tmp_path))
    index = (tmp_path / "index.html").read_text()
    assert 'href="group_a.html"' in index
    assert 'href="group_b.html"' in index


def test_build_group_page_embeds_figure_html(tmp_path):
    report = WebReport(catalogue_version="v1.0.0")
    report.add_figure("mf", _fig([1, 2, 3]), title="Mass Function",
                      description="dn/dlogM", group="halo_properties")
    report.build(str(tmp_path))

    page = (tmp_path / "halo_properties.html").read_text()
    assert "Mass Function" in page
    assert "dn/dlogM" in page
    assert "plotly" in page.lower()


def test_build_footer_stamps_catalogue_version(tmp_path):
    report = WebReport(catalogue_version="v3.2.1")
    report.add_figure("a", _fig([1]), group="g")
    report.build(str(tmp_path))

    for page_name in ("index.html", "g.html"):
        assert "v3.2.1" in (tmp_path / page_name).read_text()


def test_build_creates_out_dir_if_missing(tmp_path):
    report = WebReport(catalogue_version="v1.0.0")
    report.add_figure("a", _fig([1]), group="g")
    nested = tmp_path / "nested" / "report"

    report.build(str(nested))
    assert (nested / "index.html").exists()


def test_build_with_no_figures_writes_empty_index(tmp_path):
    report = WebReport(catalogue_version="v1.0.0")
    report.build(str(tmp_path))

    assert (tmp_path / "index.html").exists()
    assert not any(p.name != "index.html" for p in tmp_path.iterdir())


def test_build_multiple_figures_in_same_group_both_embedded(tmp_path):
    report = WebReport(catalogue_version="v1.0.0")
    report.add_figure("a", _fig([1]), title="Figure A", group="g")
    report.add_figure("b", _fig([2]), title="Figure B", group="g")
    report.build(str(tmp_path))

    page = (tmp_path / "g.html").read_text()
    assert "Figure A" in page
    assert "Figure B" in page
