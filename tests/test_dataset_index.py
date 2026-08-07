"""
Tests for api.dataset.Dataset.index (added to let backends.py map a
matched-view row back to its row in the model's own tables, e.g. for SFH
lookups keyed by galaxy row).
"""
import os

import numpy as np
import pytest

import analysistools as at

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
VRCAT = os.path.join(HERE, "data", "halos",
                     "snap_0031.VELOCIraptor.properties.0")
have_data = os.path.exists(VRCAT)
needs_data = pytest.mark.skipif(not have_data, reason="example data missing")


@needs_data
def test_full_table_index_is_arange():
    cat = at.load(VRCAT, kind="halos", fileformat="VELOCIraptor",
                  native_includes_h=False, snapnum=31)
    np.testing.assert_array_equal(cat.index, np.arange(len(cat)))
    assert cat.is_view is False


@needs_data
def test_view_index_matches_selection_rows():
    cat = at.load(VRCAT, kind="halos", fileformat="VELOCIraptor",
                  native_includes_h=False, snapnum=31)
    mass = np.asarray(cat["mass"])
    expected_rows = np.flatnonzero(mass > np.median(mass))

    view = cat.select(mask=mass > np.median(mass))
    assert view.is_view is True
    np.testing.assert_array_equal(view.index, expected_rows)
    # the index really does recover the parent's values at those rows
    np.testing.assert_array_equal(view["mass"], mass[expected_rows])


@needs_data
def test_nested_selection_composes_index():
    cat = at.load(VRCAT, kind="halos", fileformat="VELOCIraptor",
                  native_includes_h=False, snapnum=31)
    mass = np.asarray(cat["mass"])

    once = cat.select(mask=mass > np.percentile(mass, 50))
    twice = once.select(mask=np.asarray(once["mass"])
                        > np.percentile(once["mass"], 50))

    # twice.index must be indices into the ORIGINAL cat, not into `once`
    np.testing.assert_array_equal(twice["mass"], mass[twice.index])
