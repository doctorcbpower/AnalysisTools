"""
Tests for GriddingTools' projection-plotting infrastructure:
- _collapse's normalize_by_area (sum -> areal density)
- _shared_clim / plot_3d_projections defaulting to one colour scale shared
  across all three panels, rather than each auto-scaling independently
  (previously the case -- see plot_3d_projections's old implementation,
  which re-imshow'd each panel's array with no vmin/vmax at all)
- plot_ratio_projections (e.g. dust-to-gas ratio: two grids collapsed and
  divided, not divided cell-by-cell before projection)
"""
import matplotlib
matplotlib.use("Agg")

import numpy as np
import pytest

from analysistools.gridding_tools import GriddingTools

GRID_LIMITS = [-10.0, 10.0, -10.0, 10.0, -10.0, 10.0]


def _uniform_grid(shape, value):
    return np.full(shape, value, dtype=float)


def test_normalize_by_area_converts_sum_to_density():
    # A uniform grid with every voxel = 1.0, 10x10x10, over a 20x20x20 box
    # (voxel size 2x2x2). Summing along z gives 10 per column (10 voxels
    # each worth 1.0); the XY pixel area is 2*2=4, so the areal density
    # should be 10/4 = 2.5 everywhere.
    grid = _uniform_grid((10, 10, 10), 1.0)
    gt = GriddingTools()

    data_raw, *_ = gt._collapse(grid, GRID_LIMITS, slice_axis='z',
                                 mode='projection', projection='sum',
                                 normalize_by_area=False)
    data_density, *_ = gt._collapse(grid, GRID_LIMITS, slice_axis='z',
                                     mode='projection', projection='sum',
                                     normalize_by_area=True)

    np.testing.assert_allclose(data_raw, 10.0)
    np.testing.assert_allclose(data_density, 2.5)


def test_shared_clim_ignores_nonfinite_values():
    a = np.array([1.0, 10.0, 100.0])
    b = np.array([0.0, 1000.0])  # log10(0) -> -inf, must not corrupt the max
    vmin, vmax = GriddingTools._shared_clim([a, b])

    np.testing.assert_allclose(vmin, 0.0)   # log10(1.0)
    np.testing.assert_allclose(vmax, 3.0)   # log10(1000.0)


def test_plot_3d_projections_defaults_to_one_shared_scale():
    # A grid with a very different dynamic range along each axis pairing
    # would, under independent auto-scaling, give each of the 3 panels a
    # different colour range. With the new shared default, all 3 images
    # must carry the identical vmin/vmax.
    grid = np.ones((8, 8, 8))
    grid[0, 0, 0] = 1e6  # one hot voxel, visible from any projection axis

    fig, axes = GriddingTools().plot_3d_projections(grid, GRID_LIMITS, projection='sum')

    clims = [ax.images[0].get_clim() for ax in axes]
    assert clims[0] == clims[1] == clims[2]


def test_plot_3d_projections_explicit_vmin_vmax_overrides_default():
    grid = np.ones((8, 8, 8))
    fig, axes = GriddingTools().plot_3d_projections(
        grid, GRID_LIMITS, projection='sum', vmin=-5.0, vmax=5.0,
    )
    for ax in axes:
        assert ax.images[0].get_clim() == (-5.0, 5.0)


def test_plot_ratio_projections_matches_manual_ratio():
    numerator = _uniform_grid((8, 8, 8), 2.0)
    denominator = _uniform_grid((8, 8, 8), 4.0)

    fig, axes = GriddingTools().plot_ratio_projections(
        numerator, denominator, GRID_LIMITS, projection='sum',
    )

    # sum over z: numerator column = 2*8=16, denominator column = 4*8=32,
    # ratio = 0.5 everywhere -> log10(0.5)
    expected = np.log10(0.5)
    for ax in axes:
        np.testing.assert_allclose(ax.images[0].get_array(), expected)


def test_plot_ratio_projections_masks_small_denominator():
    numerator = _uniform_grid((8, 8, 8), 1.0)
    denominator = _uniform_grid((8, 8, 8), 1.0)
    denominator[0, 0, :] = 0.0  # this XY column's denominator collapses to 0

    fig, axes = GriddingTools().plot_ratio_projections(
        numerator, denominator, GRID_LIMITS, projection='sum',
    )

    xy_data = axes[0].images[0].get_array()
    mask = np.ma.getmaskarray(xy_data)
    assert mask[0, 0]
    assert not np.any(mask[1:, 1:])


def test_plot_ratio_projections_requires_matching_shapes():
    a = _uniform_grid((8, 8, 8), 1.0)
    b = _uniform_grid((4, 4, 4), 1.0)
    with pytest.raises(ValueError):
        GriddingTools().plot_ratio_projections(a, b, GRID_LIMITS)
