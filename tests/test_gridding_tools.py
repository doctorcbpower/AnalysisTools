"""
Tests for analysistools.gridding_tools.GriddingTools.sph_assign /
smooth_to_grid(method="SPH") -- per-particle SPH-kernel (cubic spline)
deposition using each particle's own smoothing length, added to replace a
single global Gaussian sigma that either over-smooths dense regions or
leaves particle-scale gaps in sparse ones (see the AGORA disc example's
gas-mass grid projection, which showed exactly that blob/gap pattern).
"""
import numpy as np
import pytest

from analysistools.gridding_tools import GriddingTools


GRID_LIMITS = [-10.0, 10.0, -10.0, 10.0, -10.0, 10.0]
GRID_SIZE = np.array([64, 64, 64], dtype=np.int64)


def _new_grid():
    return np.zeros(tuple(GRID_SIZE), dtype=float)


def test_single_particle_mass_conserved_away_from_edge():
    grid = _new_grid()
    pos = np.array([[0.0, 0.0, 0.0]])
    values = np.array([5.0])
    h = np.array([1.0])

    GriddingTools().sph_assign(grid, pos, values, h, GRID_LIMITS, GRID_SIZE)

    np.testing.assert_allclose(grid.sum(), 5.0, rtol=1e-10)
    assert np.all(grid >= 0.0)


def test_conservation_holds_across_a_range_of_smoothing_lengths():
    # h both much smaller and much larger than the grid spacing
    # (spacing = 20/64 = 0.3125) -- conservation must hold either way,
    # unlike a single global smoothing scale.
    for h_val in [0.01, 0.3125, 1.0, 3.0]:
        grid = _new_grid()
        pos = np.array([[1.3, -2.7, 0.4]])
        values = np.array([2.0])
        h = np.array([h_val])

        GriddingTools().sph_assign(grid, pos, values, h, GRID_LIMITS, GRID_SIZE)

        np.testing.assert_allclose(grid.sum(), 2.0, rtol=1e-8,
                                    err_msg=f"failed for h={h_val}")


def test_degenerate_smoothing_length_falls_back_to_ngp():
    grid = _new_grid()
    pos = np.array([[0.4, 0.4, 0.4]])
    values = np.array([7.0])
    h = np.array([0.0])

    GriddingTools().sph_assign(grid, pos, values, h, GRID_LIMITS, GRID_SIZE)

    spacing = 20.0 / 64
    i = int(np.floor((0.4 - (-10.0)) / spacing))
    assert grid[i, i, i] == pytest.approx(7.0)
    np.testing.assert_allclose(grid.sum(), 7.0)


def test_multiple_particles_independently_conserved():
    grid = _new_grid()
    pos = np.array([[-5.0, -5.0, -5.0], [3.0, 3.0, 3.0], [0.0, 0.0, 0.0]])
    values = np.array([1.0, 2.0, 3.0])
    h = np.array([0.5, 1.5, 0.2])

    GriddingTools().sph_assign(grid, pos, values, h, GRID_LIMITS, GRID_SIZE)

    np.testing.assert_allclose(grid.sum(), 6.0, rtol=1e-8)


def test_deposit_is_more_concentrated_for_smaller_smoothing_length():
    # A smaller h should spread its value over fewer nonzero cells than a
    # larger h -- the actual property that fixes the "one global sigma
    # either over- or under-smooths" artifact.
    small_h_grid = _new_grid()
    large_h_grid = _new_grid()
    pos = np.array([[0.0, 0.0, 0.0]])
    values = np.array([1.0])

    GriddingTools().sph_assign(small_h_grid, pos, values, np.array([0.2]),
                                GRID_LIMITS, GRID_SIZE)
    GriddingTools().sph_assign(large_h_grid, pos, values, np.array([2.0]),
                                GRID_LIMITS, GRID_SIZE)

    assert np.count_nonzero(small_h_grid) < np.count_nonzero(large_h_grid)


def test_symmetric_under_anisotropic_grid_spacing():
    # dz is half of dx/dy here (mirrors the AGORA example's lzgrid=lgrid/4
    # with the same ngrid per axis) -- deposition must be based on
    # physical distance, not grid-index distance, so an isotropic particle
    # kernel should still look isotropic in physical units despite the
    # grid being finer along z. Odd grid sizes so a cell centre lands
    # exactly on the particle position (0,0,0); dz = dx/4 exactly, so an
    # offset of N cells in z reaches exactly the same physical distance
    # as N/4 cells in x, with no rounding error.
    limits = [-10.0, 10.0, -10.0, 10.0, -2.5, 2.5]
    grid_size = np.array([65, 65, 65], dtype=np.int64)
    grid = np.zeros(tuple(grid_size), dtype=float)
    pos = np.array([[0.0, 0.0, 0.0]])
    values = np.array([1.0])
    h = np.array([1.0])

    GriddingTools().sph_assign(grid, pos, values, h, limits, grid_size)

    centre = 32  # cell centre == 0.0 exactly for odd grid_size, see above
    ix_off, iz_off = 3, 12  # 3*dx == 12*dz == 3*(20/65)
    w_x = grid[centre + ix_off, centre, centre]
    w_z = grid[centre, centre, centre + iz_off]
    np.testing.assert_allclose(w_x, w_z, rtol=1e-8)


def test_smooth_to_grid_sph_requires_smoothing_lengths():
    with pytest.raises(ValueError):
        GriddingTools().smooth_to_grid(
            positions=np.zeros((1, 3)), values=np.ones(1),
            grid_size=GRID_SIZE, grid_limits=GRID_LIMITS, method="SPH",
        )


def test_smooth_to_grid_sph_matches_direct_sph_assign():
    pos = np.array([[1.0, -2.0, 0.5], [-3.0, 4.0, -1.0]])
    values = np.array([1.5, 2.5])
    h = np.array([0.5, 1.2])

    via_smooth_to_grid = GriddingTools().smooth_to_grid(
        positions=pos, values=values, grid_size=GRID_SIZE,
        grid_limits=GRID_LIMITS, method="SPH", smoothing_lengths=h,
    )

    direct = _new_grid()
    GriddingTools().sph_assign(direct, pos, values, h, GRID_LIMITS, GRID_SIZE)

    np.testing.assert_array_equal(via_smooth_to_grid, direct)
