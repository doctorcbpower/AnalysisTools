"""analysistools.projection: backend-neutral SPH projection."""
import numpy as np
import pytest

from analysistools import projection as pj

BACKENDS = ["numba"]
try:
    import cosmokit.backend  # noqa: F401
    BACKENDS.append("cosmokit")
except Exception:                      # not built / not on PYTHONPATH
    pass


def _blob(n=4000, seed=0, scale=1.0):
    rng = np.random.default_rng(seed)
    return rng.normal(scale=scale, size=(n, 3))


def _disc(n=6000, seed=1, rot=(0.7, 1.1)):
    """Thin rotating exponential-ish disc, then tilted by `rot` angles."""
    rng = np.random.default_rng(seed)
    r = rng.exponential(2.0, n)
    phi = rng.uniform(0, 2 * np.pi, n)
    pos = np.c_[r * np.cos(phi), r * np.sin(phi), rng.normal(scale=0.1, size=n)]
    vel = np.c_[-np.sin(phi), np.cos(phi), np.zeros(n)] * 5.0
    a, b = rot
    Rx = np.array([[1, 0, 0], [0, np.cos(a), -np.sin(a)], [0, np.sin(a), np.cos(a)]])
    Ry = np.array([[np.cos(b), 0, np.sin(b)], [0, 1, 0], [-np.sin(b), 0, np.cos(b)]])
    T = Ry @ Rx
    return pos @ T.T, vel @ T.T, np.ones(n)


@pytest.mark.parametrize("backend", BACKENDS)
def test_mass_conserved(backend):
    pos = _blob()
    m = np.full(len(pos), 2.5)
    img = pj.project(pos, m, half_extent=8.0, npix=96, backend=backend)
    dA = (2 * 8.0 / 96) ** 2
    assert img.kind == "surface_density"
    assert img.image.sum() * dA == pytest.approx(m.sum(), rel=0.02)


@pytest.mark.parametrize("backend", BACKENDS)
def test_weighted_mean_of_constant(backend):
    pos = _blob()
    m = np.ones(len(pos))
    img = pj.project(pos, m, half_extent=6.0, npix=64, weight_by=np.full(len(pos), 3.0),
                     backend=backend)
    good = np.isfinite(img.image)
    assert good.sum() > 100
    np.testing.assert_allclose(img.image[good], 3.0, rtol=1e-6)


@pytest.mark.parametrize("backend", BACKENDS)
def test_disc_face_and_edge(backend):
    pos, vel, m = _disc()
    R = pj.orientation_from_angular_momentum(pos, vel, m, r_max=6.0)
    # new z axis is the disc normal: thin in z, wide in x/y
    xr = pos @ R.T
    assert xr[:, 2].std() < 0.2 and xr[:, 0].std() > 1.5

    def axis_ratio(img):
        n = img.image.shape[0]
        g = (np.arange(n) + 0.5) / n * (img.extent[1] - img.extent[0]) + img.extent[0]
        X, Y = np.meshgrid(g, g, indexing="ij")
        w = img.image
        sx = np.sqrt((w * X**2).sum() / w.sum())
        sy = np.sqrt((w * Y**2).sum() / w.sum())
        return sy / sx

    face = pj.project(pos, m, rotation=R, view="face", half_extent=10, npix=80, backend=backend)
    edge = pj.project(pos, m, rotation=R, view="edge", half_extent=10, npix=80, backend=backend)
    assert axis_ratio(face) == pytest.approx(1.0, abs=0.1)
    assert axis_ratio(edge) < 0.25


@pytest.mark.parametrize("backend", BACKENDS)
def test_periodic_matches_shifted(backend):
    """A blob straddling the box edge, viewed through the periodic wrap, equals
    the same blob viewed in open space."""
    L = 20.0
    blob = _blob(n=3000, seed=3, scale=0.8)
    m = np.ones(len(blob))
    open_img = pj.project(blob, m, half_extent=4.0, npix=48, backend=backend)
    wrapped = np.mod(blob + 0.0, L)                 # puts half the blob near L
    per_img = pj.project(wrapped, m, centre=np.zeros(3), boxsize=L, half_extent=4.0,
                         npix=48, backend=backend)
    np.testing.assert_allclose(per_img.image, open_img.image, rtol=1e-6, atol=1e-8)


@pytest.mark.parametrize("backend", BACKENDS)
def test_periodic_images_when_region_exceeds_half_box(backend):
    L = 10.0
    rng = np.random.default_rng(4)
    pos = rng.random((4000, 3)) * L
    m = np.ones(len(pos))
    img = pj.project(pos, m, centre=np.full(3, L / 2), boxsize=L, half_extent=L / 2,
                     half_depth=L / 2, npix=40, backend=backend)
    dA = (L / 40) ** 2
    assert img.image.sum() * dA == pytest.approx(m.sum(), rel=0.03)


def test_sph_density_uniform():
    rng = np.random.default_rng(5)
    L = 10.0
    pos = rng.random((20000, 3)) * L
    for b in BACKENDS:
        rho = pj.sph_density(pos, np.ones(len(pos)), k=48, boxsize=L, backend=b,
                             include_self=False)
        assert np.median(rho) == pytest.approx(len(pos) / L**3, rel=0.1)
    if len(BACKENDS) > 1:        # both conventions agree between backends
        a = pj.sph_density(pos, np.ones(len(pos)), k=48, boxsize=L, backend="numba")
        c = pj.sph_density(pos, np.ones(len(pos)), k=48, boxsize=L, backend="cosmokit")
        np.testing.assert_allclose(a, c, rtol=1e-6)


@pytest.mark.skipif(len(BACKENDS) < 2, reason="cosmokit backend unavailable")
def test_backends_agree():
    pos = _blob(n=8000, seed=6)
    m = np.ones(len(pos))
    a = pj.project(pos, m, half_extent=5.0, npix=64, backend="numba").image
    b = pj.project(pos, m, half_extent=5.0, npix=64, backend="cosmokit").image
    assert np.abs(a - b).sum() / a.sum() < 0.08
    ha = pj.smoothing_lengths(pos, backend="numba")
    hb = pj.smoothing_lengths(pos, backend="cosmokit")
    np.testing.assert_allclose(ha, hb, rtol=1e-10)


def test_sph_assign_drop_outside_loses_mass_off_grid():
    """A particle centred on the grid edge keeps ~half its mass with
    drop_outside=True; the default renormalises it all onto the grid."""
    from analysistools.gridding_tools import GriddingTools
    pos = np.array([[0.0, 5.0, 5.0]])
    lim, gs = (0.0, 10.0, 0.0, 10.0, 0.0, 10.0), np.array([40, 40, 40], dtype=np.int64)
    kept = {}
    for flag in (False, True):
        g = np.zeros((40, 40, 40))
        GriddingTools.sph_assign(g, pos, np.array([1.0]), np.array([1.0]), lim, gs,
                                 drop_outside=flag)
        kept[flag] = g.sum()
    assert kept[False] == pytest.approx(1.0)
    assert kept[True] == pytest.approx(0.5, abs=0.08)
