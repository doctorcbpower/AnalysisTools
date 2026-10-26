"""
projection.py -- backend-neutral SPH projection of particle/cell data.

One set of routines (smoothing lengths, SPH density, orientation, weighted
projection) shared by every deposit backend, so that figures made with the
numba gridder, cosmokit, or anything registered later are directly comparable.

Conventions
-----------
* ``H`` is the SPH **support radius**: the distance to the k-th nearest
  neighbour (self excluded). A particle's cubic-spline kernel is zero beyond
  ``H`` (``GriddingTools.sph_assign`` takes ``h = H / 2``).
* Positions, centre, extents and ``H`` share one length unit.
* Maps are indexed ``[ix, iy]`` of the in-plane axes (plot ``image.T`` with
  ``origin="lower"``); ``Map.plot`` does this for you.
* The deposit is non-periodic. Periodic boxes are handled here, before the
  backend sees the data: positions are wrapped about ``centre`` (minimum
  image) and periodic images are added when the view region is larger than
  half the box.

Backends implement three methods (see ``NumbaBackend``): ``knn_h``,
``sph_density`` and ``deposit_projected``. ``get_backend`` resolves
``"numba"`` (default), ``"cosmokit"`` (optional; ``import cosmokit.backend``),
or any object with those methods.
"""
from __future__ import annotations

import itertools
import os
from dataclasses import dataclass
from typing import Optional

import numpy as np

# Proper rotation taking a disc-aligned frame (L along z) to an edge-on view:
# new x = x, new y = z, new z (line of sight) = -y.
_EDGE_ON = np.array([[1.0, 0.0, 0.0],
                     [0.0, 0.0, 1.0],
                     [0.0, -1.0, 0.0]])


# ---------------------------------------------------------------------------
# result container
# ---------------------------------------------------------------------------
@dataclass
class Map:
    """A projected image. ``image[ix, iy]``; ``extent = (xmin, xmax, ymin, ymax)``."""
    image: np.ndarray
    extent: tuple
    kind: str            # "surface_density" or "weighted_mean"

    def log10(self):
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.log10(self.image)

    def plot(self, ax=None, log=True, **imshow_kw):
        """imshow with the right orientation; returns (ax, image artist)."""
        import matplotlib.pyplot as plt
        ax = ax or plt.gca()
        data = self.log10() if log else self.image
        imshow_kw.setdefault("origin", "lower")
        im = ax.imshow(data.T, extent=self.extent, **imshow_kw)
        return ax, im


# ---------------------------------------------------------------------------
# backends
# ---------------------------------------------------------------------------
def _wrap(pos, boxsize):
    return np.mod(pos, boxsize) if boxsize else pos


def _cubic_w(r, H):
    """Cubic spline with support radius H (normalised, 3D)."""
    q = r / H
    w = np.where(q < 0.5, 1 - 6 * q**2 + 6 * q**3, 2 * (1 - q) ** 3)
    return np.where(q < 1.0, w, 0.0) * (8.0 / (np.pi * H**3))


class NumbaBackend:
    """Default backend: scipy k-d tree + ``GriddingTools.sph_assign`` (numba)."""
    name = "numba"

    @staticmethod
    def knn_h(pos, query, k, boxsize=None):
        from scipy.spatial import cKDTree
        tree = cKDTree(_wrap(pos, boxsize), boxsize=boxsize or None)
        # query points are members of the tree, so column 0 is the point itself
        return tree.query(_wrap(query, boxsize), k=k + 1, workers=-1)[0][:, -1]

    @staticmethod
    def sph_density(pos, mass, k, boxsize=None, include_self=True):
        from scipy.spatial import cKDTree
        p = _wrap(pos, boxsize)
        tree = cKDTree(p, boxsize=boxsize or None)
        dist, idx = tree.query(p, k=k + 1, workers=-1)
        H = dist[:, -1]
        w = mass[idx] * _cubic_w(dist, H[:, None])
        return w.sum(axis=1) if include_self else w[:, 1:].sum(axis=1)

    @staticmethod
    def deposit_projected(pos, values, H, lo, hi, shape):
        """Deposit each array in ``values`` onto a 3D grid over [lo, hi] with
        shape ``shape`` and sum along axis 2. Returns a list of 2D arrays."""
        from .gridding_tools import GriddingTools
        limits = (lo[0], hi[0], lo[1], hi[1], lo[2], hi[2])
        gs = np.asarray(shape, dtype=np.int64)
        pos = np.ascontiguousarray(pos, dtype=np.float64)
        h = np.ascontiguousarray(H, dtype=np.float64) / 2.0   # kernel support is 2h
        out = []
        for v in values:
            grid = np.zeros(tuple(int(n) for n in shape))
            GriddingTools.sph_assign(grid, pos, np.ascontiguousarray(v, dtype=np.float64),
                                     h, limits, gs, drop_outside=True)
            out.append(grid.sum(axis=2))
        return out


def get_backend(backend="auto"):
    """Resolve a backend name or object. ``"auto"`` is ``"numba"`` unless the
    ``ANALYSISTOOLS_PROJECTION_BACKEND`` environment variable says otherwise."""
    if not isinstance(backend, str):
        return backend
    if backend == "auto":
        backend = os.environ.get("ANALYSISTOOLS_PROJECTION_BACKEND", "numba")
    if backend == "numba":
        return NumbaBackend()
    if backend == "cosmokit":
        try:
            from cosmokit.backend import CosmokitBackend
        except ImportError as exc:
            raise ImportError("backend='cosmokit' needs cosmokit on PYTHONPATH "
                              "(built library + python/ directory)") from exc
        return CosmokitBackend()
    raise ValueError(f"unknown projection backend {backend!r}")


# ---------------------------------------------------------------------------
# public helpers
# ---------------------------------------------------------------------------
def _min_image(d, boxsize):
    return d - boxsize * np.round(d / boxsize) if boxsize else d


def smoothing_lengths(pos, k=32, boxsize=None, backend="auto"):
    """Support radius H (distance to the k-th neighbour) for every particle."""
    pos = np.ascontiguousarray(pos, dtype=np.float64)
    return get_backend(backend).knn_h(pos, pos, k, boxsize)


def sph_density(pos, mass, k=32, boxsize=None, backend="auto", include_self=True):
    """SPH density estimate from the k nearest neighbours.
    Use this to get a local density for particles that carry none (stars).
    ``include_self=True`` (the SPH convention, and cosmokit's default) adds the
    particle's own W(0, H), biasing a Poisson-sampled field high by roughly
    W(0, H) / rho (~20% for k=48); ``False`` removes it."""
    pos = np.ascontiguousarray(pos, dtype=np.float64)
    mass = np.ascontiguousarray(mass, dtype=np.float64)
    return get_backend(backend).sph_density(pos, mass, k, boxsize, include_self)


def orientation_from_angular_momentum(pos, vel, mass, r_max, centre=None, boxsize=None):
    """Rotation matrix ``R`` (rows = new x, y, z axes) with new z along the
    angular momentum of the particles within ``r_max`` of ``centre``.
    Apply as ``(pos - centre) @ R.T``; pass it to ``project(rotation=R)``."""
    centre = np.zeros(3) if centre is None else np.asarray(centre, dtype=float)
    d = _min_image(np.asarray(pos, dtype=float) - centre, boxsize)
    sel = np.linalg.norm(d, axis=1) < r_max
    if not np.any(sel):
        raise ValueError("no particles within r_max of centre")
    m = np.asarray(mass, dtype=float)[sel]
    v = np.asarray(vel, dtype=float)[sel]
    v = v - np.average(v, weights=m, axis=0)           # remove bulk motion
    L = np.sum(np.cross(d[sel], v) * m[:, None], axis=0)
    zhat = L / np.linalg.norm(L)
    ref = np.array([0.0, 1.0, 0.0]) if abs(zhat[1]) < 0.9 else np.array([1.0, 0.0, 0.0])
    xhat = np.cross(ref, zhat)
    xhat /= np.linalg.norm(xhat)
    return np.vstack([xhat, np.cross(zhat, xhat), zhat])


# ---------------------------------------------------------------------------
# projection
# ---------------------------------------------------------------------------
def _shifts(half3, boxsize):
    """Periodic image shifts (in box units) needed to cover a region with
    half-widths ``half3`` about the centre."""
    if not boxsize:
        return [(0, 0, 0)]
    radius = float(np.linalg.norm(half3))
    if radius <= 0.5 * boxsize:
        return [(0, 0, 0)]
    if radius > 1.5 * boxsize:
        raise ValueError("view region is more than 1.5 box sizes across; "
                         "periodic images beyond one box are not supported")
    return list(itertools.product((-1, 0, 1), repeat=3))


def _select(d, R, half3, margin, boxsize):
    """Indices and rotated coordinates of the particles (and periodic images)
    inside the region expanded by ``margin``. ``d`` is relative to the centre
    and already minimum-imaged."""
    lim = np.asarray(half3) + margin
    idx_parts, x_parts = [], []
    for s in _shifts(lim, boxsize):
        ds = d + boxsize * np.asarray(s) if boxsize else d
        x = ds @ R.T
        m = np.all(np.abs(x) < lim, axis=1)
        if np.any(m):
            idx_parts.append(np.nonzero(m)[0])
            x_parts.append(x[m])
    if not idx_parts:
        return np.empty(0, dtype=np.int64), np.empty((0, 3))
    return np.concatenate(idx_parts), np.vstack(x_parts)


def project(pos, weights=None, h=None, *, view="face", rotation=None, centre=None,
            half_extent, half_depth=None, npix=512, nz=None, weight_by=None,
            boxsize=None, k=32, backend="auto") -> Map:
    """
    Project SPH-smoothed particle data onto an image.

    Parameters
    ----------
    pos : (N, 3)
    weights : (N,), optional
        Quantity to project (e.g. mass). Default 1 per particle.
    h : (N,), optional
        Per-particle support radius. Default: distance to the ``k``-th neighbour.
    view : "face" | "edge"
        Line of sight along z ("face") or along -y ("edge") of the frame set by
        ``rotation`` (identity if None). With ``rotation`` from
        ``orientation_from_angular_momentum`` these are disc face-on/edge-on.
    rotation : (3, 3), optional
        Applied to ``pos - centre`` before the view.
    centre : (3,), optional
        View centre (default origin).
    half_extent, half_depth : float
        Half-width of the square image and half-depth of the slab (default =
        ``half_extent``), in the rotated frame.
    npix : int
        Pixels per side. ``nz`` is the depth cell count (default ``npix // 4``);
        backends with isotropic cells may ignore it.
    weight_by : (N,), optional
        If given, return the ``weights``-weighted mean of this quantity along the
        line of sight, Σ w f / Σ w. Omitted: the surface density Σ w / dA.
        E.g. ``weights=mass, weight_by=rho`` is the ρ-weighted density
        (∫ρ² dz / ∫ρ dz), since m = ρV.
    boxsize : float, optional
        Periodic box size; positions are wrapped about ``centre``.
    """
    be = get_backend(backend)
    pos = np.ascontiguousarray(pos, dtype=np.float64)
    n = len(pos)
    w = np.ones(n) if weights is None else np.asarray(weights, dtype=np.float64)
    centre = np.zeros(3) if centre is None else np.asarray(centre, dtype=np.float64)
    half_depth = half_extent if half_depth is None else half_depth
    nz = max(8, npix // 4) if nz is None else int(nz)

    R = np.eye(3) if rotation is None else np.asarray(rotation, dtype=np.float64)
    if view == "edge":
        R = _EDGE_ON @ R
    elif view != "face":
        raise ValueError("view must be 'face' or 'edge'")

    half3 = (half_extent, half_extent, half_depth)
    d = _min_image(pos - centre, boxsize)

    # Pass 1: particles in the region -> their H sets the margin of particles
    # outside it whose kernels still reach in.
    idx, _ = _select(d, R, half3, 0.0, boxsize)
    if len(idx) == 0:
        raise ValueError("no particles inside the view region")
    h_arr = None if h is None else np.asarray(h, dtype=np.float64)

    def support(ix):
        return h_arr[ix] if h_arr is not None else be.knn_h(pos, pos[ix], k, boxsize)

    margin = min(float(support(np.unique(idx)).max()), 0.5 * max(half_extent, half_depth))

    # Pass 2: expanded selection, with H for every selected particle
    idx, x = _select(d, R, half3, margin, boxsize)
    uniq, inv = np.unique(idx, return_inverse=True)
    H = support(uniq)[inv]

    lo = (-half_extent, -half_extent, -half_depth)
    hi = (half_extent, half_extent, half_depth)
    ws = w[idx]
    vals = [ws] if weight_by is None else [ws * np.asarray(weight_by, dtype=np.float64)[idx], ws]
    maps = be.deposit_projected(x, vals, H, lo, hi, (npix, npix, nz))

    extent = (-half_extent, half_extent, -half_extent, half_extent)
    if weight_by is None:
        dA = (2.0 * half_extent / npix) ** 2
        return Map(maps[0] / dA, extent, "surface_density")
    num, den = maps
    with np.errstate(invalid="ignore", divide="ignore"):
        img = num / np.where(den > 0, den, np.nan)
    return Map(img, extent, "weighted_mean")
