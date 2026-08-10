#!/usr/bin/env python3
"""
analysistools.api.simulation
----------------------------
Simulation: binds the matched data products of one simulation (snapshots,
halo catalogues, merger trees, SHARK galaxies) and does the cross-matching
that otherwise has to be done by hand.

    sim = at.Simulation(
        snapshots = {0.0: "snap_0031.hdf5"},
        halos     = {0.0: "snap_0031.VELOCIraptor.properties.0"},
        trees     = "VELOCIraptor.walkabletree.hdf5",
        galaxies  = {0.0: "199/0/galaxies.hdf5"},
        snapnums  = {0.0: 31},
        label     = "CDM",
    )
    epoch = sim.at(0.0)                    # Epoch of matched Dataset views
    parts = epoch.particles_in_halo(index=0)
    gals  = epoch.galaxies_in_halo(index=0)          # position match
    track = sim.track_of(index=0, redshift=0.0)

Design decisions (DEVELOPMENT.md section 7): SHARK<->halo matching is by
position (periodic KD-tree) by default because id_halo consistency with
the N-body catalogue is not guaranteed; pass match_by="id" to use IDs.
Every component is optional -- a Simulation with only snapshots is fine.
"""
from __future__ import annotations

from typing import Any, Dict, Optional, Sequence, Union

import logging

import numpy as np

from .dataset import Dataset
from .galaxies import GalaxyCatalogue
from .halos import HaloCatalogue
from .snapshot import SnapshotDataset
from .trees import MergerTree, TrackDataset

logger = logging.getLogger(__name__)

# component spec: prebuilt object | path | {redshift: path} | EpochModel
Spec = Union[str, Dict[float, str], Dataset, Any, None]


# ---------------------------------------------------------------------------
# Position cross-matching (periodic KD-tree)
# ---------------------------------------------------------------------------

def match_positions(points: np.ndarray, centre: Sequence[float],
                    radius: float,
                    boxsize: Optional[float] = None) -> np.ndarray:
    """Indices of `points` (N, 3) within `radius` of `centre`, with
    periodic wrapping when `boxsize` is given."""
    from scipy.spatial import cKDTree
    points = np.asarray(points, dtype=float)
    centre = np.asarray(centre, dtype=float)
    if boxsize:
        tree = cKDTree(np.mod(points, boxsize), boxsize=boxsize)
        centre = np.mod(centre, boxsize)
    else:
        tree = cKDTree(points)
    return np.asarray(sorted(tree.query_ball_point(centre, radius)),
                      dtype=np.int64)


# ---------------------------------------------------------------------------
# Epoch: matched single-redshift views
# ---------------------------------------------------------------------------

class Epoch:
    """The components of a Simulation frozen at (nearest to) one redshift.

    Components load lazily on first attribute access and are cached.
    """

    def __init__(self, sim: "Simulation", redshift: float):
        self.sim = sim
        self.redshift = float(redshift)
        self._cache: Dict[str, Any] = {}

    # -- component access ------------------------------------------------

    def _component(self, kind: str):
        if kind not in self._cache:
            self._cache[kind] = self.sim._resolve(kind, self.redshift)
        return self._cache[kind]

    @property
    def snapshot(self) -> Optional[SnapshotDataset]:
        """Snapshot Dataset at this epoch's redshift, or None if unset."""
        return self._component("snapshots")

    @property
    def halos(self) -> Optional[HaloCatalogue]:
        """Halo catalogue at this epoch's redshift, or None if unset."""
        return self._component("halos")

    @property
    def galaxies(self) -> Optional[GalaxyCatalogue]:
        """Galaxy catalogue at this epoch's redshift, or None if unset."""
        return self._component("galaxies")

    @property
    def tree(self) -> Optional[MergerTree]:
        """The parent Simulation's merger tree (not epoch-specific)."""
        return self.sim.tree

    @property
    def snapnum(self) -> Optional[int]:
        """
        Snapshot number at this epoch's redshift.

        Returns
        -------
        int or None
            From `sim.snapnums` if given, else from the halo catalogue's
            metadata; None if neither is available.
        """
        sn = self.sim._snapnum_at(self.redshift)
        if sn is not None:
            return sn
        h = self.halos
        return h.meta.get("snapnum") if h is not None else None

    @property
    def boxsize(self) -> Optional[float]:
        """
        Simulation box size, read from whichever component's metadata has
        it, preferring already-loaded components and cheaper ones (halos
        before snapshot/galaxies) before forcing loads.

        Returns
        -------
        float or None
        """
        comps = [self.halos, self.snapshot, self.galaxies]
        # first from anything already loaded, then force loads in order of
        # increasing cost (halo catalogues are far lighter than snapshots)
        for force in (False, True):
            for ds in comps:
                if ds is None:
                    continue
                if force:
                    ds.preload()
                if ds.meta.get("boxsize"):
                    return float(ds.meta["boxsize"])
        return None

    # -- halo lookup helper ----------------------------------------------

    def _halo_row(self, halo_id: Optional[int] = None,
                  index: Optional[int] = None) -> int:
        if (halo_id is None) == (index is None):
            raise ValueError("Pass exactly one of halo_id= or index=.")
        cat = self.halos
        if cat is None:
            raise ValueError("This Simulation has no halo catalogue.")
        if index is not None:
            return int(index)
        ids = np.asarray(cat["halo_id"])
        rows = np.flatnonzero(ids == halo_id)
        if not rows.size:
            raise KeyError(f"halo_id={halo_id} not found in catalogue "
                           f"at z={self.redshift:g}.")
        return int(rows[0])

    # -- cross-matching ----------------------------------------------------

    def particles_in_halo(self, halo_id: Optional[int] = None, *,
                          index: Optional[int] = None,
                          r_scale: float = 1.0,
                          species: Optional[str] = None,
                          match_by: str = "auto") -> SnapshotDataset:
        """
        Snapshot particles belonging to one halo, as a Dataset view.

        Parameters
        ----------
        halo_id / index :
            The halo, by catalogue ID or by row index (exactly one).
        r_scale : float
            Radius multiplier for the position query (r < r_scale * radius).
        species : str, optional
            Restrict to one species ("dm", "gas", ...).
        match_by : {"auto", "groupid", "radius"}
            "groupid" uses the snapshot's stored group membership;
            "radius" does a periodic radius query around the halo centre;
            "auto" (default) prefers groupid when the snapshot carries it.
        """
        snap = self.snapshot
        if snap is None:
            raise ValueError("This Simulation has no snapshot.")
        row = self._halo_row(halo_id, index)
        cat = self.halos
        base = snap.select(species=species) if species else snap

        has_groupid = "groupid" in base
        if match_by == "auto":
            match_by = "groupid" if has_groupid else "radius"
        if match_by == "groupid":
            if not has_groupid:
                raise ValueError("Snapshot has no groupid block; use "
                                 "match_by='radius'.")
            gid = halo_id if halo_id is not None \
                else int(np.asarray(cat["halo_id"])[row])
            return base.select(mask=np.asarray(base["groupid"]) == gid)
        if match_by != "radius":
            raise ValueError("match_by must be 'auto', 'groupid' or "
                             "'radius'.")

        centre = np.asarray(cat["pos"])[row]
        radius = float(np.asarray(cat["radius"])[row]) * r_scale
        idx = match_positions(base["pos"], centre, radius, self.boxsize)
        view = base.select(mask=idx)
        view.label = f"{snap.label}:halo{cat['halo_id'][row]}"
        return view

    def galaxies_in_halo(self, halo_id: Optional[int] = None, *,
                         index: Optional[int] = None,
                         r_scale: float = 1.0,
                         match_by: str = "position") -> GalaxyCatalogue:
        """
        SHARK galaxies belonging to one halo, as a Dataset view.

        match_by : {"position", "id"}
            "position" (default): periodic radius query around the halo
            centre (id_halo consistency with the N-body catalogue is not
            guaranteed). "id": match galaxies["id_halo"] == halo_id.

        Position matching converts the halo catalogue's centre/radius/
        boxsize into SHARK's own fixed native convention (comoving
        Mpc/h -- see GalaxyCatalogue's meta["units"]) before the search,
        using the halo catalogue's own meta["comoving"]/["little_h"]/
        ["h0"]/["scale_factor"] -- the halo catalogue can be loaded in
        *any* comoving/little_h state (the default is comoving, h-free,
        which does *not* match SHARK's h-included convention) and this
        still matches correctly. If the halo catalogue has no HubbleParam
        available, no h correction can be applied and a warning is
        logged -- position matching may then silently find zero (or the
        wrong) galaxies whenever the two catalogues' conventions differ.
        """
        gals = self.galaxies
        if gals is None:
            raise ValueError("This Simulation has no galaxy catalogue.")
        cat = self.halos

        if match_by == "id":
            if halo_id is None:
                if cat is None:
                    raise ValueError("index= requires a halo catalogue.")
                halo_id = int(np.asarray(cat["halo_id"])
                              [self._halo_row(None, index)])
            return gals.select(mask=np.asarray(gals["id_halo"]) == halo_id)
        if match_by != "position":
            raise ValueError("match_by must be 'position' or 'id'.")

        if cat is None:
            raise ValueError("Position matching requires a halo catalogue.")
        row = self._halo_row(halo_id, index)
        length_factor = self._length_factor_to_shark_native(cat)
        centre = np.asarray(cat["pos"])[row] * length_factor
        radius = float(np.asarray(cat["radius"])[row]) * r_scale * length_factor
        boxsize = cat.meta.get("boxsize")
        boxsize = boxsize * length_factor if boxsize is not None \
            else self.boxsize
        idx = match_positions(gals["pos"], centre, radius, boxsize)
        view = gals.select(mask=idx)
        view.label = f"{gals.label}:halo{cat['halo_id'][row]}"
        return view

    @staticmethod
    def _length_factor_to_shark_native(cat: HaloCatalogue) -> float:
        """Multiplicative factor converting a length (position/radius/
        boxsize) from `cat`'s own comoving/little_h convention into
        SHARK's fixed native convention (comoving Mpc/h). ``1.0`` (no
        correction, with a warning) if `cat` has no HubbleParam to
        convert with -- see ``galaxies_in_halo``'s docstring."""
        cat.preload()  # meta isn't populated until first load
        h = cat.meta.get("h0")
        if not h:
            logger.warning(
                "galaxies_in_halo(): halo catalogue has no HubbleParam "
                "(meta['h0']) -- cannot convert to SHARK's native "
                "comoving Mpc/h convention. Position matching will be "
                "wrong (likely finding zero galaxies) if the halo "
                "catalogue's actual comoving/little_h state differs "
                "from SHARK's.")
            return 1.0
        factor = 1.0
        if not cat.meta.get("little_h"):
            factor *= h                        # h-free -> h-included
        if not cat.meta.get("comoving"):
            a = cat.meta.get("scale_factor")
            if a:
                factor /= a                     # physical -> comoving
            else:
                logger.warning(
                    "galaxies_in_halo(): halo catalogue is physical "
                    "(comoving=False) but has no scale_factor -- cannot "
                    "convert to SHARK's comoving convention; only the "
                    "little_h correction is applied.")
        return factor

    def track_of(self, halo_id: Optional[int] = None, *,
                 index: Optional[int] = None,
                 object_type: str = "Group") -> TrackDataset:
        """Merger-tree main branch of one catalogue halo (needs trees=)."""
        tree = self.tree
        if tree is None:
            raise ValueError("This Simulation has no merger tree.")
        row = self._halo_row(halo_id, index)
        if self.snapnum is not None:
            tree.link_halos(self.halos, snapnum=self.snapnum)
        return tree.from_halo(self.halos, row, object_type=object_type,
                              snapnum=self.snapnum)

    def summary(self) -> None:
        """
        Print this epoch's redshift/snapnum, the load state of each
        component (snapshots/halos/galaxies), and whether a tree is set.
        """
        print(f"Epoch z={self.redshift:g} of Simulation "
              f"'{self.sim.label}'"
              + (f" (snapnum {self.snapnum})" if self.snapnum is not None
                 else ""))
        for kind in ("snapshots", "halos", "galaxies"):
            spec = self.sim._specs.get(kind)
            state = "-" if spec is None else (
                repr(self._cache[kind]) if kind in self._cache
                else "available (not loaded)")
            print(f"  {kind:10s}: {state}")
        print(f"  tree      : "
              f"{self.sim.tree if self.sim._specs.get('trees') else '-'}")

    def __repr__(self) -> str:
        have = [k for k in ("snapshots", "halos", "galaxies", "trees")
                if self.sim._specs.get(k) is not None]
        return (f"<Epoch z={self.redshift:g} of '{self.sim.label}' "
                f"[{', '.join(have)}]>")


# ---------------------------------------------------------------------------
# Simulation
# ---------------------------------------------------------------------------

class Simulation:
    """
    Matched data products of one simulation. All components optional.

    Parameters
    ----------
    snapshots, halos, galaxies :
        Each may be: a path (single epoch); a dict {redshift: path};
        a prebuilt Dataset (single epoch); or a multi-epoch model
        (HaloModel / SharkModel -- anything with .at(z)).
    trees :
        Tree file path or a prebuilt MergerTree (trees span all epochs).
    fileformats : dict, optional
        Per-kind format overrides, e.g. dict(halos="VELOCIraptor").
    snapnums : dict {redshift: int}, optional
        Snapshot numbers per redshift (needed for tree lookups when the
        halo catalogue path doesn't record one).
    load_kwargs : dict, optional
        Per-kind extra kwargs passed to load(), e.g.
        dict(snapshots=dict(convention="SWIFT")).
    label : str
    """

    def __init__(self, snapshots: Spec = None, halos: Spec = None,
                 trees: Spec = None, galaxies: Spec = None,
                 fileformats: Optional[Dict[str, str]] = None,
                 snapnums: Optional[Dict[float, int]] = None,
                 load_kwargs: Optional[Dict[str, dict]] = None,
                 label: str = "simulation"):
        self._specs: Dict[str, Spec] = {
            "snapshots": snapshots, "halos": halos,
            "trees": trees, "galaxies": galaxies,
        }
        self.fileformats = fileformats or {}
        self.snapnums = {float(k): int(v)
                         for k, v in (snapnums or {}).items()}
        self.load_kwargs = load_kwargs or {}
        self.label = label
        self._tree: Optional[MergerTree] = None
        self._epochs: Dict[float, Epoch] = {}

    # ------------------------------------------------------------------

    @property
    def tree(self) -> Optional[MergerTree]:
        """
        The Simulation's MergerTree, built lazily from the `trees` spec
        and cached.

        Returns
        -------
        MergerTree or None
            None if no `trees` spec was given.
        """
        spec = self._specs.get("trees")
        if spec is None:
            return None
        if self._tree is None:
            if isinstance(spec, MergerTree):
                self._tree = spec
            else:
                from . import load
                kwargs = dict(self.load_kwargs.get("trees", {}))
                if "trees" in self.fileformats:
                    kwargs["fileformat"] = self.fileformats["trees"]
                self._tree = load(spec, kind="tree", **kwargs)
        return self._tree

    def _snapnum_at(self, redshift: float) -> Optional[int]:
        if not self.snapnums:
            return None
        zs = np.array(sorted(self.snapnums))
        return self.snapnums[float(zs[np.argmin(np.abs(zs - redshift))])]

    def _kind_name(self, component: str) -> str:
        return {"snapshots": "snapshot", "halos": "halos",
                "galaxies": "galaxies"}[component]

    def _resolve(self, component: str, redshift: float):
        """Build the Dataset for one component at (nearest to) redshift."""
        spec = self._specs.get(component)
        if spec is None:
            return None
        if isinstance(spec, Dataset):
            return spec
        if hasattr(spec, "at") and not isinstance(spec, (str, dict)):
            return spec.at(redshift)              # EpochModel

        if isinstance(spec, dict):
            zs = np.array(sorted(float(z) for z in spec))
            znear = float(zs[np.argmin(np.abs(zs - redshift))])
            path = spec[znear] if znear in spec else spec[
                [z for z in spec if float(z) == znear][0]]
        else:
            path = spec

        from . import load
        kwargs = dict(self.load_kwargs.get(component, {}))
        if component in self.fileformats:
            kwargs["fileformat"] = self.fileformats[component]
        if component == "halos":
            kwargs.setdefault("snapnum", self._snapnum_at(redshift))
        ds = load(path, kind=self._kind_name(component), **kwargs)
        if ds.label and self.label != "simulation":
            ds.label = f"{self.label}:{ds.label}"
        return ds

    # ------------------------------------------------------------------

    def at(self, redshift: float) -> Epoch:
        """Epoch view at (nearest available to) `redshift` (cached)."""
        z = float(redshift)
        if z not in self._epochs:
            self._epochs[z] = Epoch(self, z)
        return self._epochs[z]

    def track_of(self, halo_id: Optional[int] = None, *,
                 index: Optional[int] = None, redshift: float = 0.0,
                 object_type: str = "Group") -> TrackDataset:
        """Convenience: sim.at(redshift).track_of(...)."""
        return self.at(redshift).track_of(halo_id, index=index,
                                          object_type=object_type)

    def summary(self) -> None:
        """
        Print each component's spec (path, epochs, or object) and the
        registered snapnums, if any.
        """
        print(f"Simulation '{self.label}'")
        for kind in ("snapshots", "halos", "trees", "galaxies"):
            spec = self._specs.get(kind)
            if spec is None:
                desc = "-"
            elif isinstance(spec, dict):
                desc = f"{len(spec)} epochs " \
                       f"(z = {', '.join(f'{float(z):g}' for z in sorted(spec))})"
            else:
                desc = getattr(spec, "label", None) or str(spec)
            print(f"  {kind:10s}: {desc}")
        if self.snapnums:
            print(f"  snapnums  : { {f'{z:g}': n for z, n in sorted(self.snapnums.items())} }")

    def __repr__(self) -> str:
        have = [k for k, v in self._specs.items() if v is not None]
        return f"<Simulation '{self.label}' [{', '.join(have)}]>"
