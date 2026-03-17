"""
shark.model
-----------
SharkModel: lazy, cached access to SHARK HDF5 galaxy catalogues.

Responsibilities
~~~~~~~~~~~~~~~~
* Read raw HDF5 fields on demand for any snapshot / subvolume set.
* Cache loaded arrays so repeated access costs nothing.
* Expose derived quantities (stellar mass, sSFR, …) as properties.

It deliberately knows nothing about binning, statistics, or plotting —
those concerns live in shark.analysis and shark.plots.
"""

from __future__ import annotations

import os
from typing import Any, Dict, Optional, Sequence, Set, Tuple

import numpy as np

from . import common
from .common import _redshift_table, parse_subvolumes

# ---------------------------------------------------------------------------
# Field catalogue
# ---------------------------------------------------------------------------
# Maps a logical group name to the HDF5 (group, dataset) pairs needed
# Analysis code requests fields by logical name; SharkModel fetches only
# what is actually needed.

#: Every galaxy field available from the standard galaxies.hdf5 catalogue.
GALAXY_FIELDS: Dict[str, Tuple[str, str]] = {
    # Stellar masses
    "mstars_disk":      ("galaxies", "mstars_disk"),
    "mstars_bulge":     ("galaxies", "mstars_bulge"),
    "mstars_burst_mergers":   ("galaxies", "mstars_burst_mergers"),
    "mstars_burst_diskinstabilities": ("galaxies", "mstars_burst_diskinstabilities"),
    # Stellar metal masses
    "mstars_metals_disk":  ("galaxies", "mstars_metals_disk"),
    "mstars_metals_bulge": ("galaxies", "mstars_metals_bulge"),
    # Cold gas
    "mgas_disk":        ("galaxies", "mgas_disk"),
    "mgas_bulge":       ("galaxies", "mgas_bulge"),
    "matom_disk":       ("galaxies", "matom_disk"),     
    "matom_bulge":      ("galaxies", "matom_bulge"),     
    "mmol_disk":        ("galaxies", "mmol_disk"),      
    "mmol_bulge":       ("galaxies", "mmol_bulge"),    
    # ISM metal masses
    "mgas_metals_disk": ("galaxies", "mgas_metals_disk"),
    "mgas_metals_bulge":("galaxies", "mgas_metals_bulge"),
    # Hot gas
    "mhot":             ("galaxies", "mhot"),
    "mhot_metals":      ("galaxies", "mhot_metals"),
    "mhot_stripped":    ("galaxies", "mhot_stripped"),    
    # Ejected / reheated gas  
    "mreheated":        ("galaxies", "mreheated"),
    "mreheated_metals": ("galaxies", "mreheated_metals"),
    "mlost":            ("galaxies", "mlost"),         
    "mlost_metals":     ("galaxies", "mlost_metals"),
    # Ram-pressure stripped ISM  
    "mism_stripped":          ("galaxies", "mism_stripped"),
    "mism_metals_stripped":   ("galaxies", "mism_metals_stripped"),
    # Halo / subhalo masses
    "mvir_hosthalo":    ("galaxies", "mvir_hosthalo"),
    "mvir_subhalo":     ("galaxies", "mvir_subhalo"),
    "vmax_subhalo":     ("galaxies", "vmax_subhalo"),     
    "cnfw_subhalo":     ("galaxies", "cnfw_subhalo"),     
    # Sizes
    "rstar_disk":       ("galaxies", "rstar_disk"),
    "rstar_bulge":      ("galaxies", "rstar_bulge"),
    "rgas_disk":        ("galaxies", "rgas_disk"),
    "rgas_bulge":       ("galaxies", "rgas_bulge"),
    # Star formation
    "sfr_disk":         ("galaxies", "sfr_disk"),
    "sfr_bulge":        ("galaxies", "sfr_burst"),
    "sfr_burst_mergers":   ("galaxies", "sfr_burst_mergers"),
    "sfr_burst_diskins":   ("galaxies", "sfr_burst_diskins"),
    # Black holes
    "mbh":              ("galaxies", "m_bh"),
    "mbh_from_mergers": ("galaxies", "m_bh_assembly"),
    "bh_acc_hothalo":   ("galaxies", "bh_accretion_rate_hh"),
    "bh_acc_starburst": ("galaxies", "bh_accretion_rate_sb"),
    # Identifiers
    "id_galaxy":        ("galaxies", "id_galaxy"),
    "id_halo":          ("galaxies", "id_halo"),
    "descendant_id_galaxy":   ("galaxies", "descendant_id_galaxy"),
    # Positions [cMpc/h] and velocities [km/s]
    "position_x":       ("galaxies", "position_x"),
    "position_y":       ("galaxies", "position_y"),
    "position_z":       ("galaxies", "position_z"),
    "velocity_x":       ("galaxies", "velocity_x"),
    "velocity_y":       ("galaxies", "velocity_y"),
    "velocity_z":       ("galaxies", "velocity_z"),
    # Type flag  (0=central, 1=satellite, 2=orphan)
    "type":             ("galaxies", "type"),
}


# ---------------------------------------------------------------------------
# SharkModel
# ---------------------------------------------------------------------------

class SharkModel:
    """Lazy, cached interface to a single SHARK model run.

    Parameters
    ----------
    model_dir : str
        Path to the model output directory (parent of per-snapshot dirs).
    redshift_table : common._redshift_table
        Redshift → snapshot lookup built by ``common._redshift_table``.
    subvols : set or sequence of int
        Sub-volume identifiers to include when reading HDF5 files.
    label : str, optional
        Human-readable name (used in plot legends).
        Defaults to the last path component of *model_dir*.
    colour : str, optional
        Matplotlib colour string.  If omitted the plotter assigns one.

    Examples
    --------
    >>> from shark.common import _redshift_table, parse_subvolumes
    >>> from shark.model import SharkModel
    >>> rt  = _redshift_table("redshift_list.txt")
    >>> sv  = parse_subvolumes("0-63")
    >>> m   = SharkModel("./CDM/base_model", rt, sv, label="CDM")
    >>> mstar = m.get("mstars_disk", redshift=0.0) + m.get("mstars_bulge", redshift=0.0)
    """

    def __init__(
        self,
        model_dir:      str,
        redshift_table: _redshift_table,
        subvols:        Sequence[int],
        label:          Optional[str] = None,
        colour:         Optional[str] = None,
    ):
        self.model_dir      = model_dir
        self.redshift_table = redshift_table
        self.subvols        = set(subvols)
        self.label          = label or os.path.basename(model_dir.rstrip("/\\"))
        self.colour         = colour

        # Cache: (snapshot, field_name) -> ndarray
        self._cache: Dict[Tuple[int, str], np.ndarray] = {}
        # Cosmological / volume info keyed by snapshot
        self._meta:  Dict[int, Dict[str, float]]       = {}

    # ------------------------------------------------------------------
    # Core access
    # ------------------------------------------------------------------

    def get(self, field: str, redshift: float) -> np.ndarray:
        """Return the array for *field* at the snapshot nearest *redshift*.

        The result is cached; subsequent calls for the same
        (field, redshift) are free.

        Parameters
        ----------
        field : str
            Logical field name from ``GALAXY_FIELDS`` (e.g. ``"mstars_disk"``).
        redshift : float
            Target redshift; matched to nearest snapshot.

        Returns
        -------
        arr : ndarray, shape (n_galaxies,)
        """
        if field not in GALAXY_FIELDS:
            raise KeyError(
                f"Unknown field '{field}'. "
                f"Available: {sorted(GALAXY_FIELDS)}"
            )
        snapshot = int(self.redshift_table[redshift])
        key      = (snapshot, field)
        if key not in self._cache:
            self._load_snapshot(snapshot)
        arr = self._cache[key]
        if arr is None:
            raise AttributeError(
                f"Field '{field}' was not found in the HDF5 files for "
                f"snapshot {snapshot} of model '{self.label}'. "
                f"Check that this field is written by your SHARK configuration."
            )
        return arr

    def get_meta(self, redshift: float) -> Dict[str, float]:
        """Return cosmological / volume metadata for the snapshot nearest *redshift*.

        Keys: ``h0``, ``vol_h``  (volume in (Mpc/h)^3), ``vol`` (Mpc^3).
        """
        snapshot = int(self.redshift_table[redshift])
        if snapshot not in self._meta:
            self._load_snapshot(snapshot)
        return self._meta[snapshot]

    def preload(self, fields: Sequence[str], redshifts: Sequence[float]) -> "SharkModel":
        """Eagerly load *fields* for all *redshifts* into the cache.

        Useful before a long analysis loop to batch all HDF5 reads upfront.

        Returns
        -------
        self : SharkModel   (for chaining)
        """
        for z in redshifts:
            snapshot = int(self.redshift_table[z])
            if not self._snapshot_fully_cached(snapshot, fields):
                self._load_snapshot(snapshot)
        return self

    def clear_cache(self) -> None:
        """Release all cached arrays, freeing memory."""
        self._cache.clear()
        self._meta.clear()

    # ------------------------------------------------------------------
    # Derived quantities
    # ------------------------------------------------------------------

    def mstars(self, redshift: float) -> np.ndarray:
        """Total stellar mass = disk + bulge  [M_sun / h]."""
        return (
            self.get("mstars_disk",  redshift)
            + self.get("mstars_bulge", redshift)
        )

    def mgas(self, redshift: float) -> np.ndarray:
        """Total cold gas mass = disk + bulge  [M_sun / h]."""
        return (
            self.get("mgas_disk",  redshift)
            + self.get("mgas_bulge", redshift)
        )

    def sfr(self, redshift: float) -> np.ndarray:
        """Total star formation rate = disk + bulge  [M_sun / yr / h]."""
        return (
            self.get("sfr_disk",  redshift)
            + self.get("sfr_bulge", redshift)
        )

    def ssfr(self, redshift: float) -> np.ndarray:
        """Specific SFR = SFR / mstars  [yr^-1].  Zero where mstars == 0."""
        ms  = self.mstars(redshift)
        sfr = self.sfr(redshift)
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(ms > 0.0, sfr / ms, 0.0)

    def rstar(self, redshift: float) -> np.ndarray:
        """Half-stellar-mass radius.

        Returns disk size for disk-dominated galaxies
        (mstars_disk > mstars_bulge), bulge size otherwise  [kpc / h].
        """
        rd = self.get("rstar_disk",  redshift)
        rb = self.get("rstar_bulge", redshift)
        return np.where(
            self.get("mstars_disk", redshift) >= self.get("mstars_bulge", redshift),
            rd, rb
        )

    def mhalo(self, redshift: float) -> np.ndarray:
        """Host halo virial mass  [M_sun / h]."""
        return self.get("mvir_hosthalo", redshift)

    def msubhalo(self, redshift: float) -> np.ndarray:
        """Subhalo virial mass  [M_sun / h]."""
        return self.get("mvir_subhalo", redshift)

    def mbh(self, redshift: float) -> np.ndarray:
        """Central black hole mass  [M_sun / h]."""
        return self.get("mbh", redshift)

    def mbulge(self, redshift: float) -> np.ndarray:
        """Bulge stellar mass  [M_sun / h]."""
        return self.get("mstars_bulge", redshift)

    def galaxy_type(self, redshift: float) -> np.ndarray:
        """Galaxy type flag  (0=central, 1=satellite, 2=orphan)."""
        return self.get("type", redshift)

    def h0(self, redshift: float) -> float:
        """Hubble parameter h for the snapshot nearest *redshift*."""
        return self.get_meta(redshift)["h0"]

    def volume(self, redshift: float) -> float:
        """Comoving survey volume in Mpc^3 for the snapshot nearest *redshift*."""
        return self.get_meta(redshift)["vol"]

    # ------------------------------------------------------------------
    # Repr
    # ------------------------------------------------------------------

    def __repr__(self) -> str:
        n_cached = len({snap for snap, _ in self._cache})
        return (
            f"SharkModel(label={self.label!r}, "
            f"model_dir={self.model_dir!r}, "
            f"subvols={len(self.subvols)}, "
            f"snapshots_cached={n_cached})"
        )

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _load_snapshot(self, snapshot: int) -> None:
        """Read all available GALAXY_FIELDS for *snapshot* and populate the cache.

        Fields that are absent from the HDF5 file are skipped and stored as
        None in the cache.  Callers that request a missing field will receive
        a clear AttributeError rather than a raw h5py KeyError.
        """
        import os
        import h5py

        # Discover which fields actually exist in this snapshot by inspecting
        # the first subvolume file before attempting any reads.
        first_subvol = sorted(self.subvols)[0]
        fname = os.path.join(
            self.model_dir, str(snapshot), str(first_subvol), "galaxies.hdf5"
        )
        available: Dict[str, Set[str]] = {}   # hdf5_group -> {dataset, ...}
        missing_logical: list = []

        with h5py.File(fname, "r") as f:
            for logical, (grp, ds) in GALAXY_FIELDS.items():
                if grp in f and ds in f[grp]:
                    available.setdefault(grp, set()).add(ds)
                else:
                    missing_logical.append(logical)
                    self._cache[(snapshot, logical)] = None

        if missing_logical:
            print(
                f"  [{self.label}] snapshot {snapshot}: "
                f"fields not found in HDF5 and will be unavailable: "
                f"{missing_logical}"
            )

        if not available:
            # No readable fields — still need metadata
            with h5py.File(fname, "r") as f:
                h0    = float(f["cosmology/h"][()])
                vol_h = float(f["run_info/effective_volume"][()]) * len(self.subvols)
            self._meta[snapshot] = {
                "h0": h0, "vol_h": vol_h, "vol": vol_h / h0**3
            }
            return

        # Build fields dict with only what exists
        hdf5_fields = {grp: tuple(ds_set) for grp, ds_set in available.items()}

        raw = common.read_data(
            self.model_dir, snapshot, hdf5_fields, self.subvols
        )

        # common.read_data returns [h0, vol_h, field0, field1, ...]
        h0, vol_h = raw[0], raw[1]
        arrays    = raw[2:]

        self._meta[snapshot] = {
            "h0":    float(h0),
            "vol_h": float(vol_h),
            "vol":   float(vol_h) / float(h0)**3,
        }

        # Map returned arrays back to logical field names.
        # The order matches the iteration order of hdf5_fields / available,
        # which may differ from GALAXY_FIELDS order — rebuild the mapping
        # explicitly from (grp, ds) pairs.
        ds_to_logical = {(grp, ds): logical
                         for logical, (grp, ds) in GALAXY_FIELDS.items()}
        idx = 0
        for grp, ds_set in available.items():
            for ds in ds_set:
                logical = ds_to_logical[(grp, ds)]
                self._cache[(snapshot, logical)] = arrays[idx]
                idx += 1

    def _snapshot_fully_cached(
        self, snapshot: int, fields: Sequence[str]
    ) -> bool:
        return all((snapshot, f) in self._cache for f in fields)
