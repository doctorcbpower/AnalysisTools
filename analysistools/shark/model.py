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
import h5py

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

#: Fields available from star_formation_histories.hdf5.
#: Maps logical name -> (hdf5_group, dataset).
#: delta_t and lbt_mean are snapshot-level scalars, not per-galaxy arrays;
#: they are stored separately in _sfh_meta rather than _cache.
#:
#: SHARK tracks bulge growth as two separate channels -- disk-instability
#: driven (bulges_diskins) and merger-driven (bulges_mergers) -- mirroring
#: GALAXY_FIELDS' own mstars_burst_diskinstabilities/mstars_burst_mergers
#: split. mstars_bulge (used for StellarMass elsewhere) is already their
#: sum, so the *_diskins-only entries below are individually incomplete;
#: sfh_bulge()/sfh_metals_bulge() sum both channels -- use those, not
#: these raw keys directly, unless you specifically want one channel.
SFH_FIELDS: Dict[str, Tuple[str, str]] = {
    "sfh_disk":                ("disks", "star_formation_rate_histories"),
    "sfh_bulge_diskins":       ("bulges_diskins", "star_formation_rate_histories"),
    "sfh_bulge_mergers":       ("bulges_mergers", "star_formation_rate_histories"),
    "sfh_metals_disk":         ("disks", "metallicity_histories"),
    "sfh_metals_bulge_diskins": ("bulges_diskins", "metallicity_histories"),
    "sfh_metals_bulge_mergers": ("bulges_mergers", "metallicity_histories"),
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
        # SFH time arrays keyed by snapshot: {"delta_t": ndarray, "lbt_mean": ndarray}
        self._sfh_meta: Dict[int, Dict[str, np.ndarray]] = {}

    # ------------------------------------------------------------------
    # Core access
    # ------------------------------------------------------------------

    def get(self, field: str, redshift: float) -> np.ndarray:
        """Return the array for *field* at the snapshot nearest *redshift*.

        Covers both ``GALAXY_FIELDS`` (from galaxies.hdf5) and
        ``SFH_FIELDS`` (from star_formation_histories.hdf5).

        The result is cached; subsequent calls for the same
        (field, redshift) are free.

        Parameters
        ----------
        field : str
            Logical field name from ``GALAXY_FIELDS`` or ``SFH_FIELDS``.
        redshift : float
            Target redshift; matched to nearest snapshot.

        Returns
        -------
        arr : ndarray
            Shape (n_galaxies,) for galaxy fields.
            Shape (n_galaxies, n_sfh_bins) for SFH fields.
        """
        all_fields = {**GALAXY_FIELDS, **SFH_FIELDS}
        if field not in all_fields:
            raise KeyError(
                f"Unknown field '{field}'. "
                f"Available galaxy fields: {sorted(GALAXY_FIELDS)}. "
                f"Available SFH fields: {sorted(SFH_FIELDS)}."
            )
        snapshot = int(self.redshift_table[redshift])
        key      = (snapshot, field)
        if key not in self._cache:
            if field in SFH_FIELDS:
                self._load_sfh_snapshot(snapshot)
            else:
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

        Covers both galaxy fields and SFH fields.  The two file types are
        loaded independently so requesting only galaxy fields does not
        trigger an SFH read, and vice versa.

        Returns
        -------
        self : SharkModel   (for chaining)
        """
        galaxy_fields = [f for f in fields if f in GALAXY_FIELDS]
        sfh_fields    = [f for f in fields if f in SFH_FIELDS]
        unknown       = [f for f in fields if f not in GALAXY_FIELDS and f not in SFH_FIELDS]
        if unknown:
            raise KeyError(f"Unknown fields in preload: {unknown}")

        for z in redshifts:
            snapshot = int(self.redshift_table[z])
            if galaxy_fields and not self._snapshot_fully_cached(snapshot, galaxy_fields):
                self._load_snapshot(snapshot)
            if sfh_fields and not self._snapshot_fully_cached(snapshot, sfh_fields):
                self._load_sfh_snapshot(snapshot)
        return self

    def clear_cache(self) -> None:
        """Release all cached arrays, freeing memory."""
        self._cache.clear()
        self._meta.clear()
        self._sfh_meta.clear()

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

    def age_at_z(self, z: float) -> float:
        """
        Age of the Universe [Gyr] at redshift *z*, using this model's
        cosmology (H0, OmegaM, OmegaB read from the HDF5 cosmology group).

        Builds an astropy FlatLambdaCDM lazily from cached meta values —
        triggers an HDF5 read on first call for whichever snapshot is
        nearest *z* if not already cached.
        """
        from astropy.cosmology import FlatLambdaCDM
        meta = self.get_meta(z)
        cosmo = FlatLambdaCDM(
            H0=meta["h0"] * 100.0, Om0=meta["omega_m"], Ob0=meta["omega_b"]
        )
        return float(cosmo.age(z).value)

    def lookback_at_z(self, z: float) -> float:
        """Lookback time [Gyr] to redshift *z*, using this model's cosmology."""
        from astropy.cosmology import FlatLambdaCDM
        meta = self.get_meta(z)
        cosmo = FlatLambdaCDM(
            H0=meta["h0"] * 100.0, Om0=meta["omega_m"], Ob0=meta["omega_b"]
        )
        return float(cosmo.lookback_time(z).value)

    # ------------------------------------------------------------------
    # SFH-specific accessors
    # ------------------------------------------------------------------

    def get_sfh_meta(self, redshift: float) -> Dict[str, np.ndarray]:
        """Return the SFH time arrays for the snapshot nearest *redshift*.

        Keys
        ----
        ``delta_t``  : ndarray, shape (n_sfh_bins,)
            Width of each SFH time bin in Myr.
        ``lbt_mean`` : ndarray, shape (n_sfh_bins,)
            Mean lookback time of each bin in Gyr.

        These are snapshot-level scalars shared by all galaxies and are
        read from star_formation_histories.hdf5 on first access.
        """
        snapshot = int(self.redshift_table[redshift])
        if snapshot not in self._sfh_meta:
            self._load_sfh_snapshot(snapshot)
        return self._sfh_meta[snapshot]

    def sfh_disk(self, redshift: float) -> np.ndarray:
        """Disk SFH array  [M_sun / yr], shape (n_galaxies, n_sfh_bins)."""
        return self.get("sfh_disk", redshift)

    def sfh_bulge(self, redshift: float) -> np.ndarray:
        """Bulge SFH array  [M_sun / yr], shape (n_galaxies, n_sfh_bins) --
        sum of the disk-instability-driven and merger-driven channels
        (see SFH_FIELDS' docstring for why both are needed: mstars_bulge,
        used for StellarMass elsewhere, is already their sum, so using
        only bulges_diskins here would make the SFH-integrated formed
        mass look smaller than the current bulge mass for any galaxy
        with real merger-driven bulge growth)."""
        return (self.get("sfh_bulge_diskins", redshift)
               + self.get("sfh_bulge_mergers", redshift))

    def sfh_metals_disk(self, redshift: float) -> np.ndarray:
        """Disk stellar metallicity history  [M_sun], shape (n_galaxies, n_sfh_bins)."""
        return self.get("sfh_metals_disk", redshift)

    def sfh_metals_bulge(self, redshift: float) -> np.ndarray:
        """Bulge stellar metal mass history  [M_sun], shape (n_galaxies,
        n_sfh_bins) -- sum of the disk-instability-driven and
        merger-driven channels, see ``sfh_bulge``."""
        return (self.get("sfh_metals_bulge_diskins", redshift)
               + self.get("sfh_metals_bulge_mergers", redshift))

    def Z_disk_history(self, redshift: float) -> np.ndarray:
        """
        Mass-weighted metallicity of disk stars per SFH bin.

        Derived as mz_disk / mstar_disk per bin, clipped to [1e-4, 0.03]
        to stay within FSPS SSP grid bounds.  Where the stellar mass in a
        bin is zero the metallicity falls back to solar (0.02).

        Returns
        -------
        Z : ndarray, shape (n_galaxies, n_sfh_bins)
        """
        return _sfh_metallicity(self.sfh_metals_disk(redshift),
                                self.sfh_disk(redshift))

    def Z_bulge_history(self, redshift: float) -> np.ndarray:
        """
        Mass-weighted metallicity of bulge stars per SFH bin.
        See ``Z_disk_history`` for details.

        Returns
        -------
        Z : ndarray, shape (n_galaxies, n_sfh_bins)
        """
        return _sfh_metallicity(self.sfh_metals_bulge(redshift),
                                self.sfh_bulge(redshift))

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
                h0      = float(f["cosmology/h"][()])
                omega_m = float(f["cosmology/omega_m"][()])
                omega_b = float(f["cosmology/omega_b"][()])
                vol_h   = float(f["run_info/effective_volume"][()]) * len(self.subvols)
            self._meta[snapshot] = {
                "h0": h0, "vol_h": vol_h, "vol": vol_h / h0**3,
                "omega_m": omega_m, "omega_b": omega_b,
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

        # OmegaM/OmegaB aren't returned by common.read_data; one cheap
        # direct read of two scalars from the file already inspected above.
        with h5py.File(fname, "r") as f:
            omega_m = float(f["cosmology/omega_m"][()])
            omega_b = float(f["cosmology/omega_b"][()])

        self._meta[snapshot] = {
            "h0":      float(h0),
            "vol_h":   float(vol_h),
            "vol":     float(vol_h) / float(h0)**3,
            "omega_m": omega_m,
            "omega_b": omega_b,
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

    def _load_sfh_snapshot(self, snapshot: int) -> None:
        """Read all SFH_FIELDS for *snapshot* from star_formation_histories.hdf5.

        Uses ``common.read_sfh``, which already handles subvolume
        concatenation.  Its ``fields`` parameter is a plain ``{group: dataset}``
        dict (one dataset per group per call) — since SFH_FIELDS has two
        datasets per group (disk/burst under both star_formation_history and
        metallicity_history), each logical field is fetched with its own
        call rather than batched, to avoid dict-key collisions overwriting
        one of the two entries sharing a group.

        delta_t and lbt_mean are stored in ``_sfh_meta`` (read once, on the
        first successful field call); per-galaxy arrays go into ``_cache``.

        Fields absent from the HDF5 file are stored as None and will raise
        a clear AttributeError on access.

        **Row alignment**: ``star_formation_histories.hdf5`` commonly
        covers only a *subset* of the galaxies in ``galaxies.hdf5`` for
        the same snapshot (e.g. centrals only) — its arrays' row `i` is
        *not* the same galaxy as ``galaxies.hdf5``'s row `i`. Confirmed on
        a real run: one subvolume had 124 SFH-file rows against >28,000
        rows in the corresponding galaxies.hdf5, and naively indexing the
        SFH array with a galaxies.hdf5 row either raised ``IndexError``
        or (for smaller indices) silently returned a *different* galaxy's
        SFH entirely. Both files carry a ``galaxies/id_galaxy`` dataset;
        every SFH array here is reindexed onto ``galaxies.hdf5``'s own row
        order via that ID (NaN rows for galaxies absent from the SFH
        file), so ``model.sfh_disk(z)[row]`` uses the same `row`
        convention as every ``GALAXY_FIELDS`` accessor. If either file is
        missing ``id_galaxy`` (older SHARK output), there is no safe way
        to establish this correspondence -- every SFH field is left
        unavailable (``None``) rather than risk silently mismatched rows.
        """
        first_subvol = sorted(self.subvols)[0]

        fname = os.path.join(self.model_dir, str(snapshot), str(first_subvol), "star_formation_histories.hdf5")

        with h5py.File(fname, "r") as f:
            self._sfh_meta[snapshot] = {
                "delta_t": f["delta_t"][()],
                "lbt_mean": f["lbt_mean"][()],
            }
            available = set(f.keys())
            sfh_has_id_galaxy = "galaxies" in f and "id_galaxy" in f["galaxies"]

        row_map = None
        if sfh_has_id_galaxy:
            row_map = self._sfh_row_map(snapshot)
        if row_map is None:
            print(
                f"  [{self.label}] snapshot {snapshot}: could not "
                f"establish a galaxies.hdf5 <-> star_formation_histories"
                f".hdf5 row correspondence (missing/unreadable "
                f"'galaxies/id_galaxy' in one or both files) -- every SFH "
                f"field will be unavailable rather than risk silently "
                f"mismatched rows."
            )
            for logical in SFH_FIELDS:
                self._cache[(snapshot, logical)] = None
            return

        for logical, (grp, ds) in SFH_FIELDS.items():
            # No galaxies/components stored in this SFH file
            if grp not in available:
                print(
                    f"  [{self.label}] snapshot {snapshot}: "
                    f"missing SFH group '{grp}'. "
                    f"Field '{logical}' will be unavailable."
                )
                self._cache[(snapshot, logical)] = None
                continue
            try:
                arrays, _, _ = common.read_sfh(
                    self.model_dir, snapshot, {grp: ds}, self.subvols
                )
                raw = arrays[0]
            except (KeyError, OSError) as exc:
                print(
                    f"  [{self.label}] snapshot {snapshot}: "
                    f"could not read '{grp}/{ds}' from "
                    f"star_formation_histories.hdf5: {exc}. "
                    f"Field '{logical}' will be unavailable."
                )
                self._cache[(snapshot, logical)] = None
                continue

            self._cache[(snapshot, logical)] = self._reindex_sfh_array(
                raw, row_map)

    def _sfh_row_map(self, snapshot: int) -> Optional[np.ndarray]:
        """``galaxies.hdf5`` row -> ``star_formation_histories.hdf5`` row,
        joined via each file's own ``galaxies/id_galaxy``. ``-1`` for a
        galaxy with no SFH entry. ``None`` if either side's ID array
        couldn't be read (caller then disables SFH entirely for this
        snapshot -- see ``_load_sfh_snapshot``)."""
        try:
            sfh_ids_list, _, _ = common.read_sfh(
                self.model_dir, snapshot, {"galaxies": "id_galaxy"},
                self.subvols)
            full_ids_list = common.read_data(
                self.model_dir, snapshot, {"galaxies": ("id_galaxy",)},
                self.subvols, include_h0_volh=False)
        except (KeyError, OSError) as exc:
            print(f"  [{self.label}] snapshot {snapshot}: could not read "
                 f"id_galaxy for SFH row mapping: {exc}.")
            return None
        if not sfh_ids_list or not full_ids_list:
            return None

        sfh_ids = np.asarray(sfh_ids_list[0])
        full_ids = np.asarray(full_ids_list[0])
        id_to_sfh_row = {int(gid): i for i, gid in enumerate(sfh_ids)}
        return np.fromiter(
            (id_to_sfh_row.get(int(gid), -1) for gid in full_ids),
            dtype=np.int64, count=len(full_ids))

    @staticmethod
    def _reindex_sfh_array(raw: np.ndarray, row_map: np.ndarray) -> np.ndarray:
        """Reindex an (n_sfh_galaxies, n_bins) array onto
        (len(row_map), n_bins) using `row_map` (see `_sfh_row_map`), NaN
        for entries with no corresponding SFH row."""
        n_bins = raw.shape[1] if raw.ndim > 1 else 1
        out = np.full((len(row_map),) + raw.shape[1:], np.nan)
        valid = row_map >= 0
        out[valid] = raw[row_map[valid]]
        return out


# ---------------------------------------------------------------------------
# Module-level helpers
# ---------------------------------------------------------------------------

def _sfh_metallicity(
    mz_history: np.ndarray,
    sfr_history: np.ndarray,
    z_floor: float = 1e-4,
    z_ceil:  float = 0.03,
) -> np.ndarray:
    """Derive metallicity Z(t) from metal-mass and SFR histories.

    Parameters
    ----------
    mz_history  : (n_gal, n_bins)  cumulative metal stellar mass per bin [M_sun]
    sfr_history : (n_gal, n_bins)  SFR per bin [M_sun / yr]
    z_floor, z_ceil : float
        Clip range matching FSPS SSP grid bounds.

    Returns
    -------
    Z : (n_gal, n_bins), clipped to [z_floor, z_ceil].
    """
    # Stellar mass formed per bin ∝ SFR (delta_t is common to all galaxies
    # so the ratio mz/mstar is proportional to mz/sfr; we use sfr as a proxy
    # for the mass weight within each bin).
    with np.errstate(divide="ignore", invalid="ignore"):
        Z = np.where(sfr_history > 0, mz_history / sfr_history, 0.02)
    return np.clip(Z, z_floor, z_ceil).astype(np.float64)
