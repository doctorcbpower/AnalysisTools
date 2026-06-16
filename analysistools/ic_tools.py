"""
ic_tools.py — Tools for preparing initial conditions from simulation components.
Extends SnapshotTools with merge, taper, and mass-rescaling operations.
"""

import numpy as np
from typing import Optional, List
from .snapshot_tools import SnapshotTools

class ICTools(SnapshotTools):
    """
    Extends SnapshotTools with operations for building initial conditions:
    merging components, tapering the DM halo, rescaling component masses,
    and excising particle types.

    Usage
    -----
    ic = ICTools(snapfileformat="HDF5")
    ic.read_snapshot("dm_halo.hdf5")

    gas = ICTools(snapfileformat="HDF5")
    gas.read_snapshot("gas_disc.hdf5")

    combined = ic.merge([gas])
    combined.taper_halo(mass_fraction=0.95)
    combined.rescale_component_mass(ptype=combined.dm_type, target_mass=dm_target)
    combined.write_snapshot("combined.hdf5", ...)
    """

    def merge(self, others: List["ICTools"], cosmology: Optional[dict] = None, ) -> "ICTools":
        """
        Merge this snapshot with one or more others.

        Parameters
        ----------
        others : list of SnapshotTools
            Snapshots to merge into this one. All must have data loaded.

        Returns
        -------
        SnapshotTools
            A new SnapshotTools instance containing all particles.

        Notes
        -----
        - IDs are resequenced from 1 to npart to avoid collisions between
          independently generated snapshots.
        - massarr is set to 0 for any type where input snapshots disagree
          on the per-type mass, forcing a per-particle mass block.
        - Gas-only arrays (u, rho, hsml) are only merged if any snapshot
          has gas particles (ptype == gas_type).
        """
        snaps = [self] + list(others)

        # Validate all have data loaded
        for i, s in enumerate(snaps):
            if not hasattr(s, 'pos'):
                raise RuntimeError(f"Snapshot {i} has no data loaded — call read_snapshot() first.")

        merged = ICTools(snapfileformat=self.snapfileformat)

        # --- core arrays ---
        merged.pos  = np.concatenate([s.pos  for s in snaps], axis=0)
        merged.vel  = np.concatenate([s.vel  for s in snaps], axis=0)
        merged.mass = np.concatenate([s.mass for s in snaps], axis=0)

        # --- ptype: derive ngas from this ---
        if all(hasattr(s, 'ptype') for s in snaps):
            merged.ptype = np.concatenate([s.ptype for s in snaps], axis=0)
        else:
            # reconstruct ptype from num_part_total if not stored explicitly
            ptypes = []
            for s in snaps:
                for t, n in enumerate(s.num_part_total):
                    ptypes.append(np.full(n, t, dtype=np.int32))
            merged.ptype = np.concatenate(ptypes)

        # --- resequence IDs to avoid collisions ---
        npart = len(merged.pos)
        merged.pids = np.arange(1, npart + 1, dtype=np.int32)

        # --- gas-only arrays ---
        gas_mask = merged.ptype == merged.gas_type
        ngas = int(np.sum(gas_mask))
        if ngas > 0:
            # Gas-only arrays — only include snaps that actually have gas
            for attr in ('u', 'rho', 'hsml', 'smoothinglength'):
                sources = [
                    getattr(s, attr) for s in snaps
                    if hasattr(s, attr) and int(s.num_part_total[s.gas_type]) > 0
                ]
                if sources:
                    # Normalise to hsml regardless of which name the reader used
                    target_attr = 'hsml' if attr == 'smoothinglength' else attr
                    setattr(merged, target_attr, np.concatenate(sources))
        # --- num_part_total ---
        merged.num_part_total = np.zeros(self.num_part_type, dtype=np.int64)
        for t in range(self.num_part_type):
            merged.num_part_total[t] = int(np.sum(merged.ptype == t))

        # --- massarr: only set if all snaps agree for that type ---
        merged.mass_table = np.zeros(self.num_part_type, dtype=np.float64)
        for t in range(self.num_part_type):
            vals = [s.mass_table[t] for s in snaps
                    if hasattr(s, 'mass_table') and s.num_part_total[t] > 0]
            if vals and len(set(vals)) == 1 and vals[0] > 0.0:
                merged.mass_table[t] = vals[0]
            # else leave as 0 — per-particle mass block will be used

        # --- cosmology: take from self, warn if snapshots differ ---
        if cosmology is not None:
            for attr in ('scale_factor', 'box_size', 'omega_0',
                         'omega_lambda', 'hubble_param'):
                vals = [getattr(s, attr) for s in snaps if hasattr(s, attr)]
                if vals:
                    if len(set(vals)) > 1:
                        self.logger.warning(
                            f"Snapshots differ in '{attr}': {vals}. Using value from first snapshot."
                        )
                    setattr(merged, attr, vals[0])

        self.logger.info(
            f"Merged {len(snaps)} snapshots: {npart} particles total, "
            f"{ngas} gas particles."
        )
        return merged

    def taper_halo(self, mass_fraction: float, ptype: Optional[int] = None) -> "ICTools":
        """
        Clip the DM halo by retaining only particles within the radius
        enclosing mass_fraction of the total DM mass, with an exponential
        taper beyond that radius. All other particle types are untouched.

        Parameters
        ----------
        mass_fraction : float
            Fraction of total DM mass defining the taper radius (0 < f < 1).
        ptype : int, optional
            Particle type to taper. Defaults to self.dm_type.

        Returns
        -------
        self : SnapshotTools (modified in place)
        """
        if not 0.0 < mass_fraction < 1.0:
            raise ValueError("mass_fraction must be between 0 and 1")

        ptype = ptype if ptype is not None else self.dm_type
        dm_mask = self.ptype == ptype

        if not np.any(dm_mask):
            raise ValueError(f"No particles of type {ptype} found")

        r2_dm = np.sum(self.pos[dm_mask] ** 2, axis=1)
        total_dm_mass = np.sum(self.mass[dm_mask])

        # Find taper radius from DM particles only
        isort = np.argsort(r2_dm)
        cumulative = np.cumsum(self.mass[dm_mask][isort])
        itaper = np.searchsorted(cumulative, mass_fraction * total_dm_mass)
        itaper = np.clip(itaper, 1, len(r2_dm) - 1)
        r_taper = np.sqrt(r2_dm[isort[itaper]])

        # Exponential taper beyond r_taper, applied to DM only
        r_dm = np.sqrt(r2_dm)
        keep_prob = np.where(r_dm <= r_taper, 1.0, np.exp(-(r_dm - r_taper) / r_taper))

        rng = np.random.default_rng(seed=187619)
        keep_dm = rng.random(int(np.sum(dm_mask))) < keep_prob

        # Build full mask: keep all non-DM, apply taper only to DM
        full_mask = np.ones(len(self.ptype), dtype=bool)
        full_mask[dm_mask] = keep_dm

        self.logger.info(
            f"Taper: ptype={ptype}, r_taper={r_taper:.3f}, "
            f"retained {np.sum(keep_dm)}/{int(np.sum(dm_mask))} DM particles, "
            f"{np.sum(~dm_mask)} non-DM particles untouched"
        )
        return self._apply_mask(full_mask)

    def rescale_component_mass(self, ptype: int, target_mass: float) -> "ICTools":
        """
        Rescale the mass of a particle type to a target total mass,
        preserving the mass ratios between individual particles of that type.
        Used to conserve total system mass when combining components.

        Parameters
        ----------
        ptype : int
            Particle type to rescale (e.g. gas_type, star_type).
        target_mass : float
            Desired total mass for this component (same units as mass array).

        Returns
        -------
        self : SnapshotTools (modified in place)
        """
        mask = self.ptype == ptype
        current_mass = np.sum(self.mass[mask])

        if current_mass <= 0.0:
            raise ValueError(f"No particles of type {ptype} or zero total mass")

        scale = target_mass / current_mass
        self.mass[mask] *= scale

        # Keep mass_table consistent if this type uses a uniform mass
        if hasattr(self, 'mass_table') and self.mass_table[ptype] > 0.0:
            self.mass_table[ptype] *= scale

        self.logger.info(
            f"Rescaled ptype {ptype}: {current_mass:.4e} -> {target_mass:.4e} "
            f"(factor {scale:.6f})"
        )
        return self


    def excise_type(self, ptype: int) -> "ICTools":
        """
        Remove all particles of a given type.

        Parameters
        ----------
        ptype : int
            Particle type to remove.

        Returns
        -------
        self : SnapshotTools (modified in place)
        """
        mask = self.ptype != ptype
        n_removed = int(np.sum(~mask))
        self.logger.info(f"Excising ptype {ptype}: removing {n_removed} particles")
        return self._apply_mask(mask)


    def _apply_mask(self, mask: np.ndarray) -> "ICTools":
        """Apply a boolean mask to all particle arrays in place."""
        # Identify gas particles before masking
        gas_before = self.ptype == self.gas_type

        for attr in ('pos', 'vel', 'mass', 'pids', 'ptype', 'potential'):
            if hasattr(self, attr):
                setattr(self, attr, getattr(self, attr)[mask])

        # Gas-only arrays need a gas-specific submask
        gas_mask = mask[gas_before]
        for attr in ('u', 'rho', 'hsml'):
            if hasattr(self, attr):
                setattr(self, attr, getattr(self, attr)[gas_mask])

        # Update counts
        self.num_part_total = np.zeros(self.num_part_type, dtype=np.int64)
        for t in range(self.num_part_type):
            self.num_part_total[t] = int(np.sum(self.ptype == t))

        return self
