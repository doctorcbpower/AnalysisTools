#!/usr/bin/env python3
"""
analysistools.catalogue.backends
----------------------------------
GalaxyBackend: the one piece of the catalogue pipeline that genuinely
differs between projects, made pluggable so the schema and every other
stage stay simulation-agnostic.

Dorcha gets galaxy properties from SHARK, a semi-analytic model bolted onto
a DMO merger tree. A hydrodynamic rerun of the same halo (Arepo-Solas,
SWIFT) has stars *in the snapshot* -- there is no SAM step, but the target
fields (StellarMass, MetallicityStellar, StarFormationRate, ...) are the
same physical quantities, just measured a different way. A GalaxyBackend
takes an Epoch and a halo row and returns one row of
``Satellites/GalaxyProperties`` fields; ``pipeline.CatalogueBuilder`` never
knows or cares which backend produced them.

See docs/master_catalogue.md ("The one piece that genuinely differs
between projects") and DEVELOPMENT.md Phase 6.
"""
from __future__ import annotations

from typing import Any, Dict, Optional, Protocol

import logging

import numpy as np

logger = logging.getLogger(__name__)


class GalaxyBackend(Protocol):
    """Interface every galaxy-property source implements."""

    def galaxy_properties(self, epoch: "analysistools.api.simulation.Epoch",
                          halo_row: int) -> Dict[str, Any]:
        """Return one row of Satellites/GalaxyProperties fields (native
        field names, unconverted) for the halo at ``halo_row`` in
        ``epoch.halos``. Missing quantities should be omitted (not NaN-
        filled) so ``validation.IntegrityValidator`` can distinguish
        "not applicable" from "failed to compute"."""
        ...

    def star_formation_history(
            self, epoch: "analysistools.api.simulation.Epoch",
            halo_row: int, time_bin_edges: np.ndarray) -> Optional[np.ndarray]:
        """Star formation history for the halo at ``halo_row``, rebinned
        onto ``time_bin_edges`` (ascending lookback time in Gyr, length
        n_bins+1, e.g. ``Snapshots/time_bin_edges_sfh``), as an (n_bins,)
        Msun/yr array -- or ``None`` if not computable (no match, or no
        native SFH source available for this backend/catalogue). Used by
        ``derived.StarFormationHistoryStage``, which is also where the
        rebinned result becomes ``Satellites/GalaxyProperties/SFH`` and
        feeds ``MeanStellarAge``/``QuenchingTime``/``IsQuenched_z0``."""
        ...


def rebin_sfh(native_edges: np.ndarray, native_sfr: np.ndarray,
              output_edges: np.ndarray) -> np.ndarray:
    """Rebin a piecewise-constant SFR history (constant within each native
    bin) onto ``output_edges``, conserving total mass formed via
    fractional bin-overlap -- a mechanical numerical operation, not a
    modelling choice, so shared by every backend rather than each
    reimplementing it.

    Parameters
    ----------
    native_edges : ndarray (n_native+1,)
        Ascending bin edges of the source history (same units as
        `output_edges`, e.g. Gyr lookback time).
    native_sfr : ndarray (n_native,)
        SFR in each native bin (Msun/yr).
    output_edges : ndarray (n_output+1,)
        Ascending bin edges to rebin onto.

    Returns
    -------
    ndarray (n_output,)
        SFR in each output bin (Msun/yr).
    """
    native_edges = np.asarray(native_edges, dtype=float)
    native_sfr = np.asarray(native_sfr, dtype=float)
    output_edges = np.asarray(output_edges, dtype=float)

    native_widths = np.diff(native_edges)
    native_mass = native_sfr * native_widths  # Msun formed per native bin

    n_output = len(output_edges) - 1
    output_mass = np.zeros(n_output)
    for i in range(n_output):
        lo, hi = output_edges[i], output_edges[i + 1]
        overlap_lo = np.maximum(native_edges[:-1], lo)
        overlap_hi = np.minimum(native_edges[1:], hi)
        overlap = np.clip(overlap_hi - overlap_lo, 0.0, None)
        frac = np.divide(overlap, native_widths,
                         out=np.zeros_like(overlap), where=native_widths > 0)
        output_mass[i] = np.sum(frac * native_mass)

    output_widths = np.diff(output_edges)
    return output_mass / output_widths


class SharkGalaxyBackend:
    """Galaxy properties from a SHARK semi-analytic catalogue, matched to
    the halo catalogue via ``Epoch.galaxies_in_halo`` (position matching by
    default -- see api/simulation.py).

    When more than one SHARK galaxy falls within the match aperture
    (``r_scale`` * halo radius), the most massive one (by stellar mass) is
    treated as "the" galaxy hosted by this halo -- the central, in SAM
    convention. Returns ``{}`` (per the ``GalaxyBackend`` protocol) when
    nothing matches, or when no stellar mass is resolvable at all (so no
    central can be identified).

    Fields returned, and where they come from (see
    ``analysistools.shark.model.GALAXY_FIELDS`` for the full native-name
    registry):

    - ``StellarMass``, ``GasMass_Cold``, ``StarFormationRate``,
      ``BlackHoleMass``: the ``GalaxyCatalogue``'s own standardised
      aliases (``mass``, ``mgas``, ``sfr``, ``mbh``).
    - ``GasMass_Hot``: native ``"mhot"`` (not part of the standardised
      alias set).
    - ``MetallicityStellar``/``MetallicityGas``: metal mass / total mass,
      from the native ``*_metals_disk``/``*_metals_bulge`` fields -- a
      genuinely unit- and h-independent ratio, unlike everything else
      here (see units caveat below).

    Units caveat (same as ``pipeline.HaloExtractStage``/
    ``derived.EnvironmentStage``): values are returned exactly as
    ``GalaxyCatalogue`` provides them -- SHARK's own native convention
    (``Msun/h`` masses, ``Msun/Gyr/h`` star formation rates -- note that's
    *also* a Gyr-vs-yr mismatch against schema.py's declared "Msun/yr",
    not just little-h), not schema.py's declared units. Unit
    reconciliation isn't implemented anywhere in this pipeline yet.

    Deliberately **not** returned, documented rather than guessed:
    ``LuminosityV``/``Luminosity_ugriz``/``AbsoluteMagnitude_ugriz`` (need
    ``shark.photometry``'s separate, optional fsps-based pipeline -- not
    wired into this backend), ``HalfLightRadius`` (depends on
    luminosity), ``HalfMassRadiusStellar`` (``rstar_disk``/``rstar_bulge``
    are exponential/profile scale lengths, not half-mass radii --
    converting one to the other needs an assumed disk/bulge profile
    shape), ``SersicIndex`` (SHARK produces no structural fit at all).
    """

    name = "shark"

    def __init__(self, match_by: str = "position", r_scale: float = 1.0):
        self.match_by = match_by
        self.r_scale = r_scale

    def _match_central(self, epoch, halo_row: int):
        """Shared by galaxy_properties() and star_formation_history(): the
        matched view and the row within it treated as "the" (central)
        galaxy, or (None, None) if nothing matches."""
        matched = epoch.galaxies_in_halo(
            index=halo_row, match_by=self.match_by, r_scale=self.r_scale)
        if len(matched) == 0 or "mass" not in matched:
            return None, None
        mass = np.asarray(matched["mass"])
        return matched, int(np.argmax(mass))

    def galaxy_properties(self, epoch, halo_row: int) -> Dict[str, Any]:
        matched, central = self._match_central(epoch, halo_row)
        if matched is None:
            return {}

        def scalar(field: str, default: Optional[float] = None):
            arr = matched.get(field)
            return (float(np.asarray(arr)[central]) if arr is not None
                    else default)

        def metal_mass(disk_field: str, bulge_field: str):
            """Sum of disk+bulge metal mass, or None if *neither* component
            is available -- unlike scalar()'s own default=0.0 use below,
            this must not silently produce 0.0 for "field never existed"
            (that would compute a fake MetallicityStellar/Gas=0.0 instead
            of omitting it, exactly what the GalaxyBackend protocol says
            not to do)."""
            disk, bulge = scalar(disk_field), scalar(bulge_field)
            if disk is None and bulge is None:
                return None
            return (disk or 0.0) + (bulge or 0.0)

        stellar_mass = scalar("mass")
        gas_mass_cold = scalar("mgas")
        stellar_metal_mass = metal_mass("mstars_metals_disk",
                                        "mstars_metals_bulge")
        gas_metal_mass = metal_mass("mgas_metals_disk", "mgas_metals_bulge")

        props: Dict[str, Any] = {}
        if stellar_mass is not None:
            props["StellarMass"] = stellar_mass
            if stellar_metal_mass is not None and stellar_mass > 0:
                props["MetallicityStellar"] = \
                    stellar_metal_mass / stellar_mass
        if gas_mass_cold is not None:
            props["GasMass_Cold"] = gas_mass_cold
            if gas_metal_mass is not None and gas_mass_cold > 0:
                props["MetallicityGas"] = gas_metal_mass / gas_mass_cold
        hot = scalar("mhot")
        if hot is not None:
            props["GasMass_Hot"] = hot
        sfr = scalar("sfr")
        if sfr is not None:
            props["StarFormationRate"] = sfr
        mbh = scalar("mbh")
        if mbh is not None:
            props["BlackHoleMass"] = mbh

        return props

    def star_formation_history(self, epoch, halo_row: int, time_bin_edges):
        """SHARK's native disk+bulge SFH (``SharkModel.sfh_disk``/
        ``sfh_bulge``), rebinned from its own native lookback-time grid
        (``SharkModel.get_sfh_meta``'s ``delta_t``/``lbt_mean``) onto
        ``time_bin_edges`` via ``rebin_sfh``.

        Only available for a *model-backed* ``GalaxyCatalogue``
        (``epoch.galaxies`` built via ``GalaxyCatalogue.from_model``,
        i.e. ``.model`` is not ``None``) -- a file-backed one only reads
        ``galaxies.hdf5``, not the separate ``star_formation_histories
        .hdf5`` this needs. Returns ``None`` if unavailable for any
        reason (no match, file-backed catalogue, no native SFH).
        """
        matched, central = self._match_central(epoch, halo_row)
        if matched is None:
            return None

        model = getattr(matched, "model", None)
        if model is None:
            return None

        redshift = getattr(epoch, "redshift", None)
        if redshift is None:
            return None

        row = int(matched.index[central])
        try:
            sfh_disk = np.asarray(model.sfh_disk(redshift))[row]
            sfh_bulge = np.asarray(model.sfh_bulge(redshift))[row]
            sfh_meta = model.get_sfh_meta(redshift)
        except (KeyError, IndexError, AttributeError):
            return None

        native_sfr = sfh_disk + sfh_bulge
        lbt_mean = np.asarray(sfh_meta["lbt_mean"], dtype=float)      # Gyr
        delta_t = np.asarray(sfh_meta["delta_t"], dtype=float) / 1.0e3  # Myr -> Gyr

        # Sort by ascending lookback time -- SHARK's own storage order
        # isn't verified against real data anywhere in this codebase, so
        # don't assume it; reconstructing edges from (possibly unsorted)
        # centres+widths requires ascending order first.
        order = np.argsort(lbt_mean)
        lbt_mean, delta_t, native_sfr = \
            lbt_mean[order], delta_t[order], native_sfr[order]
        native_edges = np.empty(len(lbt_mean) + 1)
        native_edges[:-1] = lbt_mean - delta_t / 2.0
        native_edges[-1] = lbt_mean[-1] + delta_t[-1] / 2.0

        return rebin_sfh(native_edges, native_sfr, time_bin_edges)


def _cosmic_time_gyr(a, h0: float, omega_m: float, omega_lambda: float):
    """Age of a flat LambdaCDM universe at scale factor `a`, in Gyr.

    Standard analytic closed form, t(a) = (2 / (3 H0 sqrt(OL))) *
    asinh(sqrt(OL/Om) * a^1.5). Ignores radiation -- an excellent
    approximation for a > ~0.01 (z < ~100), far beyond anything relevant to
    stellar ages. Assumes flat curvature, matching essentially every
    cosmological simulation's own initial conditions.
    """
    hubble_time_gyr = 9.7779 / h0  # 1/H0 in Gyr, for H0 = 100 h km/s/Mpc
    return (2.0 / (3.0 * np.sqrt(omega_lambda))
           * np.arcsinh(np.sqrt(omega_lambda / omega_m) * np.asarray(a) ** 1.5)
           * hubble_time_gyr)


class HydroGalaxyBackend:
    """Galaxy properties synthesised from star particles in a hydrodynamic
    snapshot, matched to the halo via
    ``Epoch.particles_in_halo(species="star", ...)``.

    Fields returned (matching ``SharkGalaxyBackend``'s field names so a
    hydro catalogue is schema-identical wherever the physics is
    comparable), and how each is computed:

    - ``StellarMass``: sum of star particle ``mass`` (current bound mass,
      not initial mass -- accounts for stellar mass loss).
    - ``MetallicityStellar``: mass-weighted mean of ``stellar_Z``.
    - ``MeanStellarAge``: mass-weighted mean of (age at this epoch - age at
      formation), both converted from ``age`` (``StellarFormationTime``,
      a scale factor, per docs/snapshots.md) to cosmic time via
      ``_cosmic_time_gyr`` -- needs ``h0``/``omega_0``/``omega_lambda``/
      ``scale_factor`` on the matched particles' ``.meta``; omitted if any
      are unavailable, not silently skipped-with-wrong-units.
    - ``StarFormationRate``: (``initmass`` if available, else current
      ``mass``, of particles younger than ``young_star_age_threshold``,
      summed) / (``young_star_age_threshold`` converted Gyr -> yr).
      ``young_star_age_threshold`` is in Gyr (default 0.1 = 100 Myr, a
      standard window for this kind of "instantaneous" SFR estimate) --
      genuinely 0.0 (not omitted) when there are resolved ages but none
      young enough.
    - ``HalfMassRadiusStellar``: radius from the halo centre
      (``epoch.halos["pos"][halo_row]``) enclosing half the total star
      particle mass -- a direct geometric measurement, no assumed profile.

    Returns ``{}`` (per the ``GalaxyBackend`` protocol) when no star
    particles are found (e.g. a DMO snapshot, or an empty halo) or when
    ``mass`` isn't resolvable on them at all.

    Units caveat (same as ``SharkGalaxyBackend``/``pipeline.
    HaloExtractStage``): values are returned exactly as the matched
    ``SnapshotDataset`` view provides them (whatever comoving/little_h
    state the Epoch's snapshot was configured with), not necessarily
    schema.py's declared units. Unit reconciliation isn't implemented
    anywhere in this pipeline yet.

    Deliberately **not** returned, documented rather than guessed:
    ``HalfLightRadius``/luminosity/magnitude fields -- need a stellar
    population synthesis model (mass + age + metallicity -> luminosity),
    which nothing in this codebase provides for particle data (unlike
    SHARK, which computes luminosities itself when its own photometry
    pipeline is run). ``GasMass_Cold``/``GasMass_Hot``/``BlackHoleMass``
    would need gas-temperature and black-hole particle aggregation
    respectively -- out of scope for a star-particle-only backend as
    originally scoped here; add a similar aggregation over
    ``species="gas"``/``"bh"`` if you need them.
    """

    name = "hydro"

    def __init__(self, r_scale: float = 1.0,
                 young_star_age_threshold: float = 0.1):
        self.r_scale = r_scale
        self.young_star_age_threshold = young_star_age_threshold

    @staticmethod
    def _star_ages_gyr(stars):
        """(age_gyr, formation_mass) per star particle, or (None, None) if
        either the ``age``/``initmass``-or-``mass`` inputs or the cosmology
        needed to convert formation scale factor -> Gyr aren't available.
        Shared by galaxy_properties() and star_formation_history() so both
        use the exact same age definition."""
        if "age" not in stars:
            return None, None
        h0 = stars.meta.get("h0")
        omega_m = stars.meta.get("omega_0")
        omega_lambda = stars.meta.get("omega_lambda")
        a_now = stars.meta.get("scale_factor")
        if None in (h0, omega_m, omega_lambda, a_now):
            return None, None

        a_form = np.asarray(stars["age"])
        t_form = _cosmic_time_gyr(a_form, h0, omega_m, omega_lambda)
        t_now = _cosmic_time_gyr(a_now, h0, omega_m, omega_lambda)
        age_gyr = t_now - t_form

        formation_mass = (np.asarray(stars["initmass"]) if "initmass" in stars
                          else np.asarray(stars["mass"]))
        return age_gyr, formation_mass

    def galaxy_properties(self, epoch, halo_row: int) -> Dict[str, Any]:
        from ..merger_tree_types import periodic_delta

        stars = epoch.particles_in_halo(index=halo_row, r_scale=self.r_scale,
                                        species="star")
        if len(stars) == 0 or "mass" not in stars:
            return {}

        mass = np.asarray(stars["mass"])
        total_mass = float(np.sum(mass))

        props: Dict[str, Any] = {}
        if total_mass <= 0:
            return props
        props["StellarMass"] = total_mass

        if "stellar_Z" in stars:
            metallicity = np.asarray(stars["stellar_Z"])
            props["MetallicityStellar"] = float(
                np.average(metallicity, weights=mass))

        age_gyr, formation_mass = self._star_ages_gyr(stars)
        if age_gyr is not None:
            props["MeanStellarAge"] = float(np.average(age_gyr, weights=mass))

            young = age_gyr < self.young_star_age_threshold
            if np.any(young):
                window_yr = self.young_star_age_threshold * 1.0e9
                props["StarFormationRate"] = \
                    float(np.sum(formation_mass[young])) / window_yr
            else:
                props["StarFormationRate"] = 0.0

        centre = np.asarray(epoch.halos["pos"])[halo_row]
        pos = np.asarray(stars["pos"])
        d = periodic_delta(pos - centre, epoch.boxsize)
        distance = np.linalg.norm(d, axis=1)
        order = np.argsort(distance)
        cum_mass = np.cumsum(mass[order])
        idx = min(int(np.searchsorted(cum_mass, total_mass / 2.0)),
                  len(distance) - 1)
        props["HalfMassRadiusStellar"] = float(distance[order][idx])

        return props

    def star_formation_history(self, epoch, halo_row: int, time_bin_edges):
        """Histogram of formation mass (``initmass`` if available, else
        current ``mass``) over each particle's age (Gyr lookback,
        ``_star_ages_gyr`` -- the same age definition
        ``galaxy_properties()`` uses for ``MeanStellarAge``), binned
        directly onto ``time_bin_edges`` (no rebinning needed: unlike
        SHARK's pre-tabulated history, particle ages can be binned onto
        any grid directly). Returns ``None`` if no star particles are
        found, or age/cosmology aren't resolvable (see
        ``_star_ages_gyr``).
        """
        stars = epoch.particles_in_halo(index=halo_row, r_scale=self.r_scale,
                                        species="star")
        if len(stars) == 0:
            return None

        age_gyr, formation_mass = self._star_ages_gyr(stars)
        if age_gyr is None:
            return None

        edges = np.asarray(time_bin_edges, dtype=float)
        mass_per_bin, _ = np.histogram(age_gyr, bins=edges,
                                       weights=formation_mass)
        widths_gyr = np.diff(edges)
        return mass_per_bin / (widths_gyr * 1.0e9)


BACKENDS = {
    "shark": SharkGalaxyBackend,
    "hydro": HydroGalaxyBackend,
}


def get_backend(name: str, **kwargs) -> GalaxyBackend:
    """Resolve a project config's ``galaxy_backend: shark|hydro`` string to
    an instantiated backend."""
    try:
        cls = BACKENDS[name]
    except KeyError:
        raise ValueError(
            f"Unknown galaxy_backend '{name}'; available: "
            f"{sorted(BACKENDS)}")
    return cls(**kwargs)
