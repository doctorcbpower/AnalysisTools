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

from typing import Any, Dict, Optional, Protocol, Sequence

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

    def native_comoving_little_h(self, epoch=None):
        """Optional: (comoving, little_h) of the values
        ``galaxy_properties()`` returns, or ``(None, None)``/absent if
        unknown. Not part of the required protocol (checked with
        ``hasattr`` by ``pipeline.CrossMatchStage``) since a backend with
        no fixed or epoch-derived native state has nothing meaningful to
        report here. Used by ``validation.SchemaValidator``'s little-h/
        comoving cross-check against schema.py's declared units."""
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

    Native mass at lookback times outside ``[output_edges[0],
    output_edges[-1]]`` cannot be placed in any output bin and is
    dropped, not folded into the nearest edge bin -- e.g. a satellite
    whose star formation genuinely extends further back than the
    requested grid (a common case: ``time_bin_edges``'s upper edge needs
    to reach at least the age of the universe at z=0 for *this*
    cosmology, not a plausible-looking round number -- see
    ``backends._cosmic_time_gyr``) silently loses that early-formed
    mass. Logs a warning (once per call, not per satellite) whenever
    that drop is non-negligible, since a caller comparing the rebinned
    SFH's integral against an independently-computed total stellar mass
    (e.g. ``validation.PhysicalValidator``'s ``stellar_mass_exceeds_
    formed_mass`` check) would otherwise see a silent, confusing
    mass-conservation violation with no indication of the cause.
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

    total_native = float(np.sum(native_mass))
    if total_native > 0:
        lost_frac = (total_native - float(np.sum(output_mass))) / total_native
        if lost_frac > 1e-6:
            logger.warning(
                "rebin_sfh: %.2f%% of the native-grid formed mass falls "
                "outside output_edges=[%.3g, %.3g] (native range "
                "[%.3g, %.3g]) and was dropped, not folded into the "
                "nearest bin. If this is lookback time in Gyr, widen "
                "time_bin_edges' upper bound to at least the age of the "
                "universe at z=0 for this cosmology (see "
                "backends._cosmic_time_gyr) if satellites' true star "
                "formation extends further back than the requested grid.",
                lost_frac * 100.0, output_edges[0], output_edges[-1],
                native_edges[0], native_edges[-1])

    output_widths = np.diff(output_edges)
    return output_mass / output_widths


#: FSPS bands SharkGalaxyBackend requests when compute_photometry=True --
#: "v" (Buser V) plus fsps's standard SDSS filter names. Not validated
#: here (a wrong name only surfaces as an fsps error at first photometry
#: call, which needs a real python-fsps install to hit regardless).
DEFAULT_PHOTOMETRY_BANDS = ("v", "sdss_u", "sdss_g", "sdss_r", "sdss_i",
                           "sdss_z")

#: AB absolute solar magnitude per band -- the reference needed to convert
#: an FSPS absolute AB magnitude into a solar-luminosity ratio (L/Lsun =
#: 10^(-0.4*(M - M_sun))). "v" (4.83) matches the constant
#: ``shark.photometry.PhotometryPipeline.mass_to_light`` already uses, for
#: consistency within this package; SDSS ugriz values are the AB solar
#: magnitudes from Willmer (2018, ApJS 236, 47), Table 3.
SOLAR_ABSOLUTE_MAGNITUDE = {
    "v": 4.83,
    "sdss_u": 6.39, "sdss_g": 5.11, "sdss_r": 4.65,
    "sdss_i": 4.53, "sdss_z": 4.50,
}
_UGRIZ_BANDS = ("sdss_u", "sdss_g", "sdss_r", "sdss_i", "sdss_z")


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
    - ``SharkGalaxyID``: native ``id_galaxy``, the foreign key into the
      SHARK output the schema field is documented as.
    - ``LuminosityV``, ``Luminosity_ugriz``, ``AbsoluteMagnitude_ugriz``:
      only when constructed with ``compute_photometry=True`` -- see below.
      Off by default: FSPS SSP convolution is orders of magnitude more
      expensive per galaxy than every other field here, needs
      python-fsps installed with its SSP template data (``SPS_HOME``) --
      a heavier, optional dependency this backend doesn't otherwise need
      -- and is only meaningful for a *model-backed* ``GalaxyCatalogue``
      (same restriction as ``star_formation_history()``, since it needs
      the native per-snapshot SFH ``shark.photometry`` convolves).

    Units caveat (same as ``pipeline.HaloExtractStage``/
    ``derived.EnvironmentStage``): values are returned exactly as
    ``GalaxyCatalogue`` provides them -- SHARK's own native convention
    (``Msun/h`` masses, ``Msun/Gyr/h`` star formation rates -- note that's
    *also* a Gyr-vs-yr mismatch against schema.py's declared "Msun/yr",
    not just little-h), not schema.py's declared units. Unit
    reconciliation isn't implemented anywhere in this pipeline yet.
    ``LuminosityV``/``Luminosity_ugriz`` are the one exception: FSPS
    absolute magnitudes are h-independent, so the Lsun values derived from
    them are too.

    Deliberately **not** returned, documented rather than guessed:
    ``HalfLightRadius`` (depends on a light *profile*, not just total
    luminosity -- FSPS gives an SED, not a spatial distribution),
    ``HalfMassRadiusStellar`` (``rstar_disk``/``rstar_bulge`` are
    exponential/profile scale lengths, not half-mass radii -- converting
    one to the other needs an assumed disk/bulge profile shape),
    ``SersicIndex`` (SHARK produces no structural fit at all).
    """

    name = "shark"

    def __init__(self, match_by: str = "position", r_scale: float = 1.0,
                 compute_photometry: bool = False,
                 photometry_bands: Sequence[str] = DEFAULT_PHOTOMETRY_BANDS,
                 photometry_options: Optional[Dict[str, Any]] = None,
                 photometry_pipeline_factory=None):
        """
        Parameters
        ----------
        compute_photometry : bool
            If True, ``galaxy_properties()`` also computes
            ``LuminosityV``/``Luminosity_ugriz``/``AbsoluteMagnitude_ugriz``
            via ``shark.photometry.PhotometryPipeline`` (real python-fsps
            SSP convolution -- requires it installed with ``SPS_HOME``
            set). Default False.
        photometry_bands : sequence of str
            FSPS filter names to request; must include ``"v"`` and every
            name in ``("sdss_u", "sdss_g", "sdss_r", "sdss_i", "sdss_z")``
            for the corresponding output field to be populated (a band
            missing from this list simply leaves that field/entry
            unpopulated, same "omit, don't fabricate" contract as every
            other field). Default ``DEFAULT_PHOTOMETRY_BANDS`` (all six).
        photometry_options : dict, optional
            Extra kwargs forwarded to ``PhotometryPipeline`` (e.g.
            ``{"imf_type": 1, "add_dust": True}``).
        photometry_pipeline_factory : callable, optional
            ``(model, z_obs, bands=..., progress=..., **photometry_options)
            -> PhotometryPipeline``-like object exposing ``.abs_mags(idx)``.
            Defaults to the real ``shark.photometry.PhotometryPipeline``;
            override for testing (avoids a real fsps/SPS_HOME dependency)
            or to reuse a pipeline built elsewhere.
        """
        self.match_by = match_by
        self.r_scale = r_scale
        self.compute_photometry = compute_photometry
        self.photometry_bands = list(photometry_bands)
        self.photometry_options = dict(photometry_options or {})
        self._photometry_pipeline_factory = photometry_pipeline_factory
        self._photometry_cache: Dict[Any, Any] = {}

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
        galaxy_id = scalar("id_galaxy")
        if galaxy_id is not None:
            props["SharkGalaxyID"] = int(galaxy_id)

        if self.compute_photometry:
            self._add_photometry(props, matched, central, epoch)

        return props

    def _photometry_pipeline_for(self, model, z_obs: float):
        """One ``PhotometryPipeline`` per (model, z_obs), cached and
        reused across every ``galaxy_properties()`` call in a catalogue
        build -- constructing FSPS's ``StellarPopulation`` is far too
        expensive to redo per satellite (see ``shark.photometry.sps``'s
        own "reused across all galaxies" design note)."""
        key = (id(model), float(z_obs))
        if key not in self._photometry_cache:
            factory = self._photometry_pipeline_factory
            if factory is None:
                from ..shark.photometry import PhotometryPipeline
                factory = PhotometryPipeline
            self._photometry_cache[key] = factory(
                model, z_obs, bands=self.photometry_bands, progress=False,
                **self.photometry_options)
        return self._photometry_cache[key]

    def _add_photometry(self, props: Dict[str, Any], matched, central: int,
                        epoch) -> None:
        """Adds ``LuminosityV``/``Luminosity_ugriz``/
        ``AbsoluteMagnitude_ugriz`` to ``props`` in place, or leaves them
        unset -- never NaN-fills -- on any of: no model-backed catalogue
        (same restriction as ``star_formation_history()``), no epoch
        redshift, or an FSPS failure for this specific galaxy (logged as a
        warning, not raised -- one galaxy's SED failing shouldn't abort
        an entire catalogue build)."""
        model = getattr(matched, "model", None)
        if model is None:
            return
        redshift = getattr(epoch, "redshift", None)
        if redshift is None:
            return

        row = int(matched.index[central])
        try:
            pipeline = self._photometry_pipeline_for(model, redshift)
            mags = np.asarray(pipeline.abs_mags(np.array([row])))[0]
        except Exception:
            logger.warning(
                "%s: photometry failed for galaxy row %d; LuminosityV/"
                "Luminosity_ugriz/AbsoluteMagnitude_ugriz omitted.",
                self.name, row, exc_info=True)
            return

        band_index = {b: i for i, b in enumerate(self.photometry_bands)}

        v_idx = band_index.get("v")
        if v_idx is not None and not np.isnan(mags[v_idx]):
            props["LuminosityV"] = float(10.0 ** (
                -0.4 * (mags[v_idx] - SOLAR_ABSOLUTE_MAGNITUDE["v"])))

        if all(band in band_index for band in _UGRIZ_BANDS):
            ugriz_mags = np.array(
                [mags[band_index[band]] for band in _UGRIZ_BANDS])
            if not np.all(np.isnan(ugriz_mags)):
                props["AbsoluteMagnitude_ugriz"] = ugriz_mags
                sun_mags = np.array(
                    [SOLAR_ABSOLUTE_MAGNITUDE[band] for band in _UGRIZ_BANDS])
                with np.errstate(invalid="ignore"):
                    props["Luminosity_ugriz"] = \
                        10.0 ** (-0.4 * (ugriz_mags - sun_mags))

    def native_comoving_little_h(self, epoch=None):
        """(comoving, little_h) of the values ``galaxy_properties()``
        returns -- always ``(True, True)`` for SHARK: its output is always
        comoving Mpc/h positions and Msun/h masses regardless of ``epoch``
        (unlike ``HydroGalaxyBackend``, whose native state follows the
        snapshot it reads). Used by ``pipeline.CrossMatchStage`` to record
        the native unit state for ``validation.SchemaValidator``'s
        little-h/comoving cross-check -- see the units caveat above."""
        return True, True

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
        # delta_t: confirmed, against a real Dorcha SHARK run, to already
        # be in Gyr -- the *same* unit as lbt_mean, not Myr (an earlier,
        # never-empirically-verified assumption here divided it by 1e3
        # expecting a Myr->Gyr conversion, making every bin 1000x too
        # narrow and silently understating every formed-mass integral by
        # the same factor -- caught via PhysicalValidator's
        # stellar_mass_exceeds_formed_mass check against real data:
        # consecutive lbt_mean bin centres are spaced by exactly delta_t,
        # not delta_t/1000). Used as-is now; see the sanity check below
        # for a SHARK output where that's genuinely not true.
        delta_t = np.asarray(sfh_meta["delta_t"], dtype=float)

        # Sort by ascending lookback time -- SHARK's own storage order
        # isn't verified against real data anywhere in this codebase, so
        # don't assume it; reconstructing edges from (possibly unsorted)
        # centres+widths requires ascending order first.
        order = np.argsort(lbt_mean)
        lbt_mean, delta_t, native_sfr = \
            lbt_mean[order], delta_t[order], native_sfr[order]

        if len(lbt_mean) > 1:
            span_from_centres = float(lbt_mean[-1] - lbt_mean[0])
            span_from_widths = float(np.sum(delta_t))
            if span_from_widths > 0 and not (
                    0.1 < span_from_centres / span_from_widths < 10.0):
                logger.warning(
                    "SharkGalaxyBackend.star_formation_history(): SFH "
                    "bin widths (delta_t, sum=%.3g) are inconsistent "
                    "with the lookback-time span implied by bin centres "
                    "(lbt_mean, range=%.3g) by more than 10x -- delta_t "
                    "may not be in Gyr (the same unit as lbt_mean) for "
                    "this SHARK output, in which case every formed-mass "
                    "integral computed from this SFH will be wrong by "
                    "whatever factor separates them. Check "
                    "SharkModel.get_sfh_meta's units against your own "
                    "star_formation_histories.hdf5 before trusting SFH-"
                    "derived fields.", span_from_widths, span_from_centres)

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

    def native_comoving_little_h(self, epoch=None):
        """(comoving, little_h) of the values ``galaxy_properties()``
        returns -- follows the ``Epoch``'s snapshot config (star particles
        inherit whatever comoving/little_h state the ``SnapshotDataset``
        was constructed with), unlike ``SharkGalaxyBackend``'s fixed
        native convention. ``(None, None)`` if ``epoch``/``epoch.snapshot``
        isn't available. Used by ``pipeline.CrossMatchStage`` -- see
        ``SharkGalaxyBackend.native_comoving_little_h``."""
        snapshot = getattr(epoch, "snapshot", None) if epoch is not None else None
        if snapshot is None:
            return None, None
        return snapshot.meta.get("comoving"), snapshot.meta.get("little_h")

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
