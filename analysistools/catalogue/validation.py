#!/usr/bin/env python3
"""
analysistools.catalogue.validation
--------------------------------------
Validator / ValidationReport: automated checks run by
``pipeline.QualityControlStage`` before a catalogue is written. See
docs/dorcha_master_catalogue_design.md section 7 for the full recommended
check list (schema, integrity, physical plausibility, cross-match
verification, statistical sanity, reproducibility); this module provides
the framework and the first three categories as concrete classes.
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from collections import Counter
from dataclasses import asdict, dataclass, field
from typing import Any, Dict, List

import json
import logging

import numpy as np

logger = logging.getLogger(__name__)


@dataclass
class ValidationIssue:
    severity: str            # "error" | "warning"
    check: str
    message: str
    field: str = ""


@dataclass
class ValidationReport:
    issues: List[ValidationIssue] = field(default_factory=list)

    def add(self, severity: str, check: str, message: str,
            field: str = "") -> None:
        self.issues.append(ValidationIssue(severity, check, message, field))

    @property
    def errors(self) -> List[ValidationIssue]:
        return [i for i in self.issues if i.severity == "error"]

    @property
    def warnings(self) -> List[ValidationIssue]:
        return [i for i in self.issues if i.severity == "warning"]

    @property
    def passed(self) -> bool:
        return not self.errors

    def summary(self) -> None:
        print(f"ValidationReport: {len(self.errors)} error(s), "
              f"{len(self.warnings)} warning(s)")
        for issue in self.issues:
            print(f"  [{issue.severity.upper():7s}] {issue.check}: "
                  f"{issue.message}" + (f" ({issue.field})" if issue.field
                                        else ""))

    def to_json(self) -> str:
        """Serialise for ``Provenance/validation_report`` (written by
        ``pipeline.WriteStage`` when a report is attached to
        ``context.meta``)."""
        return json.dumps({
            "passed": self.passed,
            "n_errors": len(self.errors),
            "n_warnings": len(self.warnings),
            "issues": [asdict(issue) for issue in self.issues],
        }, indent=2)


class Validator(ABC):
    """One category of automated check."""

    name: str = "validator"

    @abstractmethod
    def check(self, context: "analysistools.catalogue.pipeline.PipelineContext",
              schema: "analysistools.catalogue.schema.CatalogueSchema"
              ) -> ValidationReport:
        ...


#: Substring -> unit family, checked in order (first match wins) so e.g.
#: "Msun/yr" (sfr) is classified before the bare "Msun" (mass) check.
_UNIT_FAMILY_PATTERNS = (
    (("msun", "/yr"), "sfr"), (("msun", "/gyr"), "sfr"),
    (("msun",), "mass"),
    (("kpc",), "length"), (("mpc",), "length"), (("pc",), "length"),
    (("km/s",), "velocity"),
    (("gyr",), "time"), (("yr",), "time"),
    (("deg",), "angle"),
    (("dimensionless",), "dimensionless"),
)

#: dtype string (schema.FieldSpec.dtype) -> numpy scalar kind it must be a
#: subtype of. "str" isn't included: HDF5/numpy string handling is varied
#: enough (fixed-width bytes, object dtype, h5py variable-length strings)
#: that a strict numpy-dtype check would false-positive constantly.
_DTYPE_KIND = {
    "float64": np.floating, "float32": np.floating,
    "int64": np.integer, "int32": np.integer, "int16": np.integer,
    "int8": np.integer, "bool": np.bool_,
}


def _unit_family(units: str) -> str:
    """Coarse physical-dimension family of a units string (schema.py's or
    a recorded native one), or "unknown" if no pattern matches. A
    heuristic classifier, not a unit parser -- used only to catch the
    specific, concretely-documented "little-h/comoving not reflected in
    the declared units string" mismatch below, not to validate arbitrary
    unit strings."""
    lowered = units.lower()
    for needles, family in _UNIT_FAMILY_PATTERNS:
        if all(n in lowered for n in needles):
            return family
    return "unknown"


class SchemaValidator(Validator):
    """Every field in the schema is present with declared dtype and units
    attribute; no undeclared datasets exist; units strings are recognised
    (matches the declared unit family, e.g. a "mass" field must resolve to
    a mass-like unit string).

    Three checks, run per schema group against ``context.columns``:

    1. **Presence** (via ``schema.validate_columns``): a schema field with
       no matching column is a *warning* (many fields are legitimately
       still deferred -- see docs/phase6_remaining_work.md -- this is
       informational, not a build failure); a column present under a
       group that group's schema field list doesn't declare is an
       *error* (there is no legitimate reason for this -- it means a
       typo'd field name or a schema/pipeline version mismatch).
    2. **dtype**: the column's numpy dtype must be a subtype of the
       schema's declared dtype's kind (see ``_DTYPE_KIND`` -- ``"str"``
       fields are skipped, see its docstring). Mismatch is a *warning*:
       a narrower/wider numeric dtype than declared isn't necessarily
       wrong, just worth flagging before ``WriteStage`` casts it.
    3. **little-h / comoving vs. declared units**: cross-checks
       ``context.meta["comoving_little_h"]`` (recorded by
       ``pipeline.HaloExtractStage``/``CrossMatchStage`` from the actual
       source's ``little_h``/``comoving`` state -- see their docstrings'
       "units caveat") against whether the schema's declared units string
       for that field says so (``"/h"``, ``"(comoving)"``). This is the
       concrete version of the "unit reconciliation isn't implemented
       anywhere in the pipeline yet" caveat repeated across
       ``derived.py``/``backends.py`` -- it catches exactly the silent
       mismatch those caveats warn about (e.g. SHARK's native ``Msun/h``
       written under a schema field declared plain ``"Msun"``), not a
       full physical unit *conversion* (this pipeline still writes native
       values as-is; only the *labelling* mismatch is checked). Always a
       *warning*: this validator doesn't convert or reject values, only
       flags them for a human (or a future ``WriteStage`` unit-conversion
       pass) to look at.

    Only groups actually declared in the schema are checked; working
    state under ``Satellites/_internal/*`` and ``MergerTrees/*`` (neither
    of which is part of the schema -- they're intermediate pipeline state,
    not catalogue output) is intentionally not compared against anything.
    """

    name = "schema"

    def check(self, context, schema) -> ValidationReport:
        report = ValidationReport()
        comoving_little_h = context.meta.get("comoving_little_h", {})

        columns_by_group: Dict[str, Dict[str, Any]] = {}
        for path, value in context.columns.items():
            group, _, name = path.rpartition("/")
            if group in ("Satellites/_internal", "MergerTrees") or not group:
                continue
            columns_by_group.setdefault(group, {})[name] = value

        checked_groups = set(schema.groups) | set(columns_by_group)
        for group in sorted(checked_groups):
            present = columns_by_group.get(group, {})
            specs = {s.name: s for s in schema.by_group(group)}

            if group not in schema.groups:
                report.add("error", "schema_group",
                          f"group '{group}' has columns but is not "
                          f"declared anywhere in schema v{schema.version}.",
                          field=group)
                continue

            for problem in schema.validate_columns(group, present.keys()):
                severity = "error" if problem.startswith("undeclared") \
                    else "warning"
                report.add(severity, "schema_presence", problem,
                          field=group)

            state = comoving_little_h.get(group, {})
            comoving = state.get("comoving")
            little_h = state.get("little_h")

            for name, value in present.items():
                spec = specs.get(name)
                if spec is None:
                    continue  # already reported above as "undeclared"

                arr = np.asarray(value) if not isinstance(value, list) \
                    else None
                expected_kind = _DTYPE_KIND.get(spec.dtype)
                if (arr is not None and expected_kind is not None
                        and not np.issubdtype(arr.dtype, expected_kind)):
                    report.add(
                        "warning", "schema_dtype",
                        f"{group}/{name}: dtype {arr.dtype} does not "
                        f"match schema-declared '{spec.dtype}'.",
                        field=f"{group}/{name}")

                if comoving is None and little_h is None:
                    continue
                family = _unit_family(spec.units)
                units_lower = spec.units.lower()
                if little_h and family in ("mass", "length", "sfr") \
                        and "/h" not in units_lower:
                    report.add(
                        "warning", "schema_units",
                        f"{group}/{name}: native values are little-h-"
                        f"scaled but schema declares units "
                        f"'{spec.units}' without an '/h' suffix.",
                        field=f"{group}/{name}")
                if comoving and family == "length" \
                        and "comoving" not in units_lower:
                    report.add(
                        "warning", "schema_units",
                        f"{group}/{name}: native values are comoving but "
                        f"schema declares units '{spec.units}' without "
                        f"'(comoving)'.", field=f"{group}/{name}")

        return report


class IntegrityValidator(Validator):
    """SatelliteID uniqueness and row-count consistency across every
    subgroup (the row-order invariant); no orphaned foreign keys; no
    unexpected inf (this pipeline's documented "not computed" sentinel is
    always NaN, never inf -- see every ``derived.py`` stage's own
    docstring -- so an inf anywhere is always a real bug, e.g. a
    divide-by-zero, not a legitimate missingness value; that makes "no
    unexpected inf" precise without needing a separately-specified
    allowed-missingness mask); particle-tag CSR offsets monotonically
    non-decreasing.

    ``schema`` is accepted (per the ``Validator`` interface) but unused --
    every check here is about internal consistency of ``context.columns``
    itself, not about the schema."""

    name = "integrity"

    def check(self, context, schema) -> ValidationReport:
        report = ValidationReport()
        cols = context.columns

        lengths: Dict[str, int] = {}
        for path, value in cols.items():
            if not (path.startswith("Satellites/")
                    or path.startswith("MergerTrees/")):
                continue
            try:
                lengths[path] = len(value)
            except TypeError:
                continue

        if lengths:
            canonical_path = "Satellites/Identification/SatelliteID"
            canonical_n = lengths.get(canonical_path)
            if canonical_n is None:
                canonical_n = Counter(lengths.values()).most_common(1)[0][0]
            for path, n in lengths.items():
                if n != canonical_n:
                    report.add(
                        "error", "row_count",
                        f"{path}: length {n} does not match the "
                        f"catalogue's satellite row count ({canonical_n}).",
                        field=path)

        satellite_id = cols.get("Satellites/Identification/SatelliteID")
        if satellite_id is not None:
            ids = np.asarray(satellite_id)
            unique, counts = np.unique(ids, return_counts=True)
            duplicated = unique[counts > 1]
            if duplicated.size:
                shown = sorted(duplicated.tolist())[:10]
                report.add(
                    "error", "satellite_id_uniqueness",
                    f"{duplicated.size} duplicated SatelliteID value(s): "
                    f"{shown}" + (" ..." if duplicated.size > 10 else ""),
                    field="Satellites/Identification/SatelliteID")

        host_ids = cols.get("Satellites/Identification/HostHaloID")
        known_hosts = cols.get("Haloes/HostHaloID")
        if host_ids is not None and known_hosts is not None:
            known = set(np.asarray(known_hosts).tolist())
            orphaned = sorted(set(np.asarray(host_ids).tolist()) - known)
            if orphaned:
                shown = orphaned[:10]
                report.add(
                    "error", "orphaned_foreign_key",
                    f"{len(orphaned)} satellite(s) reference HostHaloID "
                    f"value(s) not present in Haloes/HostHaloID: {shown}"
                    + (" ..." if len(orphaned) > 10 else ""),
                    field="Satellites/Identification/HostHaloID")

        for path, value in cols.items():
            if isinstance(value, list):
                continue
            arr = np.asarray(value)
            if not np.issubdtype(arr.dtype, np.floating):
                continue
            n_inf = int(np.count_nonzero(np.isinf(arr)))
            if n_inf:
                report.add(
                    "error", "unexpected_inf",
                    f"{path}: {n_inf} inf value(s) -- 'not computed' is "
                    f"always represented as NaN in this pipeline, never "
                    f"inf, so this indicates a real bug rather than a "
                    f"documented missingness sentinel.", field=path)

        offsets = cols.get("Satellites/ParticleTags/particle_id_offsets")
        if offsets is not None:
            offsets_arr = np.asarray(offsets)
            if offsets_arr.size > 1 and np.any(np.diff(offsets_arr) < 0):
                report.add(
                    "error", "csr_offsets_not_monotonic",
                    "Satellites/ParticleTags/particle_id_offsets is not "
                    "monotonically non-decreasing -- the CSR particle-tag "
                    "index is corrupt.",
                    field="Satellites/ParticleTags/particle_id_offsets")

        return report


class PhysicalValidator(Validator):
    """Mpeak >= M200c_z0; OrbitalApocentre >= OrbitalPericentre;
    HeliocentricDistance > 0; CompletenessWeight >= 1; a backsplash
    satellite has at least one recorded infall; total SFH-formed mass >=
    StellarMass (warning only, not an error -- see below).

    Two checks deliberately deviate from the one-line description
    originally sketched for this validator, for the same reason
    ``derived.StarFormationHistoryStage`` deviates from the schema's
    literal quenching definition -- the literal version isn't actually
    true of realistic physics, so checking it as written would flag
    correct catalogues as broken:

    - "RedshiftInfall monotonic along accretion history" isn't a
      well-defined check: ``RedshiftInfall`` is a single scalar (the one
      infall event ``HaloPropertiesStage`` records), not a per-snapshot
      history, so there is nothing to walk monotonically. Substituted:
      every satellite flagged ``IsBacksplash`` must have
      ``NumberOfInfalls >= 1`` -- a backsplash satellite (currently
      outside the host after having been inside) must, by definition, have
      had at least one recorded infall.
    - "SFH bins sum to StellarMass within tolerance" ignores stellar mass
      loss/return (evolved stars return gas to the ISM), so an *exact*
      match would only hold for a population with zero mass loss --
      checking equality would flag every realistically-evolved galaxy.
      Substituted: total mass formed (the SFH integral) should be >=
      ``StellarMass`` (mass loss can only reduce today's stellar mass
      below what was formed in-situ, never increase it) -- a *warning*,
      not an error (see below), since it is not actually a hard
      invariant: ``StellarMass`` (e.g. ``mstars_disk``/``mstars_bulge``
      from a SAM like SHARK) includes stars accreted *ex-situ* via
      mergers, but a galaxy's own tracked SFH typically records only the
      stars it formed *in-situ* -- a real central galaxy built mostly
      through mergers can legitimately have StellarMass far exceeding
      its own SFH integral with a near-empty recorded history (raw
      ``sfh_disk``/``sfh_bulge`` genuinely ~0 for its full history) and
      nothing wrong with the data. Confirmed against a real Dorcha
      catalogue after ruling out three genuine bugs this check
      previously caught (a bulge-SFH channel omission, a `time_bin_edges`
      truncation, and a `delta_t` unit error, all now fixed elsewhere in
      this codebase) -- once those were fixed, this specific relationship
      still didn't hold for an otherwise clean, fully-consistent
      (`type=0` central, matching mass fields, consistent time units)
      galaxy, which is what downgraded this check from error to warning.

    Every check skips NaN/absent inputs rather than flagging them --
    "not computed" is a legitimate, documented state (see e.g.
    ``StarFormationHistoryStage``'s "no computable SFH" case); a physical
    plausibility check only applies to values that were actually
    computed. ``schema`` is accepted (per the ``Validator`` interface) but
    unused -- these are physical relationships between already-identified
    fields, not schema-driven.
    """

    name = "physical"

    def check(self, context, schema) -> ValidationReport:
        report = ValidationReport()
        cols = context.columns

        def _pair(path_a, path_b):
            a, b = cols.get(path_a), cols.get(path_b)
            if a is None or b is None:
                return None
            a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
            return a, b, ~np.isnan(a) & ~np.isnan(b)

        pair = _pair("Satellites/HaloProperties/Mpeak",
                    "Satellites/HaloProperties/M200c_z0")
        if pair is not None:
            mpeak, m200c_z0, valid = pair
            bad = valid & (mpeak < m200c_z0)
            if np.any(bad):
                report.add(
                    "error", "mpeak_below_m200c_z0",
                    f"{int(np.count_nonzero(bad))} satellite(s) have "
                    f"Mpeak < M200c_z0 -- the peak historical mass cannot "
                    f"be below the present-day mass.",
                    field="Satellites/HaloProperties/Mpeak")

        pair = _pair("Satellites/HaloProperties/OrbitalApocentre",
                    "Satellites/HaloProperties/OrbitalPericentre")
        if pair is not None:
            apo, peri, valid = pair
            bad = valid & (apo < peri)
            if np.any(bad):
                report.add(
                    "error", "apocentre_below_pericentre",
                    f"{int(np.count_nonzero(bad))} satellite(s) have "
                    f"OrbitalApocentre < OrbitalPericentre.",
                    field="Satellites/HaloProperties/OrbitalApocentre")

        helio = cols.get("Satellites/Observability/HeliocentricDistance")
        if helio is not None:
            helio = np.asarray(helio, dtype=float)
            valid = ~np.isnan(helio)
            bad = valid & (helio <= 0)
            if np.any(bad):
                report.add(
                    "error", "non_positive_heliocentric_distance",
                    f"{int(np.count_nonzero(bad))} satellite(s) have "
                    f"HeliocentricDistance <= 0.",
                    field="Satellites/Observability/HeliocentricDistance")

        completeness = cols.get("Satellites/Observability/CompletenessWeight")
        if completeness is not None:
            completeness = np.asarray(completeness, dtype=float)
            valid = ~np.isnan(completeness)
            bad = valid & (completeness < 1.0)
            if np.any(bad):
                report.add(
                    "error", "completeness_weight_below_one",
                    f"{int(np.count_nonzero(bad))} satellite(s) have "
                    f"CompletenessWeight < 1.",
                    field="Satellites/Observability/CompletenessWeight")

        is_backsplash = cols.get("Satellites/HaloProperties/IsBacksplash")
        n_infalls = cols.get("Satellites/HaloProperties/NumberOfInfalls")
        if is_backsplash is not None and n_infalls is not None:
            is_backsplash = np.asarray(is_backsplash, dtype=bool)
            n_infalls = np.asarray(n_infalls)
            bad = is_backsplash & (n_infalls < 1)
            if np.any(bad):
                report.add(
                    "error", "backsplash_without_infall",
                    f"{int(np.count_nonzero(bad))} satellite(s) are "
                    f"flagged IsBacksplash but have NumberOfInfalls < 1 -- "
                    f"a backsplash satellite must have had at least one "
                    f"recorded infall.",
                    field="Satellites/HaloProperties/IsBacksplash")

        sfh = cols.get("Satellites/GalaxyProperties/SFH")
        stellar_mass = cols.get("Satellites/GalaxyProperties/StellarMass")
        time_bin_edges = context.meta.get("time_bin_edges_sfh")
        if (sfh is not None and stellar_mass is not None
                and time_bin_edges is not None):
            sfh = np.asarray(sfh, dtype=float)
            stellar_mass = np.asarray(stellar_mass, dtype=float)
            widths_yr = np.diff(
                np.asarray(time_bin_edges, dtype=float)) * 1.0e9
            with np.errstate(invalid="ignore"):
                formed_mass = np.nansum(sfh * widths_yr, axis=1)
            has_sfh = ~np.all(np.isnan(sfh), axis=1)
            valid = has_sfh & ~np.isnan(stellar_mass)
            bad = valid & (formed_mass < stellar_mass * (1.0 - 1e-6))
            if np.any(bad):
                # ratio (not just counts) is the fast way to tell a units
                # bug (e.g. a Gyr/yr mixup -- ratio near 1e9 or 1e-9) from
                # a small residual (ratio near 1) at a glance, without
                # having to dig the raw arrays back out of the context.
                with np.errstate(divide="ignore", invalid="ignore"):
                    ratio = np.where(bad, stellar_mass / formed_mass, np.nan)
                worst = int(np.nanargmax(ratio))
                report.add(
                    "warning", "stellar_mass_exceeds_formed_mass",
                    f"{int(np.count_nonzero(bad))} satellite(s) have "
                    f"StellarMass exceeding the total mass formed "
                    f"(integral of SFH). Worst case: row {worst}, "
                    f"StellarMass={stellar_mass[worst]:.3e}, "
                    f"formed_mass={formed_mass[worst]:.3e} (ratio "
                    f"{ratio[worst]:.3e}). A ratio near a round power of "
                    f"ten (e.g. ~1e9/1e-9) usually means a units mismatch "
                    f"(e.g. Msun/Gyr vs. Msun/yr) somewhere upstream -- "
                    f"worth checking. A smaller/less clean ratio is often "
                    f"legitimate: StellarMass can include mass accreted "
                    f"*ex-situ* via mergers that isn't recorded in this "
                    f"galaxy's own SFH (see this validator's class "
                    f"docstring) rather than a genuine data problem.",
                    field="Satellites/GalaxyProperties/StellarMass")

        return report


DEFAULT_VALIDATORS = [SchemaValidator(), IntegrityValidator(),
                     PhysicalValidator()]
