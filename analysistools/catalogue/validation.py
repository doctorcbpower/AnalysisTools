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
from dataclasses import dataclass, field
from typing import Any, Dict, List

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
        raise NotImplementedError(
            "Phase 6c: serialise for Provenance/validation_report.")


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
    unexpected NaN/inf outside an allowed-missingness mask; particle-tag
    CSR offsets monotonically non-decreasing."""

    name = "integrity"

    def check(self, context, schema) -> ValidationReport:
        raise NotImplementedError("Phase 6c.")


class PhysicalValidator(Validator):
    """Mpeak >= M200c_z0; RedshiftInfall monotonic along accretion
    history; OrbitalApocentre >= OrbitalPericentre; HeliocentricDistance >
    0; CompletenessWeight >= 1; SFH bins sum to StellarMass within
    tolerance."""

    name = "physical"

    def check(self, context, schema) -> ValidationReport:
        raise NotImplementedError("Phase 6c.")


DEFAULT_VALIDATORS = [SchemaValidator(), IntegrityValidator(),
                     PhysicalValidator()]
