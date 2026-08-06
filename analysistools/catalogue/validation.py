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


class SchemaValidator(Validator):
    """Every field in the schema is present with declared dtype and units
    attribute; no undeclared datasets exist; units strings are recognised
    (matches the declared unit family, e.g. a "mass" field must resolve to
    a mass-like unit string)."""

    name = "schema"

    def check(self, context, schema) -> ValidationReport:
        raise NotImplementedError("Phase 6c.")


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
