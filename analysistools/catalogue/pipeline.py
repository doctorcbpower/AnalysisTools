#!/usr/bin/env python3
"""
analysistools.catalogue.pipeline
-----------------------------------
PipelineStage / CatalogueBuilder: orchestrates Extract -> Cross-match ->
Compute derived quantities -> Quality control -> Write, on top of
``api.Simulation``/``Epoch`` (Extract + Cross-match are already provided by
those classes -- see docs/master_catalogue.md).

A project is a YAML config (``configs/<project>.yaml``) naming a
``galaxy_backend`` and an ordered list of ``derived_stages``; it is not a
subclass. ``CatalogueBuilder`` is the only object permitted to write the
master file (via ``analysistools.catalogue`` -- there is deliberately no
public "just write this array" escape hatch, so every field that ends up
in the catalogue went through schema validation).
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any, Dict, List, Optional, Sequence

import logging

from .backends import get_backend
from .schema import CatalogueSchema
from .validation import ValidationReport

logger = logging.getLogger(__name__)


class PipelineContext:
    """Shared, append-only state passed between stages. Each stage reads
    the columns it declares in ``inputs`` and adds the columns it declares
    in ``outputs``; ``CatalogueBuilder`` checks declared inputs are already
    present before running a stage, so a missing prerequisite fails at
    pipeline-assembly time rather than with a mid-run ``KeyError``.
    """

    def __init__(self):
        self.columns: Dict[str, Any] = {}     # full_field_path -> ndarray
        self.meta: Dict[str, Any] = {}
        self.provenance: List[Dict[str, Any]] = []

    def record_stage(self, name: str, **info: Any) -> None:
        self.provenance.append({"stage": name, **info})


class PipelineStage(ABC):
    """One stage of catalogue construction.

    Subclasses declare the field paths they require and produce so
    ``CatalogueBuilder`` can validate stage ordering statically.
    """

    name: str = "stage"
    inputs: Sequence[str] = ()
    outputs: Sequence[str] = ()

    @abstractmethod
    def run(self, context: PipelineContext) -> PipelineContext:
        ...

    def check_inputs(self, context: PipelineContext) -> None:
        missing = [f for f in self.inputs if f not in context.columns]
        if missing:
            raise RuntimeError(
                f"Stage '{self.name}' is missing required inputs "
                f"{missing}; check derived_stages ordering in the project "
                f"config.")


class ExtractStage(PipelineStage):
    """Pull one field group from a ``Simulation``/``Epoch`` into the
    context. One subclass per source: haloes, merger trees, particle
    tagging, selection functions, Rubin detectability tables. Galaxy
    properties are handled separately via ``backends.GalaxyBackend`` since
    that step is project-specific (SHARK vs. hydro)."""

    def run(self, context: PipelineContext) -> PipelineContext:
        raise NotImplementedError("Phase 6b.")


class CrossMatchStage(PipelineStage):
    """Resolve SatelliteID <-> HostHaloID <-> MergerTreeID (and, via the
    configured GalaxyBackend, <-> galaxy properties). Thin wrapper around
    ``Epoch``'s existing cross-matching (``particles_in_halo``,
    ``galaxies_in_halo``, ``track_of``) -- this stage does not reimplement
    matching, it only fixes the catalogue's canonical row order
    (SatelliteID-sorted) that every later stage and every subgroup must
    share."""

    def run(self, context: PipelineContext) -> PipelineContext:
        raise NotImplementedError("Phase 6b.")


class QualityControlStage(PipelineStage):
    """Runs ``validation.Validator`` subclasses; hard failures halt the
    build, warnings are attached to ``Provenance/validation_report``."""

    def __init__(self, validators: Sequence["Validator"] = ()):
        self.validators = validators

    def run(self, context: PipelineContext) -> PipelineContext:
        raise NotImplementedError("Phase 6c.")


class WriteStage(PipelineStage):
    """Writes the validated context to
    ``dorcha_catalogue_v{MAJOR.MINOR.PATCH}.h5`` (or the project's naming
    convention) via h5py, with chunking/compression, per-dataset
    unit/description/provenance/is_derived attrs, and an atomic
    write-to-temp-then-rename so a released file is never partially
    written. Refuses to overwrite an existing release version."""

    def __init__(self, out_path: str, schema: CatalogueSchema):
        self.out_path = out_path
        self.schema = schema

    def run(self, context: PipelineContext) -> PipelineContext:
        raise NotImplementedError("Phase 6c.")


class CatalogueBuilder:
    """Loads a project config, assembles the stage list, and runs it end
    to end against one or more ``Simulation`` objects (one per host halo,
    typically).

    Parameters
    ----------
    config_path : str
        Path to a project YAML (e.g. ``configs/dorcha.yaml``).
    """

    def __init__(self, config_path: str):
        self.config_path = config_path
        self.config: Dict[str, Any] = {}
        self.schema: Optional[CatalogueSchema] = None
        self.stages: List[PipelineStage] = []

    def _load_config(self) -> None:
        raise NotImplementedError(
            "Phase 6a: load YAML, resolve galaxy_backend via "
            "backends.get_backend(), build self.stages from "
            "derived_stages, load self.schema from schema_version.")

    def run(self, simulations: Sequence[Any],
            out_path: Optional[str] = None) -> ValidationReport:
        """Run every stage in order against ``simulations`` and write the
        catalogue. Returns the QualityControlStage's report; raises on any
        hard validation failure before writing.
        """
        raise NotImplementedError("Phase 6b/6c.")

    def run_from_stage(self, stage_name: str, context: PipelineContext,
                       simulations: Sequence[Any]) -> PipelineContext:
        """Resume from a cached context at ``stage_name`` -- for iterative
        development so a full rebuild isn't needed after editing one
        derived-quantity stage."""
        raise NotImplementedError("Phase 6b.")
