#!/usr/bin/env python3
"""
analysistools.catalogue
------------------------
Write-side pipeline for building a project's master science catalogue: one
flat, versioned, self-describing HDF5 file (one row per satellite/subhalo,
expensive quantities computed once and stored permanently) that becomes the
sole input to science analysis.

See DEVELOPMENT.md Phase 6, docs/master_catalogue.md (mapping onto this
codebase's conventions), and docs/dorcha_master_catalogue_design.md (the
full per-field schema, derived-quantity catalogue, and validation checks).

Status: design skeleton (2026-08-06). Class shapes and docstrings are
final-ish; method bodies raise NotImplementedError until Phase 6a/6b/6c land.

Read-side note: once written, a catalogue re-enters the unified interface
via ``at.load(path, kind="satellites")`` -> ``CatalogueDataset``
(analysistools/api/catalogue.py), so nothing in this subpackage duplicates
Dataset/Simulation/Epoch machinery -- it only ever *produces* the file that
adapter reads.
"""
from __future__ import annotations

from .backends import GalaxyBackend, HydroGalaxyBackend, SharkGalaxyBackend
from .pipeline import CatalogueBuilder, PipelineStage
from .schema import (CatalogueSchema, FieldSpec, available_schema_versions,
                     default_schema)
from .validation import ValidationReport, Validator

__all__ = [
    "CatalogueBuilder",
    "PipelineStage",
    "CatalogueSchema",
    "FieldSpec",
    "default_schema",
    "available_schema_versions",
    "GalaxyBackend",
    "SharkGalaxyBackend",
    "HydroGalaxyBackend",
    "Validator",
    "ValidationReport",
]
