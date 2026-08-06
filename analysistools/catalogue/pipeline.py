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

import numpy as np

from ..merger_tree_types import MergerTreeError
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
    that step is project-specific (SHARK vs. hydro).

    This base class stays an unimplemented stub -- ``HaloExtractStage`` and
    ``TreeExtractStage`` below are the two sources actually implemented so
    far (backed by ``Epoch.halos``/``Epoch.track_of``, which already exist
    and work); particle tagging, selection functions, and Rubin
    detectability need machinery ``Epoch`` doesn't have yet, so those stay
    unimplemented."""

    def run(self, context: PipelineContext) -> PipelineContext:
        raise NotImplementedError("Phase 6b.")


def _first_available_field(dataset, candidates: Sequence[str]):
    """Return the first of `candidates` present on `dataset` (a Dataset-like
    object supporting `in`/`[]`), or None if none are. Mirrors
    HaloModel._first_available's candidate-list approach -- native field
    names for the same physical quantity vary across halo-catalogue formats
    (e.g. Vmax isn't one of HaloCatalogue's STANDARD_KEYS)."""
    for name in candidates:
        if name in dataset:
            return np.asarray(dataset[name])
    return None


class HaloExtractStage(ExtractStage):
    """Pulls host + satellite halo properties from ``Epoch.halos`` into the
    context: host-level ``Haloes/*`` fields (one row), and the primary
    (P-type -- no computation beyond what ``HaloCatalogue``/``HaloTools``
    already applied) ``Satellites/Identification``/``Satellites/
    HaloProperties`` fields.

    Which halo is the host and which are satellites is *not* decided here
    (see docs/dorcha_master_catalogue_design.md's isolation/pairing
    algorithm -- not yet implemented as automatic selection). The caller
    passes both explicitly; this keeps host/satellite classification a
    project-level decision rather than a guess baked into the pipeline.

    Units caveat: values are stored exactly as ``epoch.halos`` returns them
    -- i.e. in whatever comoving/little_h state that HaloCatalogue was
    constructed with (see docs/unified_interface.md), which will not
    generally match schema.py's per-field declared units (e.g. `M200c:
    Msun` vs `Position: Mpc/h (comoving)` -- the original design mixes
    conventions per field). Reconciling stored values against each field's
    declared schema unit is not yet implemented anywhere in the pipeline
    (arguably WriteStage's job, Phase 6c) -- until then, construct the
    Epoch's HaloCatalogue with the comoving/little_h state you actually
    want written.
    """

    name = "halo_extract"
    outputs = (
        "Haloes/HostHaloID", "Haloes/M200c", "Haloes/R200c",
        "Haloes/Vmax_host", "Haloes/Position", "Haloes/Velocity",
        "Satellites/Identification/SubhaloID_z0",
        "Satellites/Identification/Snapshot",
        "Satellites/HaloProperties/M200c_z0",
        "Satellites/HaloProperties/R200c_z0",
        "Satellites/HaloProperties/Vmax_z0",
    )

    #: candidate native field names per format for quantities not part of
    #: HaloCatalogue's standardised schema (mass/pos/vel/radius/halo_id/
    #: num_part) -- mirrors HaloModel.vmax()'s candidate list.
    VMAX_CANDIDATES = ("VMax", "vmax", "SubhaloVmax", "Vmax")

    def __init__(self, epoch, host_row: int, satellite_rows: Sequence[int]):
        self.epoch = epoch
        self.host_row = int(host_row)
        self.satellite_rows = np.asarray(satellite_rows, dtype=np.int64)

    def run(self, context: PipelineContext) -> PipelineContext:
        halos = self.epoch.halos
        if halos is None:
            raise RuntimeError(
                f"Stage '{self.name}': this Epoch has no halo catalogue.")
        sats = self.satellite_rows
        if sats.size == 0:
            raise RuntimeError(
                f"Stage '{self.name}': satellite_rows is empty.")

        host = self.host_row
        snapnum = self.epoch.snapnum
        snapnum = -1 if snapnum is None else int(snapnum)

        halo_id = np.asarray(halos["halo_id"])
        mass = np.asarray(halos["mass"])
        radius = np.asarray(halos["radius"])
        pos = np.asarray(halos["pos"])
        vel = np.asarray(halos["vel"]) if "vel" in halos else None
        vmax = _first_available_field(halos, self.VMAX_CANDIDATES)

        context.columns["Haloes/HostHaloID"] = np.asarray([halo_id[host]])
        context.columns["Haloes/M200c"] = np.asarray([mass[host]])
        context.columns["Haloes/R200c"] = np.asarray([radius[host]])
        context.columns["Haloes/Position"] = pos[[host]]
        if vel is not None:
            context.columns["Haloes/Velocity"] = vel[[host]]
        if vmax is not None:
            context.columns["Haloes/Vmax_host"] = vmax[[host]]

        context.columns["Satellites/_internal/halo_row"] = sats
        context.columns["Satellites/_internal/pos_z0"] = pos[sats]
        if vel is not None:
            context.columns["Satellites/_internal/vel_z0"] = vel[sats]

        context.columns["Satellites/Identification/SubhaloID_z0"] = halo_id[sats]
        context.columns["Satellites/Identification/Snapshot"] = np.full(
            sats.size, snapnum, dtype=np.int32)
        context.columns["Satellites/HaloProperties/M200c_z0"] = mass[sats]
        context.columns["Satellites/HaloProperties/R200c_z0"] = radius[sats]
        if vmax is not None:
            context.columns["Satellites/HaloProperties/Vmax_z0"] = vmax[sats]
        else:
            logger.warning(
                "%s: no Vmax field found under any of %s for this %s "
                "catalogue; Satellites/HaloProperties/Vmax_z0 omitted.",
                self.name, self.VMAX_CANDIDATES, halos.fileformat)

        context.record_stage(self.name, host_row=host,
                             n_satellites=int(sats.size))
        return context


class TreeExtractStage(ExtractStage):
    """Looks up each satellite's merger-tree main branch via
    ``Epoch.track_of`` and stores it for later ``DerivedQuantityStage``
    consumption (``HaloPropertiesStage``, ``OrbitalPropertiesStage``, ...,
    which declare ``MergerTrees/main_branch`` as an input). Requires
    ``HaloExtractStage`` to have already run -- reads the halo-catalogue
    row indices it recorded.

    Satellites with no resolvable tree entry (e.g. below the tree's mass/
    particle resolution) get ``MergerTreeID=-1`` and a ``None`` entry in
    ``main_branch`` rather than aborting the whole stage -- flagging that
    as a QC warning belongs to ``IntegrityValidator`` (Phase 6c), not here.
    """

    name = "tree_extract"
    inputs = ("Satellites/_internal/halo_row",)
    outputs = ("Satellites/Identification/MergerTreeID",
              "MergerTrees/main_branch")

    def __init__(self, epoch):
        self.epoch = epoch

    def run(self, context: PipelineContext) -> PipelineContext:
        if self.epoch.tree is None:
            raise RuntimeError(
                f"Stage '{self.name}': this Epoch has no merger tree.")

        rows = context.columns["Satellites/_internal/halo_row"]
        tree_ids = np.full(len(rows), -1, dtype=np.int64)
        main_branch = []
        n_missing = 0

        for i, row in enumerate(rows):
            try:
                track = self.epoch.track_of(index=int(row))
            except MergerTreeError:
                main_branch.append(None)
                n_missing += 1
                continue
            main_branch.append(track)
            tree_ids[i] = int(track.track.halo_id)

        if n_missing:
            logger.warning(
                "%s: %d/%d satellites had no resolvable tree entry "
                "(MergerTreeID=-1, main_branch=None).",
                self.name, n_missing, len(rows))

        context.columns["Satellites/Identification/MergerTreeID"] = tree_ids
        context.columns["MergerTrees/main_branch"] = main_branch
        context.record_stage(self.name, n_satellites=len(rows),
                             n_missing=n_missing)
        return context


class CrossMatchStage(PipelineStage):
    """Resolve SatelliteID <-> HostHaloID <-> MergerTreeID (and, via the
    configured GalaxyBackend, <-> galaxy properties). Thin wrapper around
    ``Epoch``'s existing cross-matching (``particles_in_halo``,
    ``galaxies_in_halo``, ``track_of``) -- this stage does not reimplement
    matching, it only fixes the catalogue's canonical row order
    (SatelliteID-sorted) that every later stage and every subgroup must
    share.

    SatelliteID is assigned sequentially (0..N-1) in SubhaloID_z0-sorted
    order: deterministic and reproducible from the source catalogue, but
    *not* yet a globally-stable ID that survives across catalogue rebuilds
    (that would need a real ID registry -- a later phase, once satellites
    need tracking across builds rather than just within one).

    Which halo is the host is decided by ``HaloExtractStage``, not here --
    this stage only broadcasts the host's own ID to every satellite row as
    ``Satellites/Identification/HostHaloID``.

    Galaxy cross-matching (``SharkGalaxyID``, ``GalaxyProperties/*``) is
    skipped when ``galaxy_backend=None`` (the default): both
    ``GalaxyBackend`` implementations are themselves still Phase 6b stubs
    (``galaxy_properties()`` raises ``NotImplementedError``), so there is
    nothing to actually call yet. Pass one in once that's implemented.

    Row-permutation convention: every ``context.columns`` key under
    ``'Satellites/'`` or ``'MergerTrees/'`` is treated as satellite-indexed
    (length N, reordered together, whether an ndarray or a plain list --
    see ``TreeExtractStage``'s ``main_branch``); ``'Haloes/'`` keys are
    host-level and left untouched.
    """

    name = "cross_match"
    inputs = ("Satellites/Identification/SubhaloID_z0", "Haloes/HostHaloID")
    outputs = ("Satellites/Identification/SatelliteID",
              "Satellites/Identification/HostHaloID")

    def __init__(self, galaxy_backend=None):
        self.galaxy_backend = galaxy_backend

    def run(self, context: PipelineContext) -> PipelineContext:
        subhalo_id = np.asarray(
            context.columns["Satellites/Identification/SubhaloID_z0"])
        n = subhalo_id.size
        order = np.argsort(subhalo_id, kind="stable")

        for key in list(context.columns):
            if not (key.startswith("Satellites/")
                    or key.startswith("MergerTrees/")):
                continue
            value = context.columns[key]
            if isinstance(value, np.ndarray) and len(value) == n:
                context.columns[key] = value[order]
            elif isinstance(value, list) and len(value) == n:
                context.columns[key] = [value[i] for i in order]
            # else: not satellite-indexed (unexpected length) -- leave as-is

        context.columns["Satellites/Identification/SatelliteID"] = \
            np.arange(n, dtype=np.int64)

        host_id = int(np.asarray(context.columns["Haloes/HostHaloID"])[0])
        context.columns["Satellites/Identification/HostHaloID"] = \
            np.full(n, host_id, dtype=np.int64)

        if self.galaxy_backend is not None:
            logger.warning(
                "%s: galaxy_backend given but SharkGalaxyID/"
                "GalaxyProperties cross-matching isn't implemented yet "
                "(galaxy_backend.galaxy_properties() is itself still a "
                "Phase 6b stub) -- skipping.", self.name)

        context.record_stage(self.name, n_satellites=n, host_id=host_id)
        return context


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
