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
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Sequence

import json
import logging
import os

import h5py
import numpy as np

from ..merger_tree_types import MergerTreeError
from .backends import get_backend
from .schema import CatalogueSchema, default_schema
from .validation import DEFAULT_VALIDATORS, ValidationReport, Validator

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

        # native comoving/little_h state of every Haloes/* and *_z0 field
        # above -- read by validation.SchemaValidator's little-h/comoving
        # cross-check against schema.py's declared units (see this stage's
        # units caveat).
        context.meta.setdefault("comoving_little_h", {})["Haloes"] = {
            "comoving": halos.meta.get("comoving"),
            "little_h": halos.meta.get("little_h"),
        }

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

    Galaxy cross-matching: when ``galaxy_backend`` (and ``epoch``, needed
    to actually call it) are given, ``galaxy_backend.galaxy_properties()``
    is called once per satellite (post-sort, so results land in canonical
    row order already) and every field any satellite returned becomes a
    ``Satellites/GalaxyProperties/<field>`` column, ``NaN`` for satellites
    that didn't return that particular field -- matching the
    ``GalaxyBackend`` protocol's "omit, don't fabricate" contract at the
    per-satellite level, but an array output needs *some* fill value once
    stacked across satellites, hence NaN here specifically. Skipped (with
    a warning, not a crash) if ``galaxy_backend`` is given without
    ``epoch``, or before ``HaloExtractStage`` has populated
    ``Satellites/_internal/halo_row``. ``SharkGalaxyID`` (the foreign key
    into the native SHARK table) isn't produced here -- neither backend
    exposes the underlying galaxy's native ID through the
    ``GalaxyBackend`` protocol currently.

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

    def __init__(self, galaxy_backend=None, epoch=None):
        self.galaxy_backend = galaxy_backend
        self.epoch = epoch

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

        n_galaxy_matched = None
        if self.galaxy_backend is not None:
            if self.epoch is None:
                logger.warning(
                    "%s: galaxy_backend given but epoch= wasn't -- "
                    "skipping galaxy cross-matching "
                    "(galaxy_backend.galaxy_properties() needs an Epoch "
                    "to call).", self.name)
            elif "Satellites/_internal/halo_row" not in context.columns:
                logger.warning(
                    "%s: galaxy_backend given but no Satellites/_internal/"
                    "halo_row in context (run HaloExtractStage first) -- "
                    "skipping galaxy cross-matching.", self.name)
            else:
                halo_rows = context.columns["Satellites/_internal/halo_row"]
                per_satellite = [
                    self.galaxy_backend.galaxy_properties(self.epoch, int(row))
                    for row in halo_rows
                ]
                field_names = sorted(
                    {k for props in per_satellite for k in props})
                for field in field_names:
                    values = np.full(n, np.nan)
                    for i, props in enumerate(per_satellite):
                        if field in props:
                            values[i] = props[field]
                    context.columns[
                        f"Satellites/GalaxyProperties/{field}"] = values
                n_galaxy_matched = sum(1 for p in per_satellite if p)
                if n_galaxy_matched == 0:
                    logger.warning(
                        "%s: galaxy_backend matched 0/%d satellites -- no "
                        "Satellites/GalaxyProperties/* columns were "
                        "created at all, which will make any later stage "
                        "requiring them (e.g. StarFormationHistoryStage) "
                        "fail with a confusing 'missing required inputs' "
                        "error instead of this one. For SharkGalaxyBackend "
                        "with match_by='position' (the default), a common "
                        "cause is the halo catalogue's comoving/little_h "
                        "state not matching SHARK's native convention -- "
                        "see Epoch.galaxies_in_halo's docstring.",
                        self.name, n)

                if hasattr(self.galaxy_backend, "native_comoving_little_h"):
                    comoving, little_h = \
                        self.galaxy_backend.native_comoving_little_h(self.epoch)
                    context.meta.setdefault("comoving_little_h", {})[
                        "Satellites/GalaxyProperties"] = {
                            "comoving": comoving, "little_h": little_h}

        context.record_stage(self.name, n_satellites=n, host_id=host_id,
                             n_galaxy_matched=n_galaxy_matched)
        return context


class QualityControlStage(PipelineStage):
    """Runs ``validation.Validator`` subclasses (``validation.
    DEFAULT_VALIDATORS`` -- ``SchemaValidator``, ``IntegrityValidator``,
    ``PhysicalValidator`` -- by default); hard failures (any
    ``ValidationIssue`` with ``severity="error"``) halt the build by
    raising, warnings do not. The combined ``ValidationReport`` (every
    validator's issues, error or warning) is always attached to
    ``context.meta["validation_report"]`` first, even on the failing
    path, so ``WriteStage`` -- or a caller catching the raised error --
    can still inspect exactly what failed; ``WriteStage`` serialises it to
    ``Provenance/validation_report`` when present."""

    name = "quality_control"

    def __init__(self, schema: CatalogueSchema,
                 validators: Sequence["Validator"] = DEFAULT_VALIDATORS):
        self.schema = schema
        self.validators = validators

    def run(self, context: PipelineContext) -> PipelineContext:
        report = ValidationReport()
        for validator in self.validators:
            sub_report = validator.check(context, self.schema)
            report.issues.extend(sub_report.issues)

        context.meta["validation_report"] = report
        context.record_stage(self.name, n_errors=len(report.errors),
                             n_warnings=len(report.warnings))

        if not report.passed:
            preview = "; ".join(
                f"[{issue.check}] {issue.message}"
                for issue in report.errors[:5])
            more = (f" (+{len(report.errors) - 5} more)"
                   if len(report.errors) > 5 else "")
            raise RuntimeError(
                f"Stage '{self.name}': {len(report.errors)} validation "
                f"error(s), catalogue not written: {preview}{more}")
        return context


class WriteStage(PipelineStage):
    """Writes the validated context to
    ``dorcha_catalogue_v{MAJOR.MINOR.PATCH}.h5`` (or the project's naming
    convention) via h5py, with chunking/compression, per-dataset
    unit/description/provenance/is_derived attrs, and an atomic
    write-to-temp-then-rename so a released file is never partially
    written. Refuses to overwrite an existing release version.

    Only ``context.columns`` paths that are actually declared in
    ``schema`` are written -- anything else in ``context.columns`` (e.g.
    ``Satellites/_internal/*`` working state, ``MergerTrees/main_branch``)
    is intentionally not part of the catalogue file. A schema-declared
    field with no matching column is simply omitted (not NaN-filled or
    invented) -- run ``QualityControlStage`` first if you want that
    flagged before write, which is exactly what its ``SchemaValidator``
    "missing declared field" warning already reports.

    Two Metadata/ fields are populated automatically, since they are
    mechanically well-defined from the context itself rather than a
    modelling or provenance decision the pipeline would otherwise have to
    guess: ``n_satellites``/``n_hosts`` (row counts) and
    ``schema_version`` (from ``schema.version``) -- only if not already
    set in ``context.columns`` (so a caller can still override). Every
    other ``Metadata``/``Cosmology``/``SimulationConfig`` attribute
    (``catalogue_version``, ``pipeline_git_commit``, ``random_seed``,
    cosmology parameters, run config, ...) is real provenance/config the
    pipeline has no business inventing -- set it in ``context.columns``
    (e.g. ``context.columns["Metadata/catalogue_version"] = "1.0.0"``)
    before calling this stage if you want it written.

    Symbolic schema shapes (``"n_snap"``, ``"n_bins"``, ...) are not
    cross-checked against a resolved run-specific value here -- whatever
    shape the array in ``context.columns`` already has is what gets
    written; this stage does not (yet) verify it matches what the
    symbolic dimension is supposed to mean elsewhere in the file.
    """

    name = "write"

    def __init__(self, out_path: str, schema: CatalogueSchema):
        self.out_path = out_path
        self.schema = schema

    def run(self, context: PipelineContext) -> PipelineContext:
        if os.path.exists(self.out_path):
            raise FileExistsError(
                f"Stage '{self.name}': refusing to overwrite an existing "
                f"catalogue at '{self.out_path}' -- write to a new path "
                f"(e.g. bump the version) or remove it first.")

        self._fill_automatic_metadata(context)
        # bridge StarFormationHistoryStage's context.meta value into the
        # one schema-declared field that has nowhere else to live -- see
        # its own docstring for why it's meta, not a per-satellite column.
        if ("Snapshots/time_bin_edges_sfh" not in context.columns
                and context.meta.get("time_bin_edges_sfh") is not None):
            context.columns["Snapshots/time_bin_edges_sfh"] = \
                context.meta["time_bin_edges_sfh"]

        tmp_path = self.out_path + ".tmp"
        if os.path.exists(tmp_path):
            os.remove(tmp_path)

        n_written = 0
        try:
            with h5py.File(tmp_path, "w") as f:
                n_written = self._write_fields(f, context)
                self._write_documentation(f)
                self._write_provenance(f, context)
        except Exception:
            if os.path.exists(tmp_path):
                os.remove(tmp_path)
            raise

        os.replace(tmp_path, self.out_path)

        context.record_stage(self.name, out_path=self.out_path,
                             n_fields_written=n_written)
        return context

    def _fill_automatic_metadata(self, context: PipelineContext) -> None:
        satellite_id = context.columns.get(
            "Satellites/Identification/SatelliteID")
        if satellite_id is not None:
            context.columns.setdefault("Metadata/n_satellites",
                                       len(satellite_id))
        host_ids = context.columns.get("Haloes/HostHaloID")
        if host_ids is not None:
            context.columns.setdefault("Metadata/n_hosts", len(host_ids))
        context.columns.setdefault("Metadata/schema_version",
                                   self.schema.version)

    def _write_fields(self, f, context: PipelineContext) -> int:
        n_written = 0
        for path, spec in self.schema.fields.items():
            value = context.columns.get(path)
            if value is None:
                continue
            group = f.require_group(spec.group)
            if spec.is_attribute:
                group.attrs[spec.name] = self._attribute_value(value)
            else:
                self._write_dataset(group, spec, value)
            n_written += 1
        return n_written

    @staticmethod
    def _attribute_value(value):
        if isinstance(value, np.ndarray) and value.size == 1:
            return value.item()
        return value

    #: schema.py dtype string -> numpy dtype, for casting before write.
    _DTYPE = {
        "float64": np.float64, "float32": np.float32,
        "int64": np.int64, "int32": np.int32, "int16": np.int16,
        "int8": np.int8, "bool": np.bool_,
    }

    def _write_dataset(self, group, spec, value):
        if spec.dtype == "str":
            arr = np.asarray(value, dtype=object)
            dset = group.create_dataset(
                spec.name, data=arr,
                dtype=h5py.string_dtype(encoding="utf-8"))
        else:
            arr = np.asarray(value)
            target = self._DTYPE.get(spec.dtype)
            if target is not None and arr.dtype != target:
                try:
                    arr = arr.astype(target)
                except (TypeError, ValueError):
                    logger.warning(
                        "%s: could not cast %s/%s from %s to declared "
                        "dtype '%s'; writing as-is.", self.name, spec.group,
                        spec.name, arr.dtype, spec.dtype)
            kwargs = {}
            if arr.ndim >= 1 and arr.shape[0] > 1:
                kwargs = dict(compression="gzip", compression_opts=4,
                             chunks=True)
            dset = group.create_dataset(spec.name, data=arr, **kwargs)

        dset.attrs["units"] = spec.units
        dset.attrs["description"] = spec.description
        dset.attrs["provenance"] = spec.provenance
        dset.attrs["is_derived"] = spec.is_derived

    def _write_documentation(self, f) -> None:
        f.create_dataset(
            "Documentation/schema.json",
            data=json.dumps(self.schema.to_json_schema_dict(), indent=2))

    def _write_provenance(self, f, context: PipelineContext) -> None:
        f.create_dataset(
            "Provenance/pipeline_stages",
            data=json.dumps(context.provenance, default=str, indent=2))
        f.create_dataset("Provenance/schema_version", data=self.schema.version)
        report = context.meta.get("validation_report")
        if report is not None:
            f.create_dataset("Provenance/validation_report",
                            data=report.to_json())


@dataclass
class HostJob:
    """One host system to build catalogue rows for: an ``Epoch`` plus the
    explicit host/satellite row selection ``HaloExtractStage`` requires
    (this pipeline never guesses host/satellite classification -- see its
    docstring). ``CatalogueBuilder.run()`` takes a sequence of these, one
    per host halo, and concatenates their per-host ``PipelineContext``s
    into one project-wide catalogue (``Haloes/`` gains one row per job;
    ``Satellites/*`` gains one row per satellite across every job, with
    ``SatelliteID`` renumbered globally 0..N_total-1)."""

    epoch: Any
    host_row: int
    satellite_rows: Sequence[int]


#: Extract/cross-match stage names, always run (in this order) before any
#: derived stage, for every HostJob. "tree_extract" is conditionally
#: skipped per-job in CatalogueBuilder._run_job if that job's Epoch has no
#: merger tree -- so it stays out of this fixed prefix.
_CORE_STAGE_NAMES = ("halo_extract", "cross_match")


def _build_derived_stage(name: str, cls: type, job: "HostJob",
                         galaxy_backend, options: Dict[str, Any]
                         ) -> PipelineStage:
    """Construct a ``derived.py`` stage instance for one ``HostJob``.
    Each stage's constructor takes a different mix of (epoch, host_row,
    galaxy_backend, project-config kwargs) -- explicit per-name dispatch
    here rather than reflection/kwarg-injection, matching this codebase's
    preference for explicit code over inspecting constructor signatures.
    Adding a new ``derived.STAGES`` entry needs one new branch here."""
    if name == "halo_properties":
        return cls()
    if name == "orbital_properties":
        return cls(job.epoch, job.host_row)
    if name == "star_formation_history":
        return cls(job.epoch, galaxy_backend, **options)
    if name in ("host_environment", "environment"):
        return cls(job.epoch, job.host_row, **options)
    if name == "observability":
        return cls(job.epoch, **options)
    if name == "dorcha_specific":
        return cls(job.epoch, **options)
    raise ValueError(
        f"CatalogueBuilder doesn't know how to construct derived stage "
        f"'{name}' -- add a branch to _build_derived_stage.")


def _concatenate_contexts(contexts: Sequence[PipelineContext]
                          ) -> PipelineContext:
    """Merge one ``PipelineContext`` per ``HostJob`` into one project-wide
    context: ``Haloes/*``/``Satellites/*``/``MergerTrees/*`` columns are
    row-stacked (ndarray -> ``np.concatenate``, list -> concatenation),
    and ``Satellites/Identification/SatelliteID`` is renumbered globally
    (each job's own local 0..n-1 numbering, from that job's
    ``CrossMatchStage``, would otherwise collide across jobs).

    A column present in some jobs' contexts but not others is backfilled
    (NaN for arrays, ``None`` for list columns) for the jobs missing it,
    rather than raising -- this is expected, not a bug, for
    ``Satellites/GalaxyProperties/*``: ``CrossMatchStage`` only creates a
    column for a field if *some* satellite at that host returned it (see
    its docstring), so one host having no satellite with e.g.
    ``BlackHoleMass`` resolved is a perfectly normal per-host difference,
    not an inconsistency. Any other missing column (e.g. under
    ``Haloes/``, which every job's ``HaloExtractStage`` should always
    populate identically) is still treated as a real bug and raises.
    """
    if not contexts:
        raise ValueError(
            "CatalogueBuilder: no host contexts to concatenate -- "
            "`simulations` was empty.")

    sid_path = "Satellites/Identification/SatelliteID"
    lengths = []
    for i, c in enumerate(contexts):
        if sid_path not in c.columns:
            raise RuntimeError(
                f"CatalogueBuilder: host job {i}'s context has no "
                f"{sid_path} -- CrossMatchStage must run for every job.")
        lengths.append(len(c.columns[sid_path]))

    all_keys = set()
    for c in contexts:
        all_keys |= set(c.columns)

    combined = PipelineContext()
    for key in sorted(all_keys):
        per_job_values = []
        for c, n in zip(contexts, lengths):
            if key in c.columns:
                per_job_values.append(c.columns[key])
            elif key.startswith("Satellites/") or key.startswith("MergerTrees/"):
                per_job_values.append(
                    [None] * n if key.startswith("MergerTrees/")
                    else np.full(n, np.nan))
            else:
                raise RuntimeError(
                    f"CatalogueBuilder: '{key}' is present in some host "
                    f"job's context but not all, and isn't a per-satellite "
                    f"field safe to NaN-backfill -- every job must "
                    f"populate this field identically (a Haloes/*-style "
                    f"stage produced inconsistent output across jobs).")

        if isinstance(per_job_values[0], list):
            combined.columns[key] = sum(
                (list(v) for v in per_job_values), [])
        else:
            combined.columns[key] = np.concatenate(
                [np.asarray(v) for v in per_job_values], axis=0)

    for c in contexts:
        combined.provenance.extend(c.provenance)
        combined.meta.update(c.meta)

    combined.columns[sid_path] = np.arange(sum(lengths), dtype=np.int64)
    return combined


class CatalogueBuilder:
    """Loads a project config, assembles the stage list, and runs it end
    to end against one or more host systems (one ``HostJob`` per host
    halo).

    Parameters
    ----------
    config_path : str
        Path to a project YAML (e.g. ``configs/dorcha.yaml``), shaped
        like that file: ``schema_version``, ``galaxy_backend`` (+
        ``galaxy_backend_options``), ``derived_stages`` (an ordered list
        of names resolving via ``derived.STAGES``), an optional
        ``stage_options`` mapping (``{stage_name: {kwarg: value, ...}}``
        for any derived stage needing more than ``epoch``/``host_row`` --
        see ``_build_derived_stage``), an optional ``output.
        path_template`` (``{version}`` substituted with ``schema.version``,
        used when ``run()`` isn't given an explicit ``out_path``), and an
        optional ``metadata`` mapping (``{"Group/field": value, ...}``,
        e.g. ``{"Metadata/catalogue_version": "1.0.0"}``) merged into the
        combined context before ``QualityControlStage``/``WriteStage`` --
        the only way to populate ``Metadata``/``Cosmology``/
        ``SimulationConfig`` attribute fields, since no stage computes
        those itself (see ``WriteStage``'s docstring for why).
    """

    def __init__(self, config_path: str):
        self.config_path = config_path
        self.config: Dict[str, Any] = {}
        self.schema: Optional[CatalogueSchema] = None
        self.stages: List[str] = []  # resolved derived_stages, in order

        self._galaxy_backend = None
        self._stage_options: Dict[str, Dict[str, Any]] = {}

    def _load_config(self) -> None:
        try:
            import yaml
        except ImportError as exc:
            raise ImportError(
                "CatalogueBuilder requires pyyaml; install the "
                "'catalogue' extra: pip install -e '.[catalogue]'") from exc

        with open(self.config_path, "r") as fh:
            self.config = yaml.safe_load(fh)

        required = {"schema_version", "galaxy_backend", "derived_stages"}
        missing = required - set(self.config)
        if missing:
            raise ValueError(
                f"Project config '{self.config_path}' is missing required "
                f"key(s) {sorted(missing)}.")

        from . import derived  # local: avoids a pipeline<->derived import cycle

        self.schema = default_schema(self.config["schema_version"])
        self._galaxy_backend = get_backend(
            self.config["galaxy_backend"],
            **self.config.get("galaxy_backend_options", {}))

        self.stages = list(self.config["derived_stages"])
        unknown = [n for n in self.stages if n not in derived.STAGES]
        if unknown:
            raise ValueError(
                f"Project config '{self.config_path}' lists unknown "
                f"derived stage(s) {unknown}; available: "
                f"{sorted(derived.STAGES)}.")

        self._stage_options = self.config.get("stage_options", {})

        logger.info("CatalogueBuilder: loaded '%s' (schema v%s, "
                   "galaxy_backend=%s, %d derived stage(s)).",
                   self.config_path, self.schema.version,
                   self.config["galaxy_backend"], len(self.stages))

    def _build_stage(self, name: str, job: HostJob) -> PipelineStage:
        if name == "halo_extract":
            return HaloExtractStage(job.epoch, job.host_row,
                                    job.satellite_rows)
        if name == "tree_extract":
            return TreeExtractStage(job.epoch)
        if name == "cross_match":
            return CrossMatchStage(galaxy_backend=self._galaxy_backend,
                                   epoch=job.epoch)
        from . import derived
        return _build_derived_stage(
            name, derived.STAGES[name], job, self._galaxy_backend,
            self._stage_options.get(name, {}))

    def _run_job(self, job: HostJob) -> PipelineContext:
        context = PipelineContext()
        names = list(_CORE_STAGE_NAMES)
        # HaloExtractStage must run before we can check job.epoch.tree
        # meaningfully as "this job has a tree" vs. "not yet extracted" --
        # but tree presence is an Epoch-level property, not something
        # extraction changes, so it's safe to check up front and insert
        # tree_extract into the ordered name list before running any of it.
        if getattr(job.epoch, "tree", None) is not None:
            names.insert(1, "tree_extract")  # after halo_extract, before cross_match
        names += self.stages

        for name in names:
            stage = self._build_stage(name, job)
            stage.check_inputs(context)
            context = stage.run(context)
        return context

    def run(self, simulations: Sequence[HostJob],
            out_path: Optional[str] = None) -> ValidationReport:
        """Run every stage in order against ``simulations`` and write the
        catalogue. Returns the ``QualityControlStage``'s report; raises on
        any hard validation failure before writing (the report is still
        attached to the failing context's ``meta`` -- see
        ``QualityControlStage``).
        """
        if not self.schema:
            self._load_config()

        contexts = [self._run_job(job) for job in simulations]
        combined = _concatenate_contexts(contexts)

        # project-level provenance/config values (Metadata/Cosmology/
        # SimulationConfig attribute fields -- catalogue_version,
        # pipeline_git_commit, H0, ... -- see WriteStage's docstring for
        # why these aren't guessed anywhere in the pipeline itself) --
        # config["metadata"]: {"Group/field": value, ...}.
        for path, value in self.config.get("metadata", {}).items():
            combined.columns[path] = value

        combined = QualityControlStage(self.schema).run(combined)

        if out_path is None:
            template = self.config.get("output", {}).get(
                "path_template", "catalogue_v{version}.h5")
            out_path = template.format(version=self.schema.version)
        WriteStage(out_path, self.schema).run(combined)

        return combined.meta["validation_report"]

    def run_from_stage(self, stage_name: str, context: PipelineContext,
                       simulations: Sequence[HostJob]) -> PipelineContext:
        """Resume from a cached context at ``stage_name`` -- for iterative
        development so a full rebuild isn't needed after editing one
        derived-quantity stage. ``context`` must be the single-host,
        pre-``stage_name`` context from an earlier ``_run_job`` (or an
        earlier ``run_from_stage`` call) for the *one* job in
        ``simulations`` -- concatenated, multi-host contexts (as produced
        by ``run()``) can't be resumed this way, since every derived stage
        takes exactly one ``epoch``/``host_row``. Runs ``stage_name`` and
        every stage after it in this project's order (extract/cross-match
        stages, then ``derived_stages``); does **not** run
        ``QualityControlStage``/``WriteStage`` -- this is for inspecting
        one stage's output during development, not for producing a
        release file (use ``run()`` for that).
        """
        if not self.schema:
            self._load_config()
        jobs = list(simulations)
        if len(jobs) != 1:
            raise ValueError(
                "run_from_stage resumes a single-host cached context (see "
                "its docstring) -- pass a one-element `simulations` list "
                "matching the job that produced `context`.")
        job = jobs[0]

        all_names = list(_CORE_STAGE_NAMES)
        if getattr(job.epoch, "tree", None) is not None:
            all_names.insert(1, "tree_extract")
        all_names += self.stages

        if stage_name not in all_names:
            raise ValueError(
                f"'{stage_name}' is not one of this project's stages "
                f"({all_names}).")

        for name in all_names[all_names.index(stage_name):]:
            stage = self._build_stage(name, job)
            stage.check_inputs(context)
            context = stage.run(context)
        return context
