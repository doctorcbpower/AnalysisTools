#!/usr/bin/env python3
"""
analysistools.catalogue.schema
-------------------------------
FieldSpec / CatalogueSchema: the single source of truth for every dataset
in a master catalogue -- name, dtype, units, description, provenance, and
whether it is primary (extracted directly from a source) or derived
(computed by a pipeline stage from other stored fields).

This mirrors the ``FIELD_ALIASES`` pattern in ``api/dataset.py`` (a
registry, not scattered constants) but is the *write-side* schema: it is
what ``pipeline.CatalogueBuilder`` validates column sets against before
writing, and what gets serialised into ``Documentation/schema.json`` inside
the catalogue file itself (and mirrored as HDF5 attributes on every
dataset -- see docs/dorcha_master_catalogue_design.md section 2.3).

Units convention (DEVELOPMENT.md section 7.3): plain strings + explicit
conversion factors, no astropy/unyt dependency, matching ``meta["units"]``
elsewhere in the package.

The full schema is bundled as package data:
``analysistools/catalogue/configs/schema_v{X.Y}.yaml`` -- load it with
:func:`default_schema` rather than hand-building a ``CatalogueSchema``.
Project configs (``configs/<project>.yaml`` at the repo root, e.g.
``dorcha.yaml``) reference a schema *version*, not a path, so a schema
bump is a one-line change in the project config.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

import logging

logger = logging.getLogger(__name__)

#: Directory holding bundled schema_v{X.Y}.yaml files (package data).
_CONFIGS_DIR = Path(__file__).resolve().parent / "configs"


@dataclass
class FieldSpec:
    """One dataset (or attribute) in the catalogue schema.

    Parameters
    ----------
    name : str
        Dataset/attribute name, e.g. ``"Mpeak"``.
    group : str
        HDF5 group path it lives under, e.g. ``"Satellites/HaloProperties"``.
    dtype : str
        numpy dtype string, e.g. ``"float64"``, ``"int64"``, ``"bool"``, or
        ``"str"`` for HDF5 string datasets/attributes.
    units : str
        Plain-string units, e.g. ``"Msun/h"``, ``"kpc"``, ``"dimensionless"``.
    description : str
        One-line human-readable description.
    provenance : str
        Where this field comes from -- a source name (e.g. "halo
        catalogue") for primary fields, or an algorithm name/version (e.g.
        "orbit-fitting algorithm v1") for derived fields.
    is_derived : bool
        False for fields read directly from a source; True for fields
        computed by a ``derived.DerivedQuantityStage``.
    is_attribute : bool
        True for values stored as an HDF5 group attribute (Metadata/,
        Cosmology/, SimulationConfig/) rather than a dataset. Attributes
        are always scalar; ``shape`` is ignored when this is True.
    shape : tuple of (int or str), optional
        Trailing shape beyond the row axis, e.g. ``(3,)`` for a position
        vector, ``(5,)`` for ugriz bands. A string entry (``"n_snap"``,
        ``"n_bins"``, ``"n_bins_plus_1"``, ``"n_total_tagged"``,
        ``"n_plus_1"``) is a symbolic dimension resolved at write time
        against run-specific context (snapshot count, SFH binning,
        total tagged-particle count); see pipeline.WriteStage. Empty for
        scalar fields.
    """

    name: str
    group: str
    dtype: str
    units: str = "dimensionless"
    description: str = ""
    provenance: str = ""
    is_derived: bool = False
    is_attribute: bool = False
    shape: tuple = ()

    @property
    def path(self) -> str:
        """Full HDF5 path, e.g. 'Satellites/HaloProperties/Mpeak'."""
        return f"{self.group}/{self.name}"

    def to_dict(self) -> Dict[str, Any]:
        d = {
            "name": self.name, "group": self.group, "dtype": self.dtype,
            "units": self.units, "description": self.description,
            "provenance": self.provenance, "is_derived": self.is_derived,
            "is_attribute": self.is_attribute,
        }
        if self.shape:
            d["shape"] = list(self.shape)
        return d


@dataclass
class CatalogueSchema:
    """The full schema for one catalogue version: an ordered collection of
    :class:`FieldSpec`, keyed by full HDF5 path.

    Loaded from / dumped to YAML (``configs/schema_v{X.Y}.yaml``) so the
    schema is versioned in the same way as the pipeline code, and so
    ``Documentation/schema.json`` inside a built catalogue is a direct
    serialisation of this object -- not a hand-maintained duplicate.
    """

    version: str
    fields: Dict[str, FieldSpec] = field(default_factory=dict)

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------

    def add(self, spec: FieldSpec) -> None:
        if spec.path in self.fields:
            raise ValueError(f"Duplicate field path in schema: {spec.path}")
        self.fields[spec.path] = spec

    def by_group(self, group: str) -> List[FieldSpec]:
        """All fields declared directly under ``group`` (exact match, not
        a prefix -- 'Satellites/HaloProperties' does not also return
        'Satellites/GalaxyProperties' fields)."""
        return [f for f in self.fields.values() if f.group == group]

    @property
    def groups(self) -> List[str]:
        return sorted({f.group for f in self.fields.values()})

    def __len__(self) -> int:
        return len(self.fields)

    def __contains__(self, path: str) -> bool:
        return path in self.fields

    def __getitem__(self, path: str) -> FieldSpec:
        return self.fields[path]

    # ------------------------------------------------------------------
    # YAML I/O
    # ------------------------------------------------------------------

    @classmethod
    def from_yaml(cls, path: str) -> "CatalogueSchema":
        """Load a schema from a YAML file shaped like
        ``configs/schema_v1.0.yaml`` (a ``version`` key plus a ``fields``
        list of field dicts)."""
        try:
            import yaml
        except ImportError as exc:
            raise ImportError(
                "CatalogueSchema.from_yaml requires pyyaml; install the "
                "'catalogue' extra: pip install -e '.[catalogue]'") from exc

        with open(path, "r") as fh:
            data = yaml.safe_load(fh)

        if "version" not in data:
            raise ValueError(f"Schema file '{path}' is missing a 'version' key.")

        schema = cls(version=str(data["version"]))
        for entry in data.get("fields", []):
            missing = {"name", "group", "dtype"} - set(entry)
            if missing:
                raise ValueError(
                    f"Schema file '{path}': field entry {entry!r} is "
                    f"missing required key(s) {sorted(missing)}.")
            spec = FieldSpec(
                name=entry["name"],
                group=entry["group"],
                dtype=entry["dtype"],
                units=entry.get("units", "dimensionless"),
                description=entry.get("description", ""),
                provenance=entry.get("provenance", ""),
                is_derived=bool(entry.get("is_derived", False)),
                is_attribute=bool(entry.get("is_attribute", False)),
                shape=tuple(entry.get("shape", ())),
            )
            schema.add(spec)

        logger.info("Loaded catalogue schema v%s from %s (%d fields, "
                    "%d groups).", schema.version, path, len(schema),
                    len(schema.groups))
        return schema

    def to_yaml(self, path: str) -> None:
        """Dump this schema to a YAML file in the same shape
        :meth:`from_yaml` reads. Field order follows insertion order
        (Python dict), i.e. the order fields were added."""
        try:
            import yaml
        except ImportError as exc:
            raise ImportError(
                "CatalogueSchema.to_yaml requires pyyaml; install the "
                "'catalogue' extra: pip install -e '.[catalogue]'") from exc

        data = {
            "version": self.version,
            "fields": [spec.to_dict() for spec in self.fields.values()],
        }
        with open(path, "w") as fh:
            yaml.safe_dump(data, fh, sort_keys=False, default_flow_style=False)

    def to_json_schema_dict(self) -> Dict[str, Any]:
        """The dict written verbatim (via ``json.dumps``) into
        ``Documentation/schema.json`` inside a built catalogue."""
        return {
            "version": self.version,
            "fields": [spec.to_dict() for spec in self.fields.values()],
        }

    # ------------------------------------------------------------------
    # Validation support (used by validation.SchemaValidator)
    # ------------------------------------------------------------------

    def validate_columns(self, group: str,
                          columns: Sequence[str]) -> List[str]:
        """Compare an in-memory column set produced by a pipeline stage
        for ``group`` against this schema's declared fields for that
        group.

        Returns
        -------
        list of str
            Human-readable problem descriptions: one entry per field
            declared in the schema but absent from ``columns`` ("missing"),
            and one per column present but not declared in the schema
            ("undeclared"). Empty list means the column set exactly
            matches the schema for this group.
        """
        declared = {f.name for f in self.by_group(group)}
        have = set(columns)

        problems = [
            f"missing declared field '{group}/{name}'"
            for name in sorted(declared - have)
        ]
        problems += [
            f"undeclared field '{group}/{name}' present in columns "
            f"(not in schema v{self.version})"
            for name in sorted(have - declared)
        ]
        return problems


# ---------------------------------------------------------------------------
# Bundled schema versions
# ---------------------------------------------------------------------------

def available_schema_versions() -> List[str]:
    """Versions bundled with this package
    (``analysistools/catalogue/configs/schema_v*.yaml``)."""
    if not _CONFIGS_DIR.is_dir():
        return []
    versions = []
    for p in sorted(_CONFIGS_DIR.glob("schema_v*.yaml")):
        versions.append(p.stem[len("schema_v"):])
    return versions


def default_schema(version: str = "1.0") -> CatalogueSchema:
    """Load the schema bundled with this package for the given version.

    Parameters
    ----------
    version : str
        e.g. ``"1.0"`` -> ``analysistools/catalogue/configs/schema_v1.0.yaml``.
    """
    path = _CONFIGS_DIR / f"schema_v{version}.yaml"
    if not path.exists():
        raise FileNotFoundError(
            f"No bundled schema for version '{version}' at {path}. "
            f"Available versions: {available_schema_versions()}")
    return CatalogueSchema.from_yaml(str(path))
