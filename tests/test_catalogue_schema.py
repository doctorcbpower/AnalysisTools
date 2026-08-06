"""
Tests for analysistools.catalogue.schema (Phase 6a).

Covers: loading the bundled schema_v1.0.yaml, round-tripping through
to_yaml/from_yaml, duplicate-field rejection, and validate_columns()'s
missing/undeclared detection -- the checks
validation.SchemaValidator will lean on in Phase 6c.
"""
import os

import pytest

from analysistools.catalogue.schema import (CatalogueSchema, FieldSpec,
                                            available_schema_versions,
                                            default_schema)


def test_available_schema_versions_includes_1_0():
    assert "1.0" in available_schema_versions()


def test_default_schema_loads_and_is_populated():
    schema = default_schema("1.0")
    assert schema.version == "1.0"
    # sanity floor -- the transcribed table has ~85 fields across ~13 groups
    assert len(schema) > 60
    assert len(schema.groups) >= 10


def test_default_schema_unknown_version_raises():
    with pytest.raises(FileNotFoundError):
        default_schema("99.9")


@pytest.mark.parametrize("path,expected_units,expected_derived", [
    ("Satellites/HaloProperties/Mpeak", "Msun", True),
    ("Satellites/Identification/SatelliteID", "dimensionless", False),
    ("Satellites/GalaxyProperties/StellarMass", "Msun", False),
    ("Cosmology/H0", "km/s/Mpc", False),
])
def test_known_fields_have_expected_metadata(path, expected_units,
                                              expected_derived):
    schema = default_schema("1.0")
    assert path in schema
    spec = schema[path]
    assert spec.units == expected_units
    assert spec.is_derived is expected_derived


def test_metadata_and_cosmology_fields_are_attributes():
    schema = default_schema("1.0")
    for path in ("Metadata/catalogue_version", "Cosmology/H0",
                 "SimulationConfig/box_size"):
        assert schema[path].is_attribute is True
    assert schema["Satellites/HaloProperties/Mpeak"].is_attribute is False


def test_vector_fields_carry_shape():
    schema = default_schema("1.0")
    assert schema["Haloes/Position"].shape == (3,)
    assert schema["Satellites/GalaxyProperties/Luminosity_ugriz"].shape == (5,)
    # symbolic dims are resolved at write time, not schema-definition time
    assert schema["Satellites/GalaxyProperties/SFH"].shape == ("n_bins",)


def test_by_group_is_exact_not_prefix_match():
    schema = default_schema("1.0")
    halo_fields = {f.name for f in
                   schema.by_group("Satellites/HaloProperties")}
    galaxy_fields = {f.name for f in
                     schema.by_group("Satellites/GalaxyProperties")}
    assert "Mpeak" in halo_fields
    assert "Mpeak" not in galaxy_fields
    assert halo_fields.isdisjoint(galaxy_fields)


def test_add_duplicate_path_raises():
    schema = CatalogueSchema(version="test")
    spec = FieldSpec(name="Mass", group="Satellites/Test", dtype="float64")
    schema.add(spec)
    with pytest.raises(ValueError):
        schema.add(spec)


def test_yaml_round_trip(tmp_path):
    original = default_schema("1.0")
    out_path = tmp_path / "schema_roundtrip.yaml"
    original.to_yaml(str(out_path))
    assert out_path.exists()

    reloaded = CatalogueSchema.from_yaml(str(out_path))
    assert reloaded.version == original.version
    assert len(reloaded) == len(original)
    assert set(reloaded.fields) == set(original.fields)

    # spot-check one full record survives the round trip intact
    orig_spec = original["Satellites/HaloProperties/OrbitalEccentricity"]
    new_spec = reloaded["Satellites/HaloProperties/OrbitalEccentricity"]
    assert new_spec == orig_spec


def test_validate_columns_reports_missing_and_undeclared():
    schema = default_schema("1.0")
    declared = {f.name for f in
                schema.by_group("Satellites/HaloProperties")}
    assert "Mpeak" in declared

    have = (declared - {"Mpeak"}) | {"SomeUndeclaredField"}
    problems = schema.validate_columns("Satellites/HaloProperties", have)

    assert any("missing declared field" in p and "Mpeak" in p
               for p in problems)
    assert any("undeclared field" in p and "SomeUndeclaredField" in p
               for p in problems)


def test_validate_columns_clean_set_has_no_problems():
    schema = default_schema("1.0")
    declared = {f.name for f in
                schema.by_group("Satellites/Identification")}
    assert schema.validate_columns("Satellites/Identification", declared) == []
