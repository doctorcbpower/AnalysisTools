"""
Tests for analysistools.api.catalogue.CatalogueDataset (Phase 6a).

catalogue.pipeline.WriteStage isn't implemented yet (Phase 6c), so there's
no real catalogue file to read back -- these build a minimal but realistic
synthetic HDF5 file directly (matching the on-disk layout CatalogueDataset
expects: Satellites/<Group>/<Field> datasets, optional Metadata/Cosmology
attrs) rather than depending on unfinished pipeline code.
"""
import h5py
import numpy as np
import pytest

from analysistools.api.catalogue import CatalogueDataset


def _write_catalogue(path, n=5, with_metadata=True, with_cosmology=True,
                     with_duplicate=False, empty_satellites=False):
    with h5py.File(path, "w") as f:
        sats = f.create_group("Satellites")
        if not empty_satellites:
            ident = sats.create_group("Identification")
            d = ident.create_dataset("SatelliteID", data=np.arange(n, dtype=np.int64))
            d.attrs["description"] = "Unique satellite identifier"

            halo = sats.create_group("HaloProperties")
            d = halo.create_dataset("Mpeak", data=np.linspace(1e10, 1e12, n))
            d.attrs["units"] = "Msun"
            d.attrs["is_derived"] = True
            d.attrs["description"] = "Peak halo mass"

            gal = sats.create_group("GalaxyProperties")
            d = gal.create_dataset("StellarMass", data=np.linspace(1e8, 1e10, n))
            d.attrs["units"] = "Msun"
            d.attrs["is_derived"] = False

            if with_duplicate:
                # a second field literally named "Mpeak" under a different
                # group -- exercises the duplicate-field-name branch
                other = sats.create_group("OtherGroup")
                other.create_dataset("Mpeak", data=np.full(n, -1.0))

        if with_metadata:
            meta = f.create_group("Metadata")
            meta.attrs["catalogue_version"] = "v1.2.3"
            meta.attrs["schema_version"] = "1.0"

        if with_cosmology:
            cosmo = f.create_group("Cosmology")
            cosmo.attrs["H0"] = 70.0

    return path


# ---------------------------------------------------------------------------
# Field flattening
# ---------------------------------------------------------------------------

def test_load_flattens_satellite_fields(tmp_path):
    path = _write_catalogue(tmp_path / "cat.h5")
    cat = CatalogueDataset(str(path))

    np.testing.assert_array_equal(cat["SatelliteID"], np.arange(5))
    np.testing.assert_allclose(cat["Mpeak"], np.linspace(1e10, 1e12, 5))
    np.testing.assert_allclose(cat["StellarMass"], np.linspace(1e8, 1e10, 5))


def test_load_is_lazy(tmp_path):
    path = _write_catalogue(tmp_path / "cat.h5")
    cat = CatalogueDataset(str(path))
    assert cat._loaded is False
    cat["Mpeak"]
    assert cat._loaded is True


def test_kind_and_fileformat(tmp_path):
    path = _write_catalogue(tmp_path / "cat.h5")
    cat = CatalogueDataset(str(path))
    assert cat.kind == "satellites"
    assert cat.fileformat == "HDF5-catalogue"


def test_missing_satellites_group_raises(tmp_path):
    path = tmp_path / "not_a_catalogue.h5"
    with h5py.File(path, "w") as f:
        f.create_group("SomethingElse")

    cat = CatalogueDataset(str(path))
    with pytest.raises(ValueError, match="Satellites"):
        cat.preload()


def test_duplicate_field_name_keeps_one_and_does_not_crash(tmp_path):
    path = _write_catalogue(tmp_path / "cat.h5", with_duplicate=True)
    cat = CatalogueDataset(str(path))
    # whichever of the two same-named datasets h5py visits first is kept;
    # don't over-specify traversal order, just that it's one of the two
    # and the load didn't crash or silently merge/corrupt them.
    valid_values = (np.linspace(1e10, 1e12, 5), np.full(5, -1.0))
    assert any(np.allclose(cat["Mpeak"], v) for v in valid_values)


def test_empty_satellites_group(tmp_path):
    path = _write_catalogue(tmp_path / "cat.h5", empty_satellites=True,
                            with_metadata=False, with_cosmology=False)
    cat = CatalogueDataset(str(path))
    cat.preload()
    assert cat.meta["n_satellites"] == 0
    assert cat.fields == []


# ---------------------------------------------------------------------------
# Metadata
# ---------------------------------------------------------------------------

def test_metadata_and_cosmology_populated(tmp_path):
    path = _write_catalogue(tmp_path / "cat.h5")
    cat = CatalogueDataset(str(path))
    cat.preload()

    assert cat.meta["catalogue_version"] == "v1.2.3"
    assert cat.meta["schema_version"] == "1.0"
    assert cat.meta["n_satellites"] == 5
    assert cat.meta["h0"] == pytest.approx(0.7)


def test_metadata_and_cosmology_absent_default_to_none(tmp_path):
    path = _write_catalogue(tmp_path / "cat.h5", with_metadata=False,
                            with_cosmology=False)
    cat = CatalogueDataset(str(path))
    cat.preload()

    assert cat.meta["catalogue_version"] is None
    assert cat.meta["schema_version"] is None
    assert cat.meta.get("h0") is None


def test_units_metadata_only_includes_fields_with_units(tmp_path):
    path = _write_catalogue(tmp_path / "cat.h5")
    cat = CatalogueDataset(str(path))
    cat.preload()

    units = cat.meta["units"]
    assert units["Mpeak"] == "Msun"
    assert units["StellarMass"] == "Msun"
    # SatelliteID has no "units" attr in the fixture -- must not appear
    assert "SatelliteID" not in units


# ---------------------------------------------------------------------------
# field_info
# ---------------------------------------------------------------------------

def test_field_info_returns_attrs(tmp_path):
    path = _write_catalogue(tmp_path / "cat.h5")
    cat = CatalogueDataset(str(path))

    info = cat.field_info("Mpeak")
    assert info["units"] == "Msun"
    assert info["is_derived"] == True  # noqa: E712 (h5py attr, not a Python bool)
    assert info["description"] == "Peak halo mass"


def test_field_info_unknown_field_raises(tmp_path):
    path = _write_catalogue(tmp_path / "cat.h5")
    cat = CatalogueDataset(str(path))
    with pytest.raises(KeyError, match="NoSuchField"):
        cat.field_info("NoSuchField")
