"""
Tests for analysistools.catalogue.backends.

Covers: get_backend() dispatch, kwargs passthrough, and the unknown-name
error. Both backends' galaxy_properties() are implemented (Phase 6b) and
tested separately: test_catalogue_backends_shark.py and
test_catalogue_backends_hydro.py.
"""
import pytest

from analysistools.catalogue.backends import (
    BACKENDS, HydroGalaxyBackend, SharkGalaxyBackend, get_backend,
)


def test_backends_registry_has_expected_entries():
    assert set(BACKENDS) == {"shark", "hydro"}
    assert BACKENDS["shark"] is SharkGalaxyBackend
    assert BACKENDS["hydro"] is HydroGalaxyBackend


def test_get_backend_shark_defaults():
    backend = get_backend("shark")
    assert isinstance(backend, SharkGalaxyBackend)
    assert backend.name == "shark"
    assert backend.match_by == "position"
    assert backend.r_scale == 1.0


def test_get_backend_hydro_defaults():
    backend = get_backend("hydro")
    assert isinstance(backend, HydroGalaxyBackend)
    assert backend.name == "hydro"
    assert backend.r_scale == 1.0
    assert backend.young_star_age_threshold == 0.1


def test_get_backend_passes_through_kwargs():
    shark = get_backend("shark", match_by="id", r_scale=2.0)
    assert shark.match_by == "id"
    assert shark.r_scale == 2.0

    hydro = get_backend("hydro", r_scale=0.5, young_star_age_threshold=0.2)
    assert hydro.r_scale == 0.5
    assert hydro.young_star_age_threshold == 0.2


def test_get_backend_unknown_name_raises():
    with pytest.raises(ValueError, match="Unknown galaxy_backend"):
        get_backend("not_a_real_backend")
