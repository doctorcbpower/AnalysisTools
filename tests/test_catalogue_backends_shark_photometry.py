"""
Tests for analysistools.catalogue.backends.SharkGalaxyBackend's
compute_photometry=True path (Phase 6b, FSPS integration).

Uses a fake PhotometryPipeline factory (injected via
photometry_pipeline_factory) instead of the real shark.photometry.
PhotometryPipeline -- avoids depending on a real python-fsps install with
SPS_HOME configured, which this dev environment doesn't have. Reuses the
_FakeModel/_FakeGalaxyView/_FakeEpoch fixtures' shape from
test_catalogue_backends_shark_sfh.py.
"""
import numpy as np
import pytest

from analysistools.catalogue.backends import (
    DEFAULT_PHOTOMETRY_BANDS, SOLAR_ABSOLUTE_MAGNITUDE, SharkGalaxyBackend,
)


class _FakeModel:
    pass  # only identity (id()) matters for the photometry cache key


class _FakeGalaxyView:
    def __init__(self, columns, model=None, index=None):
        self.columns = {k: np.asarray(v) for k, v in columns.items()}
        n = len(next(iter(self.columns.values()))) if self.columns else 0
        self.model = model
        self.index = (np.asarray(index) if index is not None
                      else np.arange(n))

    def __len__(self):
        return len(next(iter(self.columns.values()))) if self.columns else 0

    def __contains__(self, field):
        return field in self.columns

    def __getitem__(self, field):
        return self.columns[field]

    def get(self, field, default=None):
        return self.columns.get(field, default)


class _FakeEpoch:
    def __init__(self, matched, redshift=0.0):
        self._matched = matched
        self.redshift = redshift

    def galaxies_in_halo(self, index, match_by, r_scale):
        return self._matched


class _FakePipeline:
    """Stands in for shark.photometry.PhotometryPipeline: records every
    (model, z_obs) it was constructed with and every gal_indices array
    abs_mags() was called with, and returns a fixed or row-dependent
    magnitude table."""

    instances = []

    def __init__(self, model, z_obs, bands=None, progress=True,
                mags_by_row=None, **kwargs):
        self.model = model
        self.z_obs = z_obs
        self.bands = bands
        self.kwargs = kwargs
        self.calls = []
        self.mags_by_row = mags_by_row or {}
        _FakePipeline.instances.append(self)

    def abs_mags(self, gal_indices):
        self.calls.append(np.asarray(gal_indices))
        return np.array([self.mags_by_row.get(int(i),
                                              [np.nan] * len(self.bands))
                        for i in gal_indices])


@pytest.fixture(autouse=True)
def _clear_pipeline_instances():
    _FakePipeline.instances = []
    yield
    _FakePipeline.instances = []


def _factory(mags_by_row):
    def make(model, z_obs, bands=None, progress=True, **kwargs):
        return _FakePipeline(model, z_obs, bands=bands, progress=progress,
                            mags_by_row=mags_by_row, **kwargs)
    return make


# ---------------------------------------------------------------------------
# Off by default
# ---------------------------------------------------------------------------

def test_photometry_off_by_default_no_fields_added():
    model = _FakeModel()
    matched = _FakeGalaxyView({"mass": [1e9]}, model=model, index=[0])
    epoch = _FakeEpoch(matched)

    backend = SharkGalaxyBackend(
        photometry_pipeline_factory=_factory({0: [-20.0] * 6}))
    props = backend.galaxy_properties(epoch, halo_row=0)

    assert "LuminosityV" not in props
    assert "Luminosity_ugriz" not in props
    assert "AbsoluteMagnitude_ugriz" not in props
    assert _FakePipeline.instances == []  # never even constructed


# ---------------------------------------------------------------------------
# Basic wiring
# ---------------------------------------------------------------------------

def test_luminosity_v_from_absolute_magnitude():
    model = _FakeModel()
    matched = _FakeGalaxyView({"mass": [1e9]}, model=model, index=[0])
    epoch = _FakeEpoch(matched, redshift=0.1)
    m_v = -19.0

    backend = SharkGalaxyBackend(
        compute_photometry=True,
        photometry_pipeline_factory=_factory(
            {0: [m_v, -18, -19, -20, -20.5, -20.7]}))
    props = backend.galaxy_properties(epoch, halo_row=0)

    expected = 10.0 ** (-0.4 * (m_v - SOLAR_ABSOLUTE_MAGNITUDE["v"]))
    assert props["LuminosityV"] == pytest.approx(expected)


def test_ugriz_magnitudes_and_luminosities():
    model = _FakeModel()
    matched = _FakeGalaxyView({"mass": [1e9]}, model=model, index=[0])
    epoch = _FakeEpoch(matched, redshift=0.0)
    ugriz = [-17.0, -18.0, -18.5, -18.7, -18.8]

    backend = SharkGalaxyBackend(
        compute_photometry=True,
        photometry_pipeline_factory=_factory({0: [-19.0] + ugriz}))
    props = backend.galaxy_properties(epoch, halo_row=0)

    np.testing.assert_allclose(props["AbsoluteMagnitude_ugriz"], ugriz)
    sun = [SOLAR_ABSOLUTE_MAGNITUDE[b] for b in
          ("sdss_u", "sdss_g", "sdss_r", "sdss_i", "sdss_z")]
    expected_lum = 10.0 ** (-0.4 * (np.array(ugriz) - np.array(sun)))
    np.testing.assert_allclose(props["Luminosity_ugriz"], expected_lum)


def test_pipeline_called_with_matched_galaxys_model_row_not_halo_row():
    # matched.index maps the local (central) row to a different row in
    # the model's own tables -- abs_mags must be called with that.
    model = _FakeModel()
    matched = _FakeGalaxyView({"mass": [1e9]}, model=model, index=[42])
    epoch = _FakeEpoch(matched, redshift=0.0)

    backend = SharkGalaxyBackend(
        compute_photometry=True,
        photometry_pipeline_factory=_factory({42: [-19.0] * 6}))
    backend.galaxy_properties(epoch, halo_row=5)

    pipeline = _FakePipeline.instances[0]
    np.testing.assert_array_equal(pipeline.calls[0], [42])


# ---------------------------------------------------------------------------
# Caching across satellites
# ---------------------------------------------------------------------------

def test_pipeline_cached_and_reused_across_calls_same_model_and_redshift():
    model = _FakeModel()
    factory = _factory({0: [-19.0] * 6, 1: [-18.0] * 6})
    backend = SharkGalaxyBackend(compute_photometry=True,
                                 photometry_pipeline_factory=factory)

    matched0 = _FakeGalaxyView({"mass": [1e9]}, model=model, index=[0])
    matched1 = _FakeGalaxyView({"mass": [1e9]}, model=model, index=[1])
    backend.galaxy_properties(_FakeEpoch(matched0, redshift=0.0), halo_row=0)
    backend.galaxy_properties(_FakeEpoch(matched1, redshift=0.0), halo_row=1)

    assert len(_FakePipeline.instances) == 1  # built once, reused


def test_pipeline_rebuilt_for_different_redshift():
    model = _FakeModel()
    factory = _factory({0: [-19.0] * 6})
    backend = SharkGalaxyBackend(compute_photometry=True,
                                 photometry_pipeline_factory=factory)

    matched = _FakeGalaxyView({"mass": [1e9]}, model=model, index=[0])
    backend.galaxy_properties(_FakeEpoch(matched, redshift=0.0), halo_row=0)
    backend.galaxy_properties(_FakeEpoch(matched, redshift=0.5), halo_row=0)

    assert len(_FakePipeline.instances) == 2


def test_pipeline_rebuilt_for_different_model():
    factory = _factory({0: [-19.0] * 6})
    backend = SharkGalaxyBackend(compute_photometry=True,
                                 photometry_pipeline_factory=factory)

    matched_a = _FakeGalaxyView({"mass": [1e9]}, model=_FakeModel(), index=[0])
    matched_b = _FakeGalaxyView({"mass": [1e9]}, model=_FakeModel(), index=[0])
    backend.galaxy_properties(_FakeEpoch(matched_a, redshift=0.0), halo_row=0)
    backend.galaxy_properties(_FakeEpoch(matched_b, redshift=0.0), halo_row=0)

    assert len(_FakePipeline.instances) == 2


# ---------------------------------------------------------------------------
# Restrictions / graceful omission
# ---------------------------------------------------------------------------

def test_no_photometry_for_file_backed_catalogue():
    matched = _FakeGalaxyView({"mass": [1e9]}, model=None, index=[0])
    epoch = _FakeEpoch(matched, redshift=0.0)

    backend = SharkGalaxyBackend(
        compute_photometry=True,
        photometry_pipeline_factory=_factory({0: [-19.0] * 6}))
    props = backend.galaxy_properties(epoch, halo_row=0)

    assert "LuminosityV" not in props
    assert _FakePipeline.instances == []


def test_no_photometry_without_epoch_redshift():
    model = _FakeModel()
    matched = _FakeGalaxyView({"mass": [1e9]}, model=model, index=[0])
    epoch = _FakeEpoch(matched, redshift=None)

    backend = SharkGalaxyBackend(
        compute_photometry=True,
        photometry_pipeline_factory=_factory({0: [-19.0] * 6}))
    props = backend.galaxy_properties(epoch, halo_row=0)

    assert "LuminosityV" not in props


def test_nan_v_band_leaves_luminosity_v_unset():
    model = _FakeModel()
    matched = _FakeGalaxyView({"mass": [1e9]}, model=model, index=[0])
    epoch = _FakeEpoch(matched, redshift=0.0)

    backend = SharkGalaxyBackend(
        compute_photometry=True,
        photometry_pipeline_factory=_factory(
            {0: [np.nan, -18, -19, -20, -20.5, -20.7]}))
    props = backend.galaxy_properties(epoch, halo_row=0)

    assert "LuminosityV" not in props
    # ugriz bands are all finite -> still populated independently
    assert "AbsoluteMagnitude_ugriz" in props


def test_all_nan_ugriz_leaves_ugriz_fields_unset():
    model = _FakeModel()
    matched = _FakeGalaxyView({"mass": [1e9]}, model=model, index=[0])
    epoch = _FakeEpoch(matched, redshift=0.0)

    backend = SharkGalaxyBackend(
        compute_photometry=True,
        photometry_pipeline_factory=_factory(
            {0: [-19.0] + [np.nan] * 5}))
    props = backend.galaxy_properties(epoch, halo_row=0)

    assert "LuminosityV" in props  # v-band alone was fine
    assert "AbsoluteMagnitude_ugriz" not in props
    assert "Luminosity_ugriz" not in props


def test_missing_ugriz_band_from_configured_bands_omits_ugriz_fields():
    model = _FakeModel()
    matched = _FakeGalaxyView({"mass": [1e9]}, model=model, index=[0])
    epoch = _FakeEpoch(matched, redshift=0.0)

    backend = SharkGalaxyBackend(
        compute_photometry=True,
        photometry_bands=["v", "sdss_u", "sdss_g"],  # r/i/z missing
        photometry_pipeline_factory=_factory({0: [-19.0, -17.0, -18.0]}))
    props = backend.galaxy_properties(epoch, halo_row=0)

    assert "LuminosityV" in props
    assert "AbsoluteMagnitude_ugriz" not in props


def test_photometry_failure_is_caught_and_logged_not_raised(caplog):
    model = _FakeModel()
    matched = _FakeGalaxyView({"mass": [1e9]}, model=model, index=[0])
    epoch = _FakeEpoch(matched, redshift=0.0)

    def _boom(model, z_obs, bands=None, progress=True, **kwargs):
        raise RuntimeError("fsps blew up")

    backend = SharkGalaxyBackend(compute_photometry=True,
                                 photometry_pipeline_factory=_boom)
    with caplog.at_level("WARNING"):
        props = backend.galaxy_properties(epoch, halo_row=0)

    assert "LuminosityV" not in props
    assert props["StellarMass"] == pytest.approx(1e9)  # other fields unaffected
    assert any("photometry failed" in r.message for r in caplog.records)


def test_default_photometry_bands_include_v_and_full_ugriz():
    assert "v" in DEFAULT_PHOTOMETRY_BANDS
    for band in ("sdss_u", "sdss_g", "sdss_r", "sdss_i", "sdss_z"):
        assert band in DEFAULT_PHOTOMETRY_BANDS


def test_photometry_options_forwarded_to_factory():
    model = _FakeModel()
    matched = _FakeGalaxyView({"mass": [1e9]}, model=model, index=[0])
    epoch = _FakeEpoch(matched, redshift=0.0)

    backend = SharkGalaxyBackend(
        compute_photometry=True,
        photometry_options={"imf_type": 2, "add_dust": True},
        photometry_pipeline_factory=_factory({0: [-19.0] * 6}))
    backend.galaxy_properties(epoch, halo_row=0)

    pipeline = _FakePipeline.instances[0]
    assert pipeline.kwargs == {"imf_type": 2, "add_dust": True}
