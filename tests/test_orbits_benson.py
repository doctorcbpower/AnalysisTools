"""Tests for analysistools.orbits.benson against closed-form limits and the
values readable off Benson (2005) figures 2 and 10."""
import numpy as np
import pytest

from analysistools.orbits import FlatLCDM
from analysistools.orbits import benson as B
from analysistools.orbits.cosmology import G_KPC_KMS2_MSUN as G


def test_bryan_norman_limits():
    eds = FlatLCDM(0.7, 1.0, 0.0) if False else None
    c = FlatLCDM(0.7, 0.3)
    assert B.bryan_norman_delta(c, 0.0) == pytest.approx(101.2, abs=0.6)   # LCDM z=0
    assert B.bryan_norman_delta(c, 20.0) == pytest.approx(18 * np.pi ** 2, rel=0.01)  # -> EdS
    assert B.omega_m_of_z(c, 0.0) == pytest.approx(0.3)


def test_vvir_identity_matches_gm_over_r():
    c = FlatLCDM(0.678, 0.308)
    a, r = 0.8, 250.0
    delta = float(B.bryan_norman_delta(c, 1 / a - 1))
    v = float(B.v_vir_from_r_vir(c, a, r, delta))
    rho_c = 3 * (100 * c.h * float(c.E(a)) / 1e3) ** 2 / (8 * np.pi * G)
    m = delta * 4 * np.pi / 3 * rho_c * r ** 3
    assert v == pytest.approx(np.sqrt(G * m / r), rel=1e-9)


def test_nfw_virial_radius_exceeds_r200c_for_lcdm_z0():
    from analysistools.orbits import NFWPotential
    c = FlatLCDM(0.678, 0.308)
    p = NFWPotential(2.0e12, 250.0, 25.0)
    rv = float(B.nfw_virial_radius(p.m200, p.r200, p.rs, float(B.bryan_norman_delta(c, 0.0)), c, 1.0)[0])
    assert 1.15 < rv / 250.0 < 1.45
    # enclosed mean density at the solution equals delta x rho_crit
    rho_c = 3 * (100 * c.h / 1e3) ** 2 / (8 * np.pi * G)
    assert float(p.menc(rv)) / (4 * np.pi / 3 * rv ** 3) == pytest.approx(
        float(B.bryan_norman_delta(c, 0.0)) * rho_c, rel=1e-6)


def test_two_body_closed_forms():
    # circular orbit at r_vir: Vr=0, Vt=1, f2=1 -> e=0, peri=apo=1, E=-1/2
    assert float(B.eccentricity(0.0, 1.0)) == pytest.approx(0.0, abs=1e-12)
    assert float(B.pericentre(0.0, 1.0)) == pytest.approx(1.0)
    assert float(B.apocentre(0.0, 1.0)) == pytest.approx(1.0)
    assert float(B.specific_energy(0.0, 1.0)) == pytest.approx(-0.5)
    # parabolic: Vr^2 + Vt^2 = 2 -> E=0, e=1 for any split
    for vt in (0.3, 0.9, 1.3):
        vr = np.sqrt(2 - vt ** 2)
        assert float(B.specific_energy(vr, vt)) == pytest.approx(0.0, abs=1e-12)
        assert float(B.eccentricity(vr, vt)) == pytest.approx(1.0)
    # radial infall: pericentre -> 0
    assert float(B.pericentre(1.0, 1e-4)) < 1e-6
    # e, peri, apo consistent: peri = a(1-e), apo = a(1+e) with a from E
    vr, vt = 0.8, 0.6
    e, rp, ra = (float(f(vr, vt)) for f in (B.eccentricity, B.pericentre, B.apocentre))
    assert (ra - rp) / (ra + rp) == pytest.approx(e, rel=1e-9)


def test_fit_reproduces_the_shape_read_from_figure_2():
    """Benson fig. 2 (z=0): dF/dVr peaks near Vr~1, dF/dVtheta near ~0.7,
    infall angles 0-60 deg with mean in the 30s; tabulated redshifts only."""
    v, dvr, dvt, phi, dphi = B.benson05_marginals(0.0)
    assert v[np.argmax(dvr)] == pytest.approx(1.0, abs=0.2)
    assert v[np.argmax(dvt)] == pytest.approx(0.7, abs=0.2)
    assert np.trapezoid(dvr, v) == pytest.approx(1.0, abs=1e-6)
    assert np.trapezoid(dvt, v) == pytest.approx(1.0, abs=1e-6)
    assert np.trapezoid(dphi, phi) == pytest.approx(1.0, abs=0.02)
    mean_phi = np.trapezoid(phi * dphi, phi)
    assert 25 < mean_phi < 50
    # the paper's own f is a (near-)normalised density
    _, _, _, norm = B.benson05_grid(0.0)
    assert 0.8 < norm < 1.2
    with pytest.raises(ValueError):
        B.benson05_fit(1.0, 0.5, z=0.3)


def test_fit_orbits_are_mostly_near_parabolic_as_in_figure_10():
    """Fig. 10/text: roughly half of all orbits have e = 1 +/- 0.1 and ~18%
    are unbound (E>0) for f2=1."""
    v, _, pdf, _ = B.benson05_grid(0.0, vmax=3.0, n=601)
    VR, VT = np.meshgrid(v, v, indexing="ij")
    w = pdf * (v[1] - v[0]) ** 2
    e = np.asarray(B.eccentricity(VR, np.maximum(VT, 1e-3)))
    E = np.asarray(B.specific_energy(VR, VT))
    assert 0.35 < w[np.abs(e - 1) < 0.1].sum() < 0.65
    assert 0.08 < w[E > 0].sum() < 0.30
