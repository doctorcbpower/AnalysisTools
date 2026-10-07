"""
Tests for analysistools.orbits.

Every test checks against an *independent* answer (an analytic result, or a
brute-force fine-timestep integration), not against the module's own output:
the point is to distinguish a correct implementation from one that merely
runs. The cosmology is made almost static (h = 1e-3 -> a changes by ~5e-4 over
5 Gyr) so that a physical orbit integrated in a static potential is also a
valid cosmological track to within that tolerance.
"""
import numpy as np
import pytest

from analysistools.merger_tree_types import HaloTrack
from analysistools.orbits import (
    FlatLCDM, NFWPotential, OrbitCatalogue, OrbitConfig, OrbitTracker,
    TabulatedSpherical, TimeDependentPotential, TrackUnits,
    integrate_orphans, G_KPC_KMS2_MSUN as G, KMS_TO_KPC_PER_GYR as KMS,
)

COSMO = FlatLCDM(h=1e-3, omega_m=0.3)
T0 = float(COSMO.age(1.0))
PHYS = TrackUnits(comoving=False, little_h=False, length_kpc=1.0,
                  mass_msun=1.0)
NFW = NFWPotential(m200=1.0e12, r200=210.0, rs=25.0)
R_REF = 210.0


def _host(snap_t, r200=R_REF):
    n = len(snap_t)
    a = COSMO.a_of_t(T0 + np.asarray(snap_t))
    return HaloTrack(
        halo_id=0, query_snapnum=n - 1, treefileformat="test",
        SnapNum=np.arange(n, dtype=np.int32), Redshift=1 / a - 1, Time=a,
        Mass=np.full(n, 1e12), Pos=np.zeros((n, 3)), Vel=np.zeros((n, 3)),
        IsSubhalo=np.zeros(n, dtype=bool), HostID=np.zeros(n, dtype=np.int64),
        extra={"GroupR200": np.full(n, r200)})


def _sat(host, pos, vel, n=None, is_subhalo=None, halo_id=1):
    n = len(pos) if n is None else n
    return HaloTrack(
        halo_id=halo_id, query_snapnum=n - 1, treefileformat="test",
        SnapNum=host.SnapNum[:n].copy(), Redshift=host.Redshift[:n].copy(),
        Time=host.Time[:n].copy(), Mass=np.full(n, 1e9),
        Pos=np.asarray(pos[:n], dtype=float),
        Vel=np.asarray(vel[:n], dtype=float),
        IsSubhalo=(np.zeros(n, dtype=bool) if is_subhalo is None
                   else np.asarray(is_subhalo)[:n]),
        HostID=np.full(n, 0, dtype=np.int64), extra={})


def _truth_orbit(x0, v0, t_end, dt_out):
    """Brute-force reference: leapfrog in the static NFW, tiny step."""
    t_out = np.arange(0.0, t_end + 1e-9, dt_out)
    xs, vs = integrate_orphans(NFW, [x0], [v0], t_out, dt_max=2e-4)
    return t_out, xs[0], vs[0]


def _true_pericentre(x0, v0, t_end):
    t = np.arange(0.0, t_end, 5e-4)
    xs, _ = integrate_orphans(NFW, [x0], [v0], t, dt_max=2e-4)
    r = np.linalg.norm(xs[0], axis=1)
    return r.min()


# ---------------------------------------------------------------------
# potentials
# ---------------------------------------------------------------------

def test_nfw_normalisation_and_circular_speed_peak():
    assert NFW.menc(210.0) == pytest.approx(1e12, rel=1e-10)
    r = np.geomspace(1, 300, 4000)
    rpk = r[np.argmax(NFW.vcirc(r))]
    assert rpk == pytest.approx(2.163 * 25.0, rel=5e-3)
    # Phi -> 0 at infinity, Phi'(r) = G M / r^2
    eps = 1e-4
    dphi = (NFW.phi(50 * (1 + eps)) - NFW.phi(50 * (1 - eps))) / (100 * eps)
    assert dphi == pytest.approx(NFW.accel(50.0), rel=1e-5)


def test_tabulated_from_particles_matches_uniform_sphere():
    rng = np.random.default_rng(1)
    n, R, M = 400_000, 100.0, 1.0e11
    u = rng.random(n)
    pts = (rng.normal(size=(n, 3)))
    pts = pts / np.linalg.norm(pts, axis=1)[:, None] * (R * u ** (1 / 3))[:, None]
    pot = TabulatedSpherical.from_particles(pts, M / n, [0, 0, 0],
                                            rmin=5.0, rmax=150.0, n_grid=120)
    for r in (20.0, 60.0, 100.0):
        assert pot.menc(r) == pytest.approx(M * (r / R) ** 3, rel=0.03)
    # analytic: inside Phi = -GM(3R^2 - r^2)/(2R^3), outside -GM/r
    for r in (30.0, 80.0):
        want = -G * M * (3 * R ** 2 - r ** 2) / (2 * R ** 3)
        assert pot.phi(r) == pytest.approx(want, rel=0.02)
    assert pot.phi(140.0) == pytest.approx(-G * M / 140.0, rel=0.02)


def test_tabulated_from_nfw_reproduces_phi():
    # NFW mass diverges logarithmically, so the tabulation's "point mass
    # beyond the grid" closure needs a far outer edge to be accurate
    r = np.geomspace(0.5, 1.0e6, 900)
    tab = TabulatedSpherical.from_potential(NFW, r)
    # Phi from tabulated M must reproduce analytic Phi, even though it is
    # built only by integrating the enclosed mass inwards
    for rr in (5.0, 25.0, 100.0, 400.0):
        assert tab.phi(rr) == pytest.approx(NFW.phi(rr), rel=2e-3)


def test_turning_points_and_period_roundtrip():
    rp, ra = 20.0, 100.0
    Pp, Pa = float(NFW.phi(rp)), float(NFW.phi(ra))
    E = (ra ** 2 * Pa - rp ** 2 * Pp) / (ra ** 2 - rp ** 2)
    L = np.sqrt(2 * (Pa - Pp) / (1 / rp ** 2 - 1 / ra ** 2))
    p, a = NFW.turning_points(E, L, r0=50.0)
    assert (p, a) == (pytest.approx(rp, rel=1e-6), pytest.approx(ra, rel=1e-6))
    # period vs. a brute-force integration started at apocentre
    T = NFW.radial_period(p, a, E, L)
    vt = L / ra
    t = np.arange(0, 2.5 * T, 1e-3)
    xs, _ = integrate_orphans(NFW, [[ra, 0, 0]], [[0, vt, 0]], t, dt_max=5e-4)
    r = np.linalg.norm(xs[0], axis=1)
    peri_idx = [i for i in range(1, r.size - 1) if r[i] < r[i - 1] and r[i] < r[i + 1]]
    assert len(peri_idx) >= 2
    T_meas = t[peri_idx[1]] - t[peri_idx[0]]
    assert T == pytest.approx(T_meas, rel=5e-3)


def test_integrator_conserves_energy_and_friction_removes_it():
    x0, v0 = np.array([[100.0, 0, 0]]), np.array([[0, 80.0, 0]])
    t = np.linspace(0, 6.0, 61)
    xs, vs = integrate_orphans(NFW, x0, v0, t, dt_max=2e-3)
    E = 0.5 * np.sum(vs[0] ** 2, axis=1) + NFW.phi(np.linalg.norm(xs[0], axis=1))
    assert np.max(np.abs(E / E[0] - 1)) < 1e-4
    L = np.linalg.norm(np.cross(xs[0], vs[0]), axis=1)
    assert np.max(np.abs(L / L[0] - 1)) < 1e-6        # central force

    xf, vf = integrate_orphans(NFW, x0, v0, t, m_sat=[5e10],
                               dynamical_friction=True, dt_max=2e-3)
    Ef = 0.5 * np.sum(vf[0] ** 2, axis=1) + NFW.phi(np.linalg.norm(xf[0], axis=1))
    assert Ef[-1] < E[0] - 1.0                          # genuinely loses energy
    assert np.linalg.norm(xf[0, -1]) < np.linalg.norm(xs[0, -1])


def test_time_dependent_potential_interpolates_mass():
    p1 = NFWPotential(1e12, 210.0, 25.0)
    p2 = NFWPotential(2e12, 260.0, 25.0)
    td = TimeDependentPotential([1.0, 3.0], [p1, p2])
    mid = float(td.menc(100.0, 2.0))
    assert mid == pytest.approx(0.5 * (float(p1.menc(100.0)) + float(p2.menc(100.0))),
                                rel=5e-3)
    assert float(td.menc(100.0, 0.0)) == pytest.approx(float(p1.menc(100.0)), rel=5e-3)
    assert float(td.menc(100.0, 9.0)) == pytest.approx(float(p2.menc(100.0)), rel=5e-3)


# ---------------------------------------------------------------------
# kinematics / events
# ---------------------------------------------------------------------

def _eccentric_case(dt_snap=0.17, t_end=7.0):
    x0, v0 = np.array([150.0, 0, 0]), np.array([0, 60.0, 0])
    t, xs, vs = _truth_orbit(x0, v0, t_end, dt_snap)
    return t, xs, vs, _true_pericentre(x0, v0, t_end)


def test_hermite_pericentre_beats_snapshot_minimum():
    """Averaged over snapshot phases (a single phase can land on pericentre
    by luck)."""
    x0, v0 = np.array([150.0, 0, 0]), np.array([0, 60.0, 0])
    t_fine, xs_f, vs_f = _truth_orbit(x0, v0, 7.0, 0.01)
    r_true = _true_pericentre(x0, v0, 7.0)
    err_trk, err_naive = [], []
    for off in (0, 3, 6, 9, 12, 15):
        sel = np.arange(off, t_fine.size, 17)
        t, xs, vs = t_fine[sel] - t_fine[sel][0], xs_f[sel], vs_f[sel]
        host = _host(t)
        rec = OrbitTracker(COSMO, PHYS, OrbitConfig(orphans=False)) \
            .run([_sat(host, xs, vs)], host).record(0)
        err_trk.append(abs(rec.scalars["r_peri_min"] - r_true) / r_true)
        err_naive.append(abs(np.linalg.norm(xs, axis=1).min() - r_true) / r_true)
        assert rec.scalars["n_pericentres"] >= 2
    assert max(err_trk) < 0.01
    assert np.mean(err_naive) > 3 * np.mean(err_trk)


def test_apocentre_is_not_history_maximum_regression():
    """Old stage returned max(Distance) over the whole overlap; a satellite
    falling in from 5 Mpc away must not report 5 Mpc as its apocentre."""
    n = 40
    t = np.arange(n) * 0.17
    host = _host(t)
    x0, v0 = np.array([150.0, 0, 0]), np.array([0, 60.0, 0])
    _, xs, vs = _truth_orbit(x0, v0, t[-1], 0.17)
    xs = xs.copy()
    xs[:3] = [[5000.0, 0, 0], [3000.0, 0, 0], [1200.0, 0, 0]]   # far, early
    rec = OrbitTracker(COSMO, PHYS, OrbitConfig(orphans=False)) \
        .run([_sat(host, xs, vs)], host).record(0)
    assert rec.scalars["r_apo_last"] < 400.0
    assert rec.scalars["ecc_last"] < 1.0


def test_straight_flyby_crossings_and_backsplash():
    n, v, b = 40, 200.0, 50.0                      # km/s, kpc impact param
    t = np.arange(n) * 0.17
    host = _host(t)
    x = np.zeros((n, 3))
    x[:, 0] = -700.0 + v * KMS * t
    x[:, 1] = b
    vel = np.tile([v, 0.0, 0.0], (n, 1))
    rec = OrbitTracker(COSMO, PHYS, OrbitConfig(orphans=False)) \
        .run([_sat(host, x, vel)], host).record(0)
    s = rec.scalars
    x_cross = np.sqrt(R_REF ** 2 - b ** 2)
    t_in_true = (-x_cross + 700.0) / (v * KMS)
    t_out_true = (x_cross + 700.0) / (v * KMS)
    assert s["t_infall"] - T0 == pytest.approx(t_in_true, abs=2e-3)
    assert s["n_infalls"] == 1 and s["n_outfalls"] == 1
    assert s["is_backsplash"] and not s["inside_now"]
    assert s["r_peri_min"] == pytest.approx(b, rel=0.01)
    assert not s["born_inside"]
    # infall state: |v| = 200, L = r x v with |L| = b * v
    assert np.linalg.norm(s["L_infall"]) == pytest.approx(b * v, rel=1e-3)
    assert s["vr_infall"] ** 2 + s["vt_infall"] ** 2 == pytest.approx(v ** 2, rel=1e-6)


def test_reference_radius_is_adaptable():
    n, v, b = 40, 200.0, 50.0
    t = np.arange(n) * 0.17
    host = _host(t)
    x = np.zeros((n, 3)); x[:, 0] = -700.0 + v * KMS * t; x[:, 1] = b
    vel = np.tile([v, 0.0, 0.0], (n, 1))
    base = OrbitTracker(COSMO, PHYS, OrbitConfig(orphans=False)) \
        .run([_sat(host, x, vel)], host).record(0).scalars
    wide = OrbitTracker(COSMO, PHYS, OrbitConfig(orphans=False, radius_factor=2.0)) \
        .run([_sat(host, x, vel)], host).record(0).scalars
    fixed = OrbitTracker(COSMO, PHYS, OrbitConfig(orphans=False, fixed_radius_kpc=100.0)) \
        .run([_sat(host, x, vel)], host).record(0).scalars
    assert wide["t_infall"] < base["t_infall"] < fixed["t_infall"]
    assert wide["r_ref_infall"] == pytest.approx(2 * R_REF, rel=1e-6)
    assert fixed["r_ref_infall"] == pytest.approx(100.0, rel=1e-6)


def test_hysteresis_ignores_boundary_skimming():
    n = 30
    t = np.arange(n) * 0.17
    host = _host(t)
    r = 1.0 + 0.01 * np.cos(np.arange(n) * 1.3)      # 0.99..1.01 x R_ref
    x = np.zeros((n, 3)); x[:, 0] = R_REF * r
    rec = OrbitTracker(COSMO, PHYS, OrbitConfig(orphans=False)) \
        .run([_sat(host, x, np.zeros((n, 3)))], host).record(0)
    assert rec.scalars["n_infalls"] == 0 and rec.scalars["n_outfalls"] == 0


def test_circular_orbit_has_no_pericentre_events():
    r0 = 80.0
    vc = float(NFW.vcirc(r0))
    t, xs, vs = _truth_orbit(np.array([r0, 0, 0]), np.array([0, vc, 0]), 5.0, 0.17)
    host = _host(t)
    s = OrbitTracker(COSMO, PHYS, OrbitConfig(orphans=False)) \
        .run([_sat(host, xs, vs)], host).record(0).scalars
    assert s["n_pericentres"] == 0 and np.isnan(s["ecc_last"])
    assert s["born_inside"]


def test_comoving_little_h_representation_gives_same_orbit():
    t, xs, vs, _ = _eccentric_case(t_end=5.0)
    host_p = _host(t)
    sat_p = _sat(host_p, xs, vs)
    a = host_p.Time
    h = 0.7
    # express the *same* physical orbit as comoving Mpc/h, sqrt(a)-scaled velocity
    pos_c = xs / a[:, None] / 1000.0 * h
    vel_c = vs / np.sqrt(a)[:, None]
    host_c = _host(t)
    host_c.extra["GroupR200"] = np.full(len(t), R_REF) / a / 1000.0 * h
    sat_c = _sat(host_c, pos_c, vel_c)
    unit_c = TrackUnits(comoving=True, little_h=True, h=h, length_kpc=1000.0,
                        mass_msun=1.0)
    s_p = OrbitTracker(COSMO, PHYS, OrbitConfig(orphans=False)).run([sat_p], host_p).record(0).scalars
    s_c = OrbitTracker(COSMO, unit_c, OrbitConfig(orphans=False)).run([sat_c], host_c).record(0).scalars
    for k in ("r_peri_min", "r_apo_last", "t_infall", "n_pericentres"):
        assert s_c[k] == pytest.approx(s_p[k], rel=1e-6) or \
            (np.isnan(s_c[k]) and np.isnan(s_p[k]))


def test_periodic_box_wrap_is_applied():
    n = 20
    t = np.arange(n) * 0.17
    host = _host(t)
    box = 1.0e5                                      # comoving kpc
    com = TrackUnits(comoving=True, little_h=False, length_kpc=1.0, mass_msun=1.0)
    x = np.zeros((n, 3)); x[:, 0] = box - 100.0       # just across the boundary
    rec = OrbitTracker(COSMO, com, OrbitConfig(orphans=False), boxsize=box) \
        .run([_sat(host, x, np.zeros((n, 3)))], host).record(0)
    assert rec.series["r"] == pytest.approx(100.0 * host.Time, rel=1e-6)


# ---------------------------------------------------------------------
# E, L estimators and orphans
# ---------------------------------------------------------------------

def test_infall_EL_turning_points_match_measured_orbit():
    # start outside R_ref so there is a real infall: a bound orbit that
    # crosses in from beyond R200; run long enough that its first apocentre
    # is *confirmed* (an extremum is only reported once r has moved off it)
    x0, v0 = np.array([300.0, 0, 0]), np.array([-20.0, 70.0, 0])
    t, xs, vs = _truth_orbit(x0, v0, 9.0, 0.17)
    host = _host(t)
    pot = TimeDependentPotential([T0 + t[0], T0 + t[-1]], [NFW, NFW])
    tr = OrbitTracker(COSMO, PHYS, OrbitConfig(orphans=False), potential=pot)
    s = tr.run([_sat(host, xs, vs)], host).record(0).scalars
    assert not s["born_inside"]
    assert s["r_peri_EL"] == pytest.approx(s["r_peri_first"], rel=0.03)
    assert s["r_apo_EL"] == pytest.approx(s["r_apo_first"], rel=0.03)


def test_orphan_continuation_recovers_true_orbit():
    x0, v0 = np.array([150.0, 0, 0]), np.array([0, 60.0, 0])
    t, xs, vs = _truth_orbit(x0, v0, 7.0, 0.17)
    host = _host(t)
    n_keep = 14
    sat = _sat(host, xs, vs, n=n_keep)               # tree loses it at snap 13
    pot = TimeDependentPotential([T0 + t[0], T0 + t[-1]], [NFW, NFW])
    rec = OrbitTracker(COSMO, PHYS, OrbitConfig(), potential=pot) \
        .run([sat], host).record(0)
    s = rec.scalars
    assert s["is_orphan_final"] and s["last_resolved_snap"] == n_keep - 1
    assert rec.series["snapnum"][-1] == len(t) - 1
    assert rec.series["extrapolated"].sum() == len(t) - n_keep
    r_true = np.linalg.norm(xs, axis=1)
    assert rec.series["r"] == pytest.approx(r_true, rel=0.03)
    assert s["r_final"] == pytest.approx(r_true[-1], rel=0.03)
    # events continue into the orphan phase: pericentres beyond the cut
    assert s["t_peri_last"] - T0 > t[n_keep - 1]


def test_fringe_orphans_are_not_continued():
    """A subhalo lost outside the reference radius is not orbiting the host;
    it must not be integrated in a bound spherical potential."""
    t = np.arange(20) * 0.17
    host = _host(t)
    x = np.zeros((20, 3)); x[:, 0] = 400.0           # 400 kpc > R_ref = 210
    pot = TimeDependentPotential([T0, T0 + t[-1]], [NFW, NFW])
    sat = _sat(host, x, np.zeros((20, 3)), n=10)
    rec = OrbitTracker(COSMO, PHYS, potential=pot).run([sat], host).record(0)
    assert rec.scalars["is_orphan_final"] and not rec.scalars["orphan_continued"]
    assert not rec.series["extrapolated"].any() and rec.series["snapnum"][-1] == 9
    # ... but inside it, it is, and the flag says so
    x[:, 0] = 100.0
    rec = OrbitTracker(COSMO, PHYS, potential=pot).run(
        [_sat(host, x, np.zeros((20, 3)), n=10)], host).record(0)
    assert rec.scalars["orphan_continued"] and rec.series["extrapolated"].any()
    # opting in to continue everything
    x[:, 0] = 400.0
    rec = OrbitTracker(COSMO, PHYS, OrbitConfig(orphan_max_r_over_ref=np.inf),
                       potential=pot).run(
        [_sat(host, x, np.zeros((20, 3)), n=10)], host).record(0)
    assert rec.scalars["orphan_continued"]


def test_orphans_without_potential_stop_and_warn(caplog):
    t = np.arange(20) * 0.17
    host = _host(t)
    x = np.zeros((20, 3)); x[:, 0] = 100.0
    sat = _sat(host, x, np.zeros((20, 3)), n=10)
    with caplog.at_level("WARNING"):
        rec = OrbitTracker(COSMO, PHYS).run([sat], host).record(0)
    assert "no potential" in caplog.text
    assert rec.series["snapnum"][-1] == 9 and rec.scalars["is_orphan_final"]


def test_unresolved_satellite_gets_fill_record():
    host = _host(np.arange(5) * 0.17)
    rec = OrbitTracker(COSMO, PHYS).run([None], host).record(0)
    assert rec.scalars["n_infalls"] == 0 and np.isnan(rec.scalars["t_infall"])


def test_preprocessing_flag():
    n, v, b = 40, 200.0, 50.0
    t = np.arange(n) * 0.17
    host = _host(t)
    x = np.zeros((n, 3)); x[:, 0] = -700.0 + v * KMS * t; x[:, 1] = b
    vel = np.tile([v, 0.0, 0.0], (n, 1))
    sub = np.zeros(n, dtype=bool); sub[3:6] = True        # satellite of someone else
    sat = _sat(host, x, vel, is_subhalo=sub)
    sat.HostID[3:6] = 77
    s = OrbitTracker(COSMO, PHYS, OrbitConfig(orphans=False)).run([sat], host).record(0).scalars
    assert s["preprocessed"] and s["snap_fof_infall"] == 3


# ---------------------------------------------------------------------
# sidecar
# ---------------------------------------------------------------------

def test_sidecar_roundtrip(tmp_path):
    t, xs, vs, _ = _eccentric_case(t_end=5.0)
    host = _host(t)
    pot = TimeDependentPotential([T0 + t[0], T0 + t[-1]], [NFW, NFW])
    sats = [_sat(host, xs, vs, halo_id=11), None,
            _sat(host, xs, vs, n=10, halo_id=13)]
    cat = OrbitTracker(COSMO, PHYS, potential=pot).run(sats, host, ids=[11, 12, 13])
    path = str(tmp_path / "orbits.h5")
    cat.write_hdf5(path)
    back = OrbitCatalogue.read_hdf5(path)
    assert list(back.ids) == [11, 12, 13]
    for name in ("r_peri_min", "r_final", "n_pericentres", "t_infall"):
        np.testing.assert_array_equal(back.scalar(name), cat.scalar(name))
    np.testing.assert_array_equal(back.record(2).series["r"], cat.record(2).series["r"])
    assert len(back.record(0).events) == len(cat.record(0).events)
    assert back.record(1).series == {}
    assert back.metadata["config_radius_key"] == "GroupR200"


def test_tree_nfw_potential_normalisations():
    from analysistools.orbits import tree_nfw_potential
    n = 10
    t = np.arange(n) * 0.5
    host = _host(t)
    r200 = np.full(n, R_REF)
    host.extra.update({"GroupM200": np.full(n, 1e12), "VmaxRad": np.full(n, 54.0),
                       "Vmax": np.full(n, 200.0), "GroupR200": r200})
    p_m = tree_nfw_potential(host, COSMO, PHYS, normalise="m200")
    p_v = tree_nfw_potential(host, COSMO, PHYS, normalise="vmax")
    tt = T0 + 1.0
    # "vmax": circular speed at Rmax reproduces the tree's Vmax exactly
    assert float(p_v.vcirc(54.0, tt)) == pytest.approx(200.0, rel=2e-2)
    # "m200": M(<R200) reproduces the group mass
    assert float(p_m.menc(R_REF, tt)) == pytest.approx(1e12, rel=2e-2)
    with pytest.raises(ValueError):
        tree_nfw_potential(host, COSMO, PHYS, normalise="nonsense")


# ---------------------------------------------------------------------
# satellite selection: inside R200c, irrespective of FOF group
# ---------------------------------------------------------------------

def test_select_satellites_inside_reference_radius_wraps_and_ignores_fof():
    from analysistools.orbits import select_satellites
    com = TrackUnits(comoving=True, little_h=False, length_kpc=1.0, mass_msun=1.0)
    box, a = 1.0e5, 0.5                         # comoving kpc; a = 0.5
    host = np.array([box - 50.0, 10.0, 10.0])   # near the +x box edge
    r_ref = 100.0                               # physical kpc -> 200 comoving
    # comoving offsets from host along x; positions wrap across the boundary
    off = np.array([0.0, 150.0, 199.0, 201.0, -190.0, 5000.0])
    pos = np.tile(host, (off.size, 1)); pos[:, 0] = (host[0] + off) % box
    ids = np.arange(100, 106)
    grp = np.array([0, 0, 1, 0, 1, 0])           # FOF group is *not* the criterion
    sel = select_satellites(pos, ids, host, r_ref, a, com, boxsize=box,
                            group_nr=grp, host_group=0, exclude_ids=(100,))
    # 100 is the host central (excluded); 201 and 5000 are outside r_ref
    assert list(sel.subhalo_id) == [101, 102, 104]
    np.testing.assert_allclose(sel.r, np.array([150.0, 199.0, 190.0]) * a)
    assert list(sel.in_host_fof) == [True, False, False]   # 102,104 are other-FOF
    assert sel.r_over_ref.max() < 1.0
    # tightening the reference radius (radius definition stays adaptable)
    tight = select_satellites(pos, ids, host, 0.8 * r_ref, a, com, boxsize=box,
                              exclude_ids=(100,))
    assert list(tight.subhalo_id) == [101]


def test_selection_matches_tracker_inside_now():
    from analysistools.orbits import select_satellites
    n, host = 12, None
    t = np.arange(n) * 0.17
    host = _host(t)
    rs = np.array([50.0, 150.0, 209.0, 211.0, 400.0])
    sats = []
    for k, r in enumerate(rs):
        x = np.zeros((n, 3)); x[:, 0] = r
        sats.append(_sat(host, x, np.zeros((n, 3)), halo_id=10 + k))
    cat = OrbitTracker(COSMO, PHYS, OrbitConfig(orphans=False)).run(
        sats, host, ids=[10 + k for k in range(len(rs))])
    pos = np.array([[r, 0, 0] for r in rs])
    a_last = float(host.Time[-1])
    # physical-unit tracks: comoving position = physical / a
    sel = select_satellites(pos / a_last, np.arange(len(rs)), [0, 0, 0],
                            R_REF, a_last,
                            TrackUnits(comoving=True, length_kpc=1.0, mass_msun=1.0))
    inside = np.zeros(len(rs), dtype=bool); inside[sel.index] = True
    np.testing.assert_array_equal(inside, cat.scalar("inside_now"))


# ---------------------------------------------------------------------
# user-defined mass / particle-count thresholds
# ---------------------------------------------------------------------

def test_select_satellites_thresholds_and_missing_data_raises():
    from analysistools.orbits import select_satellites
    u = TrackUnits(comoving=True, length_kpc=1.0, mass_msun=10.0)  # 1 unit = 10 Msun
    pos = np.array([[10.0, 0, 0]] * 4)
    ids = np.arange(4)
    mass = np.array([1.0, 5.0, 50.0, 500.0])          # -> 10, 50, 500, 5000 Msun
    npart = np.array([10, 20, 300, 4000])
    kw = dict(pos=pos, subhalo_id=ids, host_pos=[0, 0, 0], r_ref_phys=1e3,
              a=1.0, units=u, mass=mass, nparticles=npart)
    assert list(select_satellites(**kw).subhalo_id) == [0, 1, 2, 3]
    assert list(select_satellites(**kw, min_mass=100.0).subhalo_id) == [2, 3]
    assert list(select_satellites(**kw, min_particles=300).subhalo_id) == [2, 3]
    # both given: AND
    assert list(select_satellites(**kw, min_mass=100.0, min_particles=1000).subhalo_id) == [3]
    kw2 = {k: v for k, v in kw.items() if k not in ("mass", "nparticles")}
    with pytest.raises(ValueError, match="min_mass"):
        select_satellites(**kw2, min_mass=1.0)
    with pytest.raises(ValueError, match="min_particles"):
        select_satellites(**kw2, min_particles=10)


def test_catalogue_select_peak_vs_final_basis(tmp_path):
    """A satellite that was big and is now stripped passes a peak cut but
    fails a final cut -- the reason the basis is explicit."""
    n = 12
    t = np.arange(n) * 0.17
    host = _host(t)
    x = np.zeros((n, 3)); x[:, 0] = 100.0
    stripped = _sat(host, x, np.zeros((n, 3)), halo_id=1)
    stripped.Mass[:] = np.linspace(1e10, 1e8, n)
    stripped.extra["Len"] = np.linspace(1000, 10, n).astype(int)
    small = _sat(host, x, np.zeros((n, 3)), halo_id=2)
    small.Mass[:] = 5e7
    small.extra["Len"] = np.full(n, 5)
    cat = OrbitTracker(COSMO, PHYS, OrbitConfig(orphans=False)).run(
        [stripped, small, None], host, ids=[1, 2, 3])
    assert list(cat.select(min_mass=1e9, basis="peak")) == [True, False, False]
    assert list(cat.select(min_mass=1e9, basis="final")) == [False, False, False]
    assert list(cat.select(min_particles=100, basis="peak")) == [True, False, False]
    assert list(cat.select(min_particles=5, basis="final")) == [True, True, False]
    assert list(cat.select(min_mass=1e8, min_particles=100, basis="peak")) == [True, False, False]
    assert list(cat.select(inside_only=True)) == [True, True, False]
    with pytest.raises(ValueError):
        cat.select(basis="nonsense")
    # survives the sidecar round trip, and old sidecars (missing scalars) load
    path = str(tmp_path / "o.h5"); cat.write_hdf5(path)
    back = OrbitCatalogue.read_hdf5(path)
    assert list(back.select(min_particles=100)) == [True, False, False]
    import h5py
    with h5py.File(path, "a") as f:
        del f["Satellites/mass_peak"], f["Satellites/npart_peak"]
    old = OrbitCatalogue.read_hdf5(path)
    assert np.all(np.isnan(old.scalar("mass_peak")))
    with pytest.raises(ValueError, match="particle counts"):
        old.select(min_particles=10)


def test_catalogue_select_without_particle_counts_raises():
    n = 6
    host = _host(np.arange(n) * 0.17)
    x = np.zeros((n, 3)); x[:, 0] = 100.0
    cat = OrbitTracker(COSMO, PHYS, OrbitConfig(orphans=False)).run(
        [_sat(host, x, np.zeros((n, 3)))], host)
    assert cat.scalar("npart_peak")[0] == -1
    with pytest.raises(ValueError, match="particle counts"):
        cat.select(min_particles=10)


# ---------------------------------------------------------------------
# configured samples: counts (peak basis) vs orbits (final basis)
# ---------------------------------------------------------------------

def _sample_fixture():
    n = 12
    host = _host(np.arange(n) * 0.17)
    host.extra["GroupM200"] = np.full(n, 1e12)          # host M200c = 1e12 Msun
    x = np.zeros((n, 3)); x[:, 0] = 100.0               # inside R_ref (210)
    xo = x.copy(); xo[:, 0] = 400.0                     # outside
    mk = lambda hid, m0, m1, l0, l1, pos=x: _sat_with(host, pos, hid, m0, m1, l0, l1)
    sats = [
        mk(1, 1e10, 1e10, 1000, 1000),      # big, stays big       -> both
        mk(2, 1e10, 1e7, 1000, 40),         # big, stripped        -> counts only
        mk(3, 1e8, 1e8, 400, 400),          # below 1e-4*M (=1e8)? 1e8>=1e8 passes; use 9e7 below
        mk(4, 5e7, 5e7, 30, 30),            # always small         -> neither
        mk(5, 1e10, 1e10, 1000, 1000, xo),  # big but OUTSIDE R200 -> flags yes, satellites no
    ]
    sats[2].Mass[:] = 9e7                   # just under 1e-4 * 1e12
    return host, sats


def _sat_with(host, pos, hid, m0, m1, l0, l1):
    n = len(host.SnapNum)
    s = _sat(host, pos, np.zeros((n, 3)), halo_id=hid)
    s.Mass[:] = np.linspace(m0, m1, n)
    s.extra["Len"] = np.linspace(l0, l1, n).astype(int)
    return s


def test_samples_use_their_own_basis_and_combine_with_inside():
    host, sats = _sample_fixture()
    cfg = OrbitConfig(orphans=False, min_mass_fraction=1e-4, min_particles=50)
    cat = OrbitTracker(COSMO, PHYS, cfg).run(sats, host, ids=[1, 2, 3, 4, 5])
    assert cat.metadata["min_mass_effective_msun"] == pytest.approx(1e8)
    assert cat.metadata["host_m200_final_msun"] == pytest.approx(1e12)
    # counts (peak): 2 passes because its PEAK was big; 3 and 4 fail the mass cut
    assert list(cat.scalar("in_count_sample")) == [True, True, False, False, True]
    # orbits (final): the stripped satellite 2 no longer qualifies
    assert list(cat.scalar("in_orbit_sample")) == [True, False, False, False, True]
    # satellite samples additionally require being inside R_ref
    assert list(cat.satellites("counts")) == [True, True, False, False, False]
    assert list(cat.satellites("orbits")) == [True, False, False, False, False]
    with pytest.raises(ValueError):
        cat.satellites("nonsense")
    # select() accepts the same fraction, from the recorded host mass
    assert list(cat.select(min_mass_fraction=1e-4, min_particles=50, basis="peak")) \
        == list(cat.scalar("in_count_sample"))


def test_mass_cut_is_larger_of_absolute_and_fraction_and_bases_are_configurable():
    host, sats = _sample_fixture()
    cat = OrbitTracker(COSMO, PHYS, OrbitConfig(
        orphans=False, min_mass=5e9, min_mass_fraction=1e-4)).run(sats, host)
    assert cat.metadata["min_mass_effective_msun"] == pytest.approx(5e9)
    assert list(cat.scalar("in_count_sample")) == [True, True, False, False, True]
    # swap the bases: counts on final, orbits on peak
    cat2 = OrbitTracker(COSMO, PHYS, OrbitConfig(
        orphans=False, min_mass_fraction=1e-4, counts_basis="final",
        orbits_basis="peak")).run(sats, host)
    assert list(cat2.scalar("in_count_sample")) == [True, False, False, False, True]
    assert list(cat2.scalar("in_orbit_sample")) == [True, True, False, False, True]
    with pytest.raises(ValueError):
        OrbitTracker(COSMO, PHYS, OrbitConfig(counts_basis="x")).run(sats, host)
    # no cuts at all -> everyone resolved is in both samples
    cat3 = OrbitTracker(COSMO, PHYS, OrbitConfig(orphans=False)).run(sats, host)
    assert cat3.scalar("in_count_sample").all() and cat3.scalar("in_orbit_sample").all()


def test_fraction_needs_host_m200():
    host, sats = _sample_fixture()
    del host.extra["GroupM200"]
    with pytest.raises(ValueError, match="GroupM200"):
        OrbitTracker(COSMO, PHYS, OrbitConfig(min_mass_fraction=1e-4)).run(sats, host)


def test_sample_flags_survive_sidecar_roundtrip(tmp_path):
    host, sats = _sample_fixture()
    cfg = OrbitConfig(orphans=False, min_mass_fraction=1e-4, min_particles=50)
    cat = OrbitTracker(COSMO, PHYS, cfg).run(sats, host, ids=[1, 2, 3, 4, 5])
    path = str(tmp_path / "s.h5"); cat.write_hdf5(path)
    back = OrbitCatalogue.read_hdf5(path)
    assert list(back.satellites("counts")) == list(cat.satellites("counts"))
    assert list(back.select(min_mass_fraction=1e-4, basis="final")) == \
        list(cat.select(min_mass_fraction=1e-4, basis="final"))
