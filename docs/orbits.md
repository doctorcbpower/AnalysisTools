# `analysistools.orbits` — standalone orbit tracking

Status: phase 1 (tree-based) implemented, unit-tested against analytic and
brute-force answers, and run on the real Dorcha SubFind-HBT tree (see
"Real-data findings"). Particle-level tracking is phase 2 (not started).

## What it is

A tool that turns merger-tree tracks into host-frame orbits, events and
per-satellite orbital scalars, writes them to an HDF5 *sidecar*
(`OrbitCatalogue`), and feeds the master catalogue through
`OrbitalPropertiesStage`. The sidecar is built once and reused; the catalogue
stage only reads it (or runs the tracker, equivalently).

```python
from analysistools.orbits import (FlatLCDM, TrackUnits, OrbitConfig,
                                  OrbitTracker, tree_nfw_potential)

cosmo = FlatLCDM.from_tree(tree_backend)            # or FlatLCDM(h, Om)
units = TrackUnits(comoving=True, little_h=True, h=0.6774,
                   length_kpc=1.0, mass_msun=1.0e10)   # CHECK against your tree
pot   = tree_nfw_potential(host_track, cosmo, units)   # or any Potential
trk   = OrbitTracker(cosmo, units, OrbitConfig(radius_factor=1.0),
                     potential=pot, boxsize=box)
cat   = trk.run(satellite_tracks, host_track, ids=halo_ids)
cat.write_hdf5("dorcha_orbits.h5")
```

## Definitions (each recorded in the sidecar attributes)

| Quantity | Definition |
|---|---|
| host frame | physical kpc, centred on the host's central-subhalo position |
| reference radius | `radius_factor * host.extra[radius_key]` (default 1 × `GroupR200`), or a fixed physical `fixed_radius_kpc` |
| infall | first inward crossing of the reference radius on the interpolated orbit; a crossing counts only after the orbit moves `hysteresis` (5 %) beyond it |
| `born_inside` | track starts inside the radius: infall time is then only an upper limit |
| pericentre / apocentre | confirmed extrema of interpolated `r(t)` after infall; confirmed only once `r` has moved `prominence` (5 %) away. Unconfirmed trailing extremum is never reported |
| eccentricity | `(r_apo − r_peri)/(r_apo + r_peri)` from the two latest consecutive turning points |
| `*_EL` | turning points / radial period from the infall (E, L) in the host potential frozen at infall: a separate estimator from the measured events |
| backsplash | crossed inward at least once and now outside |
| preprocessed | subhalo of a different FOF group before first infall |
| orphan | track ends before the host's last snapshot; position continued by orbit integration |

Interpolation: cubic Hermite on the **comoving** separation (slope
`v_pec_rel/a`, so no Hubble-flow term is needed), `a(t)` exact from the
cosmology. Sampled-snapshot minima over-estimate pericentre at ~170 Myr
cadence; `tests/test_orbits.py::test_hermite_pericentre_beats_snapshot_minimum`
quantifies it.

## Satellite definition

**A satellite is any subhalo whose centre lies within R200c of the host
centre at z=0**, from any FOF group, excluding the host's own central
(`select_satellites`, `select_satellites_subfind`; the radius definition is
adaptable exactly as for infall). `OrbitInsideReference` in the catalogue is
the same test on the tracker's final host-frame radius. On the Dorcha host
(R200c = 273 kpc): 2159 satellites; 2139 are in the host's FOF group and 20
are in neighbouring groups (r/R200c = 0.90-0.99, M < 8e7 Msun); 727 FOF
members lie outside R200c and are *not* satellites under this definition
(for M_sub > 1e9 Msun: 12 satellites, vs 26 FOF members). Some selected objects
have `n_infalls == 0` because they only just crossed and have not yet moved
`hysteresis` (5 %) inside; they are still satellites by this definition.

### Resolution thresholds (user-defined)

Two levels, both optional and both explicit about what they measure:

* **At selection time, z=0 bound values**:
  `select_satellites_subfind(..., min_mass=<Msun>, min_particles=<N>)`
  (AND if both). Asking for a threshold the data cannot support (no
  `SubhaloLen` in the tree) raises instead of selecting everything.
* **On the orbit catalogue, peak or final basis**:
  `OrbitCatalogue.select(min_mass=, min_particles=, basis="peak"|"final",
  inside_only=)` returns a boolean mask. `peak` is the object's best-resolved
  moment (the usual completeness criterion for a stripped satellite); `final`
  is its last resolved snapshot. Stored per satellite as `mass_peak/final`,
  `npart_peak/final` (schema 1.1: `OrbitPeakMass`, `OrbitFinalMass`,
  `OrbitPeakParticles`, `OrbitFinalParticles`).

Dorcha host 0, 2159 satellites inside R200c (bound mass per particle is
2.32e5 Msun = 1.57e5 Msun/h; nothing is below 32 particles):

| threshold | z=0 / final | peak |
|---|---|---|
| >= 50 particles | 1469 | 2126 |
| >= 100 particles | 790 | 1804 |
| >= 300 particles | 270 | 904 |
| M >= 1e8 Msun | 188 | 656 |
| M >= 1e9 Msun | 12 | 58 |

The basis matters: 70 % of satellites have a peak particle count at least
twice their current one, so a "100 particle" cut is 790 or 1804 objects
depending on whether it is applied now or at peak.

## Assumptions to know about

* Spherical host potential about the host centre; no host triaxiality, disc or
  non-inertial-frame motion.
* Orphans integrated in the (time-dependent) spherical potential in physical
  coordinates; cosmological expansion terms neglected; mass frozen at last
  resolved value. **Dynamical friction is off by default** (Chandrasekhar,
  σ = Vc/√2, lnΛ = ln(1+M(<r)/M_sat) available via
  `OrbitConfig(dynamical_friction=True)`): a modelling choice with a material
  effect on orphan radii.
* `tree_nfw_potential`: NFW with (M200, R200) from the host track and
  `r_s = R_max/2.163`; snapshots where the implied Vmax disagrees with the
  tree's by > 30 % are skipped, never patched.
* Velocity convention follows `MergerTreeTools`: with `comoving=True`,
  `v_stored = v_pec/√a`. Override with `TrackUnits(velocity=...)`.

## Catalogue integration

`OrbitalPropertiesStage` (no default source: pass `orbit_catalogue=`,
`tracker=` or `tracker_options=`; see `configs/dorcha.yaml`).

* Schema 1.0 fields it fills: `OrbitalPericentre` (smallest confirmed
  pericentre since infall; NaN until the first one), `OrbitalApocentre`
  (latest confirmed apocentre), `OrbitalEccentricity`, `OrbitalPeriod`
  (measured), `OrbitalEnergy`, `OrbitalAngularMomentum` (at infall).
* Schema **1.1** (additive) adds 23 `Orbit*` columns (`extended_fields=True`),
  generated from `orbits/catalogue_fields.py`.
* `HaloPropertiesStage`'s `Minfall`/`RedshiftInfall`/`NumberOfInfalls`/
  `IsBacksplash` remain **FOF-membership** based (first `IsSubhalo`), although
  the 1.0 schema text describes them as "first crossing of host R200c". The
  R200c-based versions are `OrbitInfall*` / `OrbitIsBacksplash`. Decide
  whether to redefine the old ones (needs a unit-consistent mass conversion
  first: `HaloPropertiesStage` reports tree-native mass units).

## What changed from the old stage

Old: pericentre/apocentre = min/max of snapshot-sampled host distance over the
*entire* tree overlap, in whatever comoving/h unit the reader produced. That
(i) over-estimated pericentres, (ii) reported the pre-infall maximum distance
as "apocentre" (so eccentricity → 1 for anything that fell in from afar),
(iii) mixed scale factors. The old unit test pinned (ii); it was replaced.

## Real-data findings (Dorcha SubFind-HBT tree, 2026-10-06)

Tree: 123 snapshots, 73 below z=5, 6.1 M nodes, `BoxSize=100`, h=0.678.

* **Units verified from the file itself**: `UnitLength` = 1 Mpc, `UnitMass` =
  1e10 Msun, km/s; positions are **comoving Mpc/h**, velocities **physical
  peculiar km/s** (finite-difference dx/dt over 40 snapshots matches v/a to
  0.1 %; the GADGET-internal v/sqrt(a) hypothesis is off by 13 %). Load with
  `native_is_comoving=True, native_includes_h=True`. `OrbitTracker.from_epoch`
  derives the velocity mode from what the reader did, not from `comoving` alone.
* **No `Group_R_Crit200` in the tree** (only `Group_M_Crit200`): use
  `radius_key="R200c_from_M200"` (R200c from M200 and the cosmology; agrees
  with a real `Group_R_Crit200` to 0.01 %).
* **Group 1 is not a halo**: it has M200 = 1.55e12 (the largest in the box) but
  550 subhaloes, none above 1.1e9 Msun/h. Selecting hosts by group M200 picks
  it; select by *central-subhalo* mass. Group 0 (central 1.26e12 Msun/h,
  Vmax 198 km/s, R200c = 273 kpc) is the real MW-mass host. Worth checking
  what group 1 is (fragmented FOF? zoom-boundary contamination?) before
  building any MWA selection from `Group_M_Crit200`.
* Host 0, 2866 z=0 FOF members, tracked in ~10 s: 17 % of those that ever
  crossed R200c are backsplash, median infall z = 0.28, median first
  pericentre 60 kpc. 62 % flagged `preprocessed` (any other-FOF membership
  before infall, no mass-ratio or duration cut): treat as an upper bound.
* **Tree-NFW potential**: (M200, Rmax) implies a Vmax 12-30 % above the tree's
  bound-subhalo Vmax at z<2.7 (worse earlier), so the old 30 % cut discarded
  most of the history; tolerance is now explicit (default 50 %).
* **Orphans (truncation test)**: for satellites the tree follows to z=0,
  truncated at z~0.1-1.5 and continued as orphans, then compared with the
  tree's own z=0 radius:

  | normalisation | median dlog10(r) | 16-84 % | within 30 % | pred<R200c (true 0.88) |
  |---|---|---|---|---|
  | `m200` (cut z~0.3) | -0.16 | [-0.46, +0.14] | 0.32 | 0.95 |
  | `vmax` (cut z~0.3) | -0.05 | [-0.42, +0.30] | 0.31 | 0.78 |

  i.e. orphan radii are right only **statistically, at the factor ~2-3 level**
  (orbital-phase sensitivity), with a potential-normalisation-dependent bias in
  radial membership. Do not use individual orphan radii; use distributions, and
  quote the spread between normalisations as a systematic. A particle-measured
  potential is the way to reduce this.
* Orphan continuation is restricted to objects inside `orphan_max_r_over_ref`
  x R_ref when lost (default 1): without that, 68 % of "lost" fringe subhaloes
  were pushed beyond R200c, out to ~1 Mpc, because a bound spherical potential
  without Hubble flow is the wrong model for them. They are flagged
  (`is_orphan_final`, `orphan_continued=False`) with `r_final` the last
  resolved distance.

## Comparison with Benson (2005)

`analysistools.orbits.benson` implements Benson's (2005, MNRAS 358, 551)
conventions: r_vir/V_vir from the Bryan & Norman top-hat overdensity (not
R200c; r_vir/R200c = 1.31 on Dorcha at z=0, insensitive to concentration:
1.30-1.34 for c = 5-20), **Hubble flow added** to the relative peculiar
velocity (shifts V_r by sqrt(2/Delta_vir) ~ 0.12-0.14 V_vir), V_r positive
inward, phi = arctan(V_theta/V_r), two-body E/e/r_peri/r_apo (eqs 5, 7, 10,
11) and the fitting function eqs (2)-(4) with Table 2. Tests check the fit
against statements in the paper (about half of orbits e = 1 +/- 0.1, ~18 %
unbound).

**Correction (2026-10-07):** the "central at infall" subset first used above
required `not IsSubhalo` before crossing, which wrongly excluded objects already
in the *host's own* FOF group. The right Benson-like test is "not a subhalo of a
*different* group" (92-95 % of events pass, not 51 %). With it, at z~0: unbound
fraction 0.30 (not 0.17), vs Benson 0.20. The numbers in the next block use the
corrected flag only where stated; see the variants table below.

**Sample matters.** Benson measures *every* halo merging at each output. The
earlier "satellites inside R200c at z=0" plots are survivors, in R200c units,
peculiar velocities only -- not comparable. Dorcha host 0, all first r_vir
crossings with >= 50 particles (3002 events, z < 1.45):

* V_r: Dorcha <V_r> 0.86 / 0.84 / 0.83 at z ~ 0 / 0.5 / 1 vs fit 0.90 / 0.89 / 0.85
* V_theta: Dorcha <V_theta> 0.81 / 0.87 / 0.82 vs fit 0.75 / 0.71 / 0.70, i.e.
  **systematically more tangential motion** at every redshift
* unbound (E>0) fraction: 0.30 (all) / 0.17 (central at infall) at z~0 vs 0.20; at
  z ~ 0.5 0.34 vs 0.15 (this window contains the host's 1:2 merger)
* survivors (in r_vir at z=0) vs the rest: survivors are *less* energetic
  (median V_r 0.76 vs 0.94; 17 % vs 37 % unbound): the survivor bias is mostly
  about unbound fly-bys that leave, not about angular momentum.

Untested candidate causes of the V_theta offset: halo-centre/velocity definition
(Benson uses most-bound-particle position and velocity; HBT subhalo centre
here), a single host vs Benson's ensemble (mass-dependent trends are reported
in his section 3.2.1), the r_vir from NFW(M200c, Rmax), and the two-body
back-propagation he applies versus our interpolated actual path.

### Do Hubble flow and the radius definition explain the V_theta offset? No

Four variants (corrected sample, first crossings not pre-processed, >= 50
particles), <V_r> / <V_theta> / unbound fraction, vs Benson's fit:

| z window | r_vir + Hubble | r_vir, no Hubble | R200c + Hubble | R200c, no Hubble | Benson fit |
|---|---|---|---|---|---|
| ~0 | .87/.81/.30 | 1.00/.81/.42 | .90/.85/.37 | 1.00/.85/.46 | .90/.75/.20 |
| ~0.5 | .85/.86/.33 | .97/.86/.42 | .87/.90/.38 | .97/.90/.44 | .89/.71/.15 |
| ~1 | .83/.81/.16 | .94/.81/.28 | .84/.85/.17 | .94/.85/.31 | .85/.70/.11 |

Dropping the Hubble flow raises V_r by 0.10-0.13 (away from the fit) and the
unbound fraction by ~0.12; switching r_vir -> R200c raises V_theta by
~0.03-0.04. Neither touches the V_theta excess (+8 % at z~0, ~+20 % at z~0.5-1).
The Benson-definition variant (first column) is the closest to the fit.

## Open checks / next steps

1. Host assignment for paired MWAs / group-accreted satellites: one reference
   host per run; "host at first infall" is not implemented.
2. Orphan box positions -> Shark type-2 matching is not wired (sidecar holds
   host-frame series only).
3. Particle-measured potential: `TabulatedSpherical.from_particles` exists; a
   per-snapshot `TimeDependentPotential` builder from a `Simulation`'s
   snapshots does not (the local `snap_020` is a different 50 Mpc run).
4. Phase 2: particle tracking (tagged progenitor-field particles, per-particle
   pericentres, `ProgenitorParticleFraction`).
5. Decide the `preprocessed` definition (mass ratio / minimum time in group).
