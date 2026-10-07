# Orbit tracking for the Dorcha analysis: work report

Period: 2026-10-06 to 2026-10-07. Repository: `~/MyCodes/AnalysisTools`, commit `ab71236`
("Add analysistools.orbits: standalone orbit tracker feeding the catalogue"), on `main`, not pushed.
Technical reference for the code: [orbits.md](orbits.md). This document is the narrative: what was
asked, what was decided, what was built, how it was checked, what the Dorcha data show, what went
wrong along the way, and what is still open.

---

## 1. Starting point

### 1.1 What the papers need

Both drafts in `~/Papers/DorchaPapers` (`Paper2_SAM_LSST`, `Paper3_Fossils`) are skeletons: abstract,
figure plans, simulation and halo-finding methods, MWA selection. Neither has results yet. Both cite a
section `\ref{ssec:pericentre}` that does not exist, and the snapshot-cadence text says the 73
snapshots between z=5 and z=0 are "tailored to accurately track the orbits of individual particles".

* **Paper 2** (satellite luminosity function and Rubin forecast): the radial distribution (its Fig. 3)
  asks *why* satellites are distributed as they are (tidal disruption, dynamical friction,
  environment). That needs pericentre, infall time and orbit history, and an isolated-vs-paired
  comparison needs a defined host per satellite.
* **Paper 3** (reionization fossils): "dynamical memory" of the formation site needs infall time
  relative to reionization, pericentre and a pre-processing flag.
* Neither draft yet says how orbits are measured.

### 1.2 What AnalysisTools had

`OrbitalPropertiesStage` computed pericentre and apocentre as the min and max of the
snapshot-sampled host distance over the whole tree overlap, in whatever comoving/little-h unit the
tree reader produced. Defects found by inspection:

1. **Apocentre was the pre-infall maximum distance.** A satellite that was 5 Mpc away at z=5 reported
   5 Mpc as its apocentre; eccentricity then tended to 1. The old unit test pinned this behaviour.
2. **Pericentres biased high.** At ~170 Myr snapshot spacing a satellite at a few hundred km/s moves tens
   of kpc between snapshots.
3. **Mixed scale factors.** Distances were comoving, then min/max'd across snapshots.
4. **No events** (crossings, number of pericentres, backsplash); infall defined by FOF membership
   (`IsSubhalo`), though the schema text says "first crossing of host R200c".
5. **One fixed host**, so paired MWAs and group-accreted satellites were not handled.

The unit system of the real tree was unverified (the code itself noted that the SubFind default for
`native_is_comoving` was wrong for DORCHA_01).

---

## 2. Decisions

Per the standing pre-flight rule, a seven-point audit was presented before any code was written.

| Question | Decision (user) |
|---|---|
| Infall definition | First crossing of **R200c**; the radius definition must stay adaptable |
| Orphans | **Keep orphan positions** |
| Particle tracking | Useful, but **trees first** |
| Potential | Wants an **efficient tool to compute the host potential** |
| Satellite definition | **Inside R200c at z=0** (not FOF membership) |
| Resolution cuts | **User-defined** minimum mass and/or particle count |
| Cut defaults | Mass ≥ 1e-4 × host M200c and ≥ 50 particles; **peak basis for counts, final basis for orbits** |
| Reference for infall parameters | **Benson (2005)** |

Decisions deliberately **left open** (not guessed): host assignment for paired MWAs and pre-processed
satellites; dynamical friction on orphans (off by default); how orphan positions map onto Shark
type-2 galaxies; the definition of "pre-processed".

---

## 3. What was built

New package `analysistools/orbits/`:

| Module | Purpose |
|---|---|
| `cosmology.py` | Flat ΛCDM age, H(a), exact inverse a(t); G in kpc (km/s)²/Msun |
| `kinematics.py` | `TrackUnits`, host-frame series, Hermite interpolation, crossing and turning-point detection, `R200c_from_M200` reference radius |
| `potential.py` | `NFWPotential`, `TabulatedSpherical` (from M(<r), or exactly from particles), `TimeDependentPotential`; E, L → turning points and radial period |
| `orphans.py` | Orphan integrator (KDK leapfrog, optional Chandrasekhar friction) |
| `selection.py` | `select_satellites[_subfind]`: inside R200c from any FOF group, with mass and particle thresholds |
| `tracker.py` | `OrbitTracker`, `OrbitConfig`, `OrbitCatalogue` (HDF5 sidecar), `tree_nfw_potential` |
| `catalogue_fields.py` | Mapping from sidecar scalars to catalogue columns; generates the schema 1.1 fragment |
| `benson.py` | Benson (2005) conventions, two-body formulae, fitting function, Table 2 |

Other changes: `OrbitalPropertiesStage` rewritten to read the sidecar or run a tracker (no default
source: units and definitions must be stated); schema **1.1** (additive, 31 `Orbit*` columns
generated from `catalogue_fields.py`); `treeio_subfind` now reads `SubhaloLen` into `HaloTrack.extra["Len"]`;
`configs/dorcha.yaml` and the end-to-end notebook carry the verified Dorcha settings.

### 3.1 Definitions (all recorded in the sidecar attributes)

* **Host frame:** physical kpc, centred on the host's central subhalo.
* **Interpolation:** cubic Hermite on the *comoving* separation (slope = peculiar velocity / a, so no
  Hubble-flow term is needed), a(t) exact from the cosmology.
* **Infall:** first inward crossing of `radius_factor × host[radius_key]` (default 1 × R200c) on the
  interpolated orbit, with a 5 % hysteresis band so boundary-skimming is not counted as two events.
* **Pericentre / apocentre:** confirmed extrema of r(t) after infall (confirmed once r has moved 5 % away);
  an unconfirmed trailing extremum is never reported. Eccentricity from the two latest consecutive
  turning points. `*_EL` quantities are a separate estimator from (E, L) at infall in the frozen potential.
* **Satellite:** any subhalo within R200c of the host centre at z=0, from any FOF group, excluding the host central.
* **Orphans:** a lost subhalo is continued by integration in the (time-dependent, spherical) host potential
  only if it was inside `orphan_max_r_over_ref` × R_ref (default 1) when lost.
* **Samples:** counts sample = resolution cuts on the *peak* basis; orbit sample = cuts on the *final* basis.
  Mass cut = max(`min_mass`, `min_mass_fraction` × host M200c(z=0)).

---

## 4. Verification

### 4.1 Unit tests

Each test checks an independent answer (analytic or brute-force), not the module's own output.
Full suite at the commit: **506 passed, 4 skipped**. Highlights:

* Hermite pericentre vs the sampled minimum, averaged over six snapshot phases (a single phase can land on pericentre by luck).
* Straight-line fly-by: crossing times, backsplash flag, pericentre = impact parameter, |L| = b·v.
* Regression for the old bug (satellite infalling from 5 Mpc must not report 5 Mpc as apocentre).
* Comoving/little-h representation of the same orbit gives the same result; periodic wrap.
* Integrator: energy and L conservation, friction removes energy; E, L → turning points and period vs brute force.
* Orphan continuation recovers a truncated true orbit; fringe orphans are *not* continued.
* Benson: closed-form limits (circular, parabolic, radial), Bryan & Norman limits, V_vir identity, and that the transcribed
  Table 2 reproduces statements in the paper's text (about half of orbits e = 1 ± 0.1, ~18 % unbound).

Bugs the tests caught in the first pass: an inverted kick unit conversion in the integrator (energy drift), a sign error
in the tabulated potential, and geometric instead of linear time interpolation of the enclosed mass.

### 4.2 Real-data checks (`~/Currentwork/Dorcha/trees.hdf5`)

| Check | Result |
|---|---|
| File contents | 123 snapshots (73 below z=5), 6,097,531 nodes, 87,711 trees; BoxSize 100; h = 0.678, Ωm = 0.308 |
| Units | `UnitLength` = 1 Mpc, `UnitMass` = 1e10 Msun, km/s (from `Parameters`) |
| Positions / velocities | **Comoving Mpc/h and physical peculiar km/s**: finite-difference dx/dt matches v/a to 0.1 %; the v/√a hypothesis is off by 13 % |
| `Group_R_Crit200` | **Not in the tree** (only `Group_M_Crit200`); R200c derived from M200 and the cosmology matches a real catalogue's `Group_R_Crit200` to 0.01 % |
| Host selection | Group 1 has the largest M200 (1.55e12) but 550 subhaloes, none above 0.11e10 Msun/h: **not a halo**. Host = group 0 (SubhaloNr 0; central 1.26e12 Msun/h, Vmax 198 km/s, R200c 273 kpc, M200c 2.17e12 Msun) |
| Particle mass | 2.32e5 Msun (= 1.57e5 Msun/h, the paper's value); nothing below 32 particles |
| Runtime | 2866 FOF members tracked in ~10 s |

**Tree-NFW potential.** From (M200, Rmax) the NFW implies a Vmax 12–30 % above the tree's bound-subhalo Vmax at z<2.7 (median
ratio 1.19 over 67 accepted snapshots), worse earlier. The original 30 % acceptance cut discarded most of the history; the tolerance is
now an explicit parameter (default 50 %), and `normalise="m200"|"vmax"` is selectable.

**Orphan truncation test** (satellites the tree *does* follow to z=0, cut at z≈0.1–1.5, continued as orphans, compared with the tree's own z=0 radius):

| Normalisation (cut at z≈0.3) | median Δlog10 r | 16–84 % | within 30 % | predicted inside R200c (true 0.88) |
|---|---|---|---|---|
| `m200` | −0.16 | [−0.46, +0.14] | 0.32 | 0.95 |
| `vmax` | −0.05 | [−0.42, +0.30] | 0.31 | 0.78 |

Orphan radii are therefore right only **statistically, at the factor ~2–3 level**; the normalisation changes radial membership by 7–16 points.

---

## 5. Results on Dorcha host 0

All numbers are for one host (group 0). Plots were produced from scratch scripts (see §8); they are not in the repo.

### 5.1 Satellite sample

* 2866 z=0 FOF members → **2159 satellites inside R200c** (2139 in the host FOF group, 20 in neighbouring groups at 0.90–0.99 R200c, M < 8e7 Msun). 727 FOF members lie outside R200c (out to 1.79 R200c).
* Above 1e8 Msun: 188; above 1e9 Msun: 12 (26 FOF members).

| Threshold | z=0 / final | Peak |
|---|---|---|
| ≥ 50 particles | 1469 | 2126 |
| ≥ 100 particles | 790 | 1804 |
| ≥ 300 particles | 270 | 904 |
| M ≥ 1e8 Msun | 188 | 656 |
| M ≥ 1e9 Msun | 12 | 58 |

70 % of satellites have a peak particle count at least twice their current one, so the basis matters.
**Configured samples:** 1e-4 × M200c = 2.17e8 Msun ≈ 934 particles, so the 50-particle cut is not binding for this host.
Counts sample (peak): **315**; orbit sample (final): **71**.

### 5.2 Orbits of the ten most massive satellites (a = 0.5–1.0)

Subhaloes 1, 2, 6, 8, 12, 16, 17, 18, 21, 22. Object 1 (2.7e11 Msun) falls in late (a≈0.93) and has not reached pericentre; object 12 has three pericentres, one at 11 kpc;
objects 6 and 21 reach 14 and 17 kpc. Pericentres this deep are only a few softening lengths in and are less secure. (An earlier top-ten built from FOF membership included four objects, 3, 4, 7, 9, that
sit outside R200c at z=0; 3 and 4 never crossed, 7 and 9 are backsplash.)

### 5.3 Mpeak vs Mcurrent (counts sample, N = 315)

* Median Mcurrent/Mpeak = **0.21** (16–84 %: 0.06–0.51); 28 % have lost more than 90 %.
* **77 % (244 objects) were above the mass cut at peak but are below it at z=0.**
* **SubhaloNr 104** is the outlier: peak 5.2e11 Msun (2.2 million particles) as the central of its own FOF group, merged into the host at a≈0.67, now 2.5e8 Msun (1062 particles). This is a roughly 1:2 merger remnant, not a tracking glitch. The host's M200c doubled between z≈0.3 and z=0, so **host 0 is a post-major-merger system** and would fail a "quiet merger history" MWA criterion.

### 5.4 Infall parameters and pericentres of the satellites

(In R200c units with peculiar velocities; superseded for comparison purposes by §6.)
Orbit sample, first R200c crossing: median |v|/V200c = 1.31, v_r/V200c = −0.91, v_t/V200c = 0.78, infall angle 39°, j = L/(R200c V200c) = 0.78.
Deepest pericentre vs current radius: median r_now/r_peri = 3.0 (most recent pericentre: 1.95); 20 of the 71 orbit-sample satellites (46 of 315 in the counts sample) have no confirmed pericentre yet.

---

## 6. Comparison with Benson (2005)

The paper (read in full) differs from the first set of plots in four ways: r_vir/V_vir from the top-hat overdensity (Δ≈102 ρ_crit at z=0, not 200); **Hubble flow added** to peculiar velocities (footnote 3); **every** merger counted, not z=0 survivors; V_r positive inward.
`orbits/benson.py` implements these. Dorcha host 0: r_vir = 357 kpc = 1.307 R200c (1.30–1.34 for c = 5–20, so the NFW assumption is not driving it).
Sample: all first r_vir crossings with ≥ 50 particles, z < 1.45, not pre-processed (subhalo of a different FOF group): 3002 events.

Each cell: mean V_r / mean V_θ / unbound fraction (E > 0). Velocities are in units of V_vir (V200c for the R200c columns).

| z | r_vir + Hubble (Benson's definition) | r_vir, no Hubble | R200c + Hubble | R200c, no Hubble | Benson fit |
|---|---|---|---|---|---|
| ≈0 | 0.87 / 0.81 / 0.30 | 1.00 / 0.81 / 0.42 | 0.90 / 0.85 / 0.37 | 1.00 / 0.85 / 0.46 | 0.90 / 0.75 / 0.20 |
| ≈0.5 | 0.85 / 0.86 / 0.33 | 0.97 / 0.86 / 0.42 | 0.87 / 0.90 / 0.38 | 0.97 / 0.90 / 0.44 | 0.89 / 0.71 / 0.15 |
| ≈1 | 0.83 / 0.81 / 0.16 | 0.94 / 0.81 / 0.28 | 0.84 / 0.85 / 0.17 | 0.94 / 0.85 / 0.31 | 0.85 / 0.70 / 0.11 |

* **V_r agrees** with the fit in shape and position when the Benson definitions are used.
* **V_θ is systematically higher** in every variant (+8 % at z≈0, ~+15–20 % at z≈0.5–1); the Dorcha distributions peak near 0.9–1.0 where the fit peaks near 0.65.
* **Removing the Hubble flow and switching to R200c both move the results *away* from the fit** (V_r up by 0.10–0.13, unbound fraction up by ~0.12, V_θ up by 0.03–0.04). They do not explain the V_θ offset.
* Survivors (inside r_vir at z=0) are *less* energetic than the rest (median V_r 0.76 vs 0.94; 17 % vs 37 % unbound): the survivor bias is mostly about unbound fly-bys that leave, not angular momentum.
* The z≈0.5 window contains the host's major merger; the unusually large unbound fraction there (0.33 vs 0.15) has not been tested against it.

**Untested candidate causes of the V_θ excess:** halo centre/velocity definition (Benson uses most-bound-particle position and velocity; HBT subhalo centre here); one host vs his ensemble (he reports mass trends); NFW-derived r_vir; his two-body back-propagation vs our interpolated actual path.

---

## 7. Mistakes made and corrected during the work

An honest log; each item says what was affected.

1. **Velocity convention default.** `TrackUnits(velocity="auto")` assumed `v/√a` whenever `comoving=True`. For the Dorcha tree (native comoving, velocities physical) that would have put a ~13 % error in every Hermite slope, E and L. Fixed: `from_epoch` derives the mode from what the tree reader actually did. Caught by testing against the real tree before any science use.
2. **Arbitrary 30 % NFW acceptance cut** discarded 70 of 123 host snapshots. Made explicit (default 50 %), with the median mismatch logged.
3. **Orphan model applied to fringe objects** sent 68 % of lost subhaloes beyond R200c (to ~1 Mpc). Gated by `orphan_max_r_over_ref`; discovered by checking the output against the tree itself.
4. **Wrong host selection in the first real-data run** (picked group 1 by group M200); fixed by selecting on central-subhalo mass, which exposed the group-1 anomaly.
5. **"Central at infall" flag too strict** (excluded objects already in the host's *own* FOF group). The first Benson comparison therefore used a biased subset; the unbound fraction quoted for it (0.17 at z≈0) is unreliable. Corrected flag: not a subhalo of a *different* group (92–95 % pass).
6. **Over-stated claim in a plot caption:** objects 7 and 9 were described as never entering the virial sphere; that was only true in the a = 0.8–1 window (they are backsplash).
7. **Earlier "1.57e5 Msun per particle" statement** conflated Msun/h and Msun (2.32e5 Msun).
8. **Plot artefacts:** pre-infall shallow minima plotted as "pericentres" (now post-infall only); a 45° blip in the Benson φ curve is a square-grid binning artefact.
9. **Environment trap:** the system `python3.11` carries an editable install of an *old* copy at `~/Dropbox/AnalysisTools`. Tests were run in a clean venv with `PYTHONPATH=~/MyCodes/AnalysisTools`.

---

## 8. Reproducing and where things are

* **Environment:** clean venv (`numpy scipy h5py numba matplotlib plotly pyyaml pytest`), then
  `PYTHONPATH=$PWD python -m pytest tests -q` from the repo root (~25 s).
* **Dorcha tree:** `~/Currentwork/Dorcha/trees.hdf5`; load with `MergerTreeTools(path, "SubFind", comoving=True, little_h=True, native_includes_h=True, native_is_comoving=True)`.
* **Settings:** `configs/dorcha.yaml` (`orbital_properties` block: units, `radius_key: R200c_from_M200`, resolution cuts, bases, orphan gating).
* **Test sidecar:** `~/Currentwork/Dorcha/orbits_host0_test.h5` (written before several later fields existed; the reader fills missing scalars).
* **Analysis scripts:** `scripts/dorcha_orbits/` (see its `README.md`): `01`-`11`, each runnable from the repo root as
  `PYTHONPATH=. python scripts/dorcha_orbits/<script>.py`; tree and output locations via `DORCHA_TREE` / `DORCHA_OUT`.
  All were re-run from there and reproduce the numbers quoted in this report. Figures are written to `DORCHA_OUT` and
  are not committed (images are git-ignored).
* **Memory notes** (Claude's persistent notes) record the decisions, verified facts and open items.

---

## 9. Open items

| # | Item | Needs |
|---|---|---|
| 1 | Host assignment for paired MWAs / pre-processed satellites (z=0 host's main progenitor vs host at first infall) | Decision |
| 2 | Orphan positions → Shark type-2 galaxies (box coordinates and matching rule not wired) | Decision, then code |
| 3 | Dynamical friction on orphans (off; matters for massive orphans) | Decision |
| 4 | "Pre-processed" definition: currently any other-FOF membership before infall (62 % of z=0 satellites, an upper bound); needs a mass-ratio / minimum-duration cut | Decision |
| 5 | `HaloPropertiesStage` `Minfall`/`IsBacksplash` are still FOF-based although the 1.0 schema text says "R200c crossing"; unit-consistent mass conversion first | Decision |
| 6 | Particle-measured potential per snapshot (`TabulatedSpherical.from_particles` exists; the time-dependent builder does not); the one change that would reduce the orphan systematics | Dorcha snapshots |
| 7 | Halo-centre test of the V_θ excess (most-bound-particle velocity at infall) | Dorcha snapshots |
| 8 | Phase 2: particle tracking (progenitor-field fractions, per-particle pericentres) | Design |
| 9 | What is group 1 (fragmented FOF? zoom-boundary contamination?) and a merger-history criterion for the MWA selection (host 0 had a ~1:2 merger at z≈0.5) | Investigation |
| 10 | The `ssec:pericentre` section in both paper drafts | Writing, once the definitions above are settled |
| 11 | Everything here is for **one host**; no isolated-vs-paired statement is possible yet | More hosts |
