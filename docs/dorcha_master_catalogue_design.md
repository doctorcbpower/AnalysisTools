# Dorcha Master Science Catalogue — Design Document

**Author:** Architecture proposal for the Dorcha project
**Date:** 2026-08-06
**Status:** Design only — no implementation code included by request

---

## 0. Framing

The request is effectively: *design a long-lived scientific database, not a file*. That distinction drives every choice below. A "convenient" design optimises for the two papers in front of you. A durable design optimises for the tenth paper you haven't thought of yet, written by a collaborator who never ran your pipeline. Concretely, that means:

- The schema is organised around **physical objects and physical quantities**, not around any one analysis.
- Every expensive computation is **cached exactly once**, with enough provenance metadata that a skeptical reviewer (or your future self) can audit how a number was produced.
- The catalogue is **append-only by convention**: new columns and new groups can be added without altering the meaning or layout of existing ones, so old notebooks never break.
- Units, coordinate conventions, and cosmology are **stored, not assumed** — nothing downstream should ever hard-code $h=0.6774$ in a script.

---

## 1. Overall Architecture

### 1.1 Three-layer model

```
┌─────────────────────────────────────────────────────────────────┐
│  Layer 3 — Science Analysis                                      │
│  Notebooks / scripts that read ONLY the master catalogue.        │
│  No access to snapshots, trees, Shark, or particle files.        │
└─────────────────────────────────────────────────────────────────┘
                              ▲  read-only
┌─────────────────────────────────────────────────────────────────┐
│  Layer 2 — Master Catalogue (HDF5)                                │
│  dorcha_catalogue_v{MAJOR.MINOR}.h5                               │
│  Metadata / Cosmology / Haloes / Satellites / DerivedQuantities /│
│  Environment / Observability / Documentation / Provenance         │
└─────────────────────────────────────────────────────────────────┘
                              ▲  write-once per build
┌─────────────────────────────────────────────────────────────────┐
│  Layer 1 — Construction Pipeline (dorcha_catalogue package)       │
│  Extract → Cross-match → Compute derived quantities →             │
│  Quality control → Write                                          │
│  Reads: Dorcha halo catalogues, merger trees, Shark outputs,      │
│  particle tagging, progenitor fields, selection functions,        │
│  Rubin detectability tables                                        │
└─────────────────────────────────────────────────────────────────┘
```

The pipeline (Layer 1) is the only code permitted to touch raw simulation products. Everything above the master catalogue line is forbidden from opening a snapshot, a tree file, or a Shark HDF5 output directly — this is enforced by convention (documented in the README) and, if you want teeth, by giving science-analysis environments no filesystem path to the raw data at all.

### 1.2 Versioning philosophy

Because "at least ten future papers" will hang off this file, treat the catalogue like a public dataset with a release cycle, not a scratch product:

- **Semantic versioning on the file name and inside `Metadata`**: `MAJOR.MINOR.PATCH`. MAJOR bumps mean a breaking schema change (renamed/removed column, changed units, changed row ordering). MINOR bumps mean additive changes (new column, new group). PATCH bumps mean a bugfix rerun with no schema change.
- **Every published plot records the catalogue version it used**, ideally in a per-figure provenance stub (see §6). This lets you regenerate Paper I's Figure 3 in three years even after the schema has grown.
- Old catalogue versions are never deleted; they are moved to a `releases/` archive. Science code pins a version explicitly rather than reading "latest".

### 1.3 Why HDF5 (and where it stops being enough)

HDF5 is the right choice for this problem: self-describing (units/descriptions live next to the data), efficient columnar-ish access for large flat arrays, native chunking/compression, and effectively universal in astrophysics tooling (h5py, `unyt`, `astropy.table`, `pandas`, Julia, IDL legacy code if it ever comes to that). It is not a substitute for a real database once you need concurrent multi-user writes, transactional updates, or SQL-style joins across many catalogues at once — if Dorcha ever needs live, queryable access for a large collaboration, plan for a thin layer (e.g. exporting to Parquet + DuckDB, or loading into a proper DB) *on top of* the canonical HDF5 master file, rather than replacing it. For a single-PI-group, versioned-release workflow, HDF5 alone is appropriate.

---

## 2. HDF5 Schema

### 2.1 Top-level layout

```
dorcha_catalogue_v1.0.0.h5
│
├── Metadata/                 # file-level bookkeeping, versioning, run config
├── Provenance/                # per-dataset lineage, code versions, checksums
├── Cosmology/                 # cosmological parameters used throughout
├── SimulationConfig/           # Dorcha suite configuration (box, resolution, softening)
├── Snapshots/                  # snapshot number ↔ redshift ↔ scale factor ↔ time lookup
├── Haloes/                     # host halo (Milky-Way-analogue) table, one row per host per snapshot
├── Satellites/                 # the main object table, one row per satellite (z=0, unless noted)
│   ├── Identification/
│   ├── HaloProperties/
│   ├── GalaxyProperties/
│   ├── StarFormationHistory/    # ragged/2D arrays, indexed by SatelliteID
│   ├── DorchaProperties/
│   ├── ParticleTags/             # variable-length particle ID lists, indexed by SatelliteID
│   ├── Environment/
│   ├── Observability/
│   └── DerivedQuantities/
├── MergerTrees/                 # optional companion structural info (branch/pointer arrays)
├── SelectionFunctions/           # observational survey completeness/selection tables
├── CrossMatches/                  # cross-match tables to external catalogues (Gaia, DES, LSST sims...)
├── MLClassifications/              # model outputs (e.g. quenched/star-forming classifier), versioned
└── Documentation/                    # human-readable schema doc mirrored inside the file itself
```

Rationale for a few structural choices:

- **`Satellites/` is subdivided into subgroups**, not one flat table, because the satellite table will have 150+ columns by the time ten papers have contributed to it. Subgroups by physical domain (halo, galaxy, environment, observability, derived) keep `h5ls` output navigable and let you grant "add a column" permission to a subgroup without touching the rest.
- **`DerivedQuantities/` exists at two levels**: catalogue-wide derived quantities (e.g. host-level statistics) live at the top, and satellite-level derived quantities live inside `Satellites/DerivedQuantities/`. Anything computed from stored quantities via a documented, reproducible recipe belongs here — never as a "primary" input.
- **Ragged data (star formation histories, particle ID lists) do not live in the flat per-satellite table.** HDF5 doesn't do variable-length columns gracefully inside a table; instead each satellite's SFH is a fixed-length vector (binned to a common age/time grid defined in `Metadata/`) stored in a 2D array `[N_satellites, N_time_bins]`, indexed by row position identical to the main table. Particle tags, which are genuinely variable-length, use HDF5 variable-length dtypes or a separate index+concatenated-array pattern (`ParticleTags/offsets`, `ParticleTags/particle_ids`).
- **All tables share one canonical row order**, `SatelliteID`-sorted, fixed at write time. Every subgroup array under `Satellites/*` has the same length and the same row-to-object mapping, so joining columns across subgroups is a matter of `h5py` slicing, not a merge/join operation. This is the single most important design invariant in the file — document it loudly and enforce it in QC.

### 2.2 Dataset catalogue

Below, "P" = primary (extracted directly from a source with no computation beyond unit conversion), "D" = derived (computed by the pipeline from other stored quantities). Units follow the `unyt`/astropy convention (stored as an HDF5 attribute `units` on each dataset, plus mirrored in `Documentation/`).

#### `Metadata/` (attributes, not datasets, mostly)

| Name | Type | Units | Description | Provenance | P/D |
|---|---|---|---|---|---|
| `catalogue_version` | str attr | — | Semantic version `MAJOR.MINOR.PATCH` | pipeline config | — |
| `creation_date` | str attr (ISO 8601) | — | UTC timestamp of build | pipeline | — |
| `pipeline_git_commit` | str attr | — | Git SHA of `dorcha_catalogue` at build time | pipeline | — |
| `random_seed` | int attr | — | Seed used for any stochastic step (e.g. MC realisations of orbital params) | pipeline config | — |
| `n_satellites` | int attr | — | Row count of `Satellites/` table | pipeline | D |
| `n_hosts` | int attr | — | Number of independent Milky-Way analogues | pipeline | D |
| `schema_version` | str attr | — | Schema (not data) version, tracked separately from catalogue_version | pipeline | — |

#### `Cosmology/`

| Name | Type | Units | Description | Provenance | P/D |
|---|---|---|---|---|---|
| `H0` | float64 | km/s/Mpc | Hubble constant | Dorcha sim ICs | P |
| `Omega_m` | float64 | dimensionless | Matter density parameter | Dorcha sim ICs | P |
| `Omega_lambda` | float64 | dimensionless | Dark energy density parameter | Dorcha sim ICs | P |
| `Omega_b` | float64 | dimensionless | Baryon density parameter | Dorcha sim ICs | P |
| `sigma8` | float64 | dimensionless | Power spectrum normalisation | Dorcha sim ICs | P |
| `n_s` | float64 | dimensionless | Spectral index | Dorcha sim ICs | P |
| `little_h` | float64 | dimensionless | $h = H_0/100$, used for unit conversions | derived from H0 | D |

#### `SimulationConfig/`

| Name | Type | Units | Description | Provenance | P/D |
|---|---|---|---|---|---|
| `box_size` | float64 | Mpc/h | Simulation volume side length | Dorcha config | P |
| `particle_mass_dm` | float64 | Msun/h | Dark matter particle mass | Dorcha config | P |
| `softening_length` | float64 | kpc/h | Gravitational softening | Dorcha config | P |
| `n_particles` | int64 | — | Total DM particle count | Dorcha config | P |
| `halo_finder` | str | — | e.g. `ROCKSTAR`, `SUBFIND` | Dorcha pipeline docs | P |
| `tree_builder` | str | — | Merger tree code and version | Dorcha pipeline docs | P |
| `shark_version` | str | — | Shark git tag/commit used | Shark run config | P |
| `shark_parameter_file` | str | — | Path/checksum of the Shark parameter file used | Shark run config | P |

#### `Snapshots/`

| Name | Type | Units | Description | Provenance | P/D |
|---|---|---|---|---|---|
| `snapshot_number` | int32[N_snap] | — | Snapshot index | Dorcha sim | P |
| `redshift` | float64[N_snap] | dimensionless | Redshift at snapshot | Dorcha sim | P |
| `scale_factor` | float64[N_snap] | dimensionless | $a = 1/(1+z)$ | derived | D |
| `cosmic_time` | float64[N_snap] | Gyr | Age of universe at snapshot | derived from cosmology | D |
| `time_bin_edges_sfh` | float64[N_bins+1] | Gyr | Common lookback-time grid used for all stored SFHs | pipeline config | — |

#### `Haloes/` (host, i.e. Milky-Way-analogue, table — one row per host halo per stored snapshot, or z=0 only if history isn't needed at host level)

| Name | Type | Units | Description | Provenance | P/D |
|---|---|---|---|---|---|
| `HostHaloID` | int64[N_host] | — | Unique host identifier, stable across snapshots | Dorcha halo catalogue | P |
| `Snapshot` | int32[N_host] | — | Snapshot of this row | Dorcha halo catalogue | P |
| `M200c` | float64[N_host] | Msun | Host virial mass (200c definition) | Dorcha halo catalogue | P |
| `R200c` | float64[N_host] | kpc | Host virial radius | Dorcha halo catalogue | P
| `Vmax_host` | float64[N_host] | km/s | Peak circular velocity of host | Dorcha halo catalogue | P |
| `Position` | float64[N_host,3] | Mpc/h (comoving) | Host position | Dorcha halo catalogue | P |
| `Velocity` | float64[N_host,3] | km/s (peculiar) | Host bulk velocity | Dorcha halo catalogue | P |
| `IsIsolated` | bool[N_host] | — | Passes isolation criterion (no more massive halo within N×R200c) | derived, isolation algorithm v1 | D |
| `IsPaired` | bool[N_host] | — | Local Group–analogue pairing flag | derived, pairing algorithm v1 | D |
| `PairedHostID` | int64[N_host] | — | ID of paired companion, −1 if none | derived | D |
| `N_satellites_total` | int32[N_host] | — | Count of satellites above completeness limit | derived | D |
| `AccretionHistory_M200c` | float64[N_host,N_snap] | Msun | Host mass growth history | derived from merger tree | D |

#### `Satellites/Identification/`

| Name | Type | Units | Description | Provenance | P/D |
|---|---|---|---|---|---|
| `SatelliteID` | int64[N] | — | Unique, permanent ID; primary key of the whole table | assigned at catalogue construction | — |
| `HostHaloID` | int64[N] | — | Foreign key into `Haloes/` | cross-match stage | D |
| `MergerTreeID` | int64[N] | — | Foreign key into `MergerTrees/` | tree builder | P |
| `SubhaloID_z0` | int64[N] | — | Native subhalo finder ID at z=0 | Dorcha halo catalogue | P |
| `SharkGalaxyID` | int64[N] | — | Foreign key into Shark output | cross-match stage | D |
| `Snapshot` | int32[N] | — | Snapshot of the reported (z=0, typically) properties | pipeline | P |
| `SatelliteName` | str[N] (optional, only if matched to real MW satellites) | — | Common name if this is a mock built to represent a known dwarf | manual/observational cross-match | P |

#### `Satellites/HaloProperties/`

| Name | Type | Units | Description | Provenance | P/D |
|---|---|---|---|---|---|
| `M200c_z0` | float64[N] | Msun | Present-day virial mass | Dorcha halo catalogue | P |
| `Mpeak` | float64[N] | Msun | Peak halo mass over full history | derived from merger tree | D |
| `SnapshotAtMpeak` | int32[N] | — | Snapshot at which Mpeak occurred | derived | D |
| `R200c_z0` | float64[N] | kpc | Present-day virial radius | Dorcha halo catalogue | P |
| `Vpeak` | float64[N] | km/s | Peak circular velocity over history | derived from merger tree | D |
| `Vmax_z0` | float64[N] | km/s | Present-day peak circular velocity | Dorcha halo catalogue | P |
| `Minfall` | float64[N] | Msun | Halo mass at infall (first crossing of host R200c) | derived from merger tree | D |
| `Vinfall` | float64[N] | km/s | Vmax at infall | derived from merger tree | D |
| `RedshiftInfall` | float64[N] | dimensionless | Redshift of infall | derived from merger tree | D |
| `SnapshotInfall` | int32[N] | — | Snapshot of infall | derived from merger tree | D |
| `NumberOfInfalls` | int16[N] | — | Count of host-crossing events (backsplash tracking) | derived from merger tree | D |
| `IsBacksplash` | bool[N] | — | Currently outside R200c having previously fallen in | derived | D |
| `OrbitalPericentre` | float64[N] | kpc | Estimated pericentric distance | derived, orbit-fitting algorithm v1 | D |
| `OrbitalApocentre` | float64[N] | kpc | Estimated apocentric distance | derived, orbit-fitting algorithm v1 | D |
| `OrbitalEccentricity` | float64[N] | dimensionless | (apo−peri)/(apo+peri) | derived | D |
| `OrbitalPeriod` | float64[N] | Gyr | Estimated radial orbital period | derived, orbit-fitting algorithm v1 | D |
| `OrbitalEnergy` | float64[N] | (km/s)^2 | Specific orbital energy at infall | derived | D |
| `OrbitalAngularMomentum` | float64[N,3] | kpc km/s | Specific angular momentum vector at infall | derived | D |
| `TidalTrackClass` | int8[N] | — | Categorical tidal-stripping stage (enumerated in Documentation) | derived, tidal-tracks model | D |

#### `Satellites/GalaxyProperties/` (from Shark, plus derived)

| Name | Type | Units | Description | Provenance | P/D |
|---|---|---|---|---|---|
| `StellarMass` | float64[N] | Msun | Total stellar mass | Shark output | P |
| `LuminosityV` | float64[N] | Lsun (V-band) | V-band luminosity | Shark output / SPS model | P |
| `Luminosity_ugriz` | float64[N,5] | Lsun (per band) | Multi-band luminosities | Shark output / SPS model | P |
| `AbsoluteMagnitude_ugriz` | float64[N,5] | mag | Absolute magnitudes | derived from luminosities | D |
| `GasMass_Cold` | float64[N] | Msun | Cold gas mass | Shark output | P |
| `GasMass_Hot` | float64[N] | Msun | Hot gas mass (if resolved) | Shark output | P |
| `MetallicityStellar` | float64[N] | dimensionless (Z) or [Fe/H] | Mass-weighted stellar metallicity | Shark output | P |
| `MetallicityGas` | float64[N] | dimensionless | Cold gas-phase metallicity | Shark output | P |
| `StarFormationRate` | float64[N] | Msun/yr | Instantaneous SFR | Shark output | P |
| `SFH` | float32[N,N_bins] | Msun/yr per bin | Star formation history on common time grid | Shark output, rebinned | D |
| `MeanStellarAge` | float64[N] | Gyr | Mass-weighted mean stellar age | derived from SFH | D |
| `QuenchingTime` | float64[N] | Gyr (lookback) | Time SFR fell below threshold × sSFR_MS, permanently | derived from SFH | D |
| `IsQuenched_z0` | bool[N] | — | sSFR below threshold at z=0 | derived | D |
| `HalfLightRadius` | float64[N] | kpc | Projected half-light radius | Shark output / structural model | P |
| `HalfMassRadiusStellar` | float64[N] | kpc | 3D stellar half-mass radius | Shark output | P |
| `SersicIndex` | float64[N] | dimensionless | Structural fit parameter, if available | Shark output / image-fitting | P |
| `BlackHoleMass` | float64[N] | Msun | SMBH mass, if tracked | Shark output | P |

#### `Satellites/StarFormationHistory/` — implementation detail

Stored as `SFH` above (2D array on the common `time_bin_edges_sfh` grid from `Snapshots/`), not as a separate group, to preserve the one-row-per-object invariant. `Documentation/` records the binning scheme and interpolation method used to rebin from native Shark snapshot outputs.

#### `Satellites/DorchaProperties/`

| Name | Type | Units | Description | Provenance | P/D |
|---|---|---|---|---|---|
| `ProgenitorParticleFraction` | float64[N] | dimensionless (0–1) | Fraction of z=0 bound particles traced to a single progenitor | particle tagging + progenitor-field analysis | D |
| `EarliestProgenitorRedshift` | float64[N] | dimensionless | Highest-z snapshot at which a resolved progenitor exists | merger tree + particle tagging | D |
| `PeakOverdensity` | float64[N] | dimensionless (δ) | Peak local overdensity experienced by progenitor particles | progenitor-field analysis | D |
| `FormationEnvironmentClass` | int8[N] | — | Enumerated formation environment (e.g. filament/node/void progenitor origin) | derived, environment classifier | D |
| `FossilFraction` | float64[N] | dimensionless (0–1) | Fraction of stellar mass formed before reionisation (or other fossil definition, documented) | derived from SFH + reionisation model | D |
| `NumberOfMergers` | int16[N] | — | Count of resolved mergers in this object's own history | merger tree | D |
| `MassAccretionRateDM` | float64[N,N_snap] | Msun/Gyr | Dark matter accretion history | merger tree | D |

#### `Satellites/ParticleTags/`

| Name | Type | Units | Description | Provenance | P/D |
|---|---|---|---|---|---|
| `particle_id_offsets` | int64[N+1] | — | CSR-style offsets into `particle_ids` for each satellite | particle tagging | P |
| `particle_ids` | int64[N_total_tagged] | — | Concatenated tagged particle IDs across all satellites | particle tagging | P |
| `particle_tag_confidence` | float32[N_total_tagged] | dimensionless | Per-particle tagging confidence score, if the tagging method produces one | particle tagging | P |

This CSR (compressed sparse row) layout avoids ragged HDF5 datasets while keeping per-object particle lists cheaply sliceable: `particle_ids[offsets[i]:offsets[i+1]]` gives satellite `i`'s tagged particles.

#### `Satellites/Environment/`

| Name | Type | Units | Description | Provenance | P/D |
|---|---|---|---|---|---|
| `HostIsIsolated` | bool[N] | — | Copy-through of host's isolation flag, denormalised for query convenience | join from `Haloes/` | D |
| `HostIsPaired` | bool[N] | — | Copy-through of host's pairing flag | join from `Haloes/` | D |
| `LocalNumberDensity` | float64[N] | Mpc^-3 (comoving) | Number density of neighbouring haloes above a mass threshold within a fixed aperture | derived, environment estimator v1 | D |
| `DistanceToNearestMassiveNeighbour` | float64[N] | Mpc | Distance to nearest halo above threshold mass outside the host | derived | D |
| `TidalIndex` | float64[N] | dimensionless | Karachentsev-style tidal index from neighbouring haloes | derived | D |
| `CosmicWebClass` | int8[N] | — | Enumerated node/filament/sheet/void classification at satellite's position | derived, web classifier (e.g. T-web/V-web) | D |

#### `Satellites/Observability/`

| Name | Type | Units | Description | Provenance | P/D |
|---|---|---|---|---|---|
| `HeliocentricDistance` | float64[N] | kpc | Distance from mock solar position | derived, mock observer placement | D |
| `GalactocentricRadius` | float64[N] | kpc | 3D distance from host centre | derived | D |
| `GalacticLongitude` | float64[N] | deg | Galactic $\ell$ | derived, mock observer placement | D |
| `GalacticLatitude` | float64[N] | deg | Galactic $b$ | derived, mock observer placement | D |
| `RA`, `Dec` | float64[N] each | deg | Equatorial coordinates, if a specific mock sky orientation is fixed | derived | D |
| `ApparentMagnitude_ugriz` | float64[N,5] | mag | Apparent magnitudes given heliocentric distance | derived from luminosity + distance | D |
| `SurfaceBrightness_V` | float64[N] | mag/arcsec^2 | Mean surface brightness within half-light radius | derived from luminosity + half-light radius + distance | D |
| `RubinDetectionProbability` | float64[N] | dimensionless (0–1) | Detection probability from Rubin/LSST selection function model | Rubin detectability calculation | D |
| `RubinCoaddDepth_r` | float64[N] | mag | Assumed local coadd depth used in the detection calculation | selection function input | P |
| `CompletenessWeight` | float64[N] | dimensionless | $1/P_\mathrm{detect}$, capped, for completeness-corrected statistics | derived | D |
| `SurveyFootprintFlag` | bool[N] | — | Falls within the assumed survey footprint | selection function | D |

#### `Satellites/DerivedQuantities/` (cross-cutting quantities not naturally owned by one subgroup above)

See §3 for the full recommended list; representative entries:

| Name | Type | Units | Description | Provenance | P/D |
|---|---|---|---|---|---|
| `StellarMassHaloMassRatio` | float64[N] | dimensionless | $M_\star/M_\mathrm{peak}$ | derived | D |
| `MassLossFraction` | float64[N] | dimensionless | $1 - M_{200c,z0}/M_\mathrm{peak}$ | derived | D |
| `DynamicalMassWithinHalfLight` | float64[N] | Msun | Mass estimator within half-light radius (e.g. Wolf et al. estimator) | derived | D |
| `MeanDensityWithinHalfLight` | float64[N] | Msun/kpc^3 | Companion to above | derived | D |
| `IsUltraFaint` | bool[N] | — | $L_V < 10^5 L_\odot$ convenience flag | derived | D |
| `IsClassical` | bool[N] | — | Convenience flag matching classical-dSph luminosity range | derived | D |

#### `MergerTrees/`

| Name | Type | Units | Description | Provenance | P/D |
|---|---|---|---|---|---|
| `TreeID` | int64[N_tree] | — | Unique tree identifier | tree builder | P |
| `NodeID`, `DescendantID`, `FirstProgenitorID`, `NextProgenitorID` | int64[N_nodes] each | — | Standard tree-pointer arrays (consumer/community-format compatible) | tree builder | P |

Kept as a companion group rather than folded entirely into `Satellites/` because trees are graph-structured, not one-row-per-object; the derived scalars extracted from them (Mpeak, infall time, etc.) live in `Satellites/` while the full graph stays here for anyone who needs to walk it directly (rare, but keeps the "one input for all science" promise honest without bloating the main table).

#### `SelectionFunctions/`, `CrossMatches/`, `MLClassifications/`

Each is a **namespaced, versioned subgroup per external dataset or model**, e.g. `MLClassifications/quench_classifier_v2/`, containing its own prediction array (length N, aligned to `Satellites/` row order), a confidence/probability array, and attributes recording model version, training date, and training-set provenance. This is the extensibility hook for point 6 in the prompt (ML classifications) without touching the core schema.

#### `Documentation/`

| Name | Type | Description |
|---|---|---|
| `schema.json` or `schema.yaml` (stored as a string dataset) | str | Machine-readable mirror of every table above — name, dtype, units, description, provenance, P/D flag — so tooling can validate the file against its own spec |
| `README.md` (stored as a string dataset) | str | Human-readable version of this document, versioned alongside the schema |
| `enumerations.json` | str | Definitions of all categorical codes (`TidalTrackClass`, `FormationEnvironmentClass`, `CosmicWebClass`, etc.) |
| `changelog.md` | str | Human-readable per-version changelog |

Storing documentation *inside* the HDF5 file (in addition to the repo) means the file remains self-describing even if separated from the codebase — critical for a ten-paper, multi-year lifespan and for sharing with collaborators or a journal data repository.

### 2.3 Units convention

Every dataset carries an HDF5 attribute `units` (a string parseable by `astropy.units` or `unyt`), plus `description`, `provenance`, and `is_derived` (bool) attributes attached directly to the dataset — not only in the external `Documentation/schema.json` mirror. This means a user who opens the file with bare `h5py` and no other context can still recover full metadata: `f['Satellites/HaloProperties/Mpeak'].attrs['units']`. Redundancy between in-file attributes and `Documentation/schema.json` is intentional; the schema.json is what automated validation checks against, the attrs are what a human or ad hoc script encounters first.

---

## 3. Recommended Derived Quantities (beyond the two current papers)

Grouped by the kind of future paper likely to need them:

**Structure/abundance-matching papers:** peak mass, infall mass/redshift/Vmax, mass loss fraction, stellar-mass–halo-mass ratio, number of mergers, accretion rate histories, backsplash flag and count.

**Chemical evolution / stellar populations papers:** mean stellar age, age spread (e.g. IQR of formation time from SFH), quenching time, fossil fraction, metallicity, [Fe/H]-luminosity relation inputs, alpha-enhancement if Shark tracks it.

**Dynamics/orbits papers:** pericentre, apocentre, eccentricity, orbital period, orbital energy and angular momentum, number of pericentric passages, tidal track classification, dynamical mass and density within half-light radius (for dark-matter-content and core/cusp studies).

**Environment/large-scale-structure papers:** cosmic web classification, tidal index, local density, isolation/pairing flags, peak overdensity of progenitor particles, formation environment class.

**Observational-strategy / survey-forecasting papers:** heliocentric distance, apparent magnitude, surface brightness, Rubin detection probability, completeness weight, survey footprint flag — precisely so that luminosity functions, radial profiles, and detection-efficiency curves are *always* reconstructable from stored quantities rather than baked in as static plots (per Principle 3 in the prompt).

**Machine-learning / classification papers:** any model prediction goes into `MLClassifications/<model_name_vX>/`, never overwriting or duplicating an existing physical column — a classifier's output is data provenance-tagged like anything else.

**Cross-simulation comparison papers (future hydro runs):** design every derived quantity's *definition* (e.g. "peak mass" = "maximum $M_{200c}$ over the main branch") to be resolution- and simulation-independent, and record the *definition version* as an attribute, so a future hydro-Dorcha catalogue can reuse identical derived-quantity recipes and be safely concatenated or compared.

A practical rule: if you can imagine wanting it in a plot, and it takes more than a few seconds to compute from a merger tree or SFH, compute it once and store it. If it's a one-line combination of two already-stored columns and cheap to redo (e.g. absolute magnitude from luminosity), it's a judgment call — storing it anyway costs little and saves every downstream user from re-deriving the same constant/formula slightly differently.

---

## 4. Class Design

Principle: separate **schema knowledge** (what fields exist, their units, their validation rules) from **I/O mechanics** (how they're read/written) from **pipeline logic** (how they're computed). This is the layering that lets you extend the catalogue for ten years without a rewrite.

```
dorcha_catalogue/
│
├── schema/
│   ├── FieldSpec          — dataclass: name, dtype, units, description, provenance, is_derived, group_path
│   ├── GroupSpec          — dataclass: group path + ordered list of FieldSpec
│   └── CatalogueSchema    — the full schema as a Python object, loadable from / dumpable to
│                            Documentation/schema.json; single source of truth used by both
│                            the writer and the validator
│
├── io/
│   ├── CatalogueReader    — opens an HDF5 file read-only, exposes .satellites, .haloes, .cosmology
│   │                        as unit-aware table objects (thin wrapper, e.g. returning
│   │                        astropy.table.QTable views with units attached from attrs)
│   ├── CatalogueWriter    — accumulates in-memory arrays per FieldSpec, validates against
│   │                        CatalogueSchema, writes chunked+compressed datasets with attrs,
│   │                        enforces "write-once" (refuses to silently overwrite a released file)
│   └── ProvenanceRecorder — captures git commit, config hash, timestamps; writes Provenance/
│
├── pipeline/
│   ├── PipelineStage (ABC) — .name, .inputs (declared), .outputs (declared), .run(context) -> context
│   ├── ExtractStage         — one subclass per raw data source (HaloExtractor, TreeExtractor,
│   │                          SharkExtractor, ParticleTagExtractor, SelectionFunctionExtractor)
│   ├── CrossMatchStage      — resolves SatelliteID <-> HostHaloID <-> SharkGalaxyID <-> MergerTreeID
│   ├── DerivedQuantityStage — one subclass per physically-related family (OrbitalPropertiesStage,
│   │                          SFHDerivedStage, EnvironmentClassificationStage,
│   │                          RubinDetectabilityStage, ...), each pure-function-like: takes
│   │                          the accumulated context, returns new columns, no hidden state
│   ├── QualityControlStage — runs the validators in §5, raises or flags depending on severity
│   └── PipelineRunner      — orchestrates an ordered list of stages, passes a shared
│                              PipelineContext, logs timing/provenance per stage, supports
│                              "run from stage N" for iterative development
│
├── validation/
│   ├── Validator (ABC)      — .check(catalogue) -> ValidationReport
│   ├── SchemaValidator      — checks every dataset against CatalogueSchema (units, dtype, presence)
│   ├── IntegrityValidator   — duplicate IDs, missing values, foreign-key resolution
│   ├── PhysicalValidator    — monotonicity of accretion histories, impossible halo histories,
│   │                          Mpeak >= M200c_z0, orbital energy sign conventions, etc.
│   └── ValidationReport     — structured pass/fail/warn list, human-readable + machine-readable
│
├── units/
│   └── unit conventions built on astropy.units or unyt; every FieldSpec.units string is
│       validated to be parseable at schema-definition time, not at plot time
│
└── analysis/  (thin — most science code lives in separate paper-specific repos that
                depend on dorcha_catalogue as a library)
    └── convenience selectors, e.g. by_host(), isolated_only(), above_completeness(),
        that return filtered QTable views — no plotting, no statistics, just selection
```

Design notes:

- `PipelineStage` subclasses declare their `inputs`/`outputs` as explicit sets of field names. `PipelineRunner` can then statically check that a stage's declared inputs are already present in the context before running it — this catches "I forgot to run cross-match before derived quantities" at pipeline-assembly time instead of with a cryptic `KeyError` mid-run.
- `CatalogueWriter` is the *only* class permitted to write to the master file. It refuses to write if the target path already exists and is tagged as a release version, forcing intentional version bumps.
- `CatalogueReader` never exposes raw `h5py.Dataset` objects for satellite tables directly to users — it hands back `astropy.table.QTable` (or `unyt_array`-backed structures) with units already attached, so a scientist writing analysis code gets unit-aware quantities by default and cannot accidentally mix kpc and Mpc.
- Keep `DerivedQuantityStage` subclasses small and single-purpose (one physical family each). This is what lets ten future papers each contribute "their" derived-quantity stage as a self-contained, reviewable pull request instead of everyone patching one giant function.

---

## 5. Recommended Python Package Layout

```
dorcha_catalogue/                     # installable package, versioned independently of the data
├── pyproject.toml
├── README.md
├── CHANGELOG.md
├── src/
│   └── dorcha_catalogue/
│       ├── __init__.py
│       ├── schema/
│       │   ├── fields.py
│       │   ├── groups.py
│       │   └── catalogue_schema.py
│       ├── io/
│       │   ├── reader.py
│       │   ├── writer.py
│       │   └── provenance.py
│       ├── pipeline/
│       │   ├── base.py
│       │   ├── extract/
│       │   │   ├── haloes.py
│       │   │   ├── merger_trees.py
│       │   │   ├── shark.py
│       │   │   ├── particle_tagging.py
│       │   │   └── selection_functions.py
│       │   ├── crossmatch/
│       │   │   └── crossmatch.py
│       │   ├── derived/
│       │   │   ├── halo_properties.py
│       │   │   ├── orbital_properties.py
│       │   │   ├── star_formation_history.py
│       │   │   ├── environment.py
│       │   │   ├── observability.py
│       │   │   └── dorcha_specific.py    # fossil fraction, progenitor fraction, etc.
│       │   ├── quality_control.py
│       │   └── runner.py
│       ├── validation/
│       │   ├── base.py
│       │   ├── schema_checks.py
│       │   ├── integrity_checks.py
│       │   └── physical_checks.py
│       ├── units/
│       │   └── conventions.py
│       ├── analysis/
│       │   └── selectors.py
│       └── cli.py                     # `dorcha-build`, `dorcha-validate`, `dorcha-diff`
├── configs/
│   ├── pipeline_default.yaml           # stage order, thresholds, paths to raw data
│   └── schema_v1.0.yaml                # canonical schema definition (source for CatalogueSchema)
├── tests/
│   ├── unit/                            # one test module per module above
│   ├── integration/                     # small synthetic mock catalogue, run full pipeline on it
│   └── regression/                      # golden-file comparisons against a frozen mini catalogue
├── docs/
│   ├── schema_reference.md              # auto-generated from schema_v1.0.yaml
│   ├── pipeline_architecture.md
│   └── changelog/
└── notebooks/
    └── examples/                        # "how to open the catalogue and make a plot" — no
                                          # science conclusions, just usage patterns
```

Separate, downstream repos (`dorcha-paper1-analysis`, `dorcha-paper2-analysis`, ...) each `pip install dorcha_catalogue==X.Y` and contain only plotting/statistics code — reinforcing Principle 4: those repos should have no code path that opens a snapshot, tree, or Shark file.

---

## 6. Data Flow Diagram

```
 Raw Dorcha halo catalogues ─┐
 Merger trees ────────────────┤
 Shark galaxy catalogues ─────┼──▶  [ EXTRACT ]  ──▶  per-source tidy tables, unit-tagged,
 Particle tagging outputs ────┤                        indexed by native IDs
 Progenitor-field analysis ───┤
 Selection functions ─────────┤
 Rubin detectability calc ────┘
                                        │
                                        ▼
                              [ CROSS-MATCH ]
                     resolves SatelliteID <-> HostHaloID
                     <-> SharkGalaxyID <-> MergerTreeID
                     outputs: one unified index + per-source
                     lookup tables, cross-match QC report
                                        │
                                        ▼
                         [ COMPUTE DERIVED QUANTITIES ]
              (independent sub-stages, run in dependency order:
               halo-history-derived → orbit-derived → SFH-derived
               → environment-derived → observability-derived)
              outputs: full in-memory column set for Satellites/*
                                        │
                                        ▼
                            [ QUALITY CONTROL ]
                  schema / integrity / physical validators (§7)
                  outputs: ValidationReport; hard failures halt
                  the build, warnings are logged and attached
                  to Provenance/
                                        │
                                        ▼
                        [ WRITE MASTER CATALOGUE ]
                CatalogueWriter + ProvenanceRecorder produce
                dorcha_catalogue_vX.Y.Z.h5 (immutable release)
                                        │
                                        ▼
                          [ SCIENCE ANALYSIS ]
                  CatalogueReader only; notebooks/scripts in
                  paper-specific repos; zero raw-data access
```

Each stage's inputs/outputs are explicit Python objects (declared field sets in `PipelineStage`), not implicit shared state, so any stage can be re-run in isolation during development with a cached copy of its inputs (`PipelineRunner` supports "run from stage N").

---

## 7. Validation

Run as the `QualityControlStage`, producing a `ValidationReport` that is itself stored (as text/JSON) in `Provenance/validation_report`. Recommended checks, split by category:

**Schema-level:** every field in `CatalogueSchema` is present with the declared dtype and units attribute; no undeclared datasets exist (catches silent scope creep); units attribute is parseable and matches the declared unit family (e.g. a "mass" field must resolve to a mass dimension).

**Integrity:** `SatelliteID` values are unique and match row count across every subgroup (the row-order invariant from §2.1); no orphaned foreign keys (`HostHaloID`, `SharkGalaxyID`, `MergerTreeID` all resolve); no unexpected NaN/inf in fields that shouldn't have them, with an explicit allowed-missingness mask for fields that legitimately can be missing (e.g. `SersicIndex` when no structural fit converged) rather than silent NaNs; particle-tag CSR offsets are monotonically non-decreasing and terminate at `len(particle_ids)`.

**Physical plausibility:** `Mpeak >= M200c_z0` for every satellite; `RedshiftInfall` decreases monotonically along a satellite's own accretion history array; `OrbitalApocentre >= OrbitalPericentre`; `HeliocentricDistance > 0`; `CompletenessWeight >= 1`; stellar mass never exceeds halo mass by an unphysical margin (a soft check, flagged not failed); SFH bins sum to `StellarMass` within tolerance; no satellite's `SnapshotInfall` postdates its `SnapshotAtMpeak` by construction-inconsistent amounts.

**Cross-match verification:** every `SharkGalaxyID` referenced exists in the Shark output and vice versa (report unmatched Shark galaxies as a diagnostic, not necessarily an error — some may be centrals); spot-check a random subsample by independently recomputing the match from raw IDs; positional cross-checks where available (a matched satellite's halo position and Shark-galaxy assumed position should coincide within tolerance).

**Statistical sanity (catalogue-level, not per-object):** stellar-mass function and luminosity function of the output fall within a plausible envelope compared to abundance-matching expectations or prior Dorcha runs (flag large excursions for human review — this is a smoke test, not a science result); satellite counts per host are broadly consistent between suite realisations; the fraction of quenched satellites, isolation fraction, and other summary numbers are logged per version so a `dorcha-diff` CLI command can show what changed between catalogue versions.

**Reproducibility check:** re-running the pipeline on identical inputs with the same `random_seed` produces a byte-identical (or checksum-identical for floating point, within tolerance) catalogue — this is the test that actually protects the "computed once" promise over a multi-year project as dependencies (numpy, astropy versions) drift.

---

## 8. Extensibility

The schema accommodates the listed future needs as follows, without structural redesign:

- **New galaxy properties** — add fields to `Satellites/GalaxyProperties/` or `DerivedQuantities/`; MINOR version bump; existing readers unaffected because `CatalogueReader` exposes fields by name, not by positional column index.
- **Hydrodynamic simulations** — define a parallel catalogue (`dorcha_hydro_catalogue_v1.0.0.h5`) using the *same* `CatalogueSchema` where physically meaningful (identical field names/units/definitions for anything comparable, e.g. `Mpeak`, `RedshiftInfall`), with simulation-specific fields added under the same subgroup conventions. Comparability comes from schema reuse, not from cramming both simulations into one file.
- **Observational catalogues** — live under `CrossMatches/<survey_name>/`, row-aligned to `Satellites/` via an index array, keeping the boundary between "simulated" and "observed" data explicit rather than merging fields that aren't truly equivalent.
- **ML classifications** — namespaced, versioned subgroups under `MLClassifications/`, as noted in §2.2 — additive, never overwrite physical columns.
- **Future particle tagging** — a new tagging method becomes `ParticleTags/<method_name>_v2/` alongside the original rather than replacing it, since different tagging algorithms are legitimately different derived data products worth comparing, not drop-in replacements.

The invariant that makes all of this safe: **row order and `SatelliteID` are permanent**. Any future addition is either a new column aligned to the existing row order, or a new namespaced subgroup with its own explicit index into that row order. Nothing new is ever allowed to require renumbering or reordering existing rows.

---

## 9. Potential Pitfalls

**Silent unit mismatches.** The single most common failure mode in simulation catalogues: comoving vs. physical, $h^{-1}$ vs. not, kpc vs. Mpc. Mitigation: every dataset's units attribute is mandatory and machine-checked at write time; `CatalogueReader` returns unit-aware objects by default so arithmetic errors surface as `unyt`/`astropy.units` exceptions rather than silently wrong numbers.

**Row-order drift.** If any future pipeline run re-sorts or filters the satellite table differently, every downstream notebook that assumed positional alignment silently breaks. Mitigation: `SatelliteID` as the only safe join key; validation explicitly checks row-order consistency across subgroups; discourage (via documentation and code review) any analysis code that indexes by position rather than by ID or boolean mask.

**Derived-quantity definition drift.** "Infall" or "quenched" can be defined multiple reasonable ways; if the definition changes between catalogue versions without a version bump, old and new results silently disagree. Mitigation: every derived quantity's *algorithm version* is recorded (either as a field-level attribute or in `Provenance/`), and definitional changes force at least a MINOR version bump plus a changelog entry.

**Recomputation creep.** The temptation to "quickly recompute X in this one notebook because it's easier than checking if it's already there" undermines Principle 2 within months. Mitigation: `CatalogueReader`'s API surface should make using the stored quantity easier than recomputing it (good naming, discoverability via `Documentation/schema_reference.md`), and code review in paper-specific repos should flag any recomputation of a quantity that already exists in the catalogue.

**Plot-shaped storage creeping back in.** Someone eventually wants to cache a "binned luminosity function" for speed. Mitigation: cached *statistical summaries* for performance are acceptable as an explicitly labelled cache (e.g. under a clearly separate `Cache/` group, regenerable and excluded from the canonical schema, never read by anything except a plotting convenience function) — but the master schema itself must never gain a bespoke plot-shaped dataset. Draw this line early and document it, because the pressure to violate it will recur with every paper.

**Provenance rot.** Two years from now, "which Shark run produced this?" needs an answer. Mitigation: `Provenance/` stores enough (git commits, config file checksums, input file paths and checksums, timestamps) that any number in the catalogue can be traced to the exact code and inputs that produced it, without relying on anyone's memory or a wiki page that will go stale.

**HDF5 concurrency and partial writes.** HDF5 files can corrupt if a write is interrupted or accessed concurrently by multiple writers. Mitigation: only `CatalogueWriter` writes, it writes to a temporary file and atomically renames on success, and releases are treated as immutable (read-only permissions on disk) once published.

**Chunking/compression chosen carelessly.** Wrong chunk shapes make either write or (more commonly) read performance far worse than a flat binary format, undermining the "central analysis product" goal if every notebook load is slow. Mitigation: benchmark chunk sizes against realistic access patterns (whole-column reads for most science use, occasional whole-row reads) before finalising `SimulationConfig`/writer defaults; prefer chunk shapes aligned with columnar access since that's the dominant use case.

**Scope creep in "primary vs. derived."** It's tempting to store something "derived" as if it were primary once it's been recomputed a few times and feels canonical. Mitigation: the `is_derived` flag and provenance string are mandatory metadata, checked by `SchemaValidator`, so the distinction can't quietly disappear.

---

## 10. Best Practices for Long-Term Maintainability

- **Schema-as-code**: `CatalogueSchema` lives in version control as a YAML/JSON file (`configs/schema_v1.0.yaml`), and `Documentation/schema_reference.md` is auto-generated from it — the documentation cannot drift from the implementation because it *is* the implementation's description.
- **Immutable releases, mutable pipeline**: never edit a released catalogue file in place; a bugfix means a PATCH-version rebuild from raw data.
- **Golden regression tests**: keep a small, fast synthetic mock dataset (a handful of hosts/satellites with hand-computed expected derived quantities) that runs through the full pipeline in CI on every commit, catching silent regressions in derived-quantity logic.
- **Independent stage testability**: because `PipelineStage` subclasses declare explicit inputs/outputs, each can be unit-tested with synthetic inputs without running the full pipeline — critical once the pipeline has a dozen derived-quantity stages contributed by different people over years.
- **One PR per new derived quantity**, each required to include: the field's schema entry (name/units/description/provenance), the computing code, a unit test, and a one-line addition to `docs/changelog/`. This keeps the ten-future-papers growth path incremental and reviewable rather than a periodic large, risky rewrite.
- **Explicit ownership of ambiguous definitions**: anything with more than one reasonable physical definition (infall, quenching, isolation, pairing) gets a short docstring/`Documentation/` entry stating the chosen definition and citing the reference paper or convention it follows, so a new contributor doesn't have to reverse-engineer intent from code.
- **Keep analysis and catalogue-construction code in separate repositories/packages**, connected only via a pinned `dorcha_catalogue` version dependency — this is what actually enforces Principle 4 over time, since it removes the path of least resistance for "just quickly reopening the snapshot."
- **Treat the catalogue as a citable data product**: once mature, give it a DOI (e.g. via a Zenodo deposit per major release) so companion papers can cite it precisely, and so version pinning is externally verifiable by referees.
