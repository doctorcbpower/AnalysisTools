# Phase 6 (`analysistools.catalogue`) — remaining work

Status as of 2026-08-07 (updated same day: `HostEnvironmentStage`,
`SharkGalaxyID`, and `SchemaValidator`'s presence/dtype/unit checks landed).
This supersedes the stale "not started" language for 6b/6c in
[DEVELOPMENT.md](../DEVELOPMENT.md) and
[master_catalogue.md](master_catalogue.md) — those should be updated to point
here rather than re-describing status inline.

## Summary

| Phase | Status |
|---|---|
| 6a — schema + read-side | ✅ done |
| 6b — pipeline + derived quantities | 🟡 mostly done — see gaps below |
| 6c — validation + webreport + orchestration | 🟡 SchemaValidator done, rest not started |

All currently-implemented code is tested (unit tests with fakes/mocks for
logic, plus `@needs_data` integration tests against the bundled SWIFT/
VELOCIraptor fixtures for wiring). Nothing described as "done" below still
raises `NotImplementedError`.

---

## 6b — implemented

- **`pipeline.py`**: `HaloExtractStage`, `TreeExtractStage`, `CrossMatchStage`
  (SatelliteID assignment + sort, plus calling `galaxy_backend.galaxy_properties()`
  per satellite when both `galaxy_backend` and `epoch` are supplied).
- **`derived.py`**: `HaloPropertiesStage`, `OrbitalPropertiesStage`,
  `StarFormationHistoryStage`, `HostEnvironmentStage`, `EnvironmentStage`,
  `ObservabilityStage`, `DorchaSpecificStage` (partial — see below).
- **`backends.py`**: `SharkGalaxyBackend` and `HydroGalaxyBackend`, both with
  `galaxy_properties()` and `star_formation_history()`; shared helpers
  `rebin_sfh()` and `_cosmic_time_gyr()`.

`ExtractStage` itself (the base class) is deliberately left abstract
(`NotImplementedError`) — only its two concrete subclasses are meant to be
instantiated. This is intentional, not a gap.

## 6b — still open (base `ExtractStage`)

`analysistools/catalogue/pipeline.py:90` — the base class stub. No action
needed unless a new extraction mode is added that doesn't fit
`HaloExtractStage`/`TreeExtractStage`.

## 6b — deferred fields (per stage, with reason)

Each item below was deliberately **not** invented — implementing it requires
either an explicit modeling choice the caller must supply, or infrastructure
that doesn't exist yet. Full reasoning lives in the docstring of each stage
in `derived.py`/`backends.py`; this is a condensed index.

### `EnvironmentStage` / `HostEnvironmentStage`
- `Haloes/IsIsolated`, `IsPaired`, `PairedHostID`, `N_satellites_total` are
  **now implemented** (`HostEnvironmentStage`, a new stage) — isolation is
  "no more massive halo within `isolation_radius_factor * R200c_host`",
  pairing is "comparable mass (within `pairing_mass_ratio_min`) within
  `pairing_max_separation`", `N_satellites_total` counts neighbours above
  `completeness_mass_threshold` within the isolation radius. All four
  thresholds are required constructor arguments, not defaulted — same
  reasoning as `EnvironmentStage`'s `mass_threshold`/`aperture_radius`.
  `EnvironmentStage` copies the two boolean flags through onto every
  satellite row as `HostIsIsolated`/`HostIsPaired` when this stage has
  already run.
- Cosmic-web classification (filament/void/node membership) — no cosmic-web
  classifier exists in the codebase. Still open.

### `ObservabilityStage`
- Anything requiring a survey selection function / instrument model beyond
  the basic `observer_pos`-relative geometry already computed.

### `HaloPropertiesStage` / `OrbitalPropertiesStage`
- `MassAccretionRateDM` — needs a common time/snapshot grid across the tree
  to differentiate mass vs. time; no shared-grid machinery exists yet
  (each `Epoch` loads independently).
- Full progenitor-list-dependent quantities (e.g. anything needing
  `NextProgenitor`/full merger history, not just the main branch) — the
  tree readers (`treeio_subfind.py`, `treeio_treefrog.py`) only expose the
  main-branch walk today; parsing the full progenitor pointer graph is
  unimplemented in the underlying tree I/O, not just in `derived.py`.

### `SharkGalaxyBackend` / `HydroGalaxyBackend`
- `SharkGalaxyID` **is now populated** by `SharkGalaxyBackend` (native
  `id_galaxy`, omitted if absent, same "omit don't fabricate" contract as
  every other field). `HydroGalaxyBackend` still has no equivalent
  concept — a hydro run has no SAM-side galaxy ID at all, so this field is
  legitimately always absent for `galaxy_backend: hydro` projects.
- Any field requiring particle tagging (e.g. tracking a *specific* set of
  star particles across snapshots rather than an SPH/FOF re-aggregation
  each epoch) — no particle-tagging infrastructure exists.
- Unit reconciliation against `schema.py`'s declared units — both backends
  still return raw values in whatever units the underlying `Epoch`/
  `Dataset` natively provides (no numeric conversion is applied anywhere).
  **Partially closed**: `validation.SchemaValidator.check()` now cross-
  checks the little-h/comoving *labelling* (declared schema units string
  vs. the actual source's `little_h`/`comoving` state, recorded in
  `context.meta["comoving_little_h"]` by `HaloExtractStage`/
  `CrossMatchStage` — see both backends' `native_comoving_little_h()`).
  This catches exactly the silent mismatch the docstrings warned about
  (e.g. SHARK's native `Msun/h` written under a schema field declared
  plain `"Msun"`) but does **not** convert values — a project whose schema
  declares plain `Msun` and whose backend returns `Msun/h` will get a
  `SchemaValidator` warning, not a corrected value. Actual unit
  *conversion* before write is still unimplemented (would belong to
  `WriteStage`, 6c).

### `DorchaSpecificStage`
- Only `EarliestProgenitorRedshift` and `FossilFraction` are implemented.
  Every other `DorchaProperties` field in the schema is deferred — each
  requires either particle tagging, the cosmic-web classifier, or a
  Dorcha-project-specific convention not yet pinned down (see the class
  docstring in `derived.py` for the field-by-field list).

### Cross-cutting
- **`SatelliteID` stability**: assigned per-run by `CrossMatchStage`'s sort,
  not guaranteed stable across independent pipeline rebuilds (e.g. if the
  satellite list or sort order changes between runs). No persistent ID
  registry exists.
- **Reader/writer attribute-naming mismatch**: `api/snapshot.py`'s
  `_RAW_ATTR_ALIASES` (added this session) maps `gas_Z`→`gas_metallicity`,
  `stellar_Z`→`stellar_metallicity`, `initmass`→`stellarinitmass` on
  *read*, but nothing was fixed on the *write* side — if any writer path
  ever needs to round-trip these fields under their canonical names, that
  asymmetry will resurface. Not currently blocking anything in the
  catalogue pipeline (which only reads), but worth flagging before someone
  builds a snapshot writer against these names.

---

## 6c — not started (except `SchemaValidator`)

`SchemaValidator.check()` is implemented (pulled forward from 6c — see the
"unit reconciliation" note above). Everything else below still raises
`NotImplementedError`; only data containers (`ValidationIssue`,
`ValidationReport`, `FigureSpec`) besides `SchemaValidator` are real.

### `validation.py`
- ~~`SchemaValidator.check()`~~ — **done.** Presence (via
  `CatalogueSchema.validate_columns()`, missing→warning/undeclared→error),
  dtype-family checks (warning), and the little-h/comoving-vs-declared-
  units cross-check (warning) described above. Not wired into
  `QualityControlStage` yet (that stage itself is still unimplemented).
- `IntegrityValidator.check()` (line 97) — referential integrity (e.g.
  every `SatelliteID` referenced elsewhere actually exists, no duplicate
  IDs, tree links resolve).
- `PhysicalValidator.check()` (line 109) — sanity bounds on physical
  quantities (e.g. non-negative masses, `0 <= FossilFraction <= 1`,
  monotonic SFH bin edges).
- `ValidationReport.to_json()` (line 60) — serialization, no design
  blockers, just unwritten.

### `pipeline.py`
- `QualityControlStage.run()` (line 377) — should instantiate and run the
  three validators above and attach the resulting `ValidationReport` to
  the context/output.
- `WriteStage.run()` (line 393) — write the finished `PipelineContext`
  columns out as the versioned, self-describing HDF5 catalogue file
  (per the design doc's schema/group/attrs layout). No blockers, just
  unwritten — this is the actual file-output half of "master catalogue."
- `CatalogueBuilder._load_config()` (line 414) — YAML project config
  loading (selecting backend, stages, per-stage kwargs like
  `mass_threshold`/`quenched_ssfr_threshold`/`reionisation_lookback_time`).
  An example config already exists at `configs/dorcha.yaml` but nothing
  reads it yet.
- `CatalogueBuilder.run()` (line 425) — top-level orchestration: load
  config, construct stages in order, run the full pipeline
  Extract→CrossMatch→derived stages→QualityControl→Write.
- `CatalogueBuilder.run_from_stage()` (line 432) — resume/re-run from an
  intermediate stage (for iterating on a single derived-quantity stage
  without re-running extraction).

### `webreport.py`
- `WebReport.build()` (line 76) — assemble the added `FigureSpec`s into a
  single static self-contained HTML file (plotly, no server). No design
  blockers; `add_figure()`/`from_catalogue()` already work.

---

## Suggested order if picking this up

1. `WriteStage` + `SchemaValidator` — these two together let you produce
   and check an actual output file end-to-end for the first time, which
   is also the fastest way to surface any remaining unit/shape mismatches
   from the 6b backends.
2. `CatalogueBuilder._load_config()` + `run()` — turns the pipeline from
   "stages wired by hand in a test" into something runnable from a YAML
   config, matching `configs/dorcha.yaml`.
3. `IntegrityValidator` / `PhysicalValidator` / `ValidationReport.to_json()`.
4. `WebReport.build()`.
5. `QualityControlStage.run()` (thin glue once the validators exist).
6. `CatalogueBuilder.run_from_stage()` (convenience, lowest priority).

None of the deferred-field items above block this ordering — they can be
picked up independently whenever the relevant modeling decision or
infrastructure (particle tagging, cosmic-web classifier, host/satellite
pairing algorithm, full progenitor-list tree parsing, common-grid mass
accretion) is ready.
