# Master Science Catalogues (`analysistools.catalogue`)

**Status:** Phase 6a done (schema + read-side); Phase 6b (pipeline/derived quantities) and 6c (validation/webreport) not yet implemented. See DEVELOPMENT.md for the phase breakdown.

`analysistools.catalogue` turns the read-side unified interface (`Dataset`/`Simulation`/`Epoch`, see [unified_interface.md](unified_interface.md)) into a **write-side** pipeline that builds one flat, versioned, self-describing HDF5 catalogue per project — the single input every downstream science notebook reads from. The design originates from the Dorcha project (Milky-Way-analogue satellites + SHARK galaxies) but nothing in it is Dorcha-specific: it is a *template*, driven entirely by a project config, so a future hydrodynamic suite (Arepo-Solas, SWIFT) reuses the same package with a different config and a different galaxy backend rather than a different codebase.

## Why this sits on top of the existing `api/` layer, not the raw readers

The construction pipeline's Extract and Cross-match stages are already 80% built: `Simulation`/`Epoch` already resolve snapshots, halo catalogues, trees, and galaxies against each other (`epoch.particles_in_halo`, `epoch.galaxies_in_halo`, `sim.track_of`), already handle SubFind/AHF/VELOCIraptor/SWIFT-FOF and TreeFrog/SubFind-HBT/AHF-MergerTree transparently, and already sniff format from `at.load()`. `analysistools.catalogue` does not re-implement any of this — a `CatalogueBuilder` is handed a `Simulation` (or a dict of them, one per host) and pulls fields through `Epoch`, exactly as a notebook would. This is also what keeps the catalogue schema simulation-agnostic: anything that can be expressed as a `Simulation` — DMO + SHARK today, Arepo-Solas or SWIFT hydro tomorrow — can be catalogued with the same code.

## The one piece that genuinely differs between projects: the galaxy backend

Dorcha gets galaxy properties from SHARK (a semi-analytic model bolted onto a DMO tree). A hydro rerun of the same halo has stars *in the snapshot* — there is no SAM step. `analysistools/catalogue/backends.py` makes this a pluggable interface rather than a schema fork:

```python
class GalaxyBackend(Protocol):
    def galaxy_properties(self, epoch: Epoch, halo_row: int) -> dict[str, float]:
        """Return one row of Satellites/GalaxyProperties fields for one halo."""
```

* `SharkGalaxyBackend` — thin wrapper around `epoch.galaxies_in_halo(...)`; reads `GalaxyCatalogue` fields directly (mass, sfr, metallicity, ...).
* `HydroGalaxyBackend` — wraps `epoch.particles_in_halo(species="star", ...)` and aggregates star-particle properties (summed mass, mass-weighted metallicity/age, SFR from young-star counts) into the *same* field names. A future Arepo-Solas/SWIFT catalogue is schema-identical to the Dorcha one wherever the physics is comparable (`Mpeak`, `RedshiftInfall`, `StellarMass`, `MetallicityStellar`, ...), which is what actually makes cross-simulation comparison papers possible later without a translation layer.

Everything else in the schema — halo properties, orbital properties, environment, observability — is already simulation-agnostic because it's computed from the halo catalogue and merger tree, which `HaloCatalogue`/`MergerTree` already read uniformly across finders/trackers.

## Package layout

```
analysistools/catalogue/
├── __init__.py        # CatalogueBuilder, run_pipeline, re-exports
├── schema.py           # FieldSpec / CatalogueSchema — the Documentation/schema.json source of truth
├── backends.py         # GalaxyBackend protocol + SharkGalaxyBackend, HydroGalaxyBackend
├── derived.py          # DerivedQuantityStage subclasses (orbital, SFH, environment, observability, Dorcha-specific)
├── pipeline.py         # PipelineStage ABC, CatalogueBuilder orchestrator, extract/crossmatch/derive/write stages
├── validation.py       # Validator ABC + SchemaValidator/IntegrityValidator/PhysicalValidator
└── webreport.py         # static, self-contained HTML report generator (see below)
```

Read-side: the master catalogue is registered as a fifth `Dataset` kind (`kind="satellites"`) in `api/dataset.py`'s `FIELD_ALIASES`, with a `CatalogueDataset` adapter in `api/catalogue.py` that flattens the file's subgroups (`Satellites/HaloProperties/...`, `Satellites/GalaxyProperties/...`, ...) into one column namespace on load — so `at.load("dorcha_catalogue_v1.0.0.h5")` returns an ordinary `Dataset`, and every existing `api/plotting.py` function (`mass_function`, `overlay_points`, `profile`) works on it unchanged. This closes the loop: the same abstraction (`Dataset`) is both what the pipeline writes from and what science code reads from.

## The template mechanism: project config, not project code

The schema itself is versioned and bundled as package data
(`analysistools/catalogue/configs/schema_v{X.Y}.yaml`; load with
`analysistools.catalogue.default_schema("1.0")`) so a schema bump ships with
the package. A *project* is a separate, repo-root YAML file
(`configs/dorcha.yaml`, `configs/futurehydro.yaml`, ...) that references a
schema version by number, not a path, and is not a subclass of anything:

```yaml
project: dorcha
schema_version: "1.0"
galaxy_backend: shark          # shark | hydro
derived_stages:
  - halo_properties
  - orbital_properties
  - star_formation_history
  - environment
  - observability
  - dorcha_specific            # progenitor fraction, fossil fraction, particle tags
units:
  length: cMpc/h
  mass: Msun/h
```

`CatalogueBuilder(config_path).run(simulations)` loads the config, resolves `galaxy_backend` to a class, runs the listed `derived_stages` in order, validates, and writes. A hydro project's config drops `dorcha_specific`, sets `galaxy_backend: hydro`, and otherwise reuses every stage — no forking of `pipeline.py`.

## Units convention (unchanged from `api/`)

Per DEVELOPMENT.md §7.3, the project deliberately carries no astropy/unyt dependency: units are plain strings plus explicit conversion factors, recorded in `meta["units"]` on every `Dataset` and mirrored as HDF5 attributes on every catalogue field (`units`, `description`, `provenance`, `is_derived`). This is the same convention used throughout `api/`, just carried through to the catalogue's on-disk schema (see the original Dorcha design note, `dorcha_master_catalogue_design.md`, for the full per-field schema table — reproduce that table into `schema.py`'s `CatalogueSchema` when implementing).

## Shareable results (`webreport.py`)

For "let a collaborator or a referee look at figures without giving them the catalogue file," a live dashboard is the wrong tool — it needs a server and a reason to keep running. `webreport.py` instead builds a **static, self-contained HTML site**: `plotly` is already a core dependency (`pyproject.toml`), and `fig.to_html(full_html=False, include_plotlyjs="cdn")` produces embeddable, interactive (pan/zoom/hover) figures with no server. `WebReport(catalogue).add_figure(...).build(out_dir)` writes a small multi-page site (one `index.html` plus one page per figure group) that:

* is a plain folder of HTML/JS, deployable to GitHub Pages, attached to a Zenodo deposit as supplementary material, or emailed as a zip;
* embeds the catalogue version and provenance string used to generate it in the page footer, so a reader always knows which release the figures came from;
* optionally ships a small precomputed JSON alongside each figure so simple client-side filtering (by host, by mass cut) works via `plotly.js` restyle calls without any backend.

This deliberately does not call back into the live catalogue at view time — it is a snapshot, matching the "shareable results" use case rather than the "always-current dashboard" one. If a live, re-queryable dashboard is wanted later, `mcp__cowork__create_artifact`-style tooling or a small Streamlit/Dash app can read the same `CatalogueDataset` adapter; nothing about this design forecloses that, it's simply not what a supplementary-material link needs.

## Relationship to the original design note

[`dorcha_master_catalogue_design.md`](dorcha_master_catalogue_design.md) remains the authoritative per-field schema (HDF5 groups, dtypes, units, provenance, derived-quantity list, validation checks). This document is the mapping of that design onto `AnalysisTools`' actual conventions and existing code; implement `schema.py` by transcribing that document's tables, not by redesigning them.
