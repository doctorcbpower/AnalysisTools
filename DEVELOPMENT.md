# AnalysisTools — Development and Implementation Plan

**Status:** proposal for review — nothing below is implemented yet unless marked otherwise.
**Guiding decision:** a new unified layer is added *on top of* the existing classes; the current API keeps working unchanged. The `shark` package is integrated into `analysistools` as a first-class data source.

---

## 1. Motivation

AnalysisTools has grown organically around four kinds of data — snapshots, halo catalogues, merger trees, and SHARK galaxy catalogues — plus two toolkits for turning them into science (`GriddingTools`, `ProfileTools`). Each data source currently has its own idiom:

| Source | Class | Construction | Load | Field access |
|---|---|---|---|---|
| Snapshot | `SnapshotTools` | `SnapshotTools(snapfileformat="HDF5")` | `read_snapshot(path)` | attributes: `snap.pos`, `snap.dm.pos` after `LoadParticlesByType` |
| Halo catalogue | `HaloTools` | `HaloTools(comoving_units=True)` | `read_catalogue(path, fileformat=...)` | returned dict / `standardise_names()` |
| Halo catalogue (multi-z) | `HaloModel` | paths given to constructor | lazy, on first access | `get(field, redshift)` + accessors (`mvir`, `rvir`, …) |
| Merger tree | `MergerTreeTools` | `MergerTreeTools(path, treefileformat=...)` | in constructor (`read_tree`) | `HaloTrack` objects via `get_track` / `from_halo` |
| SHARK | `SharkModel` (in separate `shark` package) | `SharkModel(model_dir, redshift_table, subvols)` | lazy | `get(field, redshift)` + accessors (`mstars`, `sfr`, …) |
| SHARK (legacy) | `GalaxyTools` | `GalaxyTools(path, format)` | `ReadGalaxyCatalogue()` | attributes |

Consequences: every source needs different boilerplate; combining sources (e.g. gridding snapshot particles then over-plotting halo positions and SHARK galaxies) means manually extracting arrays with three different syntaxes; naming conventions mix `snake_case` and `CamelCase` methods; and `shark` lives outside the package with its own plotting stack.

The good news: two convergent patterns already exist in the codebase and the design below simply generalises them.

1. **The standardised field schema** (`halo_tools_standardise_names.py`): `mass`, `pos`, `vel`, `radius`, `halo_id`, `num_part`.
2. **The lazy epoch-model pattern** (`HaloModel` / `SharkModel`): `get(field, redshift)`, `get_meta`, per-quantity accessors, caching, `label`/`colour` for plotting. `HaloModel` was explicitly written to "mirror the SharkModel API" — this proposal formalises that mirror as an interface.

---

## 2. Design goals

1. **One mental model.** Open any data source the same way, ask for fields the same way, regardless of on-disk format.
2. **One field vocabulary.** A common schema so `ds["mass"]` and `ds["pos"]` mean the same thing for halos, galaxies, tree nodes, and (per particle species) snapshots.
3. **Composition without plumbing.** Plot functions accept dataset objects directly; a `Simulation` object ties matched snapshot + catalogue + tree + SHARK outputs together and handles cross-matching.
4. **No breakage.** Existing classes, method names, and the notebook keep working. The new layer wraps; it does not rewrite.
5. **HPC-friendly.** Lazy loading and caching everywhere (already the `HaloModel`/`SharkModel` behaviour); nothing reads from disk until asked.

---

## 3. Architecture

New subpackage `analysistools/api/` (working name), plus `shark` relocated to `analysistools/shark/`:

```
analysistools/
├── api/
│   ├── __init__.py        # load(), Simulation, re-exports
│   ├── dataset.py         # Dataset base class + field schema registry
│   ├── snapshot.py        # SnapshotDataset  (wraps SnapshotTools/SnapshotData)
│   ├── halos.py           # HaloCatalogue    (wraps HaloTools)
│   ├── trees.py           # MergerTree       (wraps MergerTreeTools)
│   ├── galaxies.py        # GalaxyCatalogue  (wraps SharkModel)
│   ├── simulation.py      # Simulation umbrella + cross-matching
│   └── plotting.py        # dataset-aware wrappers over GriddingTools/ProfileTools
├── shark/                 # moved from top-level shark/ (thin top-level shim retained)
└── (existing modules unchanged)
```

### 3.1 The `Dataset` interface

Every adapter implements the same small contract:

```python
class Dataset:
    # identity / provenance
    kind: str                  # "snapshot" | "halos" | "tree" | "galaxies"
    path: str
    fileformat: str
    label: str                 # for plot legends (as HaloModel/SharkModel today)

    # metadata (uniform keys: redshift, boxsize, h0, volume, units, ...)
    meta: dict

    # field access — standardised names, native names accepted as fallback
    def __getitem__(self, field: str) -> np.ndarray: ...
    def get(self, field: str, default=None): ...
    @property
    def fields(self) -> list[str]: ...       # standardised + native
    def __contains__(self, field) -> bool: ...
    def __len__(self) -> int: ...            # number of objects/particles

    # selection — returns a lightweight view, same interface, mask applied
    def select(self, mask=None, *, centre=None, size=None,
               geometry="spherical", **cuts) -> "Dataset": ...

    def summary(self) -> None: ...           # unified replacement for the
                                             # various summary() methods
```

Key points:

- `ds["pos"]` always returns `(N, 3)`; `ds["mass"]`, `ds["vel"]`, etc. follow the standardised schema. Native names (`"Group_M_Crit200"`, `"mstars_disk"`, `"SubhaloVel"`) remain accessible as a fallback, so nothing is hidden.
- `select()` subsumes `snapshot_tools.select_particles` (geometric cuts, periodic wrapping) and adds attribute cuts (`ds.select(mass=(1e11, 1e13))`). It returns a *view* carrying an index array — cheap, chainable, and every view is itself a `Dataset`, so anything that accepts a dataset accepts a selection.
- Uniform `meta` keys end the current situation where redshift lives in `snap.redshift`, `get_meta(z)["..."]`, or tree metadata depending on source.

### 3.2 The field schema

Extend `STANDARD_KEYS` from `halo_tools_standardise_names.py` into a single registry (module `api/dataset.py`), organised by kind:

| Field | snapshot (per species) | halos | tree node | galaxies |
|---|---|---|---|---|
| `pos`, `vel` | ✓ | ✓ | ✓ | ✓ (stacked from `position_x/y/z`) |
| `mass` | ✓ | ✓ (M200/Mvir) | ✓ | `mstars_disk+mstars_bulge` by default; native fields for components |
| `id` | `pids` | `halo_id` | `halo_id` | `id_galaxy` |
| `radius` | – | ✓ | – | `rstar_disk` (half-mass) |
| `num_part` | – | ✓ | – | – |
| `sfr` | ✓ (gas) | – | – | `sfr_disk+sfr_bulge` |
| `type` | species name | Group/Subhalo | IsSubhalo | 0/1/2 central/satellite/orphan |

Aliases are one-line dictionary entries, so extending the schema (e.g. `vmax`, `concentration`, `mhalo` for galaxies via `mvir_hosthalo`) is trivial and format-mapping stays where it belongs: in the per-format tables, exactly as `CATALOGUE_MAPPINGS` does now.

### 3.3 The `load()` factory

Single entry point with format sniffing:

```python
import analysistools as at

snap  = at.load("snap_0031.hdf5")                        # kind + format sniffed
halos = at.load("fof_output_0031.hdf5", kind="halos", fileformat="SWIFT_FOF")
tree  = at.load("VELOCIraptor.walkabletree.hdf5", kind="tree")
gals  = at.load("./CDM/base_model", kind="galaxies", redshift_table="redshift_list.txt")
```

Sniffing rules (HDF5 group names are distinctive: `PartType0` ⇒ snapshot, `Subhalo`/`Groups` ⇒ SubFind, `galaxies` ⇒ SHARK, …) with explicit `kind=`/`fileformat=` always available as an override. Internally each adapter just instantiates the existing class (`SnapshotTools`, `HaloTools`, …) and holds it as `._backend`, which remains reachable for anything the unified API doesn't cover.

### 3.4 Multi-epoch access: `EpochModel`

`HaloModel` and `SharkModel` already share `get(field, redshift)`, `get_meta`, `preload`, `clear_cache`, `label`, `colour`. Formalise this as the `EpochModel` interface (an abstract base or Protocol) that both inherit/register with, and give it the same field-schema resolution as `Dataset`, plus:

```python
ds_z0 = model.at(0.0)      # -> Dataset view frozen at one redshift
```

`at(z)` is the bridge between the multi-epoch world and the single-epoch `Dataset` world: everything downstream (selection, plotting, cross-matching) only ever needs to understand `Dataset`.

### 3.5 The `Simulation` umbrella

Binds matched products of one simulation and does the cross-matching that currently has to be done by hand:

```python
sim = at.Simulation(
    snapshots  = "output/snap_{snapnum:04d}.hdf5",      # template or dict {z: path}
    halos      = {0.0: "catalogues/snap_200/fof.hdf5", 1.0: "..."},
    trees      = "trees/VELOCIraptor.walkabletree.hdf5",
    galaxies   = "shark/CDM/base_model",
    fileformats = dict(snapshots="HDF5", halos="SWIFT_FOF", trees="TreeFrog"),
    label      = "CDM",
)

epoch = sim.at(z=0.0)                 # Epoch: matched Dataset views
epoch.snapshot["pos"]
epoch.halos["mass"]

# cross-matching helpers
parts = epoch.particles_in_halo(halo_id)          # via groupid or radius query
gals  = epoch.galaxies_in_halo(halo_id)           # via id_halo / position match
track = sim.track_of(epoch.halos, index=42)       # wraps MergerTreeTools.from_halo
```

Every component is optional — a `Simulation` with only snapshots is fine. Matching strategies (particle `groupid` where present, else periodic radius query; SHARK `id_halo` where IDs align, else position match) are pluggable, defaulting to the cheapest available.

### 3.6 Dataset-aware plotting

`GriddingTools` and `ProfileTools` remain the numerical engines, untouched. A thin layer (`api/plotting.py`) accepts datasets and handles extraction, units, labels, and composition:

```python
from analysistools.api import plotting as plt2

# gridded surface density of DM, halos over-plotted, galaxies over-plotted
ax = plt2.density_map(snap.select(centre=c, size=5.0), species="dm",
                      grid_size=(512, 512), method="CIC", projection="xy")
plt2.overlay_points(ax, halos.select(mass=(1e12, None)), marker="o")
plt2.overlay_points(ax, gals, colour_by="sfr")

# radial profile straight from a dataset
prof = plt2.profile(snap, species="dm", centre=c, kind="density",
                    rmin=0.1, rmax=300)          # -> ProfileTools result dict
plt2.plot_profile(prof, label=snap.label)

# same call, different source — mass function from halos or galaxies alike
plt2.mass_function([halos_cdm, halos_wdm], volume="auto")
```

Because every argument is a `Dataset` (or a selection view of one), "combine data from snapshots, halo catalogues, and shark catalogues" reduces to passing several datasets to the same axes. The `label`/`colour` conventions already in `HaloModel`/`SharkModel`/`shark.plots.Plotter` carry over so multi-model comparison plots keep working.

### 3.7 SHARK integration

- Move `shark/` → `analysistools/shark/` (history-preserving `git mv`). Keep a top-level `shark/__init__.py` shim that re-exports from `analysistools.shark` with a `DeprecationWarning`, so `import shark` keeps working.
- `GalaxyCatalogue` (the `Dataset` adapter) wraps `SharkModel`; `galaxy_tools.GalaxyTools` is superseded and marked deprecated (it double-reads `velocity_x` — a bug — and duplicates `SharkModel` functionality).
- `shark.analysis.Analysis` / `shark.plots.Plotter` stay for existing workflows; over time their generic pieces (number densities, mass functions) migrate into `api/plotting.py` where they can serve halos too.
- The photometry subpackage (fsps-based) is untouched by this proposal.

---

## 4. Conventions

Adopted for all *new* code (existing code untouched until deprecation decisions are made):

- **Naming:** `snake_case` methods and arguments everywhere. Standard verbs: `load` (register/open), `read` (pull to memory — usually implicit/lazy), `get` (field access), `select` (subset), `summary` (human-readable report).
- **Common argument names:** `fileformat` (not `snapfileformat`/`treefileformat`/`galfileformat`), `path`, `redshift`, `centre`, `comoving`.
- **Units:** each adapter records unit metadata in `meta["units"]`; conversions are explicit (`ds.to_physical()`, `ds.to_comoving()`) rather than constructor flags acting at read time. Constructor flags remain honoured on the backends.
- **Logging:** `logging` module throughout, default `WARNING` in library code (the `HaloModel` choice), never bare `print`.
- **Laziness:** no disk I/O at construction; first field access reads and caches. `preload`/`clear_cache` on every adapter.
- **Typing and docstrings:** full type hints, NumPy-style docstrings (the newer modules — `halo_model.py`, `merger_tree_types.py`, `shark/model.py` — are the template).

---

## 5. Implementation roadmap

Phased so each stage lands independently and the package is never broken in between.

**Phase 1 — schema + Dataset core** *(foundation)* — **✅ landed 2026-07-26**
Field-schema registry; `Dataset` base with `__getitem__`/`fields`/`meta`/`select`/`summary`; `SnapshotDataset` and `HaloCatalogue` adapters; `at.load()` with sniffing for these two kinds. Unit tests against `data/snap_0031.hdf5` and `data/halos/` (`pytest tests/` — 13 tests).
Bugs fixed en route: `HaloTools.read_catalogue` rejected the documented string `"VELOCIraptor"` (case mismatch with `FORMAT_READERS`); `standardise_catalogue_names` had the same case bug; `read_velociraptor` only handled the grouped `MetaData/Groups` layout, not the standard flat `.properties` layout; pos/vel stacking was AHF-only (now generic, covers VR `Xc/VXc`); `GriddingTools` jitted kernels crashed on tuple `grid_size` and had no 2-D support (kernels split into `_ngp/_cic_assign_2d/_3d` with periodic wrap and a Python-level dispatch).

**Phase 2 — trees + galaxies + epoch models** — **✅ landed 2026-07-26**
`MergerTree` and `GalaxyCatalogue` adapters; `EpochModel` interface retro-fitted to `HaloModel`/`SharkModel`; `model.at(z)`. `shark` package moves into `analysistools` with the compatibility shim.
Notes from implementation: `treeio_treefrog` gained support for the TreeFrog *walkable tree* flavour (topology-only `Snapshots/Snap_NNN` groups with Head/Tail pointers; the previous reader only handled the full-tree flavour with baked-in properties, and could not read `data/VELOCIraptor.walkabletree.hdf5`). Walkable-tree tracks fill Mass/Pos/Vel from linked halo catalogues where available (NaN otherwise), with automatic temporal-ID ↔ catalogue-row conversion. Tree queries return `TrackDataset` (rows = snapshots), which follows the full Dataset contract. `GalaxyCatalogue` is file-backed (one `galaxies.hdf5`) or model-backed (`SharkModel.at(z)`, lazy per-field fetch through the model cache via a new `Dataset._fetch_native` hook). `GalaxyTools` now emits `DeprecationWarning`. The top-level `shark` shim aliases every submodule in `sys.modules`, so `from shark.model import SharkModel` returns the *same* module objects as `analysistools.shark.model` (class identity preserved).

**Phase 3 — Simulation + cross-matching** — **✅ landed 2026-07-26**
`Simulation`, `Epoch`, `particles_in_halo`, `galaxies_in_halo`, `track_of`. Position-matching uses a KD-tree (scipy) with periodic wrapping via `merger_tree_types.periodic_delta`.
Notes from implementation: `Simulation` component specs accept a path, a `{redshift: path}` dict, a prebuilt Dataset, or a multi-epoch model (anything with `.at(z)`); `sim.at(z)` returns a cached `Epoch` whose components load lazily at the nearest available redshift. `particles_in_halo` supports `match_by="groupid"|"radius"|"auto"` (auto prefers stored group membership) plus `r_scale=` and `species=`; `galaxies_in_halo` defaults to periodic-KD-tree position matching with `match_by="id"` optional (decision §7.4). `Epoch.track_of` auto-links the epoch's halo catalogue into the tree (`MergerTree.link_halos`, additive per epoch) so walkable-tree tracks fill properties at every epoch a catalogue is available. Periodic matching lives in `api.simulation.match_positions` (scipy `cKDTree(boxsize=L)`).

**Phase 4 — plotting layer** — **✅ landed 2026-07-26**
`api/plotting.py`: `density_map`, `overlay_points`, `profile`/`plot_profile`, `mass_function`. Port the notebook (`notebooks/PlotSimulationAnalysis.ipynb`) to the new API as the living example.
Notes from implementation: `density_map` selects a (periodic, wrapped) region when given `centre=`/`size=`, converts assigned mass to surface density, and leaves the grid reachable on `ax._at_grid`; `overlay_points` takes the same `centre=` so overlays land in the same wrapped frame. `profile` dispatches on `kind=` ("density", "surface_density", "velocity_dispersion", "vertical_velocity_dispersion") and annotates the ProfileTools result dict so `plot_profile` needs no further context. `mass_function` accepts any mix of halo/galaxy Datasets, normalises by `meta["volume"]` (else boxsize³) with `volume="auto"`, and supports `cumulative=`. The old notebook predates the unified interface and used since-renamed APIs, so the port is a fresh executed `notebooks/UnifiedInterfaceDemo.ipynb`, generated by `analysistools/api/_build_demo_notebook.py` (regenerate after API changes; requires nbformat/nbclient). Projected-mass conservation of `density_map` is unit-tested at 1e-4.

**Phase 5 — consolidation** — **✅ landed 2026-07-26**
Deprecation warnings on superseded paths (`GalaxyTools`, top-level `shark` import); README examples switched to unified API (old examples retained in an appendix); repository housekeeping (below).
Notes from implementation: deprecations landed in Phase 2; README legacy per-class usage now sits under "Appendix: Legacy Per-Class API" with the unified interface primary; `halo_tools.orig.py`, `merger_tree_subfind_tools.py` (empty stub), and the commented-out block in `halo_tools.py` removed (no references anywhere); packaging modernised to `pyproject.toml` (setup.py removed) — adds the previously missing `numba` dependency, optional extras `[test]`, `[notebooks]`, `[photometry]`, and pytest configuration. Version bumped to 0.2.0.

Each phase = one PR-sized change set with tests. Estimated sizes: P1 and P4 are the substantial ones; P2/P3/P5 are mostly glue.

---

## 6. Repository housekeeping

Issues noticed during the audit, independent of the new layer — all resolved as of Phase 5:

- ✅ `analysistools/__init__.py`: `__all__` listed `HaloCatalogue` (never imported — the class was `HaloTools`) and `analysis_tools` (no such module). *(Fixed in Phase 0/docs commit.)*
- ✅ `halo_tools.orig.py` — superseded copy; large commented-out block at the bottom of `halo_tools.py`. *(Removed in Phase 5.)*
- ✅ `merger_tree_subfind_tools.py` — empty stub (imports only), no references. *(Deleted in Phase 5.)*
- ✅ `galaxy_tools.py` — duplicated `self.vel` read (bug); superseded by `SharkModel` (see §3.7). *(DeprecationWarning added in Phase 2.)*
- `snapshot_tools.py` mixes `CamelCase` (`LoadParticlesByType`, `UnitConversion`, `ParticleOffsetsByType`) and `snake_case` methods — new layer exposes `snake_case` only. *(By design; legacy names untouched.)*
- ✅ Test suite: `tests/` (pytest, 41 tests) against the small files in `data/`. *(Phases 1–4.)*
- ✅ `pyproject.toml` replaces the minimal `setup.py`: adds the missing `numba` dependency, extras `[test]`/`[notebooks]`/`[photometry]`, pytest config. *(Phase 5.)*

---

## 7. Design decisions (resolved 2026-07-26)

1. **Subpackage naming** — both: implementation lives in `analysistools/api/`; the user-facing names (`load`, `Simulation`, adapter classes) are re-exported at the `analysistools` top level.
2. **Snapshot species handling** — `ds["pos"]` returns all particles; `ds.dm`, `ds.gas`, `ds.star`, `ds.bh` return per-species `Dataset` views (preserving `snap.dm.pos` muscle memory); `select(species=...)` is the programmatic equivalent.
3. **Units** — plain strings + explicit conversion factors in `meta["units"]`; no astropy dependency.
4. **SHARK↔halo cross-matching** — position matching (periodic KD-tree) is the default, since `id_halo` consistency with the N-body catalogue is not guaranteed; ID matching available via `match_by="id"`.
5. **`FDMFieldTools`** — to be integrated as a fifth kind (`"field"`) with grid semantics in future work; the `Dataset` contract keeps `select`/`__len__` optional-per-kind so mesh data fits without contortions.
6. **`GalaxyTools`** — confirmed deprecated; superseded by `SharkModel` / the `GalaxyCatalogue` adapter.

---

*Document version: 2026-07-26. Update this file as phases land; it is the reference for what "done" means.*
