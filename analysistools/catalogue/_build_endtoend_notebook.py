#!/usr/bin/env python3
"""
Build notebooks/DorchaCatalogueEndToEnd.ipynb -- the end-to-end demo of
the Phase 6 master-catalogue pipeline (snapshot + halo catalogue + merger
tree + SHARK galaxy catalogue -> a built .h5 catalogue -> read back
through the unified interface -> a WebReport).

Unlike analysistools/api/_build_demo_notebook.py, this notebook needs
real Dorcha data (a SHARK catalogue in particular -- nothing bundled in
data/ has one) and is never auto-executed by this script; regenerate the
*structure* after API changes with:

    python -m analysistools.catalogue._build_endtoend_notebook

then fill in the CONFIG cell's paths and run it yourself.
"""
import nbformat as nbf


def md(text):
    return nbf.v4.new_markdown_cell(text.strip())


def code(text):
    return nbf.v4.new_code_cell(text.strip())


CELLS = [
    md("""
# Dorcha master catalogue -- end-to-end pipeline demo

Builds a real master science catalogue from scratch: reads a snapshot,
halo catalogue, merger tree, and SHARK galaxy catalogue for one or more
host systems, runs every Phase 6b derived-quantity stage, validates the
result, writes a versioned `.h5` catalogue file, reads it back through
the unified interface (`CatalogueDataset`), and builds a small static
`WebReport`.

**Before running**: edit the CONFIG cell below with real paths to your
own Dorcha data. Nothing in `data/` ships a SHARK catalogue, so this
notebook can't run against the bundled example data the way
`UnifiedInterfaceDemo.ipynb` does -- see
[docs/phase6_remaining_work.md](../docs/phase6_remaining_work.md) and
[docs/master_catalogue.md](../docs/master_catalogue.md) for the design
this pipeline implements.
"""),
    code("""
# =====================================================================
# CONFIG -- edit these before running the rest of this notebook.
# =====================================================================

# Snapshot / halo catalogue / merger tree for one epoch (z=0 here; use a
# {redshift: path} dict for multiple epochs -- see api.Simulation).
SNAPSHOT_PATH = "/path/to/snapshot.hdf5"
HALO_PATH     = "/path/to/halo_catalogue.properties.0"
HALO_FILEFORMAT = "VELOCIraptor"          # VELOCIraptor | AHF | SubFind | ...
TREE_PATH     = "/path/to/tree.walkabletree.hdf5"
SNAPNUM       = 199                        # matches HALO_PATH/TREE_PATH

# SharkModel (galaxy_backend: shark) -- see analysistools/shark/model.py.
SHARK_MODEL_DIR      = "/path/to/shark/model_dir"       # e.g. ".../CDM/base_model"
SHARK_REDSHIFT_LIST  = "/path/to/redshift_list.txt"
SHARK_SUBVOLUMES     = "0-63"               # parse_subvolumes() string

# Per-format h-convention override, only if HaloCatalogue's own per-format
# default guess is wrong for this catalogue (None = don't override --
# see docs/unified_interface.md little_h/comoving discussion).
NATIVE_INCLUDES_H = None

# Host/satellite selection is demonstrated with a naive "N most massive
# halos as hosts, next few as satellites" rule below -- HaloExtractStage
# never guesses this itself (see its docstring), so replace section 1
# with your own real selection once you have one.
N_HOSTS = 2
N_SATELLITES_PER_HOST = 5

OUT_DIR = "dorcha_endtoend_output"   # catalogue file + web report land here
"""),
    code("""
import os

import numpy as np

import analysistools as at
from analysistools.catalogue.pipeline import CatalogueBuilder, HostJob
from analysistools.shark.common import _redshift_table, parse_subvolumes
from analysistools.shark.model import SharkModel

os.makedirs(OUT_DIR, exist_ok=True)
"""),
    md("""
## 0. Load everything through the unified interface

`Simulation` ties the snapshot/halo catalogue/tree/galaxy catalogue
together as one `Epoch`; the `SharkModel` is passed straight in as the
`galaxies=` spec (anything with `.at(z)` works, per `Simulation`'s
docstring) rather than a plain file path, since Dorcha's galaxy
properties are computed by SHARK, not read from a `galaxies.hdf5` file.
"""),
    code("""
redshift_table = _redshift_table(SHARK_REDSHIFT_LIST)
subvolumes = parse_subvolumes(SHARK_SUBVOLUMES)
shark_model = SharkModel(SHARK_MODEL_DIR, redshift_table, subvolumes,
                         label="dorcha")

load_kwargs = {"halos": {"fileformat": HALO_FILEFORMAT}}
if NATIVE_INCLUDES_H is not None:
    load_kwargs["halos"]["native_includes_h"] = NATIVE_INCLUDES_H

sim = at.Simulation(
    snapshots={0.0: SNAPSHOT_PATH},
    halos={0.0: HALO_PATH},
    trees=TREE_PATH,
    galaxies=shark_model,
    snapnums={0.0: SNAPNUM},
    load_kwargs=load_kwargs,
    label="dorcha",
)
epoch = sim.at(0.0)
print(epoch)
print(f"{len(epoch.halos):,} halos in the catalogue")
"""),
    md("""
## 1. Select host/satellite systems

**Demo-only selection** -- the `N_HOSTS` most massive halos, each paired
with its next `N_SATELLITES_PER_HOST` most massive (excluding halos
already claimed by an earlier host). This has no physical meaning
(it doesn't check proximity, isolation, or subhalo membership); replace
it with your project's real host/satellite identification before using
this notebook for anything beyond exercising the pipeline. `HaloExtractStage`
requires the selection explicitly for exactly this reason -- see its
docstring (`analysistools/catalogue/pipeline.py`).
"""),
    code("""
mass = np.asarray(epoch.halos["mass"])
order = np.argsort(mass)[::-1]

used = set()
jobs = []
for host_row in (int(r) for r in order[:N_HOSTS]):
    used.add(host_row)
    satellite_rows = [int(r) for r in order if int(r) not in used
                      ][:N_SATELLITES_PER_HOST]
    used.update(satellite_rows)
    jobs.append(HostJob(epoch=epoch, host_row=host_row,
                        satellite_rows=satellite_rows))
    print(f"host {host_row}: M={mass[host_row]:.3e}, "
         f"{len(satellite_rows)} satellite(s)")
"""),
    md("""
## 2. Project configuration

`CatalogueBuilder` reads a project YAML naming the schema version, galaxy
backend, ordered `derived_stages`, and each stage's required project-
specific thresholds (`stage_options` -- every one of these is a genuine
modelling choice this pipeline deliberately never picks a default for,
see each stage's docstring in `analysistools/catalogue/derived.py`).
Written out here so every value is visible and editable in the notebook
itself, rather than hidden in a separate file -- see `configs/dorcha.yaml`
for the same structure as a standalone template.
""") ,
    code("""
import yaml

# time_bin_edges' upper bound MUST reach at least the age of the universe
# at z=0 for this cosmology -- a plausible-looking round number like
# 13.5 Gyr is *not* automatically enough (this cosmology's age is
# ~13.80 Gyr; SHARK.rebin_sfh() silently drops any native-grid formed
# mass at lookback times beyond the requested upper edge, which then
# makes StellarMass appear to exceed the SFH-integrated formed mass --
# exactly the failure PhysicalValidator's stellar_mass_exceeds_formed_mass
# check exists to catch). Derive it from the model's own cosmology
# instead of hardcoding a guess, with a small safety margin.
age_of_universe = shark_model.age_at_z(0.0)
time_bin_upper_edge = age_of_universe + 0.1   # Gyr, safety margin

config = {
    "schema_version": "1.0",
    "galaxy_backend": "shark",
    "galaxy_backend_options": {
        "match_by": "position",
        "r_scale": 1.0,
        # "compute_photometry": True,   # real FSPS SED convolution --
        #                                # needs python-fsps installed
        #                                # with SPS_HOME set, and is much
        #                                # slower per satellite than
        #                                # every other field. See section
        #                                # 5 below for a standalone demo
        #                                # instead of enabling it here.
    },
    "derived_stages": [
        "halo_properties", "orbital_properties", "star_formation_history",
        "host_environment", "environment", "observability",
        "dorcha_specific",
    ],
    "stage_options": {
        "star_formation_history": {
            # yaml.safe_dump can't represent numpy scalars -- cast to
            # plain Python float, not just list(...) (which still leaves
            # each element an np.float64).
            "time_bin_edges": [float(t) for t in
                               np.linspace(0.0, time_bin_upper_edge, 14)],  # Gyr lookback
            "quenched_ssfr_threshold": 1.0e-11,  # 1/yr, absolute
        },
        "host_environment": {
            "isolation_radius_factor": 3.0,
            "pairing_mass_ratio_min": 0.3,
            "pairing_max_separation": 1.0,
            # mass_threshold/completeness_mass_threshold are compared
            # directly against epoch.halos["mass"], which for GADGET/
            # Arepo-family catalogues (SubFind, and conventionally AHF/
            # VELOCIraptor) is in units of 1e10 Msun(/h), NOT plain Msun,
            # regardless of little_h -- see EnvironmentStage's docstring.
            # 1.0e-2 here means "1e8 Msun" (1e8 / 1e10); if your halo
            # catalogue's native mass unit is different (check its own
            # docs), convert accordingly -- and if EnvironmentStage logs
            # "no neighbours found", check the mass range it reports
            # before assuming these numbers are right for your data.
            "completeness_mass_threshold": 1.0e-2,
        },
        "environment": {
            "mass_threshold": 1.0e-2,   # "1e8 Msun" -- see the comment above
            "aperture_radius": 5.0,
        },
        "observability": {
            "observer_pos": [0.0, 0.0, 0.0],   # box coordinates of the "mock Sun"
        },
        "dorcha_specific": {
            "reionisation_lookback_time": 13.0,
        },
    },
    "metadata": {
        "Metadata/catalogue_version": "0.1.0-demo",
    },
}

config_path = os.path.join(OUT_DIR, "project.yaml")
with open(config_path, "w") as fh:
    yaml.safe_dump(config, fh, sort_keys=False)
print(f"wrote {config_path}")
"""),
    md("""
## 3. Build the catalogue

One call: extract + cross-match + every derived stage for each host job,
concatenate into one project-wide table, validate
(`QualityControlStage`), and write (`WriteStage`) -- see
`CatalogueBuilder.run()`. Raises before writing if validation finds a
hard error; warnings (e.g. schema fields legitimately left unpopulated,
see `docs/phase6_remaining_work.md`) don't block the write.
"""),
    code("""
builder = CatalogueBuilder(config_path)
out_path = os.path.join(OUT_DIR, "dorcha_catalogue_v1.0.0-demo.h5")

report = builder.run(jobs, out_path=out_path)

print(f"validation: {'PASSED' if report.passed else 'FAILED'} "
     f"({len(report.errors)} error(s), {len(report.warnings)} warning(s))")
report.summary()
print(f"\\nwrote {out_path}")
"""),
    md("""
## 4. Read it back through the unified interface

`CatalogueDataset` (`kind="satellites"`) flattens every `Satellites/*`
leaf dataset into one column namespace, so it drops straight into the
same `Dataset` interface (and every existing `api/plotting.py` function)
as a snapshot or halo catalogue.
"""),
    code("""
cat = at.load(out_path, kind="satellites")
cat.summary()
"""),
    code("""
mpeak = cat["Mpeak"]
mstar = cat["StellarMass"]
print(f"{len(cat)} satellites total")
print(f"Mpeak:       {np.nanmin(mpeak):.3e} -- {np.nanmax(mpeak):.3e} Msun")
print(f"StellarMass: {np.nanmin(mstar):.3e} -- {np.nanmax(mstar):.3e} Msun "
     f"({np.sum(~np.isnan(mstar))}/{len(mstar)} resolved)")
"""),
    md("""
## 5. (Optional) FSPS photometry on a single galaxy

Skipped by default in the project config above (see the commented-out
`compute_photometry` line) since it needs a real python-fsps install
with `SPS_HOME` set and is far slower per galaxy than every other field.
This cell demonstrates it standalone against the same `SharkModel`,
independent of the full catalogue build -- flip the config's
`compute_photometry: true` and re-run section 3 instead if you want it
in the written catalogue itself.
"""),
    code("""
try:
    import fsps  # noqa: F401
    HAVE_FSPS = True
except ImportError:
    HAVE_FSPS = False

if HAVE_FSPS:
    from analysistools.shark.photometry import PhotometryPipeline

    pipe = PhotometryPipeline(shark_model, z_obs=0.0, bands=["v"],
                              progress=False)
    m_v = pipe.abs_mag_V(gal_indices=np.arange(min(5, pipe.n_galaxies)))
    print(f"M_V (first {len(m_v)} galaxies): {m_v}")
else:
    print("python-fsps not installed -- skipping (pip install fsps, "
         "then set SPS_HOME per the python-fsps install docs)")
"""),
    md("""
## 6. Shareable static report

`WebReport` embeds `plotly` figures into a small static site (GitHub
Pages / Zenodo / email-attachable, no server) -- add whatever figures are
useful for your project; this notebook adds one as a demonstration.
"""),
    code("""
import plotly.graph_objects as go

from analysistools.catalogue.webreport import WebReport

fig = go.Figure(data=[go.Histogram(x=np.log10(mpeak[~np.isnan(mpeak)]))])
fig.update_layout(xaxis_title="log10(Mpeak / Msun)", yaxis_title="N")

web_report = WebReport.from_catalogue(cat, title="Dorcha end-to-end demo")
web_report.add_figure("mpeak_histogram", fig, title="Mpeak distribution",
                      description="Peak halo mass over each satellite's "
                                  "history.", group="halo_properties")

report_dir = web_report.build(os.path.join(OUT_DIR, "report"))
print(f"wrote {report_dir / 'index.html'}")
"""),
    md("""
---
*Generated by `analysistools/catalogue/_build_endtoend_notebook.py`;
regenerate the structure after API changes, then re-fill the CONFIG cell
with real paths. See `docs/phase6_remaining_work.md` for what's still
deferred in the pipeline this notebook exercises.*
"""),
]


def build(path="notebooks/DorchaCatalogueEndToEnd.ipynb"):
    nb = nbf.v4.new_notebook()
    nb.cells = CELLS
    nb.metadata["kernelspec"] = {
        "display_name": "Python 3", "language": "python", "name": "python3"}
    nbf.write(nb, path)
    print(f"wrote {path} (not executed -- fill in the CONFIG cell first)")


if __name__ == "__main__":
    build()
