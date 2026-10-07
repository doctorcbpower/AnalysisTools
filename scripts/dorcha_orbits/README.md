# Dorcha orbit analysis scripts

Analysis scripts behind `docs/orbit_tracking_report.md`, run on the Dorcha SubFind-HBT tree
(host = group 0, SubhaloNr 0). They are examples/provenance for the report's numbers and figures,
not library code; the reusable pieces live in `analysistools/orbits/`.

```bash
export DORCHA_TREE=~/Currentwork/Dorcha/trees.hdf5   # default
export DORCHA_OUT=./dorcha_orbit_output              # default: figures, pickles, sidecar
PYTHONPATH=. python scripts/dorcha_orbits/01_track_host.py
```

Run from the repository root with `PYTHONPATH=.` (or install the package). Do **not** rely on the
system Python's installed `analysistools`: it may be an old editable copy (see the report, section 7).

| Script | Produces |
|---|---|
| `01_track_host.py` | Sidecar `orbits_host0.h5`; summary statistics for all z=0 FOF members |
| `02a_validate_orphans_truncation.py`, `02b_validate_orphan_normalisation.py` | Orphan truncation test (tables in report section 4.2) |
| `03_orphan_fringe_check.py` | Why orphan continuation is gated on `orphan_max_r_over_ref` |
| `04_plot_satellite_orbits.py [a_min a_max]` | `satellite_orbits_a<min>_<max>.png` |
| `05_plot_mpeak_vs_mcur.py` | `mpeak_vs_mcur.png` |
| `06_plot_infall_and_pericentre.py` | `infall_params.png`, `peri_vs_now.png` (R200c units, peculiar velocities; *not* Benson's convention) |
| `07a_satellite_selection.py`, `07b_threshold_counts.py` | Satellite counts and threshold tables |
| `08_benson_events.py bn\|200c` | `benson_events_<mode>.pkl` (all first crossings, Benson quantities) |
| `09_benson_compare.py` | `benson_variants.png` and the variant table (needs 08 for both modes) |
| `10_diag_group_membership.py`, `11_diag_mpeak_outlier.py` | Diagnostics quoted in the report |

Selection cuts used by the plots (`min_mass_fraction=1e-4`, `min_particles=50`, peak/final bases) are
set in the scripts' `OrbitConfig(...)`; the catalogue defaults are in `configs/dorcha.yaml`.
