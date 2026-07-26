"""
shark.cli
---------
Command-line interface for the shark analysis package.

Usage
-----
python -m shark.cli \\
    -z redshift_list.txt -v '0-63' -o ./plots \\
    --model-dirs ./CDM/base_model ./WDM/base_model \\
    --redshifts 0 1 2 \\
    --labels CDM WDM \\
    --plots halo_mf stellar_mf sfr_ms size_mass bh_bulge
"""

from __future__ import annotations

import argparse
import os

import numpy as np

from . import common
from .common import _redshift_table, parse_subvolumes

from .model   import SharkModel
from .analysis import Analysis
from .plots    import Plotter


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m shark.cli",
        description="Run SHARK galaxy formation model analysis and produce plots.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    # Standard common.parse_args() flags
    parser.add_argument("-c", "--config",
                        help="SHARK configuration file")
    parser.add_argument("-m", "--model",
                        help="Model name")
    parser.add_argument("-s", "--simu",
                        help="Simulation name")
    parser.add_argument("-S", "--shark-dir",
                        help="SHARK base output directory")
    parser.add_argument("-z", "--redshift-file",
                        help="Redshift table file")
    parser.add_argument("-v", "--subvolumes", default="0",
                        help='Sub-volumes, e.g. "0-63" or "0,1,2" (default: 0)')
    parser.add_argument("-o", "--output-dir",
                        help="Output directory for plots")

    # Multi-model extensions
    parser.add_argument("--model-dirs", nargs="+", metavar="DIR", default=None,
                        help="Explicit model output directories (one or more).")
    parser.add_argument("--redshifts", nargs="+", type=float, default=[0.0],
                        metavar="Z",
                        help="Redshift values to process (default: 0).")
    parser.add_argument("--labels",  nargs="+", default=None,
                        help="Legend label per model directory.")
    parser.add_argument("--colours", nargs="+", default=None,
                        help="Hex colour per model directory.")

    # Plot selection
    _ALL_PLOTS = ["halo_mf", "stellar_mf", "sfr_ms", "size_mass", "bh_bulge"]
    parser.add_argument(
        "--plots", nargs="+", default=_ALL_PLOTS,
        choices=_ALL_PLOTS + ["all"],
        help=(
            "Which plots to produce. "
            f"Choose from: {', '.join(_ALL_PLOTS)}, all. "
            "Default: all."
        ),
    )

    return parser


def parse_args():
    parser = build_parser()
    opts   = parser.parse_args()

    # Resolve SHARK config / individual flags (mirrors common.parse_args logic)
    shark_dir = simu = model = redshift_file = None
    if opts.config:
        shark_dir, simu, model, redshift_file = common.read_configuration(opts.config)
    if opts.shark_dir:      shark_dir     = opts.shark_dir
    if opts.simu:           simu          = opts.simu
    if opts.model:          model         = opts.model
    if opts.redshift_file:  redshift_file = opts.redshift_file

    if not redshift_file:
        parser.error("A redshift file is required: -z/--redshift-file or -c/--config.")

    # Build model directory list
    if opts.model_dirs:
        model_dirs = opts.model_dirs
    else:
        if not (shark_dir and simu and model):
            parser.error(
                "Either --model-dirs or all of -S/-s/-m (or -c) must be provided."
            )
        model_dirs = [os.path.join(shark_dir, simu, model)]

    n = len(model_dirs)
    if opts.labels  and len(opts.labels)  != n:
        parser.error(f"--labels:  expected {n} value(s), got {len(opts.labels)}.")
    if opts.colours and len(opts.colours) != n:
        parser.error(f"--colours: expected {n} value(s), got {len(opts.colours)}.")

    outdir = opts.output_dir or (
        os.path.join(shark_dir, "Plots", simu, model)
        if (shark_dir and simu and model) else "."
    )

    redshift_table = _redshift_table(redshift_file)
    subvols        = parse_subvolumes(opts.subvolumes)

    # Expand "all"
    plots = opts.plots
    if "all" in plots:
        plots = ["halo_mf", "stellar_mf", "sfr_ms", "size_mass", "bh_bulge"]

    print(f"Output directory : {outdir}")
    print(f"Subvolumes       : {sorted(subvols)}")
    print(f"Plots requested  : {plots}")

    return (
        model_dirs,
        redshift_table,
        subvols,
        np.array(opts.redshifts),
        outdir,
        opts.labels,
        opts.colours,
        plots,
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    (
        model_dirs,
        redshift_table,
        subvols,
        redshifts,
        outdir,
        labels,
        colours,
        plots,
    ) = parse_args()

    labels  = labels  or [None] * len(model_dirs)
    colours = colours or [None] * len(model_dirs)

    # Build and preload models
    models = []
    for m_dir, lbl, col in zip(model_dirs, labels, colours):
        print(f"Loading model: {m_dir}")
        m = SharkModel(
            model_dir      = m_dir,
            redshift_table = redshift_table,
            subvols        = subvols,
            label          = lbl,
            colour         = col,
        )
        # Preload only the fields needed for requested plots
        needed = _fields_for_plots(plots)
        m.preload(needed, redshifts)
        models.append(m)

    # Run analysis and plot
    analysis = Analysis(models)
    plotter  = Plotter(
        outdir  = outdir,
        colours = [m.colour for m in models],
    )

    _dispatch(analysis, plotter, plots, redshifts)


def _fields_for_plots(plots: list) -> list:
    """Return the minimal set of SharkModel fields needed for *plots*."""
    fields = {"type", "mstars_disk", "mstars_bulge"}
    if "halo_mf"   in plots: fields |= {"mvir_hosthalo", "mvir_subhalo"}
    if "sfr_ms"    in plots: fields |= {"sfr_disk", "sfr_bulge"}
    if "size_mass" in plots: fields |= {"rstar_disk", "rstar_bulge"}
    if "bh_bulge"  in plots: fields |= {"mbh"}
    return list(fields)


def _dispatch(analysis: Analysis, plotter: Plotter, plots: list,
              redshifts: np.ndarray) -> None:
    """Call the appropriate compute + plot pair for each requested plot type."""
    z = list(redshifts)

    if "halo_mf" in plots:
        plotter.plot_halo_mf(analysis.compute_halo_mf(z))

    if "stellar_mf" in plots:
        plotter.plot_stellar_mf(analysis.compute_stellar_mf(z))

    if "sfr_ms" in plots:
        plotter.plot_sfr_main_sequence(analysis.compute_sfr_main_sequence(z))

    if "size_mass" in plots:
        plotter.plot_size_mass(analysis.compute_size_mass(z))

    if "bh_bulge" in plots:
        plotter.plot_bh_bulge(analysis.compute_bh_bulge(z))


if __name__ == "__main__":
    main()
