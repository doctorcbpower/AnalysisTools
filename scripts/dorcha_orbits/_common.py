"""Shared configuration for the Dorcha orbit-analysis scripts.

Environment variables:
  DORCHA_TREE  path to the SubFind-HBT tree file
               (default ~/Currentwork/Dorcha/trees.hdf5)
  DORCHA_OUT   directory for figures / pickles / sidecars
               (default ./dorcha_orbit_output, created if missing)

Tree conventions verified for the Dorcha tree (see docs/orbit_tracking_report.md):
positions comoving Mpc/h, velocities physical peculiar km/s, 1e10 Msun/h, h = 0.678;
load with native_is_comoving=True, native_includes_h=True.
"""
import os

TREE = os.environ.get("DORCHA_TREE",
                      os.path.expanduser("~/Currentwork/Dorcha/trees.hdf5"))
OUT_DIR = os.environ.get("DORCHA_OUT",
                         os.path.join(os.getcwd(), "dorcha_orbit_output"))
os.makedirs(OUT_DIR, exist_ok=True)
OUT = OUT_DIR if OUT_DIR.endswith(os.sep) else OUT_DIR + os.sep   # for string-concatenated paths
