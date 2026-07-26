"""
Deprecated import shim: the shark package now lives inside analysistools.

    import shark                    ->  from analysistools import shark
    from shark.model import ...     ->  from analysistools.shark.model import ...

Existing imports keep working through this shim (with a DeprecationWarning);
new code should import analysistools.shark directly.
"""
import importlib
import sys
import warnings

warnings.warn(
    "'import shark' is deprecated: the package has moved to "
    "'analysistools.shark'. This top-level shim will be removed in a "
    "future release.",
    DeprecationWarning,
    stacklevel=2,
)

_pkg = importlib.import_module("analysistools.shark")

# Alias every submodule under its old name so that `shark.model` and
# `analysistools.shark.model` are the *same* module object (identity of
# classes/exceptions preserved). Modules with heavy optional dependencies
# (e.g. photometry -> fsps) are skipped if they fail to import.
for _sub in ("common", "model", "analysis", "plots", "cli",
             "photometry", "photometry.io", "photometry.sps",
             "photometry.cosmology", "photometry.photometry"):
    try:
        sys.modules[f"{__name__}.{_sub}"] = importlib.import_module(
            f"analysistools.shark.{_sub}")
    except ImportError:
        pass

# Finally alias the package itself.
sys.modules[__name__] = _pkg
