# --- submodule imports ---
from . import snapshot_tools
from . import snapio_hdf5
from . import snapio_binary
from . import halo_tools
from . import haloio_subfind
from . import haloio_swiftfof
from . import haloio_ahf
from . import haloio_velociraptor
from . import ic_tools
from . import merger_tree_tools
from . import merger_tree_types
from . import treeio_subfind
from . import treeio_treefrog
from . import treeio_ahf
from . import profile_tools
from . import gridding_tools
from . import fdm_field_tools
from . import halo_model

# --- public class imports ---
from .snapshot_tools import SnapshotTools
from .halo_tools import HaloTools
from .halo_model import HaloModel

# --- unified interface (analysistools.api; see DEVELOPMENT.md) ---
from . import api
from .api import load, Dataset, SnapshotDataset, HaloCatalogue
from .ic_tools import ICTools
from .merger_tree_tools import MergerTreeTools
from .gridding_tools import GriddingTools
from .profile_tools import ProfileTools
from .fdm_field_tools import FDMFieldTools

__all__ = [
    # submodules
    "snapshot_tools",
    "snapio_hdf5",
    "snapio_binary",
    "halo_tools",
    "haloio_subfind",
    "haloio_swiftfof",
    "haloio_ahf",
    "haloio_velociraptor",
    "ic_tools",
    "merger_tree_tools",
    "merger_tree_types",
    "treeio_subfind",
    "treeio_treefrog",
    "treeio_ahf",
    "profile_tools",
    "gridding_tools",
    "fdm_field_tools",
    "halo_model",
    "api",
    # unified interface
    "load",
    "Dataset",
    "SnapshotDataset",
    "HaloCatalogue",
    # classes
    "SnapshotTools",
    "HaloTools",
    "HaloModel",
    "ICTools",
    "MergerTreeTools",
    "GriddingTools",
    "ProfileTools",
    "FDMFieldTools",
]
