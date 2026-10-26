"""write_snapshot and ICTools must work on data read via read_snapshot()."""
from pathlib import Path

import numpy as np
import pytest

from analysistools import SnapshotTools
from analysistools.ic_tools import ICTools
from analysistools.snapshot_tools import select_particles

SNAP = Path(__file__).resolve().parent.parent / "data" / "snap_0031.hdf5"
pytestmark = pytest.mark.skipif(not SNAP.exists(), reason="bundled snapshot missing")


def test_read_snapshot_mirrors_data_onto_tools():
    snap = SnapshotTools(snapfileformat="HDF5", convention="SWIFT")
    data = snap.read_snapshot(str(SNAP))
    assert snap.pos is data.pos and snap.ptype is data.ptype


def test_write_snapshot_roundtrip(tmp_path):
    snap = SnapshotTools(snapfileformat="HDF5", convention="SWIFT")
    data = snap.read_snapshot(str(SNAP))
    idx = select_particles(data.pos, data.box_size / 2 * np.ones(3), 20.0, "spherical")
    out = tmp_path / "sub"
    snap.write_snapshot(str(out), idx=idx, idx_type=np.ones(len(idx), dtype=np.int64),
                        convention="GADGET4", blocks_to_write=["pos", "vel", "pids", "mass"])
    back = SnapshotTools(snapfileformat="HDF5", convention="GADGET4").read_snapshot(str(out) + ".hdf5")
    assert back.pos.shape == (len(idx), 3)
    np.testing.assert_allclose(np.sort(back.pos[:, 0]), np.sort(data.pos[idx, 0]))


def test_ictools_merge_taper_rescale_excise(tmp_path):
    a = ICTools(snapfileformat="HDF5", convention="SWIFT"); da = a.read_snapshot(str(SNAP))
    b = ICTools(snapfileformat="HDF5", convention="SWIFT"); b.read_snapshot(str(SNAP))
    mass0 = da.mass.sum()
    comb = a.merge([b])
    assert len(comb.pos) == 2 * len(da.pos)
    comb.taper_halo(mass_fraction=0.95)
    assert 0 < len(comb.pos) <= 2 * len(da.pos)
    comb.rescale_component_mass(ptype=comb.dm_type, target_mass=1.0e3)
    assert comb.mass.sum() == pytest.approx(1.0e3)
    assert da.mass.sum() == mass0          # rescale must not alter the returned data
    assert len(comb.excise_type(comb.dm_type).pos) == 0
