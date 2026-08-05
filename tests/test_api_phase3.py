#!/usr/bin/env python3
"""
Phase-3 tests: Simulation, Epoch, cross-matching (particles_in_halo,
galaxies_in_halo, track_of), and the periodic position matcher.
"""
import os

import h5py
import numpy as np
import pytest

import analysistools as at
from analysistools.api.simulation import match_positions

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SNAP = os.path.join(HERE, "data", "snap_0031.hdf5")
VRCAT = os.path.join(HERE, "data", "halos",
                     "snap_0031.VELOCIraptor.properties.0")
TREE = os.path.join(HERE, "data", "VELOCIraptor.walkabletree.hdf5")

have_data = all(os.path.exists(p) for p in (SNAP, VRCAT, TREE))
needs_data = pytest.mark.skipif(not have_data, reason="example data missing")


def make_sim(**extra):
    # native_includes_h=False: this bundled VELOCIraptor catalogue was run
    # against the SWIFT snapshot above and inherits its h-free convention,
    # not the h-included default HaloCatalogue otherwise guesses for
    # VELOCIraptor -- see docs/unified_interface.md#units-comoving-vs-physical-and-little-h.
    extra.setdefault("load_kwargs", {"halos": {"native_includes_h": False}})
    return at.Simulation(
        snapshots={0.0: SNAP},
        halos={0.0: VRCAT},
        trees=TREE,
        snapnums={0.0: 31},
        label="test",
        **extra,
    )


# ---------------------------------------------------------------------------
# match_positions
# ---------------------------------------------------------------------------

def test_match_positions_periodic_wrap():
    box = 100.0
    pts = np.array([[1.0, 50.0, 50.0], [99.0, 50.0, 50.0],
                    [50.0, 50.0, 50.0]])
    # centre at the box edge: the first two points are within r=3 only
    # via periodic wrapping
    idx = match_positions(pts, [0.0, 50.0, 50.0], 3.0, boxsize=box)
    assert set(idx) == {0, 1}
    idx = match_positions(pts, [0.0, 50.0, 50.0], 3.0, boxsize=None)
    assert set(idx) == {0}


# ---------------------------------------------------------------------------
# Simulation / Epoch basics
# ---------------------------------------------------------------------------

@needs_data
def test_simulation_epoch_components():
    sim = make_sim()
    epoch = sim.at(0.0)
    assert sim.at(0.0) is epoch                    # cached
    assert epoch.snapnum == 31
    assert isinstance(epoch.snapshot, at.SnapshotDataset)
    assert isinstance(epoch.halos, at.HaloCatalogue)
    assert epoch.galaxies is None                  # not configured
    assert isinstance(epoch.tree, at.MergerTree)
    assert epoch.boxsize == pytest.approx(142.2475, rel=1e-4)


@needs_data
def test_simulation_nearest_redshift():
    sim = make_sim()
    epoch = sim.at(0.3)                            # nearest available: z=0
    assert len(epoch.halos) == 522


@needs_data
def test_prebuilt_components_accepted():
    cat = at.load(VRCAT, snapnum=31, native_includes_h=False)
    sim = at.Simulation(halos=cat, label="prebuilt")
    assert sim.at(0.0).halos is cat


# ---------------------------------------------------------------------------
# particles_in_halo
# ---------------------------------------------------------------------------

@needs_data
def test_particles_in_halo_radius_match():
    sim = make_sim()
    epoch = sim.at(0.0)
    cat = epoch.halos
    row = int(np.argmax(cat["mass"]))
    parts = epoch.particles_in_halo(index=row)
    assert len(parts) > 0
    # every selected particle is inside r200 (periodic frame)
    c = cat["pos"][row]
    r = cat["radius"][row]
    box = epoch.boxsize
    d = (parts["pos"] - c + 0.5 * box) % box - 0.5 * box
    assert np.all(np.linalg.norm(d, axis=1) <= r * (1 + 1e-6))
    # by halo_id gives the same particles
    parts2 = epoch.particles_in_halo(int(cat["halo_id"][row]))
    assert len(parts2) == len(parts)


@needs_data
def test_particles_in_halo_r_scale_and_species():
    sim = make_sim()
    epoch = sim.at(0.0)
    row = int(np.argmax(epoch.halos["mass"]))
    small = epoch.particles_in_halo(index=row, r_scale=0.5)
    big = epoch.particles_in_halo(index=row, r_scale=2.0, species="dm")
    assert len(small) < len(big)


@needs_data
def test_particles_in_halo_groupid_unavailable():
    sim = make_sim()
    with pytest.raises(ValueError):
        sim.at(0.0).particles_in_halo(index=0, match_by="groupid")


# ---------------------------------------------------------------------------
# galaxies_in_halo (synthetic SHARK fixture built at real halo positions)
# ---------------------------------------------------------------------------

@pytest.fixture
def sim_with_galaxies(tmp_path):
    cat = at.load(VRCAT, snapnum=31, native_includes_h=False)
    n_halos = 5
    per_halo = 4
    rng = np.random.default_rng(1)
    rows = np.argsort(np.asarray(cat["mass"]))[-n_halos:]
    centres = np.asarray(cat["pos"])[rows]
    radii = np.asarray(cat["radius"])[rows]
    ids = np.asarray(cat["halo_id"])[rows]

    pos, id_halo = [], []
    for c, r, hid in zip(centres, radii, ids):
        offsets = rng.normal(0, 0.2 * r, (per_halo, 3))
        pos.append(c + offsets)
        id_halo += [hid] * per_halo
    pos = np.vstack(pos)
    n = len(pos)

    path = tmp_path / "galaxies.hdf5"
    with h5py.File(path, "w") as f:
        g = f.create_group("galaxies")
        g.create_dataset("position_x", data=pos[:, 0])
        g.create_dataset("position_y", data=pos[:, 1])
        g.create_dataset("position_z", data=pos[:, 2])
        for name in ("velocity_x", "velocity_y", "velocity_z"):
            g.create_dataset(name, data=rng.normal(0, 100, n))
        g.create_dataset("mstars_disk", data=rng.lognormal(22, 1, n))
        g.create_dataset("mstars_bulge", data=rng.lognormal(21, 1, n))
        g.create_dataset("id_galaxy", data=np.arange(n, dtype=np.int64))
        g.create_dataset("id_halo", data=np.array(id_halo, dtype=np.int64))
        g.create_dataset("type", data=np.zeros(n, dtype=np.int32))
        f.create_dataset("cosmology/h", data=0.6751)
        f.create_dataset("run_info/redshift", data=0.0)

    sim = make_sim(galaxies={0.0: str(path)})
    return sim, rows, ids


@needs_data
def test_galaxies_in_halo_position_default(sim_with_galaxies):
    sim, rows, ids = sim_with_galaxies
    epoch = sim.at(0.0)
    found = epoch.galaxies_in_halo(index=int(rows[-1]))
    assert len(found) >= 1
    # all found galaxies carry the right id_halo (they were planted there)
    assert np.all(found["id_halo"] == ids[-1])


@needs_data
def test_galaxies_in_halo_by_id(sim_with_galaxies):
    sim, rows, ids = sim_with_galaxies
    epoch = sim.at(0.0)
    found = epoch.galaxies_in_halo(int(ids[-1]), match_by="id")
    assert len(found) == 4                          # exactly as planted
    assert np.all(found["id_halo"] == ids[-1])


# ---------------------------------------------------------------------------
# track_of
# ---------------------------------------------------------------------------

@needs_data
def test_track_of_from_epoch_and_sim():
    sim = make_sim()
    epoch = sim.at(0.0)
    tr = epoch.track_of(index=0)
    assert isinstance(tr, at.TrackDataset)
    assert len(tr) > 1
    assert tr["snapnum"][-1] == 31
    # the epoch's catalogue is auto-linked, so properties are filled at z=0
    assert np.isfinite(tr["mass"][-1])
    assert np.isclose(tr["mass"][-1], epoch.halos["mass"][0])
    tr2 = sim.track_of(index=0, redshift=0.0)
    assert len(tr2) == len(tr)
