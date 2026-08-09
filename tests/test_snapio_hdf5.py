"""
Tests for analysistools.snapio_hdf5.write_hdf5 -- specifically the
species-local-indexing bug in _write_metallicities/
_write_star_formation_rates/_write_stellar_ages/_write_stellar_init_mass
(and the SFR write being gated on the wrong particle type), found and
fixed alongside the gas_Z/stellar_Z/initmass/sfr/age attribute-naming
mismatch (see api/snapshot.py's _RAW_ATTR_ALIASES and
docs/phase6_remaining_work.md).

These four fields are each sized to only their own species' particle
count (gas-only: gas_metallicity, gas_sfr; star-only:
stellar_metallicity, stellarage, stellarinitmass -- see
_allocate_extra_block_memory), so writing them needs a species-local
index, not the global full-particle idx used for pos/vel/mass. Uses a
non-identity `idx` (reversed particle order) so a wrong index either
raises IndexError (out of the species-only array's bounds) or silently
writes the wrong particle's value -- either way this catches the bug,
unlike an identity idx which could accidentally look correct.
"""
import h5py
import numpy as np
import pytest

from analysistools.snapio_hdf5 import write_hdf5

GAS, DM, STAR = 0, 1, 4


def _make_writer(idx, idx_type, **extra_attrs):
    writer = write_hdf5(output_convention="SWIFT", idx=idx, idx_type=idx_type)
    writer.unit_length_in_cgs = 1.0
    writer.unit_mass_in_cgs = 1.0
    writer.unit_velocity_in_cgs = 1.0
    writer.unit_sfr_in_cgs = 1.0
    for name, value in extra_attrs.items():
        # "ptype" is deliberately not forwarded: write_hdf5.__init__
        # already sets self.ptype = idx_type (output-position order,
        # parallel to self.idx -- see _original_ptype's docstring for why
        # that's different from the fixture's original-row-order `ptype`
        # array); overwriting it here would silently break every
        # species-local-index test.
        if name == "ptype":
            continue
        setattr(writer, name, value)
    return writer


@pytest.fixture
def six_particle_snapshot():
    # rows: 0=gas, 1=star, 2=gas, 3=star, 4=gas, 5=dm
    ptype = np.array([GAS, STAR, GAS, STAR, GAS, DM], dtype=np.int32)
    # species-only arrays, populated in the order each species' particles
    # appear when scanning `ptype` in ascending row order (rows 0,2,4 for
    # gas; rows 1,3 for star) -- matching how the reader actually fills
    # them (_read_particle_data_single_file).
    gas_metallicity = np.array([0.01, 0.02, 0.03])       # rows 0, 2, 4
    gas_sfr = np.array([1.0, 2.0, 3.0])                  # rows 0, 2, 4
    stellar_metallicity = np.array([0.5, 0.6])           # rows 1, 3
    stellarage = np.array([0.1, 0.2])                    # rows 1, 3
    stellarinitmass = np.array([100.0, 200.0])           # rows 1, 3
    return dict(ptype=ptype, gas_metallicity=gas_metallicity,
               gas_sfr=gas_sfr, stellar_metallicity=stellar_metallicity,
               stellarage=stellarage, stellarinitmass=stellarinitmass)


def test_species_local_idx_maps_global_rows_to_species_local_positions(
        six_particle_snapshot):
    # reversed output order: idx=[5,4,3,2,1,0] -> idx_type=[dm,gas,star,gas,star,gas]
    idx = np.arange(6)[::-1]
    idx_type = six_particle_snapshot["ptype"][idx]
    writer = _make_writer(idx, idx_type, **six_particle_snapshot)

    gas_mask = idx_type == GAS      # output positions 1, 3, 5
    star_mask = idx_type == STAR    # output positions 2, 4

    gas_local = writer._species_local_idx(gas_mask, GAS)
    star_local = writer._species_local_idx(star_mask, STAR)

    # global gas rows selected (in output order): 4, 2, 0 -> local 2, 1, 0
    np.testing.assert_array_equal(gas_local, [2, 1, 0])
    # global star rows selected (in output order): 3, 1 -> local 1, 0
    np.testing.assert_array_equal(star_local, [1, 0])


def test_write_and_read_back_gas_and_stellar_metallicity(
        tmp_path, six_particle_snapshot):
    idx = np.arange(6)[::-1]
    idx_type = six_particle_snapshot["ptype"][idx]
    writer = _make_writer(idx, idx_type, **six_particle_snapshot)

    out = str(tmp_path / "snap")
    writer.write_hdf5_snapshot(out)

    with h5py.File(out + ".hdf5", "r") as f:
        gas_z = f["PartType0/Metallicity"][()]
        star_z = f["PartType4/Metallicity"][()]

    # see test_species_local_idx_... above for the derivation
    np.testing.assert_allclose(gas_z, [0.03, 0.02, 0.01])
    np.testing.assert_allclose(star_z, [0.6, 0.5])


def test_write_and_read_back_star_formation_rate_gated_on_gas_type(
        tmp_path, six_particle_snapshot):
    idx = np.arange(6)[::-1]
    idx_type = six_particle_snapshot["ptype"][idx]
    writer = _make_writer(idx, idx_type, **six_particle_snapshot)

    out = str(tmp_path / "snap")
    writer.write_hdf5_snapshot(out)

    with h5py.File(out + ".hdf5", "r") as f:
        assert "StarFormationRate" in f["PartType0"]   # gas group -- correct gating
        assert "StarFormationRate" not in f["PartType4"]  # not the star group
        sfr = f["PartType0/StarFormationRate"][()]

    np.testing.assert_allclose(sfr, [3.0, 2.0, 1.0])  # rows 4, 2, 0


def test_write_and_read_back_stellar_age_and_init_mass(
        tmp_path, six_particle_snapshot):
    idx = np.arange(6)[::-1]
    idx_type = six_particle_snapshot["ptype"][idx]
    writer = _make_writer(idx, idx_type, **six_particle_snapshot)

    out = str(tmp_path / "snap")
    writer.write_hdf5_snapshot(out)

    with h5py.File(out + ".hdf5", "r") as f:
        age = f["PartType4/StellarFormationTime"][()]
        initmass = f["PartType4/StellarInitMass"][()]

    np.testing.assert_allclose(age, [0.2, 0.1])           # rows 3, 1
    np.testing.assert_allclose(initmass, [200.0, 100.0])  # rows 3, 1


def test_identity_idx_still_correct(tmp_path, six_particle_snapshot):
    # sanity check: the simplest (identity) ordering, where the old bug's
    # global indices happened to sometimes stay in-bounds, is also correct.
    idx = np.arange(6)
    idx_type = six_particle_snapshot["ptype"]
    writer = _make_writer(idx, idx_type, **six_particle_snapshot)

    out = str(tmp_path / "snap")
    writer.write_hdf5_snapshot(out)

    with h5py.File(out + ".hdf5", "r") as f:
        gas_z = f["PartType0/Metallicity"][()]
        star_z = f["PartType4/Metallicity"][()]
        sfr = f["PartType0/StarFormationRate"][()]

    np.testing.assert_allclose(gas_z, [0.01, 0.02, 0.03])
    np.testing.assert_allclose(star_z, [0.5, 0.6])
    np.testing.assert_allclose(sfr, [1.0, 2.0, 3.0])
