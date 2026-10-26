"""Tests of the AHF readers on small synthetic files."""

import numpy as np
import pytest

from analysistools.haloio_ahf import (AHF_PTYPE_STAR, read_ahf, read_ahf_particles,
                                      read_ahf_substructure)

HEADER = ("ID(1) hostHalo(2) numSubStruct(3) Mvir(4) npart(5) Xc(6) Yc(7) Zc(8) "
          "VXc(9) VYc(10) VZc(11) Rvir(12) Rmax(13) r2(14) mbp_offset(15) "
          "com_offset(16) Vmax(17)")


def _write_halos(path, rows, header=HEADER):
    np.savetxt(path, rows, header=header, comments="#")


def _row(halo_id, host, nsub, mass, npart, x=0.0, vx=0.0, rvir=100.0, vmax=300.0):
    return [halo_id, host, nsub, mass, npart, x, 2 * x, 3 * x, vx, 0, 0, rvir, 0, 0, 0.5, 0.25, vmax]


def test_read_ahf_selects_columns_by_name(tmp_path):
    path = str(tmp_path / "a.AHF_halos")
    _write_halos(path, [_row(1, 0, 1, 1e14, 50, x=10.0, vx=7.0), _row(2, 1, 0, 1e12, 5, x=11.0)])
    meta, halos, subs = read_ahf(path)
    assert meta == {"Nhalos": 2, "Nsubhalos": 1}
    assert list(halos["HaloID"]) == [1, 2]
    assert list(halos["Pos"][0]) == [10.0, 20.0, 30.0]
    assert halos["Vel"][0, 0] == 7.0
    assert halos["Rvir"][0] == 100.0
    assert halos["com_offset"][0] == 0.25 and halos["Vmax"][0] == 300.0
    # The subhalo table is built from the halos that are inside another halo.
    assert list(subs["HaloID"]) == [2] and list(subs["ParentHaloID"]) == [1]
    assert subs["Pos"].shape == (1, 3)


def test_columns_follow_the_header_not_the_position(tmp_path):
    # Swap the order of Mvir and npart in both header and data.
    header = HEADER.replace("Mvir(4) npart(5)", "npart(4) Mvir(5)")
    row = _row(1, 0, 0, 1e14, 50)
    row[3], row[4] = row[4], row[3]
    path = str(tmp_path / "a.AHF_halos")
    _write_halos(path, [row], header=header)
    _, halos, _ = read_ahf(path)
    assert halos["Mass"][0] == 1e14 and halos["NumPart"][0] == 50


def test_single_row_file(tmp_path):
    path = str(tmp_path / "a.AHF_halos")
    _write_halos(path, [_row(1, 0, 0, 1e14, 50)])
    _, halos, subs = read_ahf(path)
    assert halos["Pos"].shape == (1, 3) and subs == {}


def test_hosthalo_minus_one_means_not_a_subhalo(tmp_path):
    path = str(tmp_path / "a.AHF_halos")
    _write_halos(path, [_row(1, -1, 1, 1e14, 50), _row(2, 1, 0, 1e12, 5)])
    meta, _, subs = read_ahf(path)
    assert meta["Nsubhalos"] == 1 and list(subs["HaloID"]) == [2]


def test_unrecognised_header_is_an_error_that_shows_the_header(tmp_path):
    path = str(tmp_path / "a.AHF_halos")
    _write_halos(path, [_row(1, 0, 0, 1e14, 50)], header="ID(1) something(2)")
    with pytest.raises(ValueError, match="something"):
        read_ahf(path)


def test_wrong_extension_is_refused(tmp_path):
    with pytest.raises(ValueError, match="AHF_halos"):
        read_ahf(str(tmp_path / "a.txt"))


def test_ragged_substructure_file(tmp_path):
    path = str(tmp_path / "a.AHF_substructure")
    with open(path, "w") as handle:
        handle.write("#header\n1 4\n2 3 4 5\n9 0\n7 1\n8\n")
    result = read_ahf_substructure(path)
    assert list(result[1]) == [2, 3, 4, 5]
    assert result[9].size == 0
    assert list(result[7]) == [8]


def test_substructure_count_mismatch_names_the_line(tmp_path):
    path = str(tmp_path / "a.AHF_substructure")
    with open(path, "w") as handle:
        handle.write("1 4\n2 3 4\n")
    with pytest.raises(ValueError, match="line 2"):
        read_ahf_substructure(path)


def test_read_halos_does_not_touch_the_substructure_file(tmp_path):
    # The old reader failed on ragged substructure files. This one never reads them.
    path = str(tmp_path / "a.AHF_halos")
    _write_halos(path, [_row(1, 0, 1, 1e14, 50), _row(2, 1, 0, 1e12, 5)])
    with open(str(tmp_path / "a.AHF_substructure"), "w") as handle:
        handle.write("1 6630\n" + " ".join(["2"] * 6630) + "\n3 1\n4\n")
    meta, _, _ = read_ahf(path)
    assert meta["Nhalos"] == 2


def _write_particles(path, blocks, with_type=True):
    with open(path, "w") as handle:
        handle.write(f"{len(blocks)}\n")
        for halo_id, particles in blocks:
            handle.write(f"{len(particles)} {halo_id}\n")
            for pid, ptype in particles:
                handle.write(f"{pid} {ptype}\n" if with_type else f"{pid}\n")


def test_particles_wanted_blocks_only(tmp_path):
    path = str(tmp_path / "a.AHF_particles")
    _write_particles(path, [(10, [(1, 1), (2, 4), (3, 4), (4, 0)]), (11, [(5, 4)]), (12, [(7, 4)])])
    result = read_ahf_particles(path, wanted=[10, 12])
    assert sorted(result) == [10, 12]
    assert result[10]["n_declared"] == 4
    stars = result[10]["ids"][result[10]["types"] == AHF_PTYPE_STAR]
    assert list(stars) == [2, 3]


def test_particles_read_everything_when_nothing_is_wanted(tmp_path):
    path = str(tmp_path / "a.AHF_particles")
    _write_particles(path, [(10, [(1, 4)]), (11, [])])
    result = read_ahf_particles(path)
    assert sorted(result) == [10, 11] and result[11]["ids"].size == 0


def test_particles_missing_halo_and_missing_type_column(tmp_path):
    path = str(tmp_path / "a.AHF_particles")
    _write_particles(path, [(10, [(1, 4)])])
    with pytest.raises(ValueError, match="not found"):
        read_ahf_particles(path, wanted=[99])
    _write_particles(path, [(10, [(1, 4)])], with_type=False)
    with pytest.raises(ValueError, match="type column"):
        read_ahf_particles(path, wanted=[10])
