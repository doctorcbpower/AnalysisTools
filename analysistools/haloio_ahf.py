#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
haloio_ahf.py

Readers for AHF (Amiga Halo Finder) output files.

    <prefix>.AHF_halos         one row per halo (hosts and subhalos alike)
    <prefix>.AHF_substructure  for each host, the IDs of its subhalos (ragged lines)
    <prefix>.AHF_particles     for each halo, the IDs and types of its particles

Functions
---------
read_ahf               the AHF_halos table, as (metadata, halos, subhalos)
read_ahf_substructure  the AHF_substructure file, as {host ID: array of subhalo IDs}
read_ahf_particles     selected blocks of the AHF_particles file

Columns of AHF_halos are selected by the names in the "#" header line, for example
"ID(1) hostHalo(2) numSubStruct(3) Mvir(4) ...", not by position, so a different AHF
version cannot silently shift the columns. If a required name is missing, the header
is printed in the error.

All values are returned exactly as stored in the file. Unit conversions (comoving and
little-h) are applied centrally by HaloTools.standardise_names(), not here.
"""

import re
from itertools import islice
from typing import Any, Dict, Iterable, Optional, Tuple

import numpy as np

# "Name(column number)", as in the AHF_halos header. Column numbers start at 1.
_HEADER_FIELD = re.compile(r"([A-Za-z_][A-Za-z_0-9]*)\((\d+)\)")

# Names that must be present in the header.
_REQUIRED = ["ID", "hostHalo", "numSubStruct", "Mvir", "npart",
             "Xc", "Yc", "Zc", "VXc", "VYc", "VZc", "Rvir"]
# Names that are kept if present.
_OPTIONAL = ["Vmax", "com_offset", "mbp_offset"]

# AHF particle types follow the GADGET convention.
AHF_PTYPE_GAS = 0
AHF_PTYPE_DM = 1
AHF_PTYPE_STAR = 4


def _read_halo_columns(filename: str) -> Dict[str, int]:
    """Return {column name: zero-based column index} from the AHF_halos header."""
    with open(filename) as handle:
        header = handle.readline()
    columns = {name: int(number) - 1 for name, number in _HEADER_FIELD.findall(header)}
    missing = [name for name in _REQUIRED if name not in columns]
    if missing:
        raise ValueError(
            f"{filename}: the header does not name the columns {missing}, so the "
            f"layout of this AHF_halos file is not recognised. Header line was:\n{header}"
        )
    return columns


def read_ahf(filename: str) -> Tuple[Dict[str, Any], Dict[str, np.ndarray], Dict[str, np.ndarray]]:
    """Read an AHF halo catalogue.

    Parameters
    ----------
    filename : str
        Path to the `.AHF_halos` file.

    Returns
    -------
    metadata : dict
        {"Nhalos": ..., "Nsubhalos": ...}. The file does not give HubbleParam or the
        redshift, so unit conversion in HaloTools cannot use them for AHF.
    halos : dict[str, np.ndarray]
        Every row of the file (hosts and subhalos), with keys HaloID, HostHaloID,
        NumSubStruct, Mass (Mvir), NumPart, Pos (n, 3), Vel (n, 3), Rvir, and
        Vmax, com_offset, mbp_offset where the header has them.
    subhalos : dict[str, np.ndarray]
        The rows of `halos` that sit inside another halo, with keys HaloID,
        ParentHaloID, SubhaloMass, NumPart, Pos, Vel, Rvir (and Vmax if present).
        AHF lists subhalos in the AHF_halos file itself, so they are taken from there.
        The .AHF_substructure file is not needed for this; see read_ahf_substructure
        if you want it.

    Notes
    -----
    In the hostHalo column a halo that is not inside another halo has 0 (some AHF
    versions write -1). Both are treated as "not a subhalo".
    """
    if not filename.endswith(".AHF_halos"):
        raise ValueError(f"read_ahf() expects a '.AHF_halos' file, got '{filename}'.")

    columns = _read_halo_columns(filename)
    # ndmin=2 keeps a one-row file two-dimensional.
    data = np.genfromtxt(filename, comments="#", ndmin=2)

    def col(name):
        return data[:, columns[name]]

    halos: Dict[str, np.ndarray] = {
        "HaloID": col("ID").astype(np.int64),
        "HostHaloID": col("hostHalo").astype(np.int64),
        "NumSubStruct": col("numSubStruct").astype(np.int64),
        "Mass": col("Mvir").astype(np.float64),
        "NumPart": col("npart").astype(np.int64),
        "Pos": np.column_stack([col("Xc"), col("Yc"), col("Zc")]).astype(np.float64),
        "Vel": np.column_stack([col("VXc"), col("VYc"), col("VZc")]).astype(np.float64),
        "Rvir": col("Rvir").astype(np.float64),
    }
    for name in _OPTIONAL:
        if name in columns:
            halos[name] = col(name).astype(np.float64)

    is_sub = halos["HostHaloID"] > 0
    subhalos: Dict[str, np.ndarray] = {}
    if np.any(is_sub):
        subhalos["HaloID"] = halos["HaloID"][is_sub]
        subhalos["ParentHaloID"] = halos["HostHaloID"][is_sub]
        subhalos["SubhaloMass"] = halos["Mass"][is_sub]
        subhalos["NumPart"] = halos["NumPart"][is_sub]
        subhalos["Pos"] = halos["Pos"][is_sub]
        subhalos["Vel"] = halos["Vel"][is_sub]
        subhalos["Rvir"] = halos["Rvir"][is_sub]
        if "Vmax" in halos:
            subhalos["Vmax"] = halos["Vmax"][is_sub]

    metadata: Dict[str, Any] = {"Nhalos": len(data), "Nsubhalos": int(np.count_nonzero(is_sub))}
    return metadata, halos, subhalos


def read_ahf_substructure(filename: str) -> Dict[int, np.ndarray]:
    """Read an AHF_substructure file.

    The file holds, for each host halo, a line with two numbers (host ID and the
    number of subhalos) followed by one line listing those subhalo IDs. Lines are
    therefore of very different lengths, so the file cannot be read as a table.

    Returns {host ID: int64 array of subhalo IDs}. A malformed file raises a
    ValueError that names the line. A host with a count of zero is assumed to have
    no list line (blank lines are ignored in any case).
    """
    lines = []
    with open(filename) as handle:
        for number, text in enumerate(handle, start=1):
            if text.strip() and not text.lstrip().startswith("#"):
                lines.append((number, text))

    result: Dict[int, np.ndarray] = {}
    i = 0
    while i < len(lines):
        number, text = lines[i]
        tokens = text.split()
        if len(tokens) != 2:
            raise ValueError(
                f"{filename}, line {number}: expected 'host ID, number of subhalos' "
                f"(2 columns), found {len(tokens)} columns."
            )
        host, count = int(tokens[0]), int(tokens[1])
        if count == 0:
            result[host] = np.array([], dtype=np.int64)
            i += 1
            continue
        if i + 1 >= len(lines):
            raise ValueError(f"{filename}: file ends after the header line of host {host}.")
        list_number, list_text = lines[i + 1]
        ids = list_text.split()
        if len(ids) != count:
            raise ValueError(
                f"{filename}, line {list_number}: host {host} declares {count} subhalos "
                f"but the line lists {len(ids)} IDs."
            )
        result[host] = np.array(ids, dtype=np.int64)
        i += 2
    return result


def read_ahf_particles(filename: str,
                       wanted: Optional[Iterable[int]] = None) -> Dict[int, Dict[str, Any]]:
    """Read blocks of an AHF_particles file.

    File layout: a first line with the number of halos, then for each halo a header
    line '<npart> <haloID> ...' followed by npart lines '<particle ID> <particle type>'.
    A file whose particle lines have no type column is refused, because the type is
    needed to tell stars from dark matter.

    Parameters
    ----------
    filename : str
    wanted : iterable of int, optional
        Halo IDs to read. Other blocks are skipped without being parsed, and reading
        stops when every wanted halo has been found. If None, all blocks are read
        (this can use a lot of memory for a large file).

    Returns
    -------
    dict {halo ID: {"n_declared": int, "ids": int64 array, "types": int64 array}}
    The ids and types arrays are in file order. n_declared is the block size in the
    header line, which should equal NumPart in AHF_halos.
    """
    remaining = None if wanted is None else {int(w) for w in wanted}
    result: Dict[int, Dict[str, Any]] = {}

    with open(filename) as handle:
        handle.readline()                              # number of halos; not needed
        while remaining is None or remaining:
            header = handle.readline()
            if not header.strip():
                break                                  # end of file
            tokens = header.split()
            npart, halo_id = int(tokens[0]), int(tokens[1])

            if remaining is not None and halo_id not in remaining:
                for _ in islice(handle, npart):        # skip the block
                    pass
                continue

            lines = list(islice(handle, npart))
            if len(lines) != npart:
                raise ValueError(f"{filename}: the block for halo {halo_id} is truncated.")
            if npart == 0:
                ids = np.array([], dtype=np.int64)
                types = np.array([], dtype=np.int64)
            else:
                n_columns = len(lines[0].split())
                if n_columns < 2:
                    raise ValueError(
                        f"{filename}: particle lines have {n_columns} column(s). A particle "
                        "type column is needed (AHF built without it?)."
                    )
                table = np.loadtxt(lines, dtype=np.int64, usecols=(0, 1), ndmin=2)
                ids, types = table[:, 0], table[:, 1]
            result[halo_id] = {"n_declared": npart, "ids": ids, "types": types}
            if remaining is not None:
                remaining.discard(halo_id)

    if remaining:
        raise ValueError(f"{filename}: halo IDs not found: {sorted(remaining)[:5]}")
    return result
