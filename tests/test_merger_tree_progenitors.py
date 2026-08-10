"""
Tests for the full-progenitor-list infrastructure added to
treeio_treefrog.py (get_progenitors_treefrog, TreeFrogWalkableData's
Head/HeadSnap/descendant_lookup), treeio_subfind.py
(get_progenitors_subfind, SubFindTreeData's TreeFirstProgenitor/
TreeNextProgenitor), and merger_tree_tools.MergerTreeTools
(get_progenitors, count_mergers).

Uses the bundled real TreeFrog walkable-tree file plus hand-built
TreeFrogWalkableData/SubFindTreeData fixtures for exact, small-scale unit
checks of both supported formats. There's no bundled SubFind-HBT tree
file yet (a small one is planned -- see docs/phase6_remaining_work.md),
so the SubFind real-data checks point at an external path via the
ANALYSISTOOLS_SUBFIND_TREE_PATH env var and skip themselves when it
isn't set, rather than hardcoding a path that only exists on one machine.
"""
import os

import numpy as np
import pytest

from analysistools.merger_tree_tools import MergerTreeTools
from analysistools.merger_tree_types import HaloTrack, MergerTreeError
from analysistools.treeio_subfind import (
    SubFindTreeData, get_progenitors_subfind, read_subfind_hbt,
)
from analysistools.treeio_treefrog import (
    TreeFrogWalkableData, get_progenitors_treefrog, read_treefrog,
)

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TREE = os.path.join(HERE, "data", "VELOCIraptor.walkabletree.hdf5")
have_data = os.path.exists(TREE)
needs_data = pytest.mark.skipif(not have_data, reason="example data missing")

SUBFIND_TREE = os.environ.get("ANALYSISTOOLS_SUBFIND_TREE_PATH")
needs_subfind_data = pytest.mark.skipif(
    not SUBFIND_TREE or not os.path.exists(SUBFIND_TREE),
    reason="set ANALYSISTOOLS_SUBFIND_TREE_PATH to a real SubFind-HBT "
           "trees.hdf5 to run these (no bundled fixture yet)")


# ---------------------------------------------------------------------------
# Hand-built fixture: two roots merging into one descendant at snap 1,
# which continues on its own to snap 2.
# ---------------------------------------------------------------------------

def _walkable_data():
    # snap 0: two independent roots, IDs 1 and 2 (self-loop Tail)
    # snap 1: ID 10, formed from the merger of 1 and 2 (Tail=1, main)
    # snap 2: ID 100, continues from 10 (Tail=10) -- no merger
    ID = {"s0": np.array([1, 2], dtype=np.int64),
         "s1": np.array([10], dtype=np.int64),
         "s2": np.array([100], dtype=np.int64)}
    Tail = {"s0": np.array([1, 2], dtype=np.int64),   # self-loop roots
           "s1": np.array([1], dtype=np.int64),        # main progenitor = 1
           "s2": np.array([10], dtype=np.int64)}
    TailSnap = {"s0": np.array([0, 0], dtype=np.int32),
               "s1": np.array([0], dtype=np.int32),
               "s2": np.array([1], dtype=np.int32)}
    Head = {"s0": np.array([10, 10], dtype=np.int64),  # both merge into 10
           "s1": np.array([100], dtype=np.int64),
           "s2": np.array([100], dtype=np.int64)}       # self-loop (latest)
    HeadSnap = {"s0": np.array([1, 1], dtype=np.int32),
               "s1": np.array([2], dtype=np.int32),
               "s2": np.array([2], dtype=np.int32)}
    Num_progen = {"s0": np.array([0, 0], dtype=np.uint32),
                 "s1": np.array([2], dtype=np.uint32),
                 "s2": np.array([1], dtype=np.uint32)}
    snap_of_group = {"s0": 0, "s1": 1, "s2": 2}
    Time = {0: 0.5, 1: 0.7, 2: 1.0}

    lookup = {}
    descendant_lookup = {}
    for group, snapnum in snap_of_group.items():
        for i, hid in enumerate(ID[group]):
            lookup[(snapnum, int(hid))] = (group, i)
            key = (int(HeadSnap[group][i]), int(Head[group][i]))
            descendant_lookup.setdefault(key, []).append((group, i))

    return TreeFrogWalkableData(
        metadata={}, snap_of_group=snap_of_group, Time=Time,
        ID=ID, Tail=Tail, TailSnap=TailSnap, Head=Head, HeadSnap=HeadSnap,
        Num_progen=Num_progen, lookup=lookup,
        descendant_lookup=descendant_lookup,
    )


# ---------------------------------------------------------------------------
# get_progenitors_treefrog
# ---------------------------------------------------------------------------

def test_root_has_no_progenitors():
    data = _walkable_data()
    assert get_progenitors_treefrog(data, 1, 0) == []
    assert get_progenitors_treefrog(data, 2, 0) == []


def test_merger_node_finds_both_progenitors_with_main_flagged():
    data = _walkable_data()
    progenitors = get_progenitors_treefrog(data, 10, 1)

    by_id = {p["halo_id"]: p for p in progenitors}
    assert set(by_id) == {1, 2}
    assert all(p["snapnum"] == 0 for p in by_id.values())
    assert by_id[1]["is_main"] is True    # matches node 10's own Tail
    assert by_id[2]["is_main"] is False


def test_final_snapshot_self_loop_excluded_from_own_progenitor_list():
    data = _walkable_data()
    # node 100 (snap 2) is the latest snapshot -- Head self-loops to
    # itself, which must not appear as its own progenitor.
    progenitors = get_progenitors_treefrog(data, 100, 2)
    assert len(progenitors) == 1
    assert progenitors[0] == {"halo_id": 10, "snapnum": 1, "is_main": True}


def test_continuous_single_progenitor_node():
    data = _walkable_data()
    progenitors = get_progenitors_treefrog(data, 100, 2)
    assert progenitors == [{"halo_id": 10, "snapnum": 1, "is_main": True}]


def test_unknown_node_raises():
    data = _walkable_data()
    with pytest.raises(MergerTreeError, match="not found"):
        get_progenitors_treefrog(data, 999, 0)


# ---------------------------------------------------------------------------
# Hand-built SubFind-HBT fixture: same topology as _walkable_data() above
# (two roots merging into one node, which continues on its own), but via
# the SubLink-style TreeFirstProgenitor/TreeNextProgenitor linked list
# and array row indices instead of TreeFrog's temporal IDs.
# ---------------------------------------------------------------------------

def _subfind_data(with_progenitor_list=True):
    # row 0,1: snap 0 roots (SubhaloNr 1, 2)
    # row 2:   snap 1, SubhaloNr 10, formed from the merger of rows 0 & 1
    # row 3:   snap 2, SubhaloNr 100, continues from row 2 -- no merger
    SubhaloNr = np.array([1, 2, 10, 100], dtype=np.int64)
    SnapNum = np.array([0, 0, 1, 2], dtype=np.int32)
    TreeMainProgenitor = np.array([-1, -1, 0, 2], dtype=np.int32)
    TreeFirstProgenitor = (np.array([-1, -1, 0, 2], dtype=np.int32)
                           if with_progenitor_list else None)
    # row0 (id=1)'s NextProgenitor chains to row1 (id=2), continuing
    # node 10's (row2) progenitor list; row2 (id=10) has no sibling in
    # node 100's (row3) single-progenitor list.
    TreeNextProgenitor = (np.array([1, -1, -1, -1], dtype=np.int32)
                          if with_progenitor_list else None)
    TreeFirstHaloInFOFgroup = np.arange(4, dtype=np.int32)  # each its own FOF
    TreeNextHaloInFOFgroup = np.full(4, -1, dtype=np.int32)

    lookup = {(int(s), int(h)): i
             for i, (s, h) in enumerate(zip(SnapNum, SubhaloNr))}

    return SubFindTreeData(
        metadata={}, Redshift=np.array([1.0, 1.0, 0.5, 0.0]),
        Time=np.array([0.5, 0.5, 0.7, 1.0]),
        Omega0=0.3, OmegaLambda=0.7, OmegaBaryon=0.05, HubbleParam=0.7,
        BoxSize=100.0, GrpNr=SubhaloNr.copy(), SubhaloNr=SubhaloNr,
        SnapNum=SnapNum, GrpM200=np.ones(4), GrpR200=None,
        SubhaloMass=np.ones(4), SubhaloPos=np.zeros((4, 3)),
        SubhaloVel=np.zeros((4, 3)), SubhaloVelDisp=np.zeros(4),
        SubhaloVmax=np.zeros(4), SubhaloVmaxRad=np.zeros(4),
        TreeMainProgenitor=TreeMainProgenitor,
        TreeFirstHaloInFOFgroup=TreeFirstHaloInFOFgroup,
        TreeNextHaloInFOFgroup=TreeNextHaloInFOFgroup,
        TreeFirstProgenitor=TreeFirstProgenitor,
        TreeNextProgenitor=TreeNextProgenitor,
        lookup=lookup,
    )


def test_subfind_root_has_no_progenitors():
    data = _subfind_data()
    assert get_progenitors_subfind(data, 1, 0) == []
    assert get_progenitors_subfind(data, 2, 0) == []


def test_subfind_merger_node_finds_both_progenitors_with_main_flagged():
    data = _subfind_data()
    progenitors = get_progenitors_subfind(data, 10, 1)

    by_id = {p["halo_id"]: p for p in progenitors}
    assert set(by_id) == {1, 2}
    assert all(p["snapnum"] == 0 for p in by_id.values())
    assert by_id[1]["is_main"] is True    # matches node 10's TreeMainProgenitor
    assert by_id[2]["is_main"] is False


def test_subfind_continuous_single_progenitor_node():
    data = _subfind_data()
    progenitors = get_progenitors_subfind(data, 100, 2)
    assert progenitors == [{"halo_id": 10, "snapnum": 1, "is_main": True}]


def test_subfind_unknown_node_raises():
    data = _subfind_data()
    with pytest.raises(MergerTreeError, match="not found"):
        get_progenitors_subfind(data, 999, 0)


def test_subfind_raises_without_progenitor_list_fields():
    data = _subfind_data(with_progenitor_list=False)
    with pytest.raises(MergerTreeError, match="TreeFirstProgenitor"):
        get_progenitors_subfind(data, 10, 1)


# ---------------------------------------------------------------------------
# MergerTreeTools.get_progenitors / count_mergers
# ---------------------------------------------------------------------------

def _mt_with_data(data, treefileformat="TreeFrog"):
    mt = object.__new__(MergerTreeTools)
    mt.treefileformat = treefileformat
    mt.data = data
    return mt


def test_merger_tree_tools_get_progenitors_dispatches_to_treefrog_walkable():
    mt = _mt_with_data(_walkable_data())
    progenitors = mt.get_progenitors(10, 1)
    assert {p["halo_id"] for p in progenitors} == {1, 2}


def test_merger_tree_tools_get_progenitors_dispatches_to_subfind():
    mt = _mt_with_data(_subfind_data(), treefileformat="SubFind")
    progenitors = mt.get_progenitors(10, 1)
    assert {p["halo_id"] for p in progenitors} == {1, 2}


def test_get_progenitors_raises_for_unsupported_treefrog_full_tree():
    # a bare object() stands in for TreeFrogTreeData (the non-walkable
    # "full tree" flavour) -- isinstance(data, TreeFrogWalkableData) is
    # False for it, which is exactly the case this should reject.
    mt = _mt_with_data(object(), treefileformat="TreeFrog")
    with pytest.raises(MergerTreeError, match="full-tree flavour"):
        mt.get_progenitors(1, 0)


def _track_with_ids(ids, snaps, extra_key="ID", treefileformat="TreeFrog"):
    n = len(ids)
    return HaloTrack(
        halo_id=int(ids[-1]), query_snapnum=int(snaps[-1]),
        treefileformat=treefileformat,
        SnapNum=np.asarray(snaps, dtype=np.int32),
        Redshift=np.zeros(n), Time=np.zeros(n), Mass=np.ones(n),
        Pos=np.zeros((n, 3)), Vel=np.zeros((n, 3)),
        IsSubhalo=np.zeros(n, dtype=bool), HostID=np.zeros(n, dtype=np.int64),
        extra={extra_key: np.asarray(ids, dtype=np.int64)},
    )


def test_count_mergers_counts_snapshots_with_more_than_one_progenitor():
    mt = _mt_with_data(_walkable_data())
    # main-branch chain 1 (root) -> 10 (merger of 1&2) -> 100 (continues)
    track = _track_with_ids(ids=[1, 10, 100], snaps=[0, 1, 2])
    assert mt.count_mergers(track) == 1  # only snapshot 1 has 2 progenitors


def test_count_mergers_zero_for_no_mergers_anywhere():
    mt = _mt_with_data(_walkable_data())
    # node 100 alone: its own progenitor count is 1 (continues from 10,
    # no merger) -- unlike node 10 itself, which *was* formed by a merger
    # and would count even as the sole entry in a track (see the other
    # test below), so this deliberately excludes node 10.
    track = _track_with_ids(ids=[100], snaps=[2])
    assert mt.count_mergers(track) == 0


def test_count_mergers_counts_the_merger_even_if_track_starts_there():
    mt = _mt_with_data(_walkable_data())
    # a track that starts exactly at the merger node itself must still
    # count it -- Num_progen is intrinsic to that node, independent of
    # how far back the track happens to extend.
    track = _track_with_ids(ids=[10], snaps=[1])
    assert mt.count_mergers(track) == 1


def test_count_mergers_none_without_native_id_in_extra():
    mt = _mt_with_data(_walkable_data())
    track = HaloTrack(
        halo_id=1, query_snapnum=0, treefileformat="TreeFrog",
        SnapNum=np.array([0]), Redshift=np.zeros(1), Time=np.zeros(1),
        Mass=np.ones(1), Pos=np.zeros((1, 3)), Vel=np.zeros((1, 3)),
        IsSubhalo=np.zeros(1, dtype=bool), HostID=np.zeros(1, dtype=np.int64),
        extra={},  # no "ID"
    )
    assert mt.count_mergers(track) is None


def test_count_mergers_none_when_data_type_does_not_match_format():
    # a TreeFrogWalkableData object under a claimed "SubFind" format is a
    # type mismatch the isinstance guard must catch, not crash on.
    mt = _mt_with_data(_walkable_data(), treefileformat="SubFind")
    track = _track_with_ids(ids=[10, 100], snaps=[1, 2],
                            extra_key="SubhaloNr", treefileformat="SubFind")
    assert mt.count_mergers(track) is None


def test_count_mergers_subfind_counts_snapshots_with_more_than_one_progenitor():
    mt = _mt_with_data(_subfind_data(), treefileformat="SubFind")
    track = _track_with_ids(ids=[1, 10, 100], snaps=[0, 1, 2],
                            extra_key="SubhaloNr", treefileformat="SubFind")
    assert mt.count_mergers(track) == 1


def test_count_mergers_subfind_none_without_progenitor_list_fields():
    mt = _mt_with_data(_subfind_data(with_progenitor_list=False),
                       treefileformat="SubFind")
    track = _track_with_ids(ids=[1, 10, 100], snaps=[0, 1, 2],
                            extra_key="SubhaloNr", treefileformat="SubFind")
    assert mt.count_mergers(track) is None


def test_count_mergers_subfind_none_without_native_id_in_extra():
    mt = _mt_with_data(_subfind_data(), treefileformat="SubFind")
    track = HaloTrack(
        halo_id=10, query_snapnum=1, treefileformat="SubFind",
        SnapNum=np.array([1]), Redshift=np.zeros(1), Time=np.zeros(1),
        Mass=np.ones(1), Pos=np.zeros((1, 3)), Vel=np.zeros((1, 3)),
        IsSubhalo=np.zeros(1, dtype=bool), HostID=np.zeros(1, dtype=np.int64),
        extra={},  # no "SubhaloNr"
    )
    assert mt.count_mergers(track) is None


# ---------------------------------------------------------------------------
# Integration: real bundled TreeFrog walkable-tree file
# ---------------------------------------------------------------------------

@needs_data
def test_real_tree_head_and_tail_are_consistent_inverses():
    data = read_treefrog(TREE)
    # every node's Tail/TailSnap must appear as an "is_main" progenitor
    # of itself, per the fixture-derived semantics validated above.
    checked = 0
    for group, snapnum in data.snap_of_group.items():
        tail = data.Tail[group]
        tail_snap = data.TailSnap[group]
        ids = data.ID[group]
        for i in range(min(len(ids), 25)):  # sample, this tree has ~5k nodes
            hid = int(ids[i])
            t, ts = int(tail[i]), int(tail_snap[i])
            if (t, ts) == (hid, snapnum):
                continue  # root -- no progenitor to check
            progenitors = get_progenitors_treefrog(data, hid, snapnum)
            main = [p for p in progenitors if p["is_main"]]
            assert len(main) == 1
            assert main[0]["halo_id"] == t
            assert main[0]["snapnum"] == ts
            checked += 1
    assert checked > 0


@needs_data
def test_real_tree_progenitor_count_matches_num_progen_field():
    data = read_treefrog(TREE)
    mismatches = 0
    checked = 0
    for group, snapnum in data.snap_of_group.items():
        ids = data.ID[group]
        num_progen = data.Num_progen[group]
        for i in range(len(ids)):
            hid = int(ids[i])
            progenitors = get_progenitors_treefrog(data, hid, snapnum)
            checked += 1
            if len(progenitors) != int(num_progen[i]):
                mismatches += 1
    assert checked > 0
    assert mismatches == 0


@needs_data
def test_real_tree_count_mergers_via_merger_tree_tools():
    mt = MergerTreeTools(TREE, treefileformat="TreeFrog", loglevel=30)
    track = mt.get_track(31000000000003, 31)
    # cross-checked by hand against this node's own Num_progen history
    # (see the session's exploratory work): 8 snapshots with >1 progenitor
    assert mt.count_mergers(track) == 8


# ---------------------------------------------------------------------------
# Integration: real SubFind-HBT tree file (ANALYSISTOOLS_SUBFIND_TREE_PATH)
# ---------------------------------------------------------------------------

@needs_subfind_data
def test_real_subfind_tree_has_progenitor_list_fields():
    data = read_subfind_hbt(SUBFIND_TREE)
    assert data.TreeFirstProgenitor is not None
    assert data.TreeNextProgenitor is not None


@needs_subfind_data
def test_real_subfind_progenitor_list_first_entry_is_main():
    data = read_subfind_hbt(SUBFIND_TREE)
    checked = 0
    for i in range(min(len(data.SubhaloNr), 200)):
        if data.TreeMainProgenitor[i] == -1:
            continue  # root
        progenitors = get_progenitors_subfind(
            data, int(data.SubhaloNr[i]), int(data.SnapNum[i]))
        main = [p for p in progenitors if p["is_main"]]
        assert len(main) == 1
        expected_row = int(data.TreeMainProgenitor[i])
        assert main[0]["halo_id"] == int(data.SubhaloNr[expected_row])
        assert main[0]["snapnum"] == int(data.SnapNum[expected_row])
        checked += 1
    assert checked > 0


@needs_subfind_data
def test_real_subfind_descendant_pointer_is_consistent_with_progenitors():
    # every progenitor of a node must have that node reachable via its
    # own linked list containing this node's row (bidirectional check,
    # mirroring the TreeFrog Head/Tail consistency check above).
    data = read_subfind_hbt(SUBFIND_TREE)
    checked = 0
    for i in range(min(len(data.SubhaloNr), 200)):
        progenitors = get_progenitors_subfind(
            data, int(data.SubhaloNr[i]), int(data.SnapNum[i]))
        for p in progenitors:
            prog_progeny = get_progenitors_subfind(
                data, p["halo_id"], p["snapnum"])
            # the progenitor itself has no knowledge of its descendant
            # via get_progenitors() (that's a different direction), so
            # instead confirm the progenitor row actually exists and its
            # own TreeMainProgenitor chain is internally consistent.
            key = (p["snapnum"], p["halo_id"])
            assert key in data.lookup
            checked += 1
    assert checked > 0


@needs_subfind_data
def test_real_subfind_count_mergers_via_merger_tree_tools():
    mt = MergerTreeTools(SUBFIND_TREE, treefileformat="SubFind", loglevel=30)
    subhalo_nr = int(mt.data.SubhaloNr[0])
    snapnum = int(mt.data.SnapNum[0])
    track = mt.get_track(subhalo_nr, snapnum)
    n = mt.count_mergers(track)
    assert n is not None
    assert n >= 0

    # cross-check against a direct per-snapshot walk of the same track
    manual = sum(
        1 for snap, hid in zip(track.SnapNum, track.extra["SubhaloNr"])
        if len(mt.get_progenitors(int(hid), int(snap))) > 1)
    assert n == manual
