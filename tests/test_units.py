import logging

import numpy as np
import pytest

from phytyl_tail_builder.mol2 import build_order, read_mol2
from phytyl_tail_builder.optimizer import CUTOFF, build_kdtree, count_clashes_for_segment, element_of
from phytyl_tail_builder.parser import leap_sessions, parse_leap_log_for_atoms, parse_tleap_in
from phytyl_tail_builder.pdb import (
    extract_all_coordinates,
    map_residues,
    read_residues,
    remove_hydrogens_from_rebuilt,
    write_coordinates,
)

MOL2 = """@<TRIPOS>MOLECULE
TST
@<TRIPOS>ATOM
  1 C1     0.0   0.0   0.0 C.3   1 TST   0.0
  2 C2     1.5   0.0   0.0 C.3   1 TST   0.0
  3 C3     3.0   0.0   0.0 C.3   1 TST   0.0
  4 O1    -1.5   0.0   0.0 O.3   1 TST   0.0
  5 H1     0.0   1.0   0.0 H     1 TST   0.0
  6 C4     3.0   1.5   0.0 C.3   1 TST   0.0
@<TRIPOS>BOND
  1  1  2  1
  2  3  2  1
  3  4  1  1
  4  5  2  1
  5  3  6  1
"""


def pdb_line(serial, name, resname, resnum, xyz, chain=" "):
    return (f"ATOM  {serial:5d} {name:<4} {resname:3} {chain}{resnum % 10000:4d}    "
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00\n")


def test_parse_tleap_in(tmp_path):
    f = tmp_path / "tleap.in"
    f.write_text("X = loadpdb in.pdb\nsavepdb X out.pdb\nsaveamberparm X a.prmtop a.inpcrd\n")
    assert parse_tleap_in(f) == ("in.pdb", "out.pdb", "a.prmtop", "a.inpcrd")


def test_parse_tleap_in_ignores_comments_and_case(tmp_path):
    f = tmp_path / "tleap.in"
    f.write_text("# X = loadpdb old.pdb\nX = loadPdb in.pdb  # the input\n"
                 "#savepdb X old_out.pdb\nsavePdb X out.pdb\n"
                 "#saveamberparm X old.prmtop old.inpcrd\nsaveAmberParm X a.prmtop a.inpcrd\n"
                 "# a comment that mentions savepdb\n")
    assert parse_tleap_in(f) == ("in.pdb", "out.pdb", "a.prmtop", "a.inpcrd")
    f.write_text("X = loadpdb in.pdb\nsavepdb X out.pdb\n#saveamberparm X a.prmtop a.inpcrd\n")
    assert parse_tleap_in(f) == ("in.pdb", "out.pdb", None, None)


SESSION = """log started: {date}

Log file: ./leap.log
>> logFile leap.log
log started: {date}

Log file: ./leap.log
> X = loadpdb in.pdb
{atoms}
> quit
\tQuit

Exiting LEaP: Errors = 0; Warnings = 0; Notes = 0.
"""


def test_leap_log_sessions_ignore_logfile_restarts():
    atoms = "  Added missing heavy atom: .R<CHL 7>.A<C1 47>"
    assert len(leap_sessions(SESSION.format(date="Mon", atoms=atoms))) == 1


def test_leap_log_atoms_and_names_come_from_last_session(tmp_path, caplog):
    f = tmp_path / "leap.log"
    f.write_text(SESSION.format(date="Mon", atoms="  Added missing heavy atom: .R<CLA 5>.A<C1 46>")
                 + SESSION.format(date="Tue", atoms="  Added missing heavy atom: .R<CHL 7>.A<C1 47>\n"
                                                     "  Added missing heavy atom: .R<CHL 7>.A<C2 48>"))
    with caplog.at_level(logging.WARNING):
        assert parse_leap_log_for_atoms(f) == [(7, "CHL", "C1"), (7, "CHL", "C2")]
    assert "2 tleap sessions" in caplog.text


def test_build_order_ignores_mol2_bond_direction(tmp_path):
    f = tmp_path / "TST.mol2"
    f.write_text(MOL2)
    coords, bonds = read_mol2(f)
    assert set(coords) == {"C1", "C2", "C3", "O1", "C4"}
    assert ("C3", "C2") in bonds  # listed child-first
    assert build_order(bonds, ["C2", "C3", "C4"], ["C1", "O1"]) == [("C2", "C1"), ("C3", "C2"), ("C4", "C3")]
    # The same atoms grown from the other end.
    assert build_order(bonds, ["C1", "C2", "O1"], ["C3", "C4"]) == [("C2", "C3"), ("C1", "C2"), ("O1", "C1")]
    with pytest.raises(ValueError):
        build_order(bonds + [("X1", "X2")], ["X2"], ["C1"])


def _system(tmp_path, post_extra=""):
    pre, post = tmp_path / "pre.pdb", tmp_path / "post.pdb"
    pre.write_text(pdb_line(1, "N", "ALA", 5, (9, 9, 9)) + "TER\n"
                   + pdb_line(2, "MG", "CLA", 40, (1, 2, 3)) + pdb_line(3, "MG", "CHL", 41, (4, 5, 6)))
    post.write_text(pdb_line(1, "N", "ALA", 1, (9, 9, 9)) + post_extra + "TER\n"
                    + pdb_line(2, "MG", "CLA", 2, (1, 2, 3)) + pdb_line(3, "C1", "CLA", 2, (0, 0, 0))
                    + pdb_line(4, "MG", "CHL", 3, (4, 5, 6)) + pdb_line(5, "C1", "CHL", 3, (7, 7, 7)))
    return read_residues(pre), read_residues(post)


def test_residues_mapped_by_sequential_index(tmp_path):
    pre, post = _system(tmp_path)
    # tleap numbers residues from the first resSeq of the loadpdb file (5 here).
    assert map_residues(pre, post, {6: "CLA", 7: "CHL"}) == {6: 1, 7: 2}


def test_unmatched_residue_raises(tmp_path):
    pre, post = _system(tmp_path)
    with pytest.raises(ValueError, match="loadpdb has CLA 40"):
        map_residues(pre, post, {6: "CHL"})
    with pytest.raises(ValueError, match="outside"):
        map_residues(pre, post, {9: "CHL"})
    # An extra residue in the savepdb file shifts everything after it: caught, not zipped.
    pre, post = _system(tmp_path, post_extra=pdb_line(9, "N", "GLY", 2, (8, 8, 8)))
    with pytest.raises(ValueError):
        map_residues(pre, post, {6: "CLA", 7: "CHL"})


def test_moved_atom_raises(tmp_path):
    pre, post = _system(tmp_path)
    post[2]['atoms']['MG'] = np.array([4.0, 5.0, 6.5])
    with pytest.raises(ValueError, match="not at the same coordinates"):
        map_residues(pre, post, {7: "CHL"})


def test_residues_split_on_records_not_resseq(tmp_path):
    f = tmp_path / "big.pdb"
    lines = [pdb_line(i + 1, "C1", "CLA" if i % 2 else "CHL", 9998 + i, (i, 0, 0)) for i in range(4)]
    lines.insert(2, pdb_line(9, "C2", "CHL", 10001, (9, 9, 9)))  # same resSeq 1 twice, wrapped
    f.write_text("".join(lines) + pdb_line(5, "C1", "CHL", 9998, (5, 0, 0)))
    residues = read_residues(f)
    assert [(r['name'], r['resseq']) for r in residues] == \
        [("CHL", "9998"), ("CLA", "9999"), ("CHL", "1"), ("CHL", "0"), ("CLA", "1"), ("CHL", "9998")]


def test_five_digit_resseq_of_savepdb(tmp_path):
    f = tmp_path / "big.pdb"
    # record, serial, name, resName, chain, five-digit resSeq (columns 23-27), coordinates from column 31
    line = "ATOM  " + "    1" + "  C1  " + "CLA" + "  " + "10234" + "   " + "   1.000   2.000   3.000\n"
    f.write_text(line + line.replace("10234", "10235"))
    assert [r['resseq'] for r in read_residues(f)] == ["10234", "10235"]


def test_remove_hydrogens_only_in_listed_residues(tmp_path):
    src, dst = tmp_path / "in.pdb", tmp_path / "out.pdb"
    src.write_text(pdb_line(1, "C1", "CLA", 1, (0, 0, 0)) + pdb_line(2, "H11", "CLA", 1, (0, 0, 1))
                   + pdb_line(3, "1H4", "CLA", 1, (0, 0, 2)) + pdb_line(4, "HA", "ALA", 2, (0, 0, 3))
                   + "TER\nEND\n")
    remove_hydrogens_from_rebuilt(str(src), str(dst), {0})
    names = [line[12:16].strip() for line in dst.read_text().splitlines() if line.startswith("ATOM")]
    assert names == ["C1", "HA"]
    assert dst.read_text().endswith("TER\nEND\n")


def test_write_coordinates_rounds_to_pdb_precision(tmp_path):
    f = tmp_path / "x.pdb"
    f.write_text(pdb_line(1, "C1", "CLA", 3, (0, 0, 0)) + pdb_line(2, "C2", "CLA", 3, (1, 1, 1)))
    write_coordinates(str(f), {(0, "C2"): np.array([1.23456, -2.0, 30.0004])})
    coords = extract_all_coordinates(str(f))
    assert np.array_equal(coords[0]["C2"], [1.235, -2.0, 30.0])
    assert np.array_equal(coords[0]["C1"], [0, 0, 0])


def test_clash_score_grows_as_atoms_approach():
    coords = {1: {"A": np.array([0.0, 0.0, 0.0])}}
    tree, labels, _ = build_kdtree(coords)
    near = count_clashes_for_segment({"X": np.array([2.0, 0, 0])}, tree, labels)
    far = count_clashes_for_segment({"X": np.array([3.3, 0, 0])}, tree, labels)
    alone = count_clashes_for_segment({"X": np.array([20.0, 0, 0])}, tree, labels)
    assert near > far > alone == 0.0
    # each pair is shifted to 0 at the cutoff, so crossing it changes nothing
    edge = [count_clashes_for_segment({"X": np.array([CUTOFF + d, 0, 0])}, tree, labels)
            for d in (-1e-6, 1e-6)]
    assert abs(edge[0] - edge[1]) < 1e-6


def test_attractive_trials_are_not_floored():
    # v0.1 clamped every score to >= 0.1, so all attractive trials tied and the first (+z) won.
    tree, labels, _ = build_kdtree({1: {"H1": np.array([0.0, 0.0, 0.0])}})  # C-H sigma 3.03 Å
    a = count_clashes_for_segment({"C2": np.array([3.2, 0, 0])}, tree, labels)
    b = count_clashes_for_segment({"C2": np.array([3.35, 0, 0])}, tree, labels)
    assert a < 0 and b < 0 and a != b


def test_bonded_partners_are_not_scored():
    """The parent (1-2) and the atoms on the parent (1-3) sit where bond length and angle put them."""
    tree, labels, _ = build_kdtree({1: {"C1": np.array([0.0, 0, 0]), "C2": np.array([1.5, 0, 0])}})
    trial = {"C3": np.array([2.3, 1.2, 0])}
    assert count_clashes_for_segment(trial, tree, labels, bonded={(1, "C1"), (1, "C2")}) == 0.0
    assert count_clashes_for_segment(trial, tree, labels) > 0


def test_element_radii():
    assert [element_of(n) for n in ("C1", "1HMA", "HMA1", "OBD", "NA", "MG", "CMA")] == \
        ["C", "H", "H", "O", "N", "MG", "C"]
    trial = {"C2": np.array([2.6, 0.0, 0.0])}
    scores = {}
    for name in ("CA", "HA", "OB", "MG"):
        tree, labels, _ = build_kdtree({1: {name: np.array([0.0, 0.0, 0.0])}})
        scores[name] = count_clashes_for_segment(trial, tree, labels)
    assert scores["CA"] > scores["OB"] > scores["HA"] > scores["MG"]


def test_fragment_neighbours_use_lennard_jones():
    points = np.array([[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]])
    labels = [(1, "A"), (1, "B")]
    tree, _, _ = build_kdtree({1: {"A": points[0], "B": points[1]}})
    trial = {"X": np.array([3.0, 0.0, 0.0])}
    mie_score = count_clashes_for_segment(trial, tree, labels)
    lj_score = count_clashes_for_segment(trial, tree, labels, fragment={(1, "A")})
    assert mie_score != lj_score
    assert count_clashes_for_segment(trial, tree, labels, fragment={(2, "A")}) == mie_score
    # An atom given in `extra` scores like a fragment atom in the tree.
    empty_tree, empty_labels, _ = build_kdtree({1: {"B": points[1]}})
    assert count_clashes_for_segment(trial, empty_tree, empty_labels,
                                     extra={(1, "A"): points[0]}) == pytest.approx(lj_score)
