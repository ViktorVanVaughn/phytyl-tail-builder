import numpy as np

from phytyl_tail_builder.main import apply_corrected_resnumbers
from phytyl_tail_builder.mol2 import parse_mol2_template
from phytyl_tail_builder.optimizer import build_kdtree, count_clashes_for_segment
from phytyl_tail_builder.parser import parse_leap_log, parse_leap_log_for_atoms, parse_tleap_in
from phytyl_tail_builder.pdb import (
    correct_residue_numbering,
    extract_all_coordinates,
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
@<TRIPOS>BOND
  1  1  2  1
  2  2  3  1
  3  4  1  1
  4  5  2  1
"""


def pdb_line(serial, name, resname, resnum, xyz):
    return (f"ATOM  {serial:5d} {name:<4} {resname:3} {resnum:5d}    "
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00\n")


def test_parse_tleap_in(tmp_path):
    f = tmp_path / "tleap.in"
    f.write_text("X = loadpdb in.pdb\nsavepdb X out.pdb\nsaveamberparm X a.prmtop a.inpcrd\n")
    assert parse_tleap_in(f) == ("in.pdb", "out.pdb", "a.prmtop", "a.inpcrd")


def test_leap_log_atoms_come_from_last_session(tmp_path):
    f = tmp_path / "leap.log"
    f.write_text(
        "Log file started: old\n"
        "  Added missing heavy atom: .R<CLA 5>.A<C1 46>\n"
        "Log file started: new\n"
        "  Added missing heavy atom: .R<CHL 7>.A<C1 47>\n"
        "  Added missing heavy atom: .R<CHL 7>.A<C2 48>\n")
    assert parse_leap_log_for_atoms(f) == [(7, "CHL", "C1"), (7, "CHL", "C2")]
    assert parse_leap_log(f) == {5: "CLA", 7: "CHL"}


def test_mol2_puts_anchor_bond_first(tmp_path):
    f = tmp_path / "TST.mol2"
    f.write_text(MOL2)
    order, bonds, atoms = parse_mol2_template(str(f), ["C2", "C3"])
    assert bonds == [("C1", "C2"), ("C1", "C2"), ("C2", "C3")]
    assert order == ["C1", "C2", "C3"]
    assert atoms == {"C1", "C2", "C3"}


def test_atoms_sorted_by_template_and_residues_by_size(tmp_path):
    (tmp_path / "TST.mol2").write_text(MOL2)
    extracted = [(1, "TST", "C3"), (1, "TST", "C2"), (2, "TST", "C3")]
    result = apply_corrected_resnumbers(extracted, tmp_path, {1: "TST", 2: "TST"}, {11: "TST", 12: "TST"})
    assert result == [(12, "TST", "C3"), (11, "TST", "C2"), (11, "TST", "C3")]


def test_residue_numbers_matched_by_coordinates(tmp_path):
    pre, post = tmp_path / "pre.pdb", tmp_path / "post.pdb"
    pre.write_text(pdb_line(1, "MG", "CLA", 40, (1, 2, 3)))
    post.write_text(pdb_line(1, "N", "ALA", 1, (9, 9, 9))
                    + pdb_line(2, "MG", "CLA", 2, (1, 2, 3))
                    + pdb_line(3, "C1", "CLA", 2, (4, 5, 6)))
    assert correct_residue_numbering(str(pre), str(post), {40: "CLA"}) == ({40: "CLA"}, {2: "CLA"})


def test_remove_hydrogens_only_in_listed_residues(tmp_path):
    src, dst = tmp_path / "in.pdb", tmp_path / "out.pdb"
    src.write_text(pdb_line(1, "C1", "CLA", 1, (0, 0, 0)) + pdb_line(2, "H11", "CLA", 1, (0, 0, 1))
                   + pdb_line(3, "1H4", "CLA", 1, (0, 0, 2)) + pdb_line(4, "HA", "ALA", 2, (0, 0, 3))
                   + "TER\nEND\n")
    remove_hydrogens_from_rebuilt(str(src), str(dst), {1})
    names = [line[12:16].strip() for line in dst.read_text().splitlines() if line.startswith("ATOM")]
    assert names == ["C1", "HA"]
    assert dst.read_text().endswith("TER\nEND\n")


def test_write_coordinates_rounds_to_pdb_precision(tmp_path):
    f = tmp_path / "x.pdb"
    f.write_text(pdb_line(1, "C1", "CLA", 3, (0, 0, 0)) + pdb_line(2, "C2", "CLA", 3, (1, 1, 1)))
    write_coordinates(str(f), {(3, "C2"): np.array([1.23456, -2.0, 30.0004])})
    coords = extract_all_coordinates(str(f))
    assert np.array_equal(coords[3]["C2"], [1.235, -2.0, 30.0])
    assert np.array_equal(coords[3]["C1"], [0, 0, 0])


def test_clash_score_grows_as_atoms_approach():
    coords = {1: {"A": np.array([0.0, 0.0, 0.0])}}
    tree, labels, _ = build_kdtree(coords)
    near = count_clashes_for_segment({"X": np.array([2.0, 0, 0])}, tree, labels, 3.4)
    far = count_clashes_for_segment({"X": np.array([3.3, 0, 0])}, tree, labels, 3.4)
    alone = count_clashes_for_segment({"X": np.array([20.0, 0, 0])}, tree, labels, 3.4)
    assert near > far > alone == 0.1


def test_fragment_neighbours_use_lennard_jones():
    points = np.array([[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]])
    labels = [(1, "A"), (1, "B")]
    tree, _, _ = build_kdtree({1: {"A": points[0], "B": points[1]}})
    trial = {"X": np.array([3.0, 0.0, 0.0])}
    mie_score = count_clashes_for_segment(trial, tree, labels, 3.4)
    lj_score = count_clashes_for_segment(trial, tree, labels, 3.4, fragment={(1, "A")})
    assert mie_score != lj_score
    assert count_clashes_for_segment(trial, tree, labels, 3.4, fragment={(2, "A")}) == mie_score
