from pathlib import Path

import numpy as np
import pytest

from phytyl_tail_builder.mol2 import build_order, neighbours, read_mol2
from phytyl_tail_builder.stereo import angle, check_residue, dihedral
from phytyl_tail_builder.torsion import build_torsion, nerf, plan

MINIMAL = Path(__file__).resolve().parent.parent / "examples" / "minimal"
TAIL = [f"C{i}" for i in range(1, 21)]


@pytest.mark.parametrize("torsion", [-179.0, -60.0, 0.0, 60.0, 123.4, 180.0])
def test_nerf_round_trip(torsion):
    a, b, c = np.array([0.3, 1.2, -0.4]), np.array([0.0, 0.0, 0.0]), np.array([1.5, 0.1, 0.3])
    d = nerf(a, b, c, 1.53, 111.0, torsion)
    assert np.linalg.norm(d - c) == pytest.approx(1.53)
    assert angle(b, c, d) == pytest.approx(111.0)
    assert np.cos(np.radians(dihedral(a, b, c, d) - torsion)) == pytest.approx(1.0)


def test_nerf_rebuilds_template_from_its_internal_coordinates():
    template, bonds = read_mol2(MINIMAL / "CLA.mol2")
    present = [n for n in template if n not in TAIL]
    known = {n: template[n] for n in present}
    for atom, p in build_order(bonds, TAIL, present):
        g = next(n for n in neighbours(bonds)[p] if n in known and n != atom)
        d = next(n for n in neighbours(bonds)[g] if n in known and n != p)
        t = template
        known[atom] = nerf(known[d], known[g], known[p], np.linalg.norm(t[atom] - t[p]),
                           angle(t[g], t[p], t[atom]), dihedral(t[d], t[g], t[p], t[atom]))
        assert np.allclose(known[atom], t[atom], atol=1e-6), atom


def test_plan_roles():
    template, bonds = read_mol2(MINIMAL / "CHL.mol2")
    present = [n for n in template if n not in TAIL]
    order = build_order(bonds, TAIL, present)
    steps = {s["atom"]: s for s in plan(order, template, neighbours(bonds), present)}
    assert "dihedral" in steps["C1"] and steps["C1"]["refs"][1:] == ("CGA", "O2A")  # ester, from template
    assert steps["C2"].get("free") and steps["C3"].get("free") and steps["C6"].get("free")
    assert steps["C5"]["dihedral"] == 180.0 and steps["C5"]["refs"] == ("C1", "C2", "C3")
    assert {a: steps[a]["sibling"] for a in ("C4", "C9", "C14", "C20")} == \
        {"C4": "C5", "C9": "C10", "C14": "C15", "C20": "C19"}
    assert steps["C10"].get("free") and steps["C19"].get("free")


def _isolated(template_name, missing, mirror=False):
    """A lone residue with the template's ring and nothing around it; mirror: its mirror image (C8, C13 S)."""
    template, bonds = read_mol2(MINIMAL / f"{template_name}.mol2")
    if mirror:
        template = {n: xyz * [1.0, 1.0, -1.0] for n, xyz in template.items()}
    coords = {0: {n: xyz.copy() for n, xyz in template.items()}}
    present = [n for n in template if n not in missing]
    return coords, build_order(bonds, missing, present), template, neighbours(bonds)


@pytest.mark.parametrize("name,mirror", [("CLA", False), ("CHL", False), ("CLA", True), ("CHL", True)])
def test_full_tail_is_natural_even_from_an_s_template(name, mirror):
    coords, order, template, graph = _isolated(name, TAIL, mirror)
    build_torsion(coords, 0, order, template, graph, excluded=set())
    result = check_residue(coords[0])
    assert (result["C8"], result["C13"], result["C2=C3"]) == ("R", "R", "E")
    assert result["bad_angles"] == [] and result["bad_bonds"] == []


def test_resolved_atoms_and_their_dihedrals_are_kept():
    missing = TAIL[5:]  # C6..C20, as in CLA 669
    coords, order, template, graph = _isolated("CLA", missing)
    before = {n: xyz.copy() for n, xyz in coords[0].items() if n not in missing}
    placed, _ = build_torsion(coords, 0, order, template, graph, excluded=set())
    assert set(placed) == set(missing)
    for n, xyz in before.items():
        assert np.array_equal(coords[0][n], xyz)


def test_chain_atom_next_to_a_resolved_branch_goes_on_the_r_side():
    """C7, C8, C9 resolved from an S,S template (mirrored CHL): C10 is placed relative to C9 so C8 is R."""
    missing = TAIL[9:]  # C10..C20
    coords, order, template, graph = _isolated("CHL", missing, mirror=True)
    build_torsion(coords, 0, order, template, graph, excluded=set())
    result = check_residue(coords[0])
    assert (result["C8"], result["C13"]) == ("R", "R") and result["bad_angles"] == []


@pytest.mark.parametrize("missing", [TAIL, TAIL[5:], TAIL[9:]])
def test_template_bond_order_does_not_matter(missing):
    """The same build from a CLA.mol2 whose bonds are listed backwards, each one reversed."""
    template, bonds = read_mol2(MINIMAL / "CLA.mol2")
    bonds = [(b, a) for a, b in reversed(bonds)]
    coords = {0: {n: xyz.copy() for n, xyz in template.items()}}
    order = build_order(bonds, missing, [n for n in template if n not in missing])
    build_torsion(coords, 0, order, template, neighbours(bonds), excluded=set())
    result = check_residue(coords[0])
    assert (result["C8"], result["C13"], result["C2=C3"]) == ("R", "R", "E") and result["bad_angles"] == []
