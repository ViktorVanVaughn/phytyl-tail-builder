from pathlib import Path

import numpy as np
import pytest

from phytyl_tail_builder.main import build
from phytyl_tail_builder.mol2 import read_mol2
from phytyl_tail_builder.pdb import read_residues
from phytyl_tail_builder.stereo import (
    CENTRES,
    R_SIGN,
    REFERENCE,
    angle,
    check_residue,
    chirality,
    dihedral,
    has_errors,
    is_ok,
    report_rows,
)

MINIMAL = Path(__file__).resolve().parent.parent / "examples" / "minimal"
V010 = MINIMAL / "expected" / "default" / "LCC_no_sc_no_Hy.pdb"


def residues_by_resseq(path, resseqs):
    residues = read_residues(path)
    return {r['resseq']: r for r in residues if r['resseq'] in resseqs and r['name'] in ("CLA", "CHL")}


def test_reference_is_the_cla_template():
    template, _ = read_mol2(MINIMAL / "CLA.mol2")
    for name, xyz in REFERENCE.items():
        assert np.allclose(template[name], xyz, atol=1e-3), name


def test_chirality_calibration_on_templates():
    """CLA.mol2 and the corrected CHL.mol2 are (2E,7R,11R)."""
    cla, _ = read_mol2(MINIMAL / "CLA.mol2")
    chl, _ = read_mol2(MINIMAL / "CHL.mol2")
    assert {k: check_residue(cla)[k] for k in ("C8", "C13", "C2=C3")} == {"C8": "R", "C13": "R", "C2=C3": "E"}
    assert {k: check_residue(chl)[k] for k in ("C8", "C13", "C2=C3")} == {"C8": "R", "C13": "R", "C2=C3": "E"}
    # Independent check: with the substituents in CIP order (C7 > C10 > C9 and C12 > C15 > C14),
    # R means 1 -> 2 -> 3 runs clockwise with H at the back, i.e. a negative signed volume.
    assert chirality(cla, "C8", "C7", "C10", "C9") == -1
    assert chirality(cla, "C13", "C12", "C15", "C14") == -1
    assert R_SIGN == {c: chirality(cla, c, *s) for c, s in CENTRES.items()}


def test_mirror_image_inverts_both_centres():
    mirrored = {n: xyz * [1, 1, -1] for n, xyz in REFERENCE.items()}
    result = check_residue(mirrored)
    assert (result["C8"], result["C13"], result["C2=C3"]) == ("S", "S", "E")


def test_planar_centre_has_no_label():
    """A centre coplanar with its substituents gets None, however its coordinates are rounded."""
    coords = dict(REFERENCE)
    c8, c7 = REFERENCE["C8"], REFERENCE["C7"]
    coords["C10"] = c8 + (c8 - c7)                              # C7-C8-C10 collinear: volume 0
    for shift in (0.0, 1e-3, -1e-3):
        coords["C9"] = REFERENCE["C9"] + shift
        assert check_residue(coords)["C8"] is None
    assert check_residue(REFERENCE)["C8"] == "R"


def test_dihedral():
    p = [np.array(v, float) for v in ([1, 0, 0], [0, 0, 0], [0, 0, 1], [0, 1, 1])]
    assert dihedral(*p) == pytest.approx(90.0)
    assert dihedral(p[0], p[1], p[2], np.array([-1.0, 0, 1])) == pytest.approx(180.0)
    assert dihedral(p[0], p[1], p[2], np.array([1.0, 0, 1])) == pytest.approx(0.0)
    assert np.isnan(dihedral(p[0], p[1], np.array([2.0, 0, 0]), p[3]))


def test_check_residue_flags_angles_bonds_and_z():
    coords = dict(REFERENCE)
    assert check_residue(coords)["bad_angles"] == [] and check_residue(coords)["bad_bonds"] == []
    coords["C5"], coords["C4"] = REFERENCE["C4"], REFERENCE["C5"]
    coords["C20"] = REFERENCE["C18"] + 2.0 * (REFERENCE["C18"] - REFERENCE["C17"]) / 1.53
    result = check_residue(coords)
    assert result["C2=C3"] == "Z"
    assert ("C17", "C18", "C20") in [a[:3] for a in result["bad_angles"]]
    assert result["worst_angle"][3] == pytest.approx(180.0)
    assert ("C18", "C20") in [b[:2] for b in result["bad_bonds"]]
    assert not is_ok(result)


def bend(coords, a, b, c, degrees):
    """Move the end atom a, in the plane of a-b-c, so that the angle a-b-c is `degrees`."""
    u = (coords[c] - coords[b]) / np.linalg.norm(coords[c] - coords[b])
    v = coords[a] - coords[b]
    w = v - (v @ u) * u
    w /= np.linalg.norm(w)
    t = np.radians(degrees)
    return dict(coords, **{a: coords[b] + np.linalg.norm(v) * (np.cos(t) * u + np.sin(t) * w)})


def test_angles_are_judged_against_natural_phytyl():
    """12 degrees off at the sp2 C3 is thermal motion and passes (the fixed 130 limit failed it);
    30 degrees off at the sp3 C18 fails, whatever other angle is further off."""
    c4 = angle(REFERENCE["C2"], REFERENCE["C3"], REFERENCE["C4"])
    warm = check_residue(bend(REFERENCE, "C4", "C3", "C2", c4 + 12))
    assert warm["bad_angles"] == [] and is_ok(warm)
    c20 = angle(REFERENCE["C17"], REFERENCE["C18"], REFERENCE["C20"])
    hot = check_residue(bend(REFERENCE, "C20", "C18", "C17", c20 + 30))
    assert ("C17", "C18", "C20") in [a[:3] for a in hot["bad_angles"]]
    assert has_errors(hot) and not is_ok(hot)


def test_overlapping_atoms_fail():
    """A natural tail with an environment atom 1.8 A from C10 fails; at 2.5 A it passes."""
    out = REFERENCE["C10"] - (REFERENCE["C8"] + REFERENCE["C11"]) / 2      # away from both neighbours
    out /= np.linalg.norm(out)
    for distance, ok in ((1.8, "no"), (2.5, "yes")):
        water = {"O": REFERENCE["C10"] + distance * out}
        residues = [{"name": "CLA", "resseq": "1", "chain": "", "atoms": dict(REFERENCE)},
                    {"name": "HOH", "resseq": "2", "chain": "", "atoms": water}]
        row = report_rows(residues, [0])[0]
        assert row["ok"] == ok and float(row["min_contact"]) == pytest.approx(distance, abs=0.01)


def test_review_findings_on_v010_output():
    """Step 0: the v0.1.0 golden output (examples/minimal/expected/default) is stereochemically wrong."""
    res = residues_by_resseq(V010, {"663", "669", "719"})
    result = {k: check_residue(r['atoms']) for k, r in res.items()}
    assert (result["669"]["C8"], result["669"]["C13"]) == ("S", "S")
    assert (result["719"]["C8"], result["719"]["C13"]) == ("S", "S")
    assert (result["663"]["C8"], result["663"]["C13"]) == ("R", "S")
    # Not E: C1-C2-C3-C5 = 100.5 and -137.2 degrees.
    assert result["663"]["dihedral"] == pytest.approx(100.5, abs=0.1)
    assert result["719"]["dihedral"] == pytest.approx(-137.2, abs=0.1)
    assert result["663"]["C2=C3"] == result["719"]["C2=C3"] == "twisted"
    worst = sorted(max(a[3] for a in r["bad_angles"]) for r in result.values())
    assert worst[0] == pytest.approx(160.0, abs=0.1) and worst[-1] == pytest.approx(180.0, abs=0.05)
    # The 180-degree angle C10-C11-C12 leaves dihedrals through it undefined.
    a = res["669"]['atoms']
    assert np.isnan(dihedral(a["C9"], a["C10"], a["C11"], a["C12"]))
    assert not any(is_ok(r) for r in result.values())


def test_report_rows_mark_incomplete_tails():
    residues = read_residues(MINIMAL / "LCC_input.pdb")
    rows = {r["resseq"]: r for r in report_rows(residues, range(len(residues))) if r["resname"] == "CHL"}
    assert rows["663"]["ok"] == "incomplete" and rows["663"]["C8"] == "-"
    assert rows["440"]["C8"] == "S" and rows["440"]["ok"] == "no"  # resolved in the cryo-EM model


def test_example_output_is_natural_phytyl(tmp_path):
    """Every rebuilt tail of the example is (2E,7R,11R) with every tail angle near natural phytyl."""
    out, rows = build(MINIMAL / "tleap.in", MINIMAL / "leap.log", output_prefix=tmp_path,
                      fix_resolved_stereo=True)
    residues = read_residues(out)
    rebuilt = [r["residue"] - 1 for r in rows]
    assert len(rebuilt) == 12
    for index in rebuilt:
        result = check_residue(residues[index]['atoms'])
        assert (result["C8"], result["C13"], result["C2=C3"]) == ("R", "R", "E"), residues[index]['resseq']
        assert result["bad_angles"] == []
        # CHL residues of the cryo-EM model have a resolved 1.45 Å CGA-O2A bond (template 1.35 Å).
        assert {b[:2] for b in result["bad_bonds"]} <= {("CGA", "O2A")}
