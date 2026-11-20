"""Output must match the reference structures to 1e-3 Å."""

import os
import shutil
from pathlib import Path

import numpy as np
import pytest

from phytyl_tail_builder.main import build
from phytyl_tail_builder.stereo import TAIL_ATOMS

MINIMAL = Path(__file__).resolve().parent.parent / "examples" / "minimal"
TOL = 1e-3


def read_atoms(path):
    labels, xyz = [], []
    with open(path) as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                labels.append(line[12:27])
                xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return labels, np.array(xyz)


def assert_same_structure(result, reference):
    labels, xyz = read_atoms(result)
    ref_labels, ref_xyz = read_atoms(reference)
    assert labels == ref_labels
    deviation = np.linalg.norm(xyz - ref_xyz, axis=1)
    assert (deviation > TOL).sum() == 0, f"max deviation {deviation.max():.4f} Å"


def run(tmp_path, name, **options):
    out, _ = build(MINIMAL / "tleap.in", MINIMAL / "leap.log", output_prefix=tmp_path / name, **options)
    return out


# examples/minimal/expected/default and all-mie hold the v0.1.0 output (3 residues, no stereo
# control); they are kept unchanged for reference and for tests/test_stereo.py.
@pytest.mark.parametrize("name,options", [
    ("torsion", {}),
    ("torsion-fix-resolved", {"fix_resolved_stereo": True}),
    ("sphere", {"builder": "sphere"}),
    ("sphere-all-mie", {"builder": "sphere", "fragment_potential": "mie", "environment_potential": "mie",
                        "keep_initial_positions": True}),
])
def test_golden_output(tmp_path, name, options):
    out = run(tmp_path, name, **options)
    assert_same_structure(out, MINIMAL / "expected" / name / "LCC_no_sc_no_Hy.pdb")
    assert (out.parent / "stereo_report.tsv").read_text() == \
        (MINIMAL / "expected" / name / "stereo_report.tsv").read_text()


@pytest.mark.parametrize("fragment,environment", [("lj", "lj"), ("mie", "lj"), ("mie", "mie")])
def test_other_combinations_run(tmp_path, fragment, environment):
    out = run(tmp_path, "x", fragment_potential=fragment, environment_potential=environment)
    labels, xyz = read_atoms(out)
    ref_labels, _ = read_atoms(MINIMAL / "expected" / "torsion" / "LCC_no_sc_no_Hy.pdb")
    assert labels == ref_labels and np.isfinite(xyz).all()


@pytest.mark.parametrize("builder", ["torsion", "sphere"])
def test_potentials_change_the_result(tmp_path, builder):
    _, a = read_atoms(run(tmp_path, "a", builder=builder))
    _, b = read_atoms(run(tmp_path, "b", builder=builder, fragment_potential="mie"))
    assert np.linalg.norm(a - b, axis=1).max() > 0.1


def test_unknown_potential_is_rejected(tmp_path):
    with pytest.raises(ValueError):
        run(tmp_path, "u", fragment_potential="morse")


def test_second_pass_separates_neighbouring_tails(tmp_path):
    """Pass 1 builds CLA 717 and CHL 719 one after the other, about 2 A apart; pass 2 sees both tails."""
    def tail_contact(**options):
        _, rows = build(MINIMAL / "tleap.in", MINIMAL / "leap.log", output_prefix=tmp_path / str(options),
                        fix_resolved_stereo=True, **options)
        return float(next(r for r in rows if r["resseq"] == "719")["min_contact"])
    assert tail_contact(passes=1) < 2.3
    assert tail_contact() > 3.0


def test_fine_dihedral_grid_does_not_end_in_a_clash(tmp_path):
    """With 5 degree steps and no diversity, the beam filled with near-identical tails and CHL 663 ended
    1.67 A from PHE 369 O; the tails must stay clear of overlaps on any grid."""
    _, rows = build(MINIMAL / "tleap.in", MINIMAL / "leap.log", output_prefix=tmp_path,
                    fix_resolved_stereo=True, dihedral_step=5.0)
    assert all(r["ok"] == "yes" for r in rows) and min(float(r["min_contact"]) for r in rows) > 2.0


@pytest.mark.parametrize("options", [{"passes": 0}, {"beam": 0}, {"dihedral_step": 0},
                                     {"dihedral_step": -10}])
def test_search_settings_are_validated(tmp_path, options):
    with pytest.raises(ValueError):
        run(tmp_path, "p", **options)


def test_fix_resolved_stereo_needs_torsion(tmp_path):
    with pytest.raises(ValueError):
        run(tmp_path, "f", builder="sphere", fix_resolved_stereo=True)


def test_minimal_example_moves_only_tail_atoms(tmp_path):
    out = run(tmp_path, "m")
    before = dict(zip(*read_atoms(MINIMAL / "LCC_no_sc.pdb")))
    after = dict(zip(*read_atoms(out)))
    moved = [k for k in after if np.linalg.norm(after[k] - before[k]) > TOL]
    assert moved and all(k[5:8] in ("CHL", "CLA") and k[:4].strip() in TAIL_ATOMS for k in moved)


def copy_example(tmp_path, names=("tleap.in", "leap.log", "LCC_input.pdb", "LCC_no_sc.pdb", "CLA.mol2",
                                  "CHL.mol2")):
    for name in names:
        shutil.copy(MINIMAL / name, tmp_path)
    return tmp_path / "tleap.in", tmp_path / "leap.log"


def test_leap_log_of_other_files_is_refused(tmp_path):
    """leap.log residue 127 shifted to 128 (also a CHL) would overwrite the resolved atoms of CHL 666."""
    tleap_in, leap_log = copy_example(tmp_path)
    leap_log.write_text(leap_log.read_text().replace("<CHL 127>", "<CHL 128>"))
    with pytest.raises(ValueError, match="already has C1"):
        build(tleap_in, leap_log, output_prefix=tmp_path / "out")


def test_resolved_atoms_are_kept_unless_discarded(tmp_path):
    inputs = dict(zip(*read_atoms(MINIMAL / "LCC_input.pdb")))
    out, rows = build(MINIMAL / "tleap.in", MINIMAL / "leap.log", output_prefix=tmp_path,
                      fix_resolved_stereo=True)
    after = read_atoms(out)
    discarded = {(r["resseq"], a) for r in rows for a in r["discarded"].split() if a != "-"}
    assert {a for _, a in discarded} == {"C8", "C9", "C10", "C11", "C12", "C13", "C14", "C15"}
    assert {res for res, _ in discarded} == {"440", "666"}
    by_name_res = {(label[5:8], label[10:14].strip(), label[:4].strip()): xyz for label, xyz in zip(*after)}
    for label, xyz in inputs.items():
        key = (label[5:8], label[10:14].strip(), label[:4].strip())
        if key[0] in ("CLA", "CHL") and not key[2].startswith("H") and (key[1], key[2]) not in discarded:
            assert np.allclose(by_name_res[key], xyz, atol=TOL), key


def test_any_residue_name_with_phytyl_atoms(tmp_path):
    """CHL renamed PHO (a pheophytin) everywhere: recognised by its atoms, built exactly as before."""
    tleap_in, leap_log = copy_example(tmp_path)
    for name in ("LCC_input.pdb", "LCC_no_sc.pdb"):
        lines = (tmp_path / name).read_text().splitlines(True)
        renamed = [f"{x[:17]}PHO{x[20:]}" if x[17:20] == "CHL" else x for x in lines]
        (tmp_path / name).write_text("".join(renamed))
    leap_log.write_text(leap_log.read_text().replace(".R<CHL ", ".R<PHO "))
    (tmp_path / "CHL.mol2").rename(tmp_path / "PHO.mol2")
    out, rows = build(tleap_in, leap_log, output_prefix=tmp_path / "out", fix_resolved_stereo=True)
    assert {r["resname"] for r in rows} == {"CLA", "PHO"} and all(r["ok"] == "yes" for r in rows)
    _, xyz = read_atoms(out)
    _, ref = read_atoms(MINIMAL / "expected" / "torsion-fix-resolved" / "LCC_no_sc_no_Hy.pdb")
    assert np.abs(xyz - ref).max() < TOL


def test_missing_chlorophyll_template_is_an_error(tmp_path):
    names = ("tleap.in", "leap.log", "LCC_input.pdb", "LCC_no_sc.pdb", "CLA.mol2")   # no CHL.mol2
    tleap_in, leap_log = copy_example(tmp_path, names)
    with pytest.raises(ValueError, match="No template for CHL"):
        build(tleap_in, leap_log, output_prefix=tmp_path / "out")


@pytest.mark.skipif("PTB_FULL_DATA" not in os.environ,
                    reason="set PTB_FULL_DATA to the full data folder")
def test_full_data(tmp_path):
    """PTB_FULL_DATA holds tleap.in, leap.log, the PDBs and the templates. Every rebuilt tail must
    pass the stereo gate (the v0.1 golden output there is not reproducible after the v0.2 fixes)."""
    data = Path(os.environ["PTB_FULL_DATA"])
    work = tmp_path / "data"
    work.mkdir()
    for f in list(data.glob("*.pdb")) + list(data.glob("*.mol2")) + [data / "tleap.in", data / "leap.log"]:
        shutil.copy(f, work)
    _, rows = build(work / "tleap.in", work / "leap.log", output_prefix=tmp_path / "out",
                    fix_resolved_stereo=True)
    assert rows and all(r["ok"] == "yes" for r in rows), [r for r in rows if r["ok"] != "yes"]
