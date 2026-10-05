"""Output must match the reference structures to 1e-3 Å."""

import os
import shutil
from pathlib import Path

import numpy as np
import pytest

from phytyl_tail_builder.main import build

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
    return build(MINIMAL / "tleap.in", MINIMAL / "leap.log", output_prefix=tmp_path / name, **options)


def test_default_potentials(tmp_path):
    out = run(tmp_path, "d")
    assert_same_structure(out, MINIMAL / "expected" / "default" / "LCC_no_sc_no_Hy.pdb")


def test_all_mie_with_tleap_positions(tmp_path):
    out = run(tmp_path, "m", fragment_potential="mie", environment_potential="mie",
              keep_initial_positions=True)
    assert_same_structure(out, MINIMAL / "expected" / "all-mie" / "LCC_no_sc_no_Hy.pdb")


@pytest.mark.parametrize("fragment,environment", [("lj", "lj"), ("mie", "lj"), ("mie", "mie")])
def test_other_combinations_run(tmp_path, fragment, environment):
    out = run(tmp_path, "x", fragment_potential=fragment, environment_potential=environment)
    labels, xyz = read_atoms(out)
    ref_labels, _ = read_atoms(MINIMAL / "expected" / "default" / "LCC_no_sc_no_Hy.pdb")
    assert labels == ref_labels and np.isfinite(xyz).all()


def test_potentials_change_the_result(tmp_path):
    _, a = read_atoms(run(tmp_path, "a"))
    _, b = read_atoms(run(tmp_path, "b", fragment_potential="mie"))
    assert np.linalg.norm(a - b, axis=1).max() > 0.1


def test_unknown_potential_is_rejected(tmp_path):
    with pytest.raises(ValueError):
        run(tmp_path, "u", fragment_potential="morse")


def test_minimal_example_moves_tail_atoms(tmp_path):
    out = build(MINIMAL / "tleap.in", MINIMAL / "leap.log", output_prefix=tmp_path)
    before = dict(zip(*read_atoms(MINIMAL / "LCC_no_sc.pdb")))
    after = dict(zip(*read_atoms(out)))
    moved = [k for k in after if np.linalg.norm(after[k] - before[k]) > TOL]
    assert moved and all(k[5:8] in ("CHL", "CLA") for k in moved)


@pytest.mark.skipif("PTB_FULL_DATA" not in os.environ,
                    reason="set PTB_FULL_DATA to the full data folder")
def test_full_data(tmp_path):
    """PTB_FULL_DATA holds tleap.in, leap.log, the PDBs, the templates and
    golden/LCC_no_sc_no_Hy.pdb written by the reference script (all-Mie, tleap positions kept)."""
    data = Path(os.environ["PTB_FULL_DATA"])
    work = tmp_path / "data"
    work.mkdir()
    for f in list(data.glob("*.pdb")) + list(data.glob("*.mol2")) + [data / "tleap.in", data / "leap.log"]:
        shutil.copy(f, work)
    out = build(work / "tleap.in", work / "leap.log", output_prefix=tmp_path / "out",
                fragment_potential="mie", environment_potential="mie", keep_initial_positions=True)
    assert_same_structure(out, data / "golden" / "LCC_no_sc_no_Hy.pdb")
