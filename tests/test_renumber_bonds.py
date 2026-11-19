import subprocess
import sys
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "renumber_bonds.py"


def pdb(path, resseqs):
    """One atom per residue, at x = position; resSeq written in columns 23-27 as tleap does."""
    path.write_text("".join(f"ATOM  {k + 1:5d}  C1  PA   {s:>5}   {float(k):8.3f}   0.000   0.000\n"
                            for k, s in enumerate(resseqs)))
    return path


@pytest.mark.parametrize("post_first,expected", [(1, "bond x.3.C21 x.5.C12"),
                                                 (9998, "bond x.10000.C21 x.10002.C12")])
def test_bond_lines_follow_tleap_sequence_numbers(tmp_path, post_first, expected):
    """tleap numbers 16, 220, 221, 440, 441 as 16..20: residue 18 is resSeq 221, residue 20 is 441."""
    pre = pdb(tmp_path / "pre.pdb", [16, 220, 221, 440, 441])
    post = pdb(tmp_path / "post.pdb", range(post_first, post_first + 5))
    script = tmp_path / "check.in"
    script.write_text("x = loadpdb post.pdb\nbond x.18.C21 x.20.C12\n# bond x.1.C1 x.2.C1\n")
    subprocess.run([sys.executable, str(SCRIPT), str(pre), str(post), str(script)], check=True)
    assert script.read_text().splitlines()[1:] == [expected, "# bond x.1.C1 x.2.C1"]


def test_unknown_residue_stops(tmp_path):
    pre = pdb(tmp_path / "pre.pdb", [16, 17])
    script = tmp_path / "check.in"
    script.write_text("bond x.30.C21 x.16.C12\n")
    result = subprocess.run([sys.executable, str(SCRIPT), str(pre), str(pre), str(script)],
                            capture_output=True, text=True)
    assert result.returncode != 0 and "residue 30" in result.stderr
