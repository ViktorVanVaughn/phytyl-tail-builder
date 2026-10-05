import subprocess
import sys
from pathlib import Path

import pytest

pytest.importorskip("matplotlib")

ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / "examples" / "minimal"


def test_plot_tail_writes_figure(tmp_path):
    out = tmp_path / "tail.png"
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts" / "plot_tail.py"),
            str(EXAMPLE / "LCC_no_sc.pdb"),
            str(EXAMPLE / "expected" / "default" / "LCC_no_sc_no_Hy.pdb"),
            str(EXAMPLE / "CLA.mol2"),
            "669",
            "-o",
            str(out),
        ],
        check=True,
    )
    assert out.stat().st_size > 10_000
