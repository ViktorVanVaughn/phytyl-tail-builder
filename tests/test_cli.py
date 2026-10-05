import subprocess
import sys
from pathlib import Path

import pytest

from phytyl_tail_builder import __version__
from phytyl_tail_builder.cli import main

MINIMAL = Path(__file__).resolve().parent.parent / "examples" / "minimal"


def test_help(capsys):
    with pytest.raises(SystemExit) as exc:
        main(["--help"])
    assert exc.value.code == 0
    out = capsys.readouterr().out
    assert "--tleap-in" in out and "--leap-log" in out and "--tail-potential" in out


def test_version(capsys):
    with pytest.raises(SystemExit):
        main(["--version"])
    assert capsys.readouterr().out.strip() == f"phytyl-tail-builder {__version__}"


def test_required_arguments():
    with pytest.raises(SystemExit) as exc:
        main([])
    assert exc.value.code == 2


def test_missing_input_returns_1(tmp_path):
    assert main(["--tleap-in", str(tmp_path / "tleap.in"),
                 "--leap-log", str(tmp_path / "leap.log")]) == 1


def test_run_on_minimal_example(tmp_path):
    result = subprocess.run(
        [sys.executable, "-m", "phytyl_tail_builder.cli",
         "--tleap-in", "tleap.in", "--leap-log", "leap.log", "-o", str(tmp_path)],
        cwd=MINIMAL, capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    assert (tmp_path / "LCC_no_sc_no_Hy.pdb").exists()
    assert result.stdout.count("Average clash score") == 3


def test_unknown_potential_choice():
    with pytest.raises(SystemExit) as exc:
        main(["--tleap-in", "a", "--leap-log", "b", "--tail-potential", "morse"])
    assert exc.value.code == 2
