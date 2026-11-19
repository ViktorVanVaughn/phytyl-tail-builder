import subprocess
import sys
from pathlib import Path

import pytest

from phytyl_tail_builder import __version__
from phytyl_tail_builder.cli import check_main, main

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


def run_cli(tmp_path, *options, module="phytyl_tail_builder.cli"):
    return subprocess.run(
        [sys.executable, "-m", module, "--tleap-in", "tleap.in", "--leap-log", "leap.log",
         "-o", str(tmp_path), *options],
        cwd=MINIMAL, capture_output=True, text=True,
    )


def test_run_on_minimal_example(tmp_path):
    result = run_cli(tmp_path, "--fix-resolved-stereo")
    assert result.returncode == 0, result.stderr
    assert (tmp_path / "LCC_no_sc_no_Hy.pdb").exists()
    assert result.stdout.count("Average clash score") == 12
    assert "CHL 440: resolved atoms have C8 = S; discarded and rebuilt" in result.stderr
    report = (tmp_path / "stereo_report.tsv").read_text().splitlines()
    assert len(report) == 13 and all(line.split("\t")[10] == "yes" for line in report[1:])


def test_stereo_gate_fails_on_resolved_inversions(tmp_path):
    result = run_cli(tmp_path)
    assert result.returncode == 3
    assert "stereo check failed for CHL 440, CHL 666" in result.stderr
    assert "CHL 440: resolved atoms have C8 = S; kept as in the input" in result.stderr
    assert "(--fix-resolved-stereo rebuilds them)" in result.stderr
    assert run_cli(tmp_path, "--allow-stereo-errors").returncode == 0


def test_sphere_builder_fails_the_gate(tmp_path):
    assert run_cli(tmp_path, "--builder", "sphere").returncode == 3
    assert run_cli(tmp_path, "--builder", "sphere", "--fix-resolved-stereo").returncode == 2


@pytest.mark.parametrize("option", [["--beam", "0"], ["--dihedral-step", "0"], ["--passes", "0"]])
def test_bad_search_setting_returns_1(tmp_path, option):
    result = run_cli(tmp_path, *option)
    assert result.returncode == 1 and "Traceback" not in result.stderr, result.stderr


def test_unknown_potential_choice():
    with pytest.raises(SystemExit) as exc:
        main(["--tleap-in", "a", "--leap-log", "b", "--tail-potential", "morse"])
    assert exc.value.code == 2


def test_check_cli(capsys):
    assert check_main([str(MINIMAL / "expected" / "torsion-fix-resolved" / "LCC_no_sc_no_Hy.pdb")]) == 0
    out = capsys.readouterr().out.splitlines()
    assert out[0].startswith("residue\tresseq") and len(out) == 13
    assert check_main([str(MINIMAL / "expected" / "default" / "LCC_no_sc_no_Hy.pdb")]) == 1
    assert check_main([str(MINIMAL / "LCC_input.pdb"), "--resnames", "CHL"]) == 1  # 440, 666: C8 = S
    assert check_main([str(MINIMAL / "LCC_input.pdb"), "--resnames", "CLA"]) == 0  # partial, nothing wrong


def test_check_cli_with_nothing_to_check(tmp_path):
    assert check_main([str(MINIMAL / "LCC_input.pdb"), "--resnames", "cla"]) == 2
    protein = tmp_path / "protein.pdb"
    protein.write_text("".join(line for line in (MINIMAL / "LCC_input.pdb").read_text().splitlines(True)
                               if line[17:20] not in ("CLA", "CHL")))
    assert check_main([str(protein)]) == 2


def test_check_cli_finds_phytyl_residues_by_atoms(tmp_path, capsys):
    renamed = tmp_path / "renamed.pdb"
    renamed.write_text((MINIMAL / "LCC_input.pdb").read_text().replace(" CLA ", " BCL "))
    check_main([str(renamed)])
    names = {line.split("\t")[2] for line in capsys.readouterr().out.splitlines()[1:]}
    assert names == {"BCL", "CHL"}


def test_check_cli_entry_point():
    result = subprocess.run(["phytyl-tail-check", "--help"], capture_output=True, text=True)
    assert result.returncode == 0 and "--resnames" in result.stdout
