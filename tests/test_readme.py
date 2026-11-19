import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
START, END = "<!-- options-table:start -->", "<!-- options-table:end -->"


def test_readme_options_table_matches_the_parsers():
    """README.md shows the options as scripts/options_table.py generates them from the parsers."""
    readme = (ROOT / "README.md").read_text()
    assert readme.count(START) == 1 and readme.count(END) == 1
    shown = readme.split(START)[1].split(END)[0].strip()
    generated = subprocess.run([sys.executable, str(ROOT / "scripts" / "options_table.py")],
                               capture_output=True, text=True, check=True).stdout.strip()
    assert shown == generated, "regenerate it: python scripts/options_table.py"
