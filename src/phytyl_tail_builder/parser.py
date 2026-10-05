"""Reading tleap input scripts and leap.log files."""

import re
from pathlib import Path


def parse_tleap_in(tleap_in_path):
    """Return (loadpdb file, savepdb file, prmtop, inpcrd) named in a tleap script.

    Entries that are not present are None; a missing script gives all None.
    """
    pdb_in = pdb_out = parm_file = crd_file = None
    path = Path(tleap_in_path)
    if not path.exists():
        return pdb_in, pdb_out, parm_file, crd_file

    with open(path) as f:
        for line in f:
            if 'loadpdb' in line:
                pdb_in = line.split()[-1]
            elif 'savepdb' in line:
                pdb_out = line.split()[-1]
            elif 'saveamberparm' in line:
                parts = line.split()
                parm_file, crd_file = parts[-2], parts[-1]
    return pdb_in, pdb_out, parm_file, crd_file


def _read(path):
    path = Path(path)
    return path.read_text() if path.exists() else ""


def parse_leap_log(leap_log_path):
    """Return {residue_number: residue_name} for residues that received heavy atoms."""
    matches = re.findall(r"Added missing heavy atom.*?(\w+?)\s+(\d+)", _read(leap_log_path))
    return {int(num): name for name, num in matches}


def parse_leap_log_for_atoms(leap_log_path):
    """Return [(residue_number, residue_name, atom_name), ...] added in the last LEaP session."""
    last_session = _read(leap_log_path).split("Log file started:")[-1]
    matches = re.findall(r"Added missing heavy atom: .R<(\w+) (\d+)>.\w+<(\w+) \d+>", last_session)
    return [(int(num), name, atom) for name, num, atom in matches]
