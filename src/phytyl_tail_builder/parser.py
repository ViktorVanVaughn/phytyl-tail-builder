"""Reading tleap input scripts and leap.log files."""

import logging
import re
from pathlib import Path

logger = logging.getLogger(__name__)


def parse_tleap_in(tleap_in_path):
    """Return (loadpdb file, savepdb file, prmtop, inpcrd) named in a tleap script.

    Comments (# to end of line) are ignored and commands match in any case.
    Entries that are not present are None; a missing script gives all None.
    """
    pdb_in = pdb_out = parm_file = crd_file = None
    path = Path(tleap_in_path)
    if not path.exists():
        return pdb_in, pdb_out, parm_file, crd_file

    with open(path) as f:
        for line in f:
            parts = [p.strip('"\'') for p in line.split('#', 1)[0].split()]
            words = [p.lower() for p in parts]
            if 'loadpdb' in words:
                pdb_in = parts[-1]
            elif 'savepdb' in words:
                pdb_out = parts[-1]
            elif 'saveamberparm' in words:
                parm_file, crd_file = parts[-2], parts[-1]
    return pdb_in, pdb_out, parm_file, crd_file


def leap_sessions(text):
    """Split a leap.log into tleap sessions.

    tleap appends to leap.log and writes "log started:" whenever a log is
    opened, also for every logFile command of a leaprc. A session starts at a
    "log started:" line that does not follow a logFile command.
    """
    sessions, previous = [[]], ""
    for line in text.splitlines():
        if re.match(r'\s*log (file )?started:', line, re.I) and not re.search(r'\blogfile\b', previous, re.I):
            sessions.append([])
        sessions[-1].append(line)
        if line.strip():
            previous = line
    return ["\n".join(s) for s in sessions if s]


def parse_leap_log_for_atoms(leap_log_path):
    """Return [(residue_number, residue_name, atom_name), ...] added in the last tleap session."""
    path = Path(leap_log_path)
    sessions = leap_sessions(path.read_text()) if path.exists() else []
    if len(sessions) > 1:
        logger.warning("%s holds %d tleap sessions; reading only the last one", path, len(sessions))
    last = sessions[-1] if sessions else ""
    matches = re.findall(r"Added missing heavy atom: .R<(\w+) (\d+)>.A<(\w+) \d+>", last)
    return [(int(num), name, atom) for name, num, atom in matches]
