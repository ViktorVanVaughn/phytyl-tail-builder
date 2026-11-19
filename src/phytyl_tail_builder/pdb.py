"""Reading and writing PDB files.

Residues are identified by their position in the file (0-based index over
contiguous ATOM/HETATM records), never by resSeq, which wraps above 9999 and
is not unique across chains.
"""

import re

import numpy as np

HYDROGEN = re.compile(r'^(H|1H|2H|3H|HH)')


def _is_atom(line):
    return line.startswith('ATOM') or line.startswith('HETATM')


def _atom_name(line):
    return line[12:16].strip().replace(" ", "")


def residue_indices(lines):
    """Yield (residue_index, line); residue_index is None for non-atom lines.

    A new residue starts when resName, chainID, resSeq or iCode changes, or
    after a TER record.
    """
    index, key = -1, None
    for line in lines:
        if _is_atom(line):
            if line[17:27] != key:
                index += 1
                key = line[17:27]
            yield index, line
        else:
            if line.startswith('TER'):
                key = None
            yield None, line


def read_residues(pdb_path):
    """Return a list of residues in file order.

    Each residue is a dict with 'name', 'resseq' (as written, with any insertion code;
    tleap's savepdb writes five digits from 10000, into the insertion-code column), 'chain' and
    'atoms' ({atom_name: np.array([x, y, z])}).
    """
    residues = []
    with open(pdb_path) as f:
        for index, line in residue_indices(f):
            if index is None:
                continue
            if index == len(residues):
                residues.append({'name': line[17:20].strip(), 'resseq': line[22:27].strip(),
                                 'chain': line[21].strip(), 'atoms': {}})
            xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            residues[index]['atoms'][_atom_name(line)] = xyz
    return residues


def extract_all_coordinates(pdb_path):
    """Return {residue_index: {atom_name: np.array([x, y, z])}}."""
    return {i: res['atoms'] for i, res in enumerate(read_residues(pdb_path))}


def residue_label(residue):
    return f"{residue['name']} {residue['resseq']}"


def first_number(residues):
    """tleap's number for the first residue of a loadpdb file: its resSeq (insertion code dropped)."""
    if not residues:
        raise ValueError("loadpdb file has no atoms")
    return int(re.match(r'-?\d+', residues[0]['resseq']).group())


def map_residues(pre_residues, post_residues, leap_residues):
    """Map leap.log residue numbers to residue indices of the savepdb file.

    tleap numbers residues sequentially, starting from the resSeq of the first
    residue of the loadpdb file. Every mapped residue must have the same name
    in leap.log, the loadpdb file and the savepdb file, and every atom of the
    loadpdb residue must sit at the same coordinates in the savepdb residue.
    Raises ValueError on any mismatch. Returns {leap_number: residue_index}.
    """
    first = first_number(pre_residues)
    mapping = {}
    for number, name in leap_residues.items():
        index = number - first
        where = f"leap.log residue {name} {number} (residue #{index + 1} of the loadpdb file)"
        if not 0 <= index < min(len(pre_residues), len(post_residues)):
            raise ValueError(f"{where} is outside the loadpdb/savepdb files "
                             f"({len(pre_residues)}/{len(post_residues)} residues)")
        pre, post = pre_residues[index], post_residues[index]
        if pre['name'] != name or post['name'] != name:
            raise ValueError(f"{where}: loadpdb has {residue_label(pre)}, "
                             f"savepdb has {residue_label(post)}")
        post_xyz = {tuple(np.round(xyz, 3)) for xyz in post['atoms'].values()}
        moved = [atom for atom, xyz in pre['atoms'].items() if tuple(np.round(xyz, 3)) not in post_xyz]
        if moved:
            raise ValueError(f"{where}: atoms {', '.join(moved)} of {residue_label(pre)} are not "
                             f"at the same coordinates in savepdb residue {residue_label(post)}")
        mapping[number] = index
    return mapping


def remove_hydrogens_from_rebuilt(input_pdb_path, output_pdb_path, residues_to_strip):
    """Copy a PDB, dropping hydrogens of the given residue indices."""
    with open(input_pdb_path) as infile, open(output_pdb_path, 'w') as outfile:
        for index, line in residue_indices(infile):
            if index in residues_to_strip and HYDROGEN.match(line[12:16].strip()):
                continue
            outfile.write(line)


def format_coordinate(position):
    """Format x, y, z as the 24 PDB coordinate columns (3 x %8.3f)."""
    return f"{position[0]:8.3f}{position[1]:8.3f}{position[2]:8.3f}"


def write_coordinates(pdb_path, updates):
    """Rewrite pdb_path in place with new positions.

    updates maps (residue_index, atom_name) to a position.
    """
    with open(pdb_path) as f:
        lines = f.readlines()
    with open(pdb_path, 'w') as f:
        for index, line in residue_indices(lines):
            key = (index, _atom_name(line))
            if index is not None and key in updates:
                line = f"{line[:30]}{format_coordinate(updates[key])}{line[54:]}"
            f.write(line)
