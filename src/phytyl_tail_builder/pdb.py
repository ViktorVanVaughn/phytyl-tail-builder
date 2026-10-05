"""Reading and writing PDB files."""

import re

import numpy as np

HYDROGEN = re.compile(r'^(H|1H|2H|3H|HH)')


def _is_atom(line):
    return line.startswith('ATOM') or line.startswith('HETATM')


def _atom_key(line):
    return int(line[22:26].strip()), line[12:16].strip().replace(" ", "")


def extract_all_atom_coords(pdb_path):
    """Return {residue_number: [(x, y, z), ...]} in file order."""
    coords = {}
    with open(pdb_path) as f:
        for line in f:
            if _is_atom(line):
                resnum = int(line[22:26].strip())
                xyz = tuple(float(x) for x in line[30:54].split())
                coords.setdefault(resnum, []).append(xyz)
    return coords


def extract_all_coordinates(pdb_path):
    """Return {residue_number: {atom_name: np.array([x, y, z])}}."""
    coords = {}
    with open(pdb_path) as f:
        for line in f:
            if _is_atom(line):
                res_num, atom_name = _atom_key(line)
                xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                coords.setdefault(res_num, {})[atom_name] = xyz
    return coords


def correct_residue_numbering(pre_pdb_path, post_pdb_path, leap_log_resnames):
    """Find the post-LEaP number of each rebuilt residue by matching coordinates.

    A post-LEaP residue matches when it contains every pre-LEaP coordinate of
    the residue. Returns (leap_log_resnames, {post_number: residue_name}).
    """
    pre_coords = extract_all_atom_coords(pre_pdb_path)
    post_coords = extract_all_atom_coords(post_pdb_path)

    corrected = {}
    for resnum, resname in leap_log_resnames.items():
        if resnum in pre_coords:
            pre_res_coords = pre_coords[resnum]
            for post_resnum, post_res_coords in post_coords.items():
                if all(c in post_res_coords for c in pre_res_coords):
                    corrected[post_resnum] = resname
                    break

    return leap_log_resnames, corrected


def remove_hydrogens_from_rebuilt(input_pdb_path, output_pdb_path, residues_to_strip):
    """Copy a PDB, dropping hydrogens of the given residue numbers."""
    with open(input_pdb_path) as infile, open(output_pdb_path, 'w') as outfile:
        for line in infile:
            if _is_atom(line):
                res_num_str = line[22:26].strip()
                res_num = int(res_num_str) if res_num_str else None
                if res_num in residues_to_strip and HYDROGEN.match(line[12:16].strip()):
                    continue
            outfile.write(line)


def format_coordinate(position):
    """Format x, y, z as the 24 PDB coordinate columns (3 x %8.3f)."""
    return f"{position[0]:8.3f}{position[1]:8.3f}{position[2]:8.3f}"


def write_coordinates(pdb_path, updates):
    """Rewrite pdb_path in place with new positions.

    updates maps (residue_number, atom_name) to a position; every matching
    ATOM/HETATM line is changed.
    """
    with open(pdb_path) as f:
        lines = f.readlines()
    with open(pdb_path, 'w') as f:
        for line in lines:
            if _is_atom(line):
                key = _atom_key(line)
                if key in updates:
                    line = f"{line[:30]}{format_coordinate(updates[key])}{line[54:]}"
            f.write(line)
