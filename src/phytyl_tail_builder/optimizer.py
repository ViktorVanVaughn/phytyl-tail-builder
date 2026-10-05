"""Placement of rebuilt tail atoms by spherical sampling and clash scoring."""

from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from .mol2 import parse_mol2_template
from .pdb import extract_all_coordinates, write_coordinates

POTENTIALS = ("lj", "mie")


def build_kdtree(coordinates, exclude_atoms=None):
    """Build a KD-tree over all atoms not in exclude_atoms.

    coordinates is {residue_number: {atom_name: xyz}}; exclude_atoms is a set of
    (residue_number, atom_name). Returns (kdtree, labels, points).
    """
    exclude_atoms = set(exclude_atoms or ())
    points = []
    labels = []
    for res_num, atoms in coordinates.items():
        for atom_name, xyz in atoms.items():
            if (res_num, atom_name) not in exclude_atoms:
                points.append(xyz)
                labels.append((res_num, atom_name))
    points = np.array(points)
    return cKDTree(points), labels, points


def extract_atom_coordinates(bonds, atom_names, res_num, coordinates):
    """Return positions of atom_names in res_num and the lengths of the given bonds."""
    residue = coordinates.get(res_num, {})
    atom_coords = {name: xyz for name, xyz in residue.items() if name in atom_names}
    bond_distances = {}
    for atom1, atom2 in bonds:
        if atom1 in atom_coords and atom2 in atom_coords:
            bond_distances[(atom1, atom2)] = np.linalg.norm(atom_coords[atom2] - atom_coords[atom1])
    return atom_coords, bond_distances


def count_clashes_for_segment(new_positions, kdtree, labels, base_distance_threshold, fragment=None,
                              fragment_potential="lj", environment_potential="mie"):
    """Clash score of trial positions against the atoms in kdtree.

    Neighbours whose label is in `fragment` (a set of (residue_number, atom_name))
    are scored with `fragment_potential`, all others with `environment_potential`;
    each is "lj" (Lennard-Jones 12-6) or "mie" (Mie 16-6). Each term is weighted
    by cutoff / distance, and the cutoff grows with the local atom density.
    """
    for name in (fragment_potential, environment_potential):
        if name not in POTENTIALS:
            raise ValueError(f"potential must be one of {POTENTIALS}, got {name!r}")
    fragment = fragment or ()
    clash_score = 0.0
    epsilon = 0.36
    sigma = 3.4
    mie_n = 16
    mie_m = 6

    def adaptive_threshold(position):
        nearby = kdtree.query_ball_point(position, base_distance_threshold)
        density = len(nearby) / (base_distance_threshold ** 3)
        return base_distance_threshold * (1 + density)

    def lennard_jones(distance):
        return 4 * epsilon * ((sigma / distance) ** 12 - (sigma / distance) ** 6)

    def mie(distance):
        return epsilon * (
            (mie_n / mie_m) ** (mie_m / (mie_n - mie_m))
            * ((sigma / distance) ** mie_n - (sigma / distance) ** mie_m)
        )

    for position in new_positions.values():
        threshold = adaptive_threshold(position)
        for idx in kdtree.query_ball_point(position, threshold):
            distance = np.linalg.norm(position - kdtree.data[idx])
            if distance < threshold:
                kind = fragment_potential if labels[idx] in fragment else environment_potential
                potential = lennard_jones(distance) if kind == "lj" else mie(distance)
                clash_score += potential * (threshold / distance)

    return max(clash_score, 1e-1)


def apply_rotation(center_atom, segment_coords, bond_distances, kdtree, labels, distance_threshold,
                   fragment=None, fragment_potential="lj", environment_potential="mie"):
    """Try 20 x 10 directions around center_atom at fixed bond length.

    Returns (best_positions, best_score). best_positions excludes the center atom.
    """
    center = segment_coords.get(center_atom)
    if center is None:
        return segment_coords, 0

    min_clashes = float('inf')
    best_positions = segment_coords.copy()

    for theta in np.linspace(0, 2 * np.pi, 20):
        for phi in np.linspace(0, np.pi, 10):
            new_positions = {}
            for atom in segment_coords:
                if atom == center_atom:
                    continue
                bond_key = (center_atom, atom)
                if bond_key not in bond_distances:
                    bond_key = (atom, center_atom)
                bond_length = bond_distances.get(bond_key, 0)
                new_positions[atom] = np.array([
                    center[0] + bond_length * np.sin(phi) * np.cos(theta),
                    center[1] + bond_length * np.sin(phi) * np.sin(theta),
                    center[2] + bond_length * np.cos(phi),
                ])

            score = count_clashes_for_segment(new_positions, kdtree, labels, distance_threshold, fragment,
                                              fragment_potential, environment_potential)
            if score < min_clashes:
                min_clashes = score
                best_positions = new_positions.copy()

    return best_positions, min_clashes


def optimize_residue_tails(pdb_path, sorted_atoms, mol2_dir, fragment_potential="lj",
                           environment_potential="mie", keep_initial_positions=False):
    """Place the rebuilt atoms of every residue in sorted_atoms, in that order.

    Residues are processed one after another, bond by bond along the template
    order, and each placement sees all earlier ones. The neighbour search uses
    coordinates rounded to the 3 decimals of the PDB format. pdb_path is
    rewritten in place at the end.

    Trial positions are scored with a pair potential, "lj" (Lennard-Jones 12-6)
    or "mie" (Mie 16-6). fragment_potential applies to the atoms of the tail
    that are already placed and to the atom the tail grows from;
    environment_potential applies to every other atom.

    By default the tleap position of the atom being placed, and of the tail
    atoms that are still to be placed, is not part of the environment. With
    keep_initial_positions=True the tleap positions of the residue stay in the
    environment throughout.
    """
    for name in (fragment_potential, environment_potential):
        if name not in POTENTIALS:
            raise ValueError(f"potential must be one of {POTENTIALS}, got {name!r}")

    mol2_dir = Path(mol2_dir)
    coordinates = extract_all_coordinates(pdb_path)
    updates = {}

    residue_atoms = {}
    for res_num, res_name, atom_name in sorted_atoms:
        residue_atoms.setdefault(res_num, {'res_name': res_name, 'atom_names': []})
        residue_atoms[res_num]['atom_names'].append(atom_name)

    # Rebuilt atoms are left out of the neighbour search until their residue is placed.
    excluded = {(res_num, atom_name) for res_num, _, atom_name in sorted_atoms}

    for res_num, details in residue_atoms.items():
        atom_names = details['atom_names']
        mol2_path = mol2_dir / f"{details['res_name']}.mol2"
        if not mol2_path.exists():
            continue

        _, bonds, _ = parse_mol2_template(str(mol2_path), atom_names)
        placed, bond_distances = extract_atom_coordinates(bonds, atom_names, res_num, coordinates)

        for atom in atom_names:
            excluded.discard((res_num, atom))

        tail_done = set()
        total_clash_score = 0
        for atom1, atom2 in bonds:
            if atom1 == atom2:
                continue

            segment = {}
            if atom2 in placed:
                segment[atom2] = placed[atom2]
            if atom1 in placed:
                segment[atom1] = placed[atom1]
            if not segment:
                continue

            if not keep_initial_positions:
                if atom1 not in segment:
                    placed.update(segment)
                    continue
                # Leave this residue's atoms that are being placed, or are still
                # to be placed, out of the environment.
                moving = {atom for atom in segment if atom != atom1}
                pending = {atom for atom in atom_names if atom not in tail_done and atom != atom1}
                skip = excluded | {(res_num, atom) for atom in moving | pending}
            else:
                skip = excluded
            fragment = {(res_num, atom) for atom in tail_done}

            kdtree, labels, _ = build_kdtree(coordinates, exclude_atoms=skip)
            rotated, score = apply_rotation(atom1, segment, bond_distances, kdtree, labels,
                                            distance_threshold=3.4, fragment=fragment,
                                            fragment_potential=fragment_potential,
                                            environment_potential=environment_potential)
            placed.update(rotated)
            total_clash_score += score
            tail_done.add(atom1)
            tail_done.update(rotated)

            for atom_name, position in rotated.items():
                residue = coordinates.get(res_num, {})
                if atom_name in residue:
                    residue[atom_name] = np.array([float(f"{v:8.3f}") for v in position])
                    updates[(res_num, atom_name)] = position

        num_atoms = len(atom_names)
        avg_clash = total_clash_score / (num_atoms - 1) if num_atoms > 0 else 0
        print(f"Residue {res_num}: Average clash score per atom = {avg_clash:.2f}")
        print(f"Progress: {100 * (1 - len(excluded) / len(sorted_atoms)):.1f}%")

    write_coordinates(pdb_path, updates)
