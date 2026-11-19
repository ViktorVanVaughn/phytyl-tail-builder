"""Clash scoring, and placement of rebuilt atoms by spherical sampling (--builder sphere)."""

import numpy as np
from scipy.spatial import cKDTree

POTENTIALS = ("lj", "mie")

# Lennard-Jones sigma per element (Å): 2 * Rmin/2 / 2**(1/6) from Amber parameters
# (C, N, O, S, P: parm99/GAFF; H: hc; MG: CLX.frcmod). Other elements use C.
SIGMA = {"C": 3.40, "N": 3.25, "O": 2.96, "H": 2.65, "MG": 1.16, "S": 3.56, "P": 3.74}


def element_of(atom_name):
    """Element of a PDB atom name: MG, else the first letter after any digits."""
    name = atom_name.upper().lstrip("0123456789")
    # shortcut: only MG is read as two letters (Na+, Cl- give N, C); add ions here if they sit near tails
    return "MG" if name == "MG" else name[:1]


def sigma_of(atom_name):
    return SIGMA.get(element_of(atom_name), SIGMA["C"])


def build_kdtree(coordinates, exclude_atoms=None):
    """Build a KD-tree over all atoms not in exclude_atoms.

    coordinates is {residue: {atom_name: xyz}}; exclude_atoms is a set of
    (residue, atom_name). Returns (kdtree, labels, points).
    """
    exclude_atoms = set(exclude_atoms or ())
    points = []
    labels = []
    for res_num, atoms in coordinates.items():
        for atom_name, xyz in atoms.items():
            if (res_num, atom_name) not in exclude_atoms:
                points.append(xyz)
                labels.append((res_num, atom_name))
    points = np.array(points).reshape(-1, 3)
    return cKDTree(points), labels, points


def count_clashes_for_segment(new_positions, kdtree, labels, base_distance_threshold, fragment=None,
                              fragment_potential="lj", environment_potential="mie", extra=None):
    """Clash score of trial positions against the atoms in kdtree (and in extra).

    new_positions is {atom_name: xyz}. Neighbours whose label is in `fragment`
    (a set of (residue, atom_name)), and every atom of `extra`
    ({(residue, atom_name): xyz}), are scored with `fragment_potential`, all
    others with `environment_potential`; each is "lj" (Lennard-Jones 12-6) or
    "mie" (Mie 16-6), with sigma = mean of the two element sigmas (SIGMA).
    Each term is weighted by cutoff / distance, and the cutoff grows with the
    local atom density.
    """
    for name in (fragment_potential, environment_potential):
        if name not in POTENTIALS:
            raise ValueError(f"potential must be one of {POTENTIALS}, got {name!r}")
    fragment = fragment or ()
    extra = extra or {}
    extra_labels = list(extra)
    extra_points = np.array([extra[k] for k in extra_labels]).reshape(-1, 3)
    clash_score = 0.0
    epsilon = 0.36
    mie_n = 16
    mie_m = 6

    def lennard_jones(distance, sigma):
        return 4 * epsilon * ((sigma / distance) ** 12 - (sigma / distance) ** 6)

    def mie(distance, sigma):
        return epsilon * (
            (mie_n / mie_m) ** (mie_m / (mie_n - mie_m))
            * ((sigma / distance) ** mie_n - (sigma / distance) ** mie_m)
        )

    for atom_name, position in new_positions.items():
        sigma_i = sigma_of(atom_name)
        extra_distances = np.linalg.norm(extra_points - position, axis=1)
        nearby = len(kdtree.query_ball_point(position, base_distance_threshold))
        nearby += int((extra_distances < base_distance_threshold).sum())
        threshold = base_distance_threshold * (1 + nearby / base_distance_threshold ** 3)

        neighbours = [(np.linalg.norm(position - kdtree.data[idx]), labels[idx], labels[idx] in fragment)
                      for idx in kdtree.query_ball_point(position, threshold)]
        neighbours += [(d, label, True) for d, label in zip(extra_distances, extra_labels) if d < threshold]
        for distance, label, in_fragment in neighbours:
            if distance < threshold:
                sigma = (sigma_i + sigma_of(label[1])) / 2
                kind = fragment_potential if in_fragment else environment_potential
                potential = lennard_jones(distance, sigma) if kind == "lj" else mie(distance, sigma)
                clash_score += potential * (threshold / distance)

    return clash_score


def apply_rotation(center, atom_name, bond_length, kdtree, labels, distance_threshold, fragment=None,
                   fragment_potential="lj", environment_potential="mie"):
    """Try 20 x 10 directions around center at fixed bond length.

    Returns (best_position, best_score). Ties keep the first direction.
    """
    best_score = float('inf')
    best_position = None
    for theta in np.linspace(0, 2 * np.pi, 20):
        for phi in np.linspace(0, np.pi, 10):
            position = center + bond_length * np.array(
                [np.sin(phi) * np.cos(theta), np.sin(phi) * np.sin(theta), np.cos(phi)])
            score = count_clashes_for_segment({atom_name: position}, kdtree, labels, distance_threshold,
                                              fragment, fragment_potential, environment_potential)
            if score < best_score:
                best_score, best_position = score, position
    return best_position, best_score


def build_sphere(coordinates, res, order, excluded, fragment_potential="lj", environment_potential="mie",
                 keep_initial_positions=False):
    """Place the atoms of one residue one by one on a sphere around their parent.

    coordinates is {residue: {atom_name: xyz}} holding tleap positions for the
    missing atoms; it is updated in place. order is build_order() output.
    Bond lengths are the tleap ones. Atoms in `excluded` ((residue, atom_name))
    are left out of the environment. Returns ({atom_name: xyz}, total_score).
    """
    atoms = coordinates[res]
    bond_length = {atom: np.linalg.norm(atoms[atom] - atoms[parent]) for atom, parent in order}
    pending = {atom for atom, _ in order}
    placed, total, tail_done = {}, 0.0, set()
    for atom, parent in order:
        skip = set(excluded)
        if not keep_initial_positions:
            skip |= {(res, a) for a in pending}
        fragment = {(res, a) for a in tail_done | {parent}}
        kdtree, labels, _ = build_kdtree(coordinates, exclude_atoms=skip)
        position, score = apply_rotation(atoms[parent], atom, bond_length[atom], kdtree, labels, 3.4,
                                         fragment, fragment_potential, environment_potential)
        atoms[atom] = np.round(position, 3)
        placed[atom] = position
        pending.discard(atom)
        tail_done |= {atom, parent}
        total += score
    return placed, total
