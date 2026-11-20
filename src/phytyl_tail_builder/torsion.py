"""Torsion-space tail building (--builder torsion).

Every missing atom is placed by NeRF from three known atoms with the
template's bond length and bond angle. Only backbone dihedrals of the phytyl
chain are free; they are scanned on a grid with a beam search scored by
optimizer.count_clashes_for_segment. Branches are placed at the template's
dihedral offset from their chain sibling, with the sign that makes C8 and C13
R. C1-C2-C3-C5 is fixed at 180 degrees (E); every other non-free dihedral
(the ester CGA-O2A, ring atoms) is copied from the template.
"""

import numpy as np

from .optimizer import build_kdtree, count_clashes_for_segment
from .stereo import CENTRES, R_SIGN, TAIL_ATOMS, angle, chirality, dihedral

FIXED_BONDS = {frozenset(("C2", "C3"))}
FIXED_DIHEDRALS = {("C1", "C2", "C3", "C5"): 180.0}
# The beam keeps at most one partial tail per DIVERSITY-degree cell of its free dihedrals. With a fine grid
# the best partial tails otherwise differ by a step or two and all run into the same clash further on.
DIVERSITY = 30.0


def nerf(a, b, c, bond, angle, torsion):
    """Position d with |cd| = bond, angle b-c-d = angle and dihedral a-b-c-d = torsion (degrees)."""
    bc = c - b
    bc /= np.linalg.norm(bc)
    n = np.cross(b - a, bc)
    n /= np.linalg.norm(n)
    m = np.cross(n, bc)
    angle, torsion = np.radians(angle), np.radians(torsion)
    return c + bond * (-np.cos(angle) * bc + np.sin(angle) * (np.cos(torsion) * m + np.sin(torsion) * n))


def plan(order, template, graph, known):
    """Internal coordinates of every missing atom.

    order is build_order() output, template {name: xyz}, graph the template
    neighbour lists, known the names with coordinates before building.
    Returns a list of dicts: atom, refs (d, g, p), bond, angle and either
    'free' (dihedral scanned), 'sibling' + 'offset' (placed relative to an
    atom on the same parent), or 'dihedral' (fixed value).
    """
    available = set(known)
    steps = []
    for atom, p in order:
        # g: a known neighbour of the parent that has a known neighbour d of its own; the first
        # neighbour in bond-list order can be a sibling with nothing behind it
        for g in (n for n in graph[p] if n in available and n != atom):
            d = next((n for n in graph.get(g, ()) if n in available and n != p), None)
            if d is not None:
                break
        else:
            raise ValueError(f"cannot place {atom}: needs three bonded known atoms behind it")
        t = template
        step = {"atom": atom, "refs": (d, g, p), "bond": np.linalg.norm(t[atom] - t[p]),
                "angle": angle(t[g], t[p], t[atom])}
        sibling = next((n for n in graph[p] if n in available and n not in (g, atom)), None)
        if sibling is not None:
            step["sibling"] = sibling
            step["offset"] = dihedral(t[d], t[g], t[p], t[atom]) - dihedral(t[d], t[g], t[p], t[sibling])
        elif (d, g, p, atom) in FIXED_DIHEDRALS:
            step["dihedral"] = FIXED_DIHEDRALS[(d, g, p, atom)]
        elif g in TAIL_ATOMS and p in TAIL_ATOMS and frozenset((g, p)) not in FIXED_BONDS:
            step["free"] = True
        else:
            step["dihedral"] = dihedral(t[d], t[g], t[p], t[atom])
        steps.append(step)
        available.add(atom)
    return steps


def _place(step, xyz, value):
    d, g, p = (xyz[n] for n in step["refs"])
    return nerf(d, g, p, step["bond"], step["angle"], value)


def _place_relative(step, xyz):
    """Place a branch at +/- the template offset from its sibling; prefer R at a stereocentre."""
    d, g, p = step["refs"]
    base = dihedral(xyz[d], xyz[g], xyz[p], xyz[step["sibling"]])
    first = _place(step, xyz, base + step["offset"])
    if p not in CENTRES:
        return first
    trial = dict(xyz, **{step["atom"]: first})
    if all(n in trial for n in CENTRES[p]) and chirality(trial, p, *CENTRES[p]) != R_SIGN[p]:
        return _place(step, xyz, base - step["offset"])
    return first


def _staggered(value):
    """Distance (degrees) from the nearest of -60, 60, 180."""
    return min(abs((value - s + 180) % 360 - 180) for s in (-60, 60, 180))


def build_torsion(coordinates, res, order, template, graph, excluded, fragment_potential="lj",
                  environment_potential="mie", keep_initial_positions=False, beam=5, dihedral_step=10.0):
    """Build the missing atoms of one residue in torsion space.

    coordinates is {residue: {atom_name: xyz}}; the missing atoms of `res`
    (those in order) are replaced in place. Atoms in `excluded` are left out
    of the environment, and so are the missing atoms of `res` unless
    keep_initial_positions. Scores are summed over atoms; before each free
    dihedral only the `beam` best partial tails are kept (ties: closest to
    staggered), at most one per DIVERSITY cell of their free dihedrals.
    Returns ({atom_name: xyz}, total_score).
    """
    atoms = coordinates[res]
    missing = [atom for atom, _ in order]
    known = {n: xyz for n, xyz in atoms.items() if n not in missing}
    steps = plan(order, template, graph, known)

    skip = set(excluded)
    if not keep_initial_positions:
        skip |= {(res, a) for a in missing}
    kdtree, labels, _ = build_kdtree(coordinates, exclude_atoms=skip)
    anchors = {(res, p) for atom, p in order if p in known}
    grid = np.arange(-180.0, 180.0, dihedral_step)

    def best(states):
        kept, cells = [], set()
        for state in sorted(states, key=lambda state: (round(state[0], 6), state[1])):
            cell = tuple(int((value + 180.0) % 360.0 // DIVERSITY) for value in state[3])
            if cell not in cells:
                cells.add(cell)
                kept.append(state)
                if len(kept) == beam:
                    break
        return kept

    # Prune just before each free dihedral, so the atoms fixed by the previous one
    # (branches, C5) are scored before its value is chosen.
    states = [(0.0, 0.0, {}, ())]  # (score, stagger penalty, placed, free dihedral values)
    for step in steps:
        if step.get("free"):
            states = best(states)
        expanded = []
        for score, penalty, placed, values in states:
            xyz = {**known, **placed}
            if step.get("free"):
                trials = [(_place(step, xyz, value), _staggered(value), (value,)) for value in grid]
            elif "sibling" in step:
                trials = [(_place_relative(step, xyz), 0.0, ())]
            else:
                trials = [(_place(step, xyz, step["dihedral"]), 0.0, ())]
            extra = {(res, n): v for n, v in placed.items()}
            p = step["refs"][2]                                 # parent: 1-2 and 1-3 partners are not scored
            bonded = {(res, p)} | {(res, n) for n in graph[p]}
            for position, stagger, value in trials:
                s = count_clashes_for_segment({step["atom"]: position}, kdtree, labels, anchors,
                                              fragment_potential, environment_potential, extra, bonded)
                expanded.append((score + s, penalty + stagger, {**placed, step["atom"]: position},
                                 values + value))
        states = expanded

    score, _, placed, _ = best(states)[0]
    for atom, position in placed.items():
        atoms[atom] = np.round(position, 3)
    return placed, score
