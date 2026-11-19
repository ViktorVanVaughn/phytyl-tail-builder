"""Reading residue templates from MOL2 files."""

import re
from collections import deque

import numpy as np

HYDROGEN = re.compile(r'^(H|1H|2H|3H|HH)')


def read_mol2(mol2_path):
    """Return (coords, bonds) of the heavy atoms of a template.

    coords is {atom_name: np.array([x, y, z])} in file order; bonds is a list
    of (atom1, atom2) pairs as listed (the order says nothing about direction).
    """
    names, coords, bonds = {}, {}, []
    section = None
    with open(mol2_path) as f:
        for line in f:
            if line.startswith("@<TRIPOS>"):
                section = line.strip()
                continue
            parts = line.split()
            if section == "@<TRIPOS>ATOM" and len(parts) > 4:
                names[int(parts[0])] = parts[1]
                if not HYDROGEN.match(parts[1]):
                    coords[parts[1]] = np.array([float(x) for x in parts[2:5]])
            elif section == "@<TRIPOS>BOND" and len(parts) >= 3:
                a, b = names.get(int(parts[1])), names.get(int(parts[2]))
                if a in coords and b in coords:
                    bonds.append((a, b))
    return coords, bonds


def neighbours(bonds):
    graph = {}
    for a, b in bonds:
        graph.setdefault(a, []).append(b)
        graph.setdefault(b, []).append(a)
    return graph


def build_order(bonds, missing, present):
    """Orient the template bonds by breadth-first search from the present atoms.

    Returns [(atom, parent), ...] for every atom of `missing` reachable from
    `present`, each after its parent. Among the children of one parent, the
    one with the most missing atoms behind it (the chain) comes first.
    Raises ValueError if a missing atom is not connected to a present one.
    """
    graph = neighbours(bonds)
    missing = set(missing)
    parent = {}
    queue = deque(a for a in graph if a in present)
    seen = set(queue)
    while queue:
        atom = queue.popleft()
        for nb in graph.get(atom, ()):
            if nb in missing and nb not in seen:
                parent[nb] = atom
                seen.add(nb)
                queue.append(nb)
    unreachable = missing - set(parent)
    if unreachable:
        raise ValueError(f"atoms {sorted(unreachable)} are not bonded to any present atom in the template")

    children = {}
    for atom, p in parent.items():
        children.setdefault(p, []).append(atom)

    def size(atom):
        return 1 + sum(size(c) for c in children.get(atom, ()))

    order = []
    queue = deque(sorted((a for a in present if a in children), key=list(graph).index))
    while queue:
        p = queue.popleft()
        for child in sorted(children.get(p, ()), key=lambda c: (-size(c), list(graph).index(c))):
            order.append((child, p))
            queue.append(child)
    return order
