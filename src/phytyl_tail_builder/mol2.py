"""Reading bond topology from MOL2 residue templates."""

import os
import re

HYDROGEN = re.compile(r'^(H|1H|2H|3H|HH)')


def parse_mol2_template(mol2_path, sorted_atom_names):
    """Heavy-atom bonds of a template that involve the rebuilt atoms.

    Returns (atom_order, bonds, atom_set):
      atom_order  atoms in the order they first appear in those bonds, with the
                  existing anchor atom first;
      bonds       (atom1, atom2) pairs in template order, with the bond from the
                  existing structure to the rebuilt atoms first;
      atom_set    set(atom_order).
    A missing file gives ([], [], set()).
    """
    atom_order = []
    bonds = []
    anchor_bond = None

    if not os.path.exists(mol2_path):
        return atom_order, bonds, set()

    names = {}
    all_names = set()
    section = None
    with open(mol2_path) as f:
        for line in f:
            if "@<TRIPOS>" in line:
                section = line.strip()
                continue
            parts = line.split()
            if section == "@<TRIPOS>ATOM" and len(parts) > 1:
                names[int(parts[0])] = parts[1]
                all_names.add(parts[1])
            elif section == "@<TRIPOS>BOND" and len(parts) >= 3:
                origin = names.get(int(parts[1]))
                target = names.get(int(parts[2]))
                if not origin or not target or HYDROGEN.match(origin) or HYDROGEN.match(target):
                    continue

                if origin in sorted_atom_names or target in sorted_atom_names:
                    bonds.append((origin, target))
                    if origin not in atom_order:
                        atom_order.append(origin)
                    if target not in atom_order:
                        atom_order.append(target)

                existing = all_names - set(sorted_atom_names)
                if (origin in sorted_atom_names and target in existing) or \
                   (target in sorted_atom_names and origin in existing):
                    anchor_bond = (origin, target)

    if anchor_bond:
        bonds.insert(0, anchor_bond)
        anchor = anchor_bond[1] if anchor_bond[1] not in sorted_atom_names else anchor_bond[0]
        if anchor not in atom_order:
            atom_order.insert(0, anchor)

    return atom_order, bonds, set(atom_order)
