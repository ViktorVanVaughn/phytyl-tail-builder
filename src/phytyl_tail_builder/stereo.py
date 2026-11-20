"""Stereochemistry of phytyl tails, in the PDB atom names of chlorophylls and pheophytins.

Natural phytyl is (2E,7R,11R): phytyl C7 and C11 are template atoms C8 and
C13, and the double bond is C2=C3 (E means C1-C2-C3-C5 near 180 degrees).
R is not hard-coded: it is the sign of the signed volume in REFERENCE, the
tail of examples/minimal/CLA.mol2, whose CIP labels are (2E,7R,11R). Tail
angles are judged against their values in REFERENCE too.
"""

import numpy as np
from scipy.spatial import cKDTree

TAIL_BONDS = [
    ("CGA", "O2A"), ("O2A", "C1"), ("C1", "C2"), ("C2", "C3"), ("C3", "C4"), ("C3", "C5"),
    ("C5", "C6"), ("C6", "C7"), ("C7", "C8"), ("C8", "C9"), ("C8", "C10"), ("C10", "C11"),
    ("C11", "C12"), ("C12", "C13"), ("C13", "C14"), ("C13", "C15"), ("C15", "C16"),
    ("C16", "C17"), ("C17", "C18"), ("C18", "C19"), ("C18", "C20"),
]
TAIL_ATOMS = ["O2A"] + [f"C{i}" for i in range(1, 21)]
# Stereocentre -> its heavy substituents (chain towards C1, methyl, chain away).
CENTRES = {"C8": ("C7", "C9", "C10"), "C13": ("C12", "C14", "C15")}
EZ_DIHEDRAL = ("C1", "C2", "C3", "C5")
# A tail angle fails more than ANGLE_TOLERANCE degrees from its value in natural phytyl (REFERENCE): thermal
# motion at 300 K spreads tail angles by about 4 degrees, while a misplaced atom gives 160-180 degrees.
ANGLE_TOLERANCE = 20.0
BOND_TOLERANCE = 0.1
# Two heavy atoms more than three bonds apart and closer than this overlap.
CONTACT_LIMIT = 2.0
# A tetrahedral tail centre has |signed volume| near 2.5 A^3. Below this its three substituents are coplanar
# with it (to under a degree), and the sign only reflects rounding of the coordinates.
PLANAR_VOLUME = 0.05

# Tail of examples/minimal/CLA.mol2, (2E,7R,11R) by CIP (checked with RDKit).
REFERENCE = {name: np.array(xyz) for name, xyz in {
    "CGA": (1.824, -0.681, -3.106), "O2A": (1.831, -2.030, -3.167),
    "C1": (2.970, -2.626, -3.848), "C2": (2.813, -4.110, -3.756), "C3": (3.718, -5.009, -3.338),
    "C4": (5.115, -4.672, -2.876), "C5": (3.363, -6.485, -3.321), "C6": (3.134, -7.053, -1.904),
    "C7": (2.752, -8.540, -1.922), "C8": (2.395, -9.144, -0.547), "C9": (3.600, -9.148, 0.408),
    "C10": (1.799, -10.556, -0.727), "C11": (1.192, -11.179, 0.539), "C12": (0.528, -12.538, 0.273),
    "C13": (-0.191, -13.170, 1.484), "C14": (0.791, -13.527, 2.612), "C15": (-1.017, -14.393, 1.031),
    "C16": (-1.951, -14.987, 2.096), "C17": (-2.825, -16.127, 1.552), "C18": (-3.855, -16.699, 2.546),
    "C19": (-3.195, -17.351, 3.771), "C20": (-4.786, -17.698, 1.842),
}.items()}


def has_phytyl(atom_names):
    """True for a residue with a phytyl ester in PDB atom names: CGA and tail atoms (O2A, C1-C20).

    The PDB chemical components CLA, CHL, CL7, F6C (chlorophylls a, b, d, f), BCL, BCB (bacteriochlorophylls
    a, b), PHO and BPH (pheophytins) all name and bond their phytyl atoms this way.
    """
    return "CGA" in atom_names and any(n in atom_names for n in TAIL_ATOMS)


def chirality(coords, center, a, b, c):
    """Sign of dot(a - center, cross(b - center, c - center)): +1, -1, or 0 below PLANAR_VOLUME (planar)."""
    o = coords[center]
    volume = float(np.dot(coords[a] - o, np.cross(coords[b] - o, coords[c] - o)))
    return 0 if abs(volume) < PLANAR_VOLUME else int(np.sign(volume))


def dihedral(p1, p2, p3, p4):
    """IUPAC dihedral angle in degrees, in (-180, 180]; nan if three points are collinear."""
    b1, b2, b3 = p2 - p1, p3 - p2, p4 - p3
    n1, n2 = np.cross(b1, b2), np.cross(b2, b3)
    if np.linalg.norm(n1) < 1e-6 or np.linalg.norm(n2) < 1e-6:
        return float("nan")
    return float(np.degrees(np.arctan2(np.linalg.norm(b2) * np.dot(b1, n2), np.dot(n1, n2))))


def angle(p1, p2, p3):
    """Angle p1-p2-p3 in degrees."""
    u, v = p1 - p2, p3 - p2
    cos = np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))
    return float(np.degrees(np.arccos(np.clip(cos, -1.0, 1.0))))


R_SIGN = {centre: chirality(REFERENCE, centre, *subst) for centre, subst in CENTRES.items()}


def centre_label(coords, centre):
    """'R' or 'S' for a stereocentre, None if an atom is missing or it is planar."""
    names = (centre,) + CENTRES[centre]
    if not all(n in coords for n in names):
        return None
    sign = chirality(coords, centre, *CENTRES[centre])
    return None if sign == 0 else ("R" if sign == R_SIGN[centre] else "S")


def ez_label(phi):
    """'E' (|phi| >= 150), 'Z' (|phi| <= 30), 'twisted' in between, None if undefined."""
    if phi is None or np.isnan(phi):
        return None
    return "E" if abs(phi) >= 150 else "Z" if abs(phi) <= 30 else "twisted"


def _neighbours(names):
    graph = {}
    for a, b in TAIL_BONDS:
        if a in names and b in names:
            graph.setdefault(a, []).append(b)
            graph.setdefault(b, []).append(a)
    return graph


def check_residue(coords_by_name):
    """Stereo and geometry of one phytyl tail.

    Returns a dict:
      C8, C13         'R', 'S' or None (atom missing / planar)
      C2=C3           'E', 'Z', 'twisted' or None (see ez_label)
      dihedral        C1-C2-C3-C5 in degrees, or None
      bad_angles      [(a, b, c, degrees, deviation)] heavy-atom tail angles more than
                      ANGLE_TOLERANCE from their value in REFERENCE
      bad_bonds       [(a, b, length, reference)] tail bonds off REFERENCE by > 0.1 Å
      worst_angle     (a, b, c, degrees, deviation) furthest from REFERENCE, or None
    """
    coords = coords_by_name
    result = {centre: centre_label(coords, centre) for centre in CENTRES}
    phi = None
    if all(n in coords for n in EZ_DIHEDRAL):
        phi = dihedral(*(coords[n] for n in EZ_DIHEDRAL))
    result["dihedral"] = phi
    result["C2=C3"] = ez_label(phi)

    angles = []
    for b, neighbours in _neighbours(coords).items():
        if b not in TAIL_ATOMS:
            continue
        for i, a in enumerate(neighbours):
            for c in neighbours[i + 1:]:
                value = angle(coords[a], coords[b], coords[c])
                angles.append((a, b, c, value, value - angle(REFERENCE[a], REFERENCE[b], REFERENCE[c])))
    result["bad_angles"] = [x for x in angles if abs(x[4]) > ANGLE_TOLERANCE]
    result["worst_angle"] = max(angles, key=lambda x: abs(x[4]), default=None)

    result["bad_bonds"] = []
    for a, b in TAIL_BONDS:
        if a in coords and b in coords:
            length = float(np.linalg.norm(coords[a] - coords[b]))
            ref = float(np.linalg.norm(REFERENCE[a] - REFERENCE[b]))
            if abs(length - ref) > BOND_TOLERANCE:
                result["bad_bonds"].append((a, b, length, ref))
    return result


def is_ok(result):
    """Gate: (R, R, E) and every tail angle within ANGLE_TOLERANCE of natural phytyl."""
    return (result["C8"] == "R" and result["C13"] == "R" and result["C2=C3"] == "E"
            and not has_errors(result))


def has_errors(result):
    """Anything wrong among what can be determined: S, not E, or a tail angle off natural phytyl."""
    return ("S" in (result["C8"], result["C13"]) or result["C2=C3"] not in (None, "E")
            or bool(result["bad_angles"]))


def _is_heavy(name):
    return not name.lstrip("0123456789").startswith("H")


def _heavy_tree(residues):
    owners = [(i, n) for i, res in enumerate(residues) for n in res['atoms'] if _is_heavy(n)]
    points = np.array([residues[i]['atoms'][n] for i, n in owners]).reshape(-1, 3)
    return cKDTree(points), owners


def min_contact(residues, index, tree=None, cutoff=6.0):
    """Shortest distance (Å) from a tail atom of residues[index] to a heavy atom.

    Pairs inside the residue that are three or fewer bonds apart are skipped.
    Bonds are TAIL_BONDS within the tail and anything under 1.9 Å elsewhere in
    the residue. Contacts beyond cutoff are ignored. Returns
    (distance, tail_atom, 'RES resSeq ATOM') or None. tree is _heavy_tree(residues),
    built here when not given.
    """
    tree, owners = tree or _heavy_tree(residues)
    atoms = residues[index]['atoms']
    tail = [n for n in TAIL_ATOMS if n in atoms]
    graph = {n: set(nb) for n, nb in _neighbours(atoms).items()}
    rest = [n for n in atoms if n not in TAIL_ATOMS and _is_heavy(n)]
    for i, a in enumerate(rest):
        for b in rest[i + 1:]:
            if np.linalg.norm(atoms[a] - atoms[b]) < 1.9:
                graph.setdefault(a, set()).add(b)
                graph.setdefault(b, set()).add(a)

    best = None
    for t in tail:
        near, frontier = {t}, {t}
        for _ in range(3):
            frontier = {nb for n in frontier for nb in graph.get(n, ())} - near
            near |= frontier
        for j in tree.query_ball_point(atoms[t], cutoff):
            i, n = owners[j]
            if i == index and n in near:
                continue
            d = float(np.linalg.norm(atoms[t] - tree.data[j]))
            if best is None or d < best[0]:
                best = (d, t, f"{residues[i]['name']} {residues[i]['resseq']} {n}")
    return best


# ok: yes/no for a complete tail (is_ok); for a partial one, no if has_errors, else incomplete;
# no for any tail with a heavy-atom contact below CONTACT_LIMIT.
REPORT_COLUMNS = ["residue", "resseq", "resname", "C8", "C13", "C2=C3", "C1-C2-C3-C5",
                  "min_contact", "contact_pair", "worst_angle", "ok"]


def report_rows(residues, indices):
    """stereo_report.tsv rows (dicts over REPORT_COLUMNS) for residues[i], i in indices."""
    tree = _heavy_tree(residues)
    rows = []
    for index in indices:
        res = residues[index]
        result = check_residue(res['atoms'])
        contact = min_contact(residues, index, tree)
        worst = result["worst_angle"]
        ok = ("yes" if is_ok(result) else "no") if all(n in res['atoms'] for n in TAIL_ATOMS) \
            else "no" if has_errors(result) else "incomplete"
        if contact is not None and contact[0] < CONTACT_LIMIT:
            ok = "no"
        rows.append({
            "residue": index + 1, "resseq": res['resseq'], "resname": res['name'],
            "C8": result["C8"] or "-", "C13": result["C13"] or "-", "C2=C3": result["C2=C3"] or "-",
            "C1-C2-C3-C5": "-" if result["dihedral"] is None else f"{result['dihedral']:.1f}",
            "min_contact": "-" if contact is None else f"{contact[0]:.2f}",
            "contact_pair": "-" if contact is None else f"{contact[1]}-{contact[2]}",
            "worst_angle": "-" if worst is None else f"{'-'.join(worst[:3])} {worst[3]:.1f}",
            "ok": ok,
        })
    return rows


def write_report(path, rows, extra_columns=()):
    columns = REPORT_COLUMNS + list(extra_columns)
    with open(path, "w") as f:
        f.write("\t".join(columns) + "\n")
        for row in rows:
            f.write("\t".join(str(row.get(c, "")) for c in columns) + "\n")
