"""End-to-end workflow."""

import logging
from pathlib import Path

from .mol2 import build_order, neighbours, read_mol2
from .optimizer import POTENTIALS, build_sphere
from .parser import parse_leap_log_for_atoms, parse_tleap_in
from .pdb import (
    extract_all_coordinates,
    first_number,
    map_residues,
    read_residues,
    remove_hydrogens_from_rebuilt,
    residue_label,
    write_coordinates,
)
from .stereo import TAIL_BONDS, check_residue, has_phytyl, report_rows, write_report
from .torsion import build_torsion

logger = logging.getLogger(__name__)

BUILDERS = ("torsion", "sphere")


def build(tleap_in_path, leap_log_path, templates_dir=None, output_prefix=".",
          fragment_potential="lj", environment_potential="mie", keep_initial_positions=False,
          builder="torsion", beam=5, dihedral_step=10.0, fix_resolved_stereo=False, passes=2):
    """Rebuild tail positions, write <savepdb name>_no_Hy.pdb and stereo_report.tsv.

    tleap_in_path  tleap script with loadpdb (before LEaP) and savepdb (after LEaP)
    leap_log_path  leap.log of that run (only its last session is read)
    templates_dir  folder with <RESNAME>.mol2 files (default: folder of tleap_in_path)
    output_prefix  output folder
    builder        "torsion" (template geometry, scanned dihedrals) or "sphere" (v0.1 algorithm)
    beam, dihedral_step
                   beam width and dihedral grid step (degrees) of the torsion builder
    passes         torsion builder: 1 builds each residue once, in order; every further pass
                   rebuilds each tail with all the other tails in place
    fix_resolved_stereo
                   discard resolved tail atoms from a wrong stereocentre (or Z/twisted
                   C2=C3) outward and rebuild them (torsion builder only)
    fragment_potential, environment_potential, keep_initial_positions
                   see optimizer.count_clashes_for_segment and build_sphere

    Returns (path, rows): the written PDB (or the savepdb file when there is
    nothing to place) and the stereo report rows of the rebuilt residues.
    """
    for name in (fragment_potential, environment_potential):
        if name not in POTENTIALS:
            raise ValueError(f"potential must be one of {POTENTIALS}, got {name!r}")
    if builder not in BUILDERS:
        raise ValueError(f"builder must be one of {BUILDERS}, got {builder!r}")
    if passes < 1:
        raise ValueError("passes must be at least 1")
    if beam < 1:
        raise ValueError("beam must be at least 1")
    if not 0 < dihedral_step <= 360:
        raise ValueError(f"dihedral_step must be in (0, 360] degrees, got {dihedral_step}")
    if fix_resolved_stereo and builder != "torsion":
        raise ValueError("fix_resolved_stereo needs the torsion builder")

    tleap_in_path = Path(tleap_in_path)
    leap_log_path = Path(leap_log_path)
    templates_dir = Path(templates_dir) if templates_dir is not None else tleap_in_path.parent

    for label, path in (("tleap input", tleap_in_path), ("leap.log", leap_log_path),
                        ("Templates directory", templates_dir)):
        if not path.exists():
            raise FileNotFoundError(f"{label} not found: {path}")

    pdb_in, pdb_out, _, _ = parse_tleap_in(tleap_in_path)
    if not pdb_in or not pdb_out:
        raise ValueError("tleap input must contain loadpdb and savepdb lines")

    pdb_in = tleap_in_path.parent / pdb_in
    pdb_out = tleap_in_path.parent / pdb_out
    for path in (pdb_in, pdb_out):
        if not path.exists():
            raise FileNotFoundError(f"PDB not found: {path}")

    extracted_atoms = parse_leap_log_for_atoms(leap_log_path)
    if not extracted_atoms:
        logger.warning("No added heavy atoms in %s", leap_log_path)
        return pdb_out, []

    templates = {}
    for _, res_name, _ in extracted_atoms:
        mol2_path = templates_dir / f"{res_name}.mol2"
        if res_name not in templates and mol2_path.exists():
            templates[res_name] = read_mol2(mol2_path)
    pre_residues, post_residues = read_residues(pdb_in), read_residues(pdb_out)
    skipped = sorted({name for _, name, _ in extracted_atoms} - set(templates))
    # a residue with a phytyl ester (any name) must not be skipped: its tail would stay unbuilt and unchecked
    first = first_number(pre_residues)
    unbuilt = sorted({name for num, name, atom in extracted_atoms if name in skipped
                      and 0 <= num - first < len(pre_residues)
                      and has_phytyl(set(pre_residues[num - first]['atoms']) | {atom})})
    if unbuilt:
        raise ValueError(f"No template for {', '.join(unbuilt)} in {templates_dir}: "
                         f"their phytyl tails would stay where tleap placed them, unchecked")
    if skipped:
        logger.warning("No template for %s in %s; their added atoms are left as placed by tleap",
                       ", ".join(skipped), templates_dir)

    mapping = map_residues(pre_residues, post_residues,
                           {num: name for num, name, _ in extracted_atoms if name in templates})
    missing = {}
    for num, _, atom in extracted_atoms:
        if num in mapping:
            label = residue_label(post_residues[mapping[num]])
            if atom not in post_residues[mapping[num]]['atoms']:
                raise ValueError(f"leap.log adds {atom} to residue {num}, but savepdb residue "
                                 f"{label} has no {atom}")
            # an atom tleap added cannot be in its input: if it is, leap.log belongs to other files
            if atom in pre_residues[mapping[num]]['atoms']:
                raise ValueError(f"leap.log adds {atom} to residue {num}, but loadpdb residue "
                                 f"{label} already has {atom}")
            missing.setdefault(mapping[num], []).append(atom)
    if not missing:
        logger.warning("No residues with a matching template")
        return pdb_out, []
    logger.info("%d added atoms in %d residues", sum(map(len, missing.values())), len(missing))

    no_h_pdb = Path(output_prefix) / f"{pdb_out.stem}_no_Hy.pdb"
    no_h_pdb.parent.mkdir(parents=True, exist_ok=True)
    remove_hydrogens_from_rebuilt(str(pdb_out), str(no_h_pdb), set(missing))
    coordinates = extract_all_coordinates(str(no_h_pdb))

    discarded = {}
    for res, atoms in missing.items():
        wrong = resolved_stereo_errors(
            {n: xyz for n, xyz in coordinates[res].items() if n not in atoms})
        if not wrong:
            continue
        label = residue_label(post_residues[res])
        if fix_resolved_stereo:
            drop = [a for a in discard_from(wrong) if a in coordinates[res] and a not in atoms]
            discarded[res] = drop
            atoms.extend(drop)
            logger.warning("%s: resolved atoms have %s; discarded and rebuilt: %s",
                           label, ", ".join(wrong.values()), " ".join(drop))
        else:
            logger.warning("%s: resolved atoms have %s; kept as in the input "
                           "(--fix-resolved-stereo rebuilds them)", label, ", ".join(wrong.values()))

    # Pass 1 builds the residues one by one, each without the tails still to come.
    # Later passes (torsion builder) rebuild every tail with all the others in place.
    excluded = {(res, atom) for res, atoms in missing.items() for atom in atoms}
    updates = {}
    passes = passes if builder == "torsion" else 1
    for n in range(1, passes + 1):
        for res in sorted(missing, key=lambda r: len(missing[r])):
            template, bonds = templates[post_residues[res]['name']]
            present = [a for a in coordinates[res] if a not in missing[res] and a in template]
            order = build_order(bonds, missing[res], present)
            excluded -= {(res, atom) for atom in missing[res]}
            options = dict(fragment_potential=fragment_potential, environment_potential=environment_potential,
                           keep_initial_positions=keep_initial_positions and n == 1)
            if builder == "sphere":
                placed, score = build_sphere(coordinates, res, order, excluded, neighbours(bonds), **options)
            else:
                placed, score = build_torsion(coordinates, res, order, template, neighbours(bonds), excluded,
                                              beam=beam, dihedral_step=dihedral_step, **options)
            updates.update({(res, atom): xyz for atom, xyz in placed.items()})
            if n == passes:
                print(f"Residue {residue_label(post_residues[res])}: "
                      f"Average clash score per atom = {score / len(placed):.2f}")

    write_coordinates(str(no_h_pdb), updates)
    logger.info("Wrote %s", no_h_pdb)

    rows = report_rows(read_residues(no_h_pdb), sorted(missing))
    for row in rows:
        row["discarded"] = " ".join(discarded.get(row["residue"] - 1, [])) or "-"
    write_report(no_h_pdb.parent / "stereo_report.tsv", rows, ["discarded"])
    return no_h_pdb, rows


def resolved_stereo_errors(coords):
    """{'C8'|'C13'|'C2=C3': description} for resolved atoms with a wrong configuration."""
    result = check_residue(coords)
    wrong = {c: f"{c} = {result[c]}" for c in ("C8", "C13") if result[c] == "S"}
    if result["C2=C3"] not in (None, "E"):
        wrong["C2=C3"] = f"C2=C3 {result['C2=C3']} ({result['dihedral']:.0f} deg)"
    return wrong


def discard_from(wrong):
    """Tail atoms to rebuild for the errors of resolved_stereo_errors().

    A wrong centre drops the centre and everything beyond it; a wrong C2=C3
    keeps C3 and drops C4, C5 and beyond.
    """
    order = build_order(TAIL_BONDS, set(neighbours(TAIL_BONDS)) - {"CGA"}, ["CGA"])
    parent = dict(order)
    roots = [c for c in ("C8", "C13") if c in wrong]
    if "C2=C3" in wrong:
        roots += ["C4", "C5"]
    drop = []
    for atom, _ in order:
        a = atom
        while a in parent and a not in roots:
            a = parent[a]
        if a in roots:
            drop.append(atom)
    return drop
