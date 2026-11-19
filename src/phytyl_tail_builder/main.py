"""End-to-end workflow."""

import logging
from pathlib import Path

from .mol2 import build_order, read_mol2
from .optimizer import POTENTIALS, build_sphere
from .parser import parse_leap_log_for_atoms, parse_tleap_in
from .pdb import (
    extract_all_coordinates,
    map_residues,
    read_residues,
    remove_hydrogens_from_rebuilt,
    residue_label,
    write_coordinates,
)

logger = logging.getLogger(__name__)


def build(tleap_in_path, leap_log_path, templates_dir=None, output_prefix=".",
          fragment_potential="lj", environment_potential="mie", keep_initial_positions=False):
    """Rebuild tail positions and write <savepdb name>_no_Hy.pdb.

    tleap_in_path  tleap script with loadpdb (before LEaP) and savepdb (after LEaP)
    leap_log_path  leap.log of that run (only its last session is read)
    templates_dir  folder with <RESNAME>.mol2 files (default: folder of tleap_in_path)
    output_prefix  output folder
    fragment_potential, environment_potential, keep_initial_positions
                   see optimizer.count_clashes_for_segment and build_sphere

    Returns the path of the written PDB, or of the savepdb file when there is
    nothing to place.
    """
    for name in (fragment_potential, environment_potential):
        if name not in POTENTIALS:
            raise ValueError(f"potential must be one of {POTENTIALS}, got {name!r}")

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
        return pdb_out

    templates = {}
    for _, res_name, _ in extracted_atoms:
        mol2_path = templates_dir / f"{res_name}.mol2"
        if res_name not in templates and mol2_path.exists():
            templates[res_name] = read_mol2(mol2_path)
    skipped = sorted({name for _, name, _ in extracted_atoms} - set(templates))
    if skipped:
        logger.warning("No template for %s in %s; their added atoms are left as placed by tleap",
                       ", ".join(skipped), templates_dir)

    post_residues = read_residues(pdb_out)
    mapping = map_residues(read_residues(pdb_in), post_residues,
                           {num: name for num, name, _ in extracted_atoms if name in templates})
    missing = {}
    for num, _, atom in extracted_atoms:
        if num in mapping:
            if atom not in post_residues[mapping[num]]['atoms']:
                raise ValueError(f"leap.log adds {atom} to residue {num}, but savepdb residue "
                                 f"{residue_label(post_residues[mapping[num]])} has no {atom}")
            missing.setdefault(mapping[num], []).append(atom)
    if not missing:
        logger.warning("No residues with a matching template")
        return pdb_out
    logger.info("%d added atoms in %d residues", sum(map(len, missing.values())), len(missing))

    no_h_pdb = Path(output_prefix) / f"{pdb_out.stem}_no_Hy.pdb"
    no_h_pdb.parent.mkdir(parents=True, exist_ok=True)
    remove_hydrogens_from_rebuilt(str(pdb_out), str(no_h_pdb), set(missing))
    coordinates = extract_all_coordinates(str(no_h_pdb))

    # Residues are built one by one, each without the tails still to come.
    excluded = {(res, atom) for res, atoms in missing.items() for atom in atoms}
    updates = {}
    for res in sorted(missing, key=lambda r: len(missing[r])):
        template, bonds = templates[post_residues[res]['name']]
        present = [a for a in coordinates[res] if a not in missing[res] and a in template]
        order = build_order(bonds, missing[res], present)
        excluded -= {(res, atom) for atom in missing[res]}
        placed, score = build_sphere(coordinates, res, order, excluded, fragment_potential,
                                     environment_potential, keep_initial_positions)
        updates.update({(res, atom): xyz for atom, xyz in placed.items()})
        print(f"Residue {residue_label(post_residues[res])}: "
              f"Average clash score per atom = {score / len(placed):.2f}")

    write_coordinates(str(no_h_pdb), updates)
    logger.info("Wrote %s", no_h_pdb)
    return no_h_pdb
