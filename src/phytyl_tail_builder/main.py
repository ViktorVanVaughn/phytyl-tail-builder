"""End-to-end workflow."""

import logging
from pathlib import Path

from .mol2 import parse_mol2_template
from .optimizer import optimize_residue_tails
from .parser import parse_leap_log, parse_leap_log_for_atoms, parse_tleap_in
from .pdb import correct_residue_numbering, remove_hydrogens_from_rebuilt

logger = logging.getLogger(__name__)


def build(tleap_in_path, leap_log_path, templates_dir=None, output_prefix=".",
          fragment_potential="lj", environment_potential="mie", keep_initial_positions=False):
    """Rebuild tail positions and write <savepdb name>_no_Hy.pdb.

    tleap_in_path  tleap script with loadpdb (before LEaP) and savepdb (after LEaP)
    leap_log_path  leap.log of that run
    templates_dir  folder with <RESNAME>.mol2 files (default: folder of tleap_in_path)
    output_prefix  output folder
    fragment_potential, environment_potential, keep_initial_positions
                   see optimize_residue_tails

    Returns the path of the written PDB, or of the savepdb file when there is
    nothing to place.
    """
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

    leap_log_resnames = parse_leap_log(leap_log_path)
    extracted_atoms = parse_leap_log_for_atoms(leap_log_path)
    if not extracted_atoms:
        logger.warning("No added heavy atoms in %s", leap_log_path)
        return pdb_out
    logger.info("%d added atoms in %d residues", len(extracted_atoms),
                len({a[0] for a in extracted_atoms}))

    pre_residues, post_residues = correct_residue_numbering(
        str(pdb_in), str(pdb_out), leap_log_resnames)
    sorted_atoms = apply_corrected_resnumbers(
        extracted_atoms, templates_dir, pre_residues, post_residues)
    if not sorted_atoms:
        logger.warning("No residues with a matching template")
        return pdb_out

    no_h_pdb = Path(output_prefix) / f"{pdb_out.stem}_no_Hy.pdb"
    no_h_pdb.parent.mkdir(parents=True, exist_ok=True)
    remove_hydrogens_from_rebuilt(str(pdb_out), str(no_h_pdb),
                                  {res_num for res_num, _, _ in sorted_atoms})

    optimize_residue_tails(str(no_h_pdb), sorted_atoms, templates_dir, fragment_potential,
                           environment_potential, keep_initial_positions)
    logger.info("Wrote %s", no_h_pdb)
    return no_h_pdb


def apply_corrected_resnumbers(extracted_atoms, templates_dir, pre_residues, post_residues):
    """Map leap.log residue numbers to post-LEaP numbers and order the atoms.

    Residues are ordered by number of added atoms (fewest first, stable), and
    atoms within a residue by template bond order. Residues without a template
    are skipped. Returns [(residue_number, residue_name, atom_name), ...].
    """
    residue_map = dict(zip(pre_residues.keys(), post_residues.keys()))
    residue_atoms = {}
    for res_num, res_name, atom_name in extracted_atoms:
        corrected = residue_map.get(res_num, res_num)
        residue_atoms.setdefault(corrected, {'res_name': res_name, 'atoms': []})
        residue_atoms[corrected]['atoms'].append(atom_name)

    sorted_atoms = []
    for res_num, details in sorted(residue_atoms.items(), key=lambda x: len(x[1]['atoms'])):
        res_name = details['res_name']
        atom_names = details['atoms']
        mol2_path = Path(templates_dir) / f"{res_name}.mol2"
        if not mol2_path.exists():
            logger.warning("No template for %s: %s", res_name, mol2_path)
            continue
        order, _, _ = parse_mol2_template(str(mol2_path), atom_names)
        atom_names = sorted(atom_names,
                            key=lambda x: order.index(x) if x in order else len(order))
        sorted_atoms.extend((res_num, res_name, atom) for atom in atom_names)
    return sorted_atoms
