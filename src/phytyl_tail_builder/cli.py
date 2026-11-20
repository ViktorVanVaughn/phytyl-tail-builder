"""Command-line interfaces."""

import argparse
import logging
import sys

from . import __version__
from .main import BUILDERS, build
from .optimizer import POTENTIALS
from .pdb import read_residues
from .stereo import REPORT_COLUMNS, has_phytyl, report_rows

STEREO_FAILED = 3


def build_parser():
    parser = argparse.ArgumentParser(
        prog="phytyl-tail-builder",
        description="Reposition tail atoms added by tleap and write a hydrogen-free PDB.",
    )
    parser.add_argument('--tleap-in', required=True, help='tleap input script that was run')
    parser.add_argument('--leap-log', required=True, help='leap.log written by that run')
    parser.add_argument('--templates', default=None,
                        help='folder with <RESNAME>.mol2 templates (default: folder of --tleap-in)')
    parser.add_argument('-o', '--output-prefix', default='.',
                        help='output folder (default: current directory)')
    parser.add_argument('--builder', choices=BUILDERS, default='torsion',
                        help="torsion: template bond lengths and angles, branches and C2=C3 fixed "
                             "so the tail is (2E,7R,11R), backbone dihedrals scanned; "
                             "sphere: the v0.1 direction search, with the v0.2 fixes (default: torsion)")
    parser.add_argument('--beam', type=int, default=5,
                        help='torsion builder: partial tails kept at each dihedral (default: 5)')
    parser.add_argument('--dihedral-step', type=float, default=10.0,
                        help='torsion builder: dihedral grid step in degrees (default: 10)')
    parser.add_argument('--passes', type=int, default=2,
                        help="torsion builder: 1 builds each tail once, in order; each further pass "
                             "rebuilds every tail with all the other tails in place (default: 2)")
    parser.add_argument('--fix-resolved-stereo', action='store_true',
                        help="torsion builder: discard resolved tail atoms from a wrong stereocentre "
                             "(or a non-E C2=C3) outward and rebuild them")
    parser.add_argument('--allow-stereo-errors', action='store_true',
                        help=f"only warn when a rebuilt tail fails the stereo check "
                             f"(default: exit with status {STEREO_FAILED})")
    parser.add_argument('--tail-potential', choices=POTENTIALS, default='lj',
                        help="potential between a new atom and the tail atoms already placed "
                             "(lj: Lennard-Jones 12-6, mie: Mie 16-6; default: lj)")
    parser.add_argument('--environment-potential', choices=POTENTIALS, default='mie',
                        help="potential between a new atom and everything else "
                             "(default: mie)")
    parser.add_argument('--keep-tleap-positions', action='store_true',
                        help="keep the tleap positions of the residue's added atoms in the "
                             "environment while the tail is built")
    parser.add_argument('-v', '--verbose', action='store_true', help='verbose logging')
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.fix_resolved_stereo and args.builder != 'torsion':
        parser.error('--fix-resolved-stereo needs --builder torsion')

    logging.basicConfig(level=logging.INFO if args.verbose else logging.WARNING,
                        format='%(levelname)s: %(message)s')
    try:
        output, rows = build(args.tleap_in, args.leap_log, args.templates, args.output_prefix,
                             fragment_potential=args.tail_potential,
                             environment_potential=args.environment_potential,
                             keep_initial_positions=args.keep_tleap_positions, builder=args.builder,
                             beam=args.beam, dihedral_step=args.dihedral_step, passes=args.passes,
                             fix_resolved_stereo=args.fix_resolved_stereo)
    except (OSError, ValueError) as e:
        logging.error("%s", e)
        return 1
    print(f"Output: {output}")
    failed = [f"{r['resname']} {r['resseq']}" for r in rows if r['ok'] != 'yes']
    if failed:
        message = f"stereo check failed for {', '.join(failed)} (see stereo_report.tsv)"
        if not args.allow_stereo_errors:
            logging.error("%s", message)
            return STEREO_FAILED
        logging.warning("%s", message)
    return 0


def build_check_parser():
    parser = argparse.ArgumentParser(
        prog="phytyl-tail-check",
        description="Report C8/C13 R/S, C2=C3 E/Z, closest contact and worst bond angle of every "
                    "phytyl tail in a PDB. Exit status 1 if any tail is not (R,R,E) or has a "
                    "bond angle more than 20 degrees from natural phytyl or a heavy atom within 2.0 A of "
                    "another (partial tails fail only on what can be determined), "
                    "2 if the PDB cannot be read or has none of the residues.",
    )
    parser.add_argument('pdb', help='PDB file to check')
    parser.add_argument('--resnames', nargs='+', default=None,
                        help='residue names to check (default: every residue with a phytyl ester, '
                             'atoms CGA and O2A/C1-C20 as in CLA, CHL, BCL, PHO)')
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')
    return parser


def check_main(argv=None):
    args = build_check_parser().parse_args(argv)
    try:
        residues = read_residues(args.pdb)
    except (OSError, ValueError) as e:
        print(f"error: {e}", file=sys.stderr)
        return 2
    indices = [i for i, r in enumerate(residues)
               if (r['name'] in args.resnames if args.resnames else has_phytyl(r['atoms']))]
    if not indices:
        wanted = " or ".join(args.resnames) if args.resnames else "phytyl"
        print(f"error: no {wanted} residue in {args.pdb}", file=sys.stderr)
        return 2
    rows = report_rows(residues, indices)
    print("\t".join(REPORT_COLUMNS))
    for row in rows:
        print("\t".join(str(row[c]) for c in REPORT_COLUMNS))
    bad = [r for r in rows if r['ok'] == 'no']
    incomplete = sum(r['ok'] == 'incomplete' for r in rows)
    print(f"# {len(rows)} residues, {incomplete} incomplete, {len(bad)} failing", file=sys.stderr)
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main())
