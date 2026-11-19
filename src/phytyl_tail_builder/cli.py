"""Command-line interfaces."""

import argparse
import logging
import sys

from . import __version__
from .main import build
from .optimizer import POTENTIALS
from .pdb import read_residues
from .stereo import REPORT_COLUMNS, report_rows

STEREO_FAILED = 3


def main(argv=None):
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
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.INFO if args.verbose else logging.WARNING,
                        format='%(levelname)s: %(message)s')
    try:
        output, rows = build(args.tleap_in, args.leap_log, args.templates, args.output_prefix,
                             fragment_potential=args.tail_potential,
                             environment_potential=args.environment_potential,
                             keep_initial_positions=args.keep_tleap_positions)
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


def check_main(argv=None):
    parser = argparse.ArgumentParser(
        prog="phytyl-tail-check",
        description="Report C8/C13 R/S, C2=C3 E/Z, closest contact and worst bond angle of every "
                    "phytyl tail in a PDB. Exit status 1 if any tail is not (R,R,E) or has a "
                    "bond angle above 130 degrees; partial tails fail only on what can be determined.",
    )
    parser.add_argument('pdb')
    parser.add_argument('--resnames', nargs='+', default=['CLA', 'CHL'],
                        help='residue names to check (default: CLA CHL)')
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')
    args = parser.parse_args(argv)

    try:
        residues = read_residues(args.pdb)
    except (OSError, ValueError) as e:
        print(f"error: {e}", file=sys.stderr)
        return 2
    rows = report_rows(residues, [i for i, r in enumerate(residues) if r['name'] in args.resnames])
    print("\t".join(REPORT_COLUMNS))
    for row in rows:
        print("\t".join(str(row[c]) for c in REPORT_COLUMNS))
    bad = [r for r in rows if r['ok'] == 'no']
    incomplete = sum(r['ok'] == 'incomplete' for r in rows)
    print(f"# {len(rows)} residues, {incomplete} incomplete, {len(bad)} failing", file=sys.stderr)
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main())
