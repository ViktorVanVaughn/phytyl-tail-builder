"""Command-line interface."""

import argparse
import logging
import sys

from . import __version__
from .main import build
from .optimizer import POTENTIALS


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
        output = build(args.tleap_in, args.leap_log, args.templates, args.output_prefix,
                       args.tail_potential, args.environment_potential,
                       args.keep_tleap_positions)
    except (OSError, ValueError) as e:
        logging.error("%s", e)
        return 1
    print(f"Output: {output}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
