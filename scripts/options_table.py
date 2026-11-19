#!/usr/bin/env python3
"""Print the README options tables from the argument parsers, so the README cannot drift from --help.

    python scripts/options_table.py

tests/test_readme.py checks that README.md holds exactly this output between its options-table markers.
"""

from phytyl_tail_builder.cli import build_check_parser, build_parser


def rows(parser):
    for action in parser._actions:          # argparse has no public list of its arguments
        name = ", ".join(action.option_strings) or action.dest
        if action.option_strings and action.nargs != 0:
            meta = action.metavar or action.dest.upper()
            if action.choices:
                meta = "{" + ",".join(action.choices) + "}"
            name += " " + meta + (" ..." if action.nargs == "+" else "")
        text = action.help.replace("|", "/").replace("<", "&lt;").replace(">", "&gt;")   # GitHub drops <tags>
        if action.required and action.option_strings:
            text += " (required)"
        yield f"| `{name}` | {text} |"


def table():
    out = []
    for parser in (build_parser(), build_check_parser()):
        out += [f"`{parser.prog}`: {parser.description}", "", "| Option | Meaning |", "|---|---|"]
        out += [*rows(parser), ""]
    return "\n".join(out).rstrip() + "\n"


if __name__ == "__main__":
    print(table(), end="")
