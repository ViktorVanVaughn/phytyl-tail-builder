#!/usr/bin/env python3
"""Rewrite the residue numbers of the bond lines of a tleap script for another PDB file.

    python scripts/renumber_bonds.py PRE.pdb POST.pdb SCRIPT

In `bond u.N.ATOM`, tleap reads N as the residue's sequence number in the loaded file: the resSeq of
its first residue plus the position. SCRIPT was written for PRE (the original loadpdb file); its bond
lines are rewritten in place for POST (written from PRE by tleap and phytyl-tail-builder), whose
residues are matched to those of PRE by atom coordinates. check_amber.sh calls it.
"""
import re
import sys


def residues(path):
    """(first resSeq, [coordinate strings of each residue]), residues split as tleap splits them."""
    first, out, key = None, [], None
    with open(path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                if line[17:27] != key:
                    key = line[17:27]
                    out.append([])
                    if first is None:
                        first = int(re.match(r"-?\d+", line[22:27].strip()).group())
                out[-1].append(line[30:54])
            elif line.startswith("TER"):
                key = None
    return first, out


def renumber(pre, post, script):
    first_pre, pre_res = residues(pre)
    first_post, post_res = residues(post)
    where = {xyz: j for j, atoms in enumerate(post_res) for xyz in atoms}

    def new_number(m):
        i = int(m.group(2)) - first_pre
        found = {where[xyz] for xyz in pre_res[i] if xyz in where} if 0 <= i < len(pre_res) else set()
        if len(found) != 1:
            sys.exit(f"bond line residue {m.group(2)}: {len(found)} matching residues in {post}")
        return f"{m.group(1)}{first_post + found.pop()}."

    with open(script) as f:
        lines = f.readlines()
    lines = [re.sub(r"(\w+\.)(\d+)\.", new_number, line) if re.match(r"\s*bond\s", line, re.I) else line
             for line in lines]
    with open(script, "w") as f:
        f.writelines(lines)


if __name__ == "__main__":
    renumber(*sys.argv[1:4])
