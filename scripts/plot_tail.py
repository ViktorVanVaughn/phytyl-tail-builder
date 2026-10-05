#!/usr/bin/env python3
"""Draw one rebuilt tail before and after phytyl-tail-builder, seen from the side.

The residue is drawn heavy atoms only, with bonds taken from its mol2 template. The
atoms that moved between the two structures are the tail. The camera looks along the
plane of the ring system, so the ring is a vertical line and the tail is seen leaving it.

    python scripts/plot_tail.py before.pdb after.pdb RESNAME.mol2 RESNUM -o tail.png

Needs numpy, scipy and matplotlib.
"""

import argparse
import re

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from scipy.spatial import cKDTree

CLASH, NEAR = 2.6, 2.9
RING = re.compile(r"^(MG|N[ABCD]|C[1-4][ABCD]|CH[ABCD])$")  # macrocycle atoms of chlorophylls
TAIL = re.compile(r"^C\d+$")  # C1-C20 of the phytyl chain
ORANGE, TEAL, GREY, RED = "#ED780D", "#007380", "#9EA3AD", "#D91A1A"


def read_pdb(path):
    atoms = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                name = line[12:16].strip()
                if name.lstrip("0123456789").startswith("H"):
                    continue
                xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                atoms.append((name, int(line[22:26]), xyz))
    return atoms


def read_bonds(path):
    names, bonds, section = {}, [], None
    with open(path) as fh:
        for line in fh:
            if line.startswith("@<TRIPOS>"):
                section = line.strip()
                continue
            p = line.split()
            if section == "@<TRIPOS>ATOM" and len(p) > 1:
                names[p[0]] = p[1]
            elif section == "@<TRIPOS>BOND" and len(p) > 2:
                bonds.append((names[p[1]], names[p[2]]))
    return bonds


def residue(atoms, resnum):
    own = {n: x for n, r, x in atoms if r == resnum}
    env_x = np.array([x for n, r, x in atoms if r != resnum])
    return own, env_x


def camera(own_b, own_a, tail):
    ring = np.array([x for n, x in own_b.items() if RING.match(n)])
    c = ring.mean(0)
    normal = np.linalg.svd(ring - c)[2][2]
    pts = np.array([own_b[n] for n in tail] + [own_a[n] for n in tail]) - c
    pts -= np.outer(pts @ normal, normal)
    u = np.linalg.svd(pts - pts.mean(0))[2][0]
    depth = np.cross(normal, u)
    return c, np.array([u, normal, depth])


def draw(ax, own, env, bonds, tail, cam, lim, label, sub):
    c, M = cam
    pos = {n: M @ (x - c) for n, x in own.items()}
    for a, b in bonds:
        if a in pos and b in pos:
            t = a in tail or b in tail
            ax.plot([pos[a][0], pos[b][0]], [pos[a][1], pos[b][1]], color=ORANGE if t else TEAL,
                    lw=3.2 if t else 1.4, solid_capstyle="round", zorder=3 if t else 2)
    tree = cKDTree(env)
    for n in tail:
        for i in tree.query_ball_point(own[n], NEAR):
            q = M @ (env[i] - c)
            close = np.linalg.norm(env[i] - own[n]) < CLASH
            ax.scatter(*q[:2], s=34, color=GREY, zorder=4, linewidths=0)
            if close:
                ax.plot([pos[n][0], q[0]], [pos[n][1], q[1]], color=RED, lw=1.3, ls=(0, (2, 1.5)), zorder=5)
    ax.set_xlim(*lim[0])
    ax.set_ylim(*lim[1])
    ax.set_aspect("equal")
    ax.axis("off")
    ax.text(0.0, 1.0, label, transform=ax.transAxes, fontsize=16, fontweight="bold", va="top")
    ax.text(0.09, 1.0, sub, transform=ax.transAxes, fontsize=11, va="top", linespacing=1.35)


def count_contacts(own, env, tail):
    tree = cKDTree(env)
    return sum(len(tree.query_ball_point(own[n], CLASH)) for n in tail)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("before")
    ap.add_argument("after")
    ap.add_argument("template")
    ap.add_argument("resnum", type=int)
    ap.add_argument("-o", "--output", default="tail.png")
    ap.add_argument("--name", default=None, help="label for the residue (default: template name and number)")
    args = ap.parse_args()

    own_b, env_b = residue(read_pdb(args.before), args.resnum)
    own_a, env_a = residue(read_pdb(args.after), args.resnum)
    tail = [n for n in own_a if n in own_b and TAIL.match(n)]
    if not tail:
        raise SystemExit("no phytyl atoms (C1-C20) found in this residue")
    bonds = read_bonds(args.template)
    cam = camera(own_b, own_a, tail)
    c, M = cam
    allp = np.array([M @ (x - c) for x in list(own_b.values()) + list(own_a.values())])
    pad = 2.0
    lim = ((allp[:, 0].min() - pad, allp[:, 0].max() + pad), (allp[:, 1].min() - pad, allp[:, 1].max() + pad))

    w = 180 / 25.4
    h = w / 2 * (lim[1][1] - lim[1][0]) / (lim[0][1] - lim[0][0])
    fig, axes = plt.subplots(1, 2, figsize=(w, h + 1.1))
    nb, na = count_contacts(own_b, env_b, tail), count_contacts(own_a, env_a, tail)
    draw(axes[0], own_b, env_b, bonds, tail, cam, lim, "a",
         f"tleap tail (as built)\n{nb} contacts < {CLASH} Å")
    draw(axes[1], own_a, env_a, bonds, tail, cam, lim, "b",
         f"phytyl-tail-builder\n{na} contacts < {CLASH} Å")
    handles = [
        Line2D([], [], color=ORANGE, lw=5),
        Line2D([], [], color=TEAL, lw=2.5),
        Line2D([], [], marker="o", ls="", color=GREY, ms=9),
        Line2D([], [], color=RED, lw=2.5, ls=(0, (2, 1.5))),
    ]
    labels = [
        f"Phytyl tail ({args.name or args.resnum})",
        "Ring system (edge-on)",
        f"Neighbouring atoms (< {NEAR} Å)",
        f"Contacts (< {CLASH} Å)",
    ]
    fig.legend(handles, labels, loc="lower center", ncol=2, fontsize=9, frameon=False, columnspacing=1.8)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.17, wspace=0.04)
    fig.savefig(args.output, dpi=300)


if __name__ == "__main__":
    main()
