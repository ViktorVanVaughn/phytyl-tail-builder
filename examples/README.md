# Examples

## minimal

A cut-out of the LHCII-CP29-CP24 model used in the paper: three chlorophylls
(CHL 663, CLA 669, CHL 719) and every residue with an atom within 10 Å of
them. Residue numbers in the PDB files are those of the full system. Twelve
chlorophylls of the cut-out have truncated tails, and all twelve are rebuilt.
In two of them (CHL 440, CHL 666) the resolved C8 is already inverted (7S).

| File | Content |
|---|---|
| `tleap.in` | tleap script run on `LCC_input.pdb`; `$PARAMS` stands for the cofactor parameter folder, with the corrected chlorophyll b template (7R,11R) and C2=C3 torsion (`X -c8-c5-X`, phase 180), `.lib` and `.frcmod` files side by side |
| `leap.log` | the log of that tleap run (AmberTools 25) from its `loadpdb` line on, with the path of the teLeap binary cut to its name; residues numbered as tleap numbers `LCC_input.pdb` (from its first resSeq, 16) |
| `LCC_input.pdb` | structure before tleap (`loadpdb`) |
| `LCC_no_sc.pdb` | structure written by tleap (`savepdb`), renumbered as in the full system, without the terminal atoms (OXT, H1-H3) tleap added to the cut protein fragments |
| `CHL.mol2`, `CLA.mol2` | residue templates of that parameter set, both (2E,7R,11R) |
| `expected/torsion/` | reference output of v0.2 with the defaults (exits with status 3: CHL 440, 666) |
| `expected/torsion-fix-resolved/` | reference output with `--fix-resolved-stereo` |
| `expected/sphere/` | reference output with `--builder sphere --allow-stereo-errors` |
| `expected/sphere-all-mie/` | the same with `--tail-potential mie --keep-tleap-positions` |
| `expected/default/`, `expected/all-mie/` | v0.1.0 output (three residues, from the earlier full-system leap.log), kept unchanged as the record of the stereo errors fixed in v0.2 |

Each v0.2 folder also holds the `stereo_report.tsv` of that run.

```bash
cd examples/minimal
phytyl-tail-builder --tleap-in tleap.in --leap-log leap.log --fix-resolved-stereo
phytyl-tail-check LCC_no_sc_no_Hy.pdb
```

This writes `LCC_no_sc_no_Hy.pdb` and `stereo_report.tsv`, which match `expected/torsion-fix-resolved/` to 1e-3 Å. Without `--fix-resolved-stereo` the run matches `expected/torsion/`, but the stereo gate stops it with status 3 because of the two resolved inversions.
Protein fragments start and end at the cut, so the amide H of each fragment's
first residue was removed and tleap adds terminal atoms (OXT, H1-H3) when the
output is loaded.

To build Amber topology from the output:

```bash
scripts/check_amber.sh -p <parameter folder> -r ASP \
    examples/minimal/tleap.in examples/minimal/LCC_no_sc_no_Hy.pdb
```

`-p` is the cofactor parameter folder (default: `$PARAMS`), and the script prints the folder it uses. The unit name is read from the `loadpdb` line (override with `-u`). `-r ASP` loads ASP as an ordinary residue, because `savepdb` writes a TER before some ASP residues that are not chain starts, and the force field would otherwise treat them as N-terminal.

## tleap_rebuild.in

tleap script for the full system: loads `LCC_no_sc_no_Hy.pdb`, rebuilds the
hydrogens and writes `LCC.prmtop` / `LCC.inpcrd`. Replace `$PARAMS` with the
corrected cofactor parameter folder before running it.

## Viewing the tails

`scripts/plot_tail.py` draws the tail of one residue before and after, from the same side view, with matplotlib:

```bash
cd examples/minimal
python ../../scripts/plot_tail.py LCC_no_sc.pdb expected/torsion-fix-resolved/LCC_no_sc_no_Hy.pdb CLA.mol2 669 -o tail.png --name "CLA 669"
```

The rebuilt residues in this example are 220, 440, 662, 663, 666, 669, 672, 675, 711, 715, 717 and 719 (CHL residues use `CHL.mol2`).
