# Examples

## minimal

A cut-out of the LHCII-CP29-CP24 model used in the paper: three chlorophylls
whose tails were rebuilt by tleap (CHL 663, CLA 669, CHL 719) and every residue
with an atom within 10 Å of them. Residue numbers are those of the full system.

| File | Content |
|---|---|
| `tleap.in` | tleap script of the full run, without the lipid `bond` lines; `$PARAMS` stands for the cofactor parameter folder |
| `leap.log` | `Added missing heavy atom` lines of the three residues |
| `LCC_input.pdb` | structure before tleap (`loadpdb`) |
| `LCC_no_sc.pdb` | structure written by tleap (`savepdb`) |
| `CHL.mol2`, `CLA.mol2` | residue templates |
| `expected/default/` | reference output with the default potentials |
| `expected/all-mie/` | reference output with `--tail-potential mie --keep-tleap-positions` |

```bash
cd examples/minimal
phytyl-tail-builder --tleap-in tleap.in --leap-log leap.log
```

This writes `LCC_no_sc_no_Hy.pdb`, which should match `expected/default/` to 1e-3 Å. With `--tail-potential mie --keep-tleap-positions` it matches `expected/all-mie/`.
Protein fragments start and end at the cut, so the amide H of each fragment's
first residue was removed and tleap adds terminal atoms (OXT, H1-H3) when the
output is loaded.

To build Amber topology from the output:

```bash
scripts/check_amber.sh examples/minimal/tleap.in examples/minimal/LCC_no_sc_no_Hy.pdb
```

The script loads ASP as an ordinary residue, because `savepdb` writes a TER before some ASP residues that are not chain starts, and the force field would otherwise treat them as N-terminal.

## tleap_rebuild.in

tleap script for the full system: loads `LCC_no_sc_no_Hy.pdb`, rebuilds the
hydrogens and writes `LCC.prmtop` / `LCC.inpcrd`. Replace `$PARAMS` with the
cofactor parameter folder before running it.

## Viewing the tails

`scripts/plot_tail.py` draws the tail of one residue before and after, from the same side view, with matplotlib:

```bash
cd examples/minimal
python ../../scripts/plot_tail.py LCC_no_sc.pdb expected/default/LCC_no_sc_no_Hy.pdb CLA.mol2 669 -o tail.png --name "CLA 669"
```

The rebuilt residues in this example are 663, 669 and 719 (CHL 663 and CHL 719 use `CHL.mol2`).
