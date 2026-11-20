# Changelog

## 0.2.1

### Fixed
- Contacts are now part of the stereo check: a tail with a heavy atom within 2.0 Å of another heavy atom fails
  the gate and `phytyl-tail-check`.
- The beam kept partial tails that differed by a grid step or two, so a fine grid could end in a clash (with
  `--dihedral-step 5` the example ended with a 1.67 Å contact). It now keeps at most one partial tail per 30°
  cell of the free dihedrals.

### Changed
- The minimal example is set up with the corrected cofactor parameters: the chlorophyll b template is
  (7R,11R) and the C2=C3 torsion `X -c8-c5-X` has phase 180. New `tleap.in` (flat parameter folder, load order
  ff14SB first), `leap.log` (the full tleap session), `LCC_no_sc.pdb` (tleap's placement with the corrected
  library), `CHL.mol2` and `CLA.mol2`; reference outputs regenerated. `examples/tleap_rebuild.in` loads the same
  set.

## 0.2.0

Tails are now built with natural stereochemistry and checked, and the clash score no longer has a step. The
v0.1.0 output of the minimal example had inverted centres (CLA 669 and CHL 719 7S,11S; CHL 663 11S), twisted
C2=C3 bonds (C1-C2-C3-C5 = 100° and −137°) and bond angles up to 180°.

### Added
- Torsion-space builder, now the default (`--builder torsion`): NeRF placement with template bond lengths and
  angles, scanned backbone dihedrals (O2A–C1, C1–C2, C3–C5 onwards) with a beam search (`--beam`,
  `--dihedral-step`), C1-C2-C3-C5 fixed at 180°, and branches on the side that makes C8 and C13 R. Resolved
  atoms are kept. A second pass (`--passes`, default 2) rebuilds every tail with all the others in place.
- `stereo` module: C8/C13 R/S from the signed volume, calibrated on the tail of CLA.mol2 (a centre within
  0.05 Å³ of planar gets no label); C2=C3 E/Z; tail angles more than 20° from their values in natural phytyl;
  tail bonds more than 0.1 Å off the template.
- Stereo gate: `stereo_report.tsv` for every run; exit status 3 if a rebuilt tail is not (R,R,E) or has an
  angle off natural phytyl (`--allow-stereo-errors` to only warn).
- `--fix-resolved-stereo`: rebuild resolved atoms from an inverted centre (or non-E C2=C3) outward.
- `phytyl-tail-check PDB [--resnames ...]`: the same check on any PDB, for use as a QC gate (exit status 1 if a
  tail fails, 2 if the PDB cannot be read or has no phytyl residue). It finds phytyl residues by their atoms
  (CGA, O2A, C1–C20 in PDB names: CLA, CHL, BCL, PHO and the other chlorophylls and pheophytins).
- `scripts/renumber_bonds.py` (used by `check_amber.sh`) and `scripts/options_table.py`, which generates the
  README options tables from the argument parsers; a test fails when the README differs.

### Changed
- Clash score: each pair is a Lennard-Jones or Mie term shifted to 0 at a fixed 4.5 Å cutoff, and the 1-2 and
  1-3 partners of the new atom are left out. v0.1.0 widened the cutoff with the number of atoms within 3.4 Å
  and weighted every term, bonded ones included, by cutoff/distance, so one atom crossing 3.4 Å could cost
  more than a real clash; it also floored scores at 0.1 and used σ = 3.4 Å for every element. σ is now
  element-specific (C, N, O, H, Mg, S, P).
- The v0.1 algorithm is kept as `--builder sphere`, with the new score and the fixes below, so its output
  differs from v0.1.0.
- `build()` returns `(path, report_rows)`.
- Residues are identified by their position in the PDB file, not by resSeq.
- `examples/minimal/leap.log` is now the log of tleap run on the cut-out itself: 12 rebuilt chlorophylls,
  numbered as tleap numbers that file. New reference outputs are in `expected/torsion*` and `expected/sphere*`.
  The v0.1.0 outputs in `expected/default` and `expected/all-mie` are unchanged.
- `scripts/check_amber.sh`: the parameter folder (`-p`, else `$PARAMS`), unit name (`-u`, else read from
  `loadpdb`) and residue maps (`-r ASP`) are arguments, and it prints the parameter folder it uses.

### Fixed
- leap.log residues were looked up by PDB resSeq and paired by order (`zip`), so one unmatched residue
  moved atoms in the wrong chlorophyll; on a log numbered as tleap writes it, nothing was moved at all.
  Residues are now mapped by sequential index (tleap counts from the first resSeq of the `loadpdb` file) and
  verified by name and coordinates, and an added atom that the `loadpdb` residue already has is refused; any
  mismatch is an error.
- leap.log: only the last tleap session is read, for atoms and residue names alike, with a warning when the
  file holds several (sessions were never split before).
- `parse_tleap_in` matched commented lines and was case sensitive.
- mol2 bonds were assumed to point parent → child; they are now oriented by breadth-first search.
- The first added atom kept its tleap position; it is now placed like the others.
- A residue with a phytyl ester but no template was skipped with a warning, its tail left where tleap put it;
  this is now an error.
- resSeq from 10000 up (five digits in tleap's savepdb) was cut to four digits in messages.
- `scripts/check_amber.sh` renumbered `bond` lines by PDB resSeq. tleap numbers residues by sequence (first
  resSeq plus position), so with gaps the bonds joined other residues. It also skipped the renumbering without
  a message when `loadpdb` named an absolute or quoted path.

### Known template and parameter issues (documented, not shipped)
- `CHL.mol2` and `CHL.lib` (PSII_cofactors/data/CLX, chlorophylls/) are (7S,11S); the torsion builder still
  builds R,R.
- `CLX.frcmod` in PSII_cofactors/data/CLX: `X -c8-c5-X` has phase 0, which puts the C2=C3 minimum at 90°; it
  should be 180, as in the chlorophylls/ copy.
