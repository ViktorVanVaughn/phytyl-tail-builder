#!/usr/bin/env bash
# Load a phytyl-tail-builder output in tleap, rebuild hydrogens and write Amber topology.
#
# usage: scripts/check_amber.sh TLEAP_IN PDB [OUTDIR]
#
# TLEAP_IN is the tleap script used for the original run. Its parameter paths
# are pointed at $PARAMS (default ~/parameters/force-fields/PSII_cofactors),
# loadpdb is pointed at PDB, and check, charge and saveamberparm are added.
# Everything is written to OUTDIR (default: a new folder next to PDB).

set -euo pipefail

if [ $# -lt 2 ]; then
    sed -n '2,10p' "$0" | sed 's/^# \{0,1\}//'
    exit 2
fi

TLEAP_IN=$1
PDB=$(realpath "$2")
OUTDIR=${3:-$(dirname "$PDB")/amber_check}
PARAMS=${PARAMS:-$HOME/parameters/force-fields/PSII_cofactors}

[ -f "$TLEAP_IN" ] || { echo "not found: $TLEAP_IN" >&2; exit 1; }
[ -f "$PDB" ] || { echo "not found: $PDB" >&2; exit 1; }
[ -d "$PARAMS" ] || { echo "parameter folder not found: $PARAMS (set PARAMS)" >&2; exit 1; }
command -v tleap >/dev/null || { echo "tleap not on PATH" >&2; exit 1; }

mkdir -p "$OUTDIR"
OUTDIR=$(realpath "$OUTDIR")
IN="$OUTDIR/check.in"
LOG="$OUTDIR/check.log"

# Drop the run's own check/charge/save/quit lines, then append ours.
sed -E \
    -e 's#^(\s*source\s+)\S*leaprc\.lipid21#\1leaprc.lipid21#' \
    -e 's#\S*/(frcmod\.ionsjc_tip3p)#\1#' \
    -e 's#\S*/PSII_cofactors/#'"$PARAMS"'/#' \
    -e 's#\$PARAMS#'"$PARAMS"'#g' \
    -e 's#^(\s*\w+\s*=\s*loadpdb\s+).*#\1'"$PDB"'#' \
    -e '/^\s*(check|charge|savepdb|saveamberparm|savemol2|quit)\b/d' \
    "$TLEAP_IN" > "$IN"

# savepdb writes a TER before ASP 442, which is not a chain start in the full
# system, but the force field maps a residue after a TER to its N-terminal form
# (NASP). Load ASP as the ordinary residue.
sed -i -E '0,/^\s*\w+\s*=\s*loadpdb\s/{s//addPdbResMap { { 0 "ASP" "ASP" } }\n&/}' "$IN"

# bond lines use the residue numbers of the original loadpdb file, while PDB
# is numbered as written by savepdb. Map each original residue to the PDB
# residue holding the coordinates of its first atom and renumber the bond lines.
PRE=$(sed -nE 's#^\s*\w+\s*=\s*loadpdb\s+(\S+).*#\1#p' "$TLEAP_IN" | head -1)
PRE="$(dirname "$TLEAP_IN")/$PRE"
if [ -f "$PRE" ] && grep -qE '^\s*bond\s' "$IN"; then
    python3 - "$PRE" "$PDB" "$IN" <<'PY'
import re
import sys

pre, post, script = sys.argv[1:]

def atoms(path):
    for line in open(path):
        if line.startswith(('ATOM', 'HETATM')):
            yield int(line[22:26]), line[30:54]

where = {xyz: res for res, xyz in atoms(post)}
new_number = {}
for res, xyz in atoms(pre):
    new_number.setdefault(res, where.get(xyz))

def renumber(m):
    res = int(m.group(2))
    if new_number.get(res) is None:
        sys.exit(f"no match for residue {res} of {pre} in {post}")
    return f"{m.group(1)}{new_number[res]}."

lines = open(script).read().splitlines(True)
lines = [re.sub(r'(\w+\.)(\d+)\.', renumber, l) if re.match(r'\s*bond\s', l) else l
         for l in lines]
open(script, 'w').writelines(lines)
PY
fi

cat >> "$IN" <<EOF

check LCC
charge LCC
saveamberparm LCC $OUTDIR/LCC.prmtop $OUTDIR/LCC.inpcrd
quit
EOF

(cd "$OUTDIR" && tleap -f "$IN" > "$LOG" 2>&1) || true

count() { grep -c -E "$1" "$LOG" || true; }
echo "tleap input:      $IN"
echo "tleap log:        $LOG"
grep -E '^Exiting LEaP' "$LOG" || echo "tleap did not exit cleanly"
echo "unknown residues: $(count 'Unknown residue')"
echo "added heavy atoms: $(count 'Added missing heavy atom')"
grep -E 'Leap added [0-9]+ missing atoms|^ +[0-9]+ (Heavy|H / lone pairs)' "$LOG" || true
grep -E 'Total (unperturbed )?charge' "$LOG" | tail -1 || true
for f in LCC.prmtop LCC.inpcrd; do
    if [ -s "$OUTDIR/$f" ]; then echo "$f: written"; else echo "$f: missing"; fi
done
[ -s "$OUTDIR/LCC.prmtop" ] && [ -s "$OUTDIR/LCC.inpcrd" ]
