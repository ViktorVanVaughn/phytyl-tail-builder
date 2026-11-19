#!/usr/bin/env bash
# Load a phytyl-tail-builder output in tleap, rebuild hydrogens and write Amber topology.
#
# usage: scripts/check_amber.sh [-p PARAMS] [-u UNIT] [-r RES]... TLEAP_IN PDB [OUTDIR]
#
#   -p PARAMS  cofactor parameter folder (default: $PARAMS); replaces $PARAMS
#              and any .../PSII_cofactors/ path in TLEAP_IN
#   -u UNIT    unit name (default: the variable of the loadpdb line in TLEAP_IN)
#   -r RES     load residue RES as its ordinary (non-terminal) form even after a
#              TER, e.g. -r ASP; may be repeated
#
# TLEAP_IN is the tleap script used for the original run. loadpdb is pointed
# at PDB, and check, charge and saveamberparm are added. Everything is written
# to OUTDIR (default: a new folder next to PDB).

set -euo pipefail

usage() { sed -n '2,15p' "$0" | sed 's/^# \{0,1\}//'; exit 2; }

UNIT=
RESMAP=()
while getopts "p:u:r:h" opt; do
    case $opt in
        p) PARAMS=$OPTARG ;;
        u) UNIT=$OPTARG ;;
        r) RESMAP+=("$OPTARG") ;;
        *) usage ;;
    esac
done
shift $((OPTIND - 1))
[ $# -ge 2 ] || usage

TLEAP_IN=$1
PDB=$(realpath "$2")
OUTDIR=${3:-$(dirname "$PDB")/amber_check}

[ -f "$TLEAP_IN" ] || { echo "not found: $TLEAP_IN" >&2; exit 1; }
[ -f "$PDB" ] || { echo "not found: $PDB" >&2; exit 1; }
[ -n "${PARAMS:-}" ] || { echo "no parameter folder: use -p PARAMS or set PARAMS" >&2; exit 1; }
[ -d "$PARAMS" ] || { echo "parameter folder not found: $PARAMS" >&2; exit 1; }
PARAMS=$(realpath "$PARAMS")
command -v tleap >/dev/null || { echo "tleap not on PATH" >&2; exit 1; }
if [ -z "$UNIT" ]; then
    UNIT=$(sed -nE 's/^\s*(\w+)\s*=\s*loadpdb\s.*/\1/ip' "$TLEAP_IN" | tail -1)
    [ -n "$UNIT" ] || { echo "no loadpdb line in $TLEAP_IN (use -u UNIT)" >&2; exit 1; }
fi
echo "parameters:       $PARAMS"
echo "unit:             $UNIT"

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
    -e 's#^(\s*\w+\s*=\s*loadpdb\s+).*#\1'"$PDB"'#I' \
    -e '/^\s*(check|charge|savepdb|saveamberparm|savemol2|quit)\b/Id' \
    "$TLEAP_IN" > "$IN"

# savepdb can write a TER before a residue that is not a chain start in the
# full system, and the force field then maps it to its N-terminal form (e.g.
# NASP). -r RES loads RES as the ordinary residue.
if [ ${#RESMAP[@]} -gt 0 ]; then
    MAP=""
    for r in "${RESMAP[@]}"; do MAP+="{ 0 \"$r\" \"$r\" } "; done
    sed -i -E '0,/^\s*\w+\s*=\s*loadpdb\s/I{s//addPdbResMap { '"$MAP"'}\n&/}' "$IN"
fi

# bond lines number residues as tleap numbered the original loadpdb file (its
# first resSeq plus the position); renumber_bonds.py rewrites them for PDB.
PRE=$(sed -nE 's#^\s*\w+\s*=\s*loadpdb\s+(\S+).*#\1#Ip' "$TLEAP_IN" | tail -1)
PRE=${PRE//\"/}
PRE=${PRE//\'/}
case $PRE in /*) ;; *) PRE="$(dirname "$TLEAP_IN")/$PRE" ;; esac
if grep -qiE '^\s*bond\s' "$IN"; then
    [ -f "$PRE" ] || { echo "bond lines need the loadpdb file of $TLEAP_IN: $PRE not found" >&2; exit 1; }
    python3 "$(dirname "$0")/renumber_bonds.py" "$PRE" "$PDB" "$IN"
fi

cat >> "$IN" <<EOF

check $UNIT
charge $UNIT
saveamberparm $UNIT $OUTDIR/$UNIT.prmtop $OUTDIR/$UNIT.inpcrd
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
for f in "$UNIT.prmtop" "$UNIT.inpcrd"; do
    if [ -s "$OUTDIR/$f" ]; then echo "$f: written"; else echo "$f: missing"; fi
done
[ -s "$OUTDIR/$UNIT.prmtop" ] && [ -s "$OUTDIR/$UNIT.inpcrd" ]
