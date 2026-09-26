#!/usr/bin/env bash
# Convert every SLD mini-DST Parquet file under $SLD_PARQUET to EDM4hep under $SLD_OUT,
# running $JOBS containers in parallel and validating each output against its source.
#
#   scripts/convert_all.sh            # all files
#   JOBS=2 scripts/convert_all.sh mc  # only the mc subdirectory
#
# Outputs: $SLD_OUT/{data,mc}/<basename>.edm4hep.root, logs in $SLD_OUT/logs/.
# Existing outputs with a matching .ok marker are skipped.
set -euo pipefail

SLD_PARQUET=${SLD_PARQUET:-/media/sqian/WD_RED/phizero_jets/data/sld_v3}
SLD_OUT=${SLD_OUT:-/media/sqian/WD_RED/sld-edm4hep/edm4hep}
JOBS=${JOBS:-3}
REPO=$(cd "$(dirname "$0")/.." && pwd)
SUBDIRS=${*:-data mc}

mkdir -p "$SLD_OUT/logs"
list=$(mktemp)
for sub in $SUBDIRS; do
  for f in "$SLD_PARQUET/$sub"/*.parquet; do
    base=$(basename "$f" .parquet)
    out="$SLD_OUT/$sub/$base.edm4hep.root"
    [ -f "$out.ok" ] && continue
    echo "$sub $f $out" >> "$list"
  done
done
n=$(wc -l < "$list")
echo "$n files to convert with $JOBS parallel jobs -> $SLD_OUT"

one() {
  sub=$1; f=$2; out=$3
  base=$(basename "$out" .edm4hep.root)
  log="$SLD_OUT/logs/$base.log"
  mkdir -p "$(dirname "$out")"
  if "$REPO/scripts/run_key4hep.sh" python -m sld_edm4hep.convert --quiet "$f" "$out" > "$log" 2>&1 \
     && "$REPO/scripts/run_key4hep.sh" python scripts/validate.py "$out" "$f" >> "$log" 2>&1; then
    touch "$out.ok"; echo "ok   $sub/$base"
  else
    echo "FAIL $sub/$base (see $log)"
  fi
}
export -f one
export REPO SLD_OUT

xargs -P "$JOBS" -L 1 bash -c 'one "$@"' _ < "$list"
rm -f "$list"
echo "done: $(ls "$SLD_OUT"/*/*.ok 2>/dev/null | wc -l) outputs validated"
