#!/usr/bin/env bash
# Confirmatory runs of H1 to H5 (PREREG Sec. 5; DEVIATIONS.md C1 to C8). Each hypothesis is run once.
# Usage: MALECNS_DATA=/abs/path/to/data bash code/confirmatory/run_confirmatory.sh
# Every attempt logs to a new logs/confirmatory_<UTC date>_attempt<k>/ ; results go to results/*.json
# and are never overwritten. A hypothesis whose scored result exists is not run again; if only its
# unscored companion is missing, the script is rerun in MALECNS_UNSCORED_ONLY mode, which must
# reproduce the saved scored result exactly (C6).
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../.." && pwd)
export MALECNS_DATA=$(realpath "${MALECNS_DATA:-$ROOT/data}")
export MALECNS_RESULTS=$(realpath -m "${MALECNS_RESULTS:-$ROOT/results}")
export PYTHONUNBUFFERED=1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1      # bit-reproducible runs (C8)
k=1
while [ -e "$ROOT/logs/confirmatory_$(date -u +%Y-%m-%d)_attempt$k" ]; do k=$((k + 1)); done
LOG=$ROOT/logs/confirmatory_$(date -u +%Y-%m-%d)_attempt$k
mkdir -p "$LOG" "$MALECNS_RESULTS"
cd "$HERE"

{
  echo "date_utc   $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "attempt    $k"
  echo "data       $MALECNS_DATA"
  echo "results    $MALECNS_RESULTS"
  echo "host       $(uname -srm)"
  echo "cpu        $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2 | xargs)  x $(nproc)"
  echo "mem        $(free -h | awk '/Mem:/{print $2}')"
  echo "threads    OPENBLAS=$OPENBLAS_NUM_THREADS OMP=$OMP_NUM_THREADS MKL=$MKL_NUM_THREADS"
} > "$LOG/00_environment.log" 2>&1

stop() { echo "STOPPED: $1" | tee -a "$LOG/progress.log"; exit 1; }

# code identity: every executed file must match the manifest committed with the frozen tag (C6)
( cd "$ROOT" && sha256sum -c code/confirmatory/MANIFEST.sha256 ) >> "$LOG/00_environment.log" 2>&1 \
  || stop "code does not match code/confirmatory/MANIFEST.sha256"
# environment identity (C8)
python3 - >> "$LOG/00_environment.log" 2>&1 <<'EOF' || stop "package versions differ from requirements-confirmatory.txt"
import sys, importlib
want = dict(l.strip().split('==') for l in open('requirements-confirmatory.txt') if '==' in l)
bad = 0
print('python    ', sys.version.split()[0], '(required 3.11.15)')
bad += sys.version.split()[0] != '3.11.15'
for name, ver in want.items():
    got = importlib.import_module(name).__version__
    print(f'{name:<10s} {got}  (required {ver})')
    bad += got != ver
sys.exit(1 if bad else 0)
EOF

run() {  # run <tag> <cmd...>
  local tag=$1; shift
  echo "=== $tag: $* ===" > "$LOG/$tag.log"
  local t0=$(date +%s)
  python3 "$ROOT/tools/peakrun.py" "$@" >> "$LOG/$tag.log" 2>&1
  local rc=$?
  echo "rc=$rc  wall=$(( $(date +%s) - t0 ))s" >> "$LOG/$tag.log"
  echo "$(date -u +%H:%M:%S)  $tag  rc=$rc" >> "$LOG/progress.log"
  [ $rc -eq 0 ] || stop "$tag exited with $rc"
}

hyp() {  # hyp <tag> <script> <scored result> <unscored companion or ->
  local tag=$1 script=$2 scored=$3 comp=$4
  if [ -e "$MALECNS_RESULTS/$scored.json" ]; then
    if [ "$comp" = "-" ] || [ -e "$MALECNS_RESULTS/$comp.json" ]; then
      echo "$(date -u +%H:%M:%S)  $tag  already run; not repeated (C6)" >> "$LOG/progress.log"
      return
    fi
    echo "$(date -u +%H:%M:%S)  $tag  scored result exists; completing the unscored part (C6)" >> "$LOG/progress.log"
    export MALECNS_UNSCORED_ONLY=1
    run "$tag" python3 "$script"
    unset MALECNS_UNSCORED_ONLY
  else
    run "$tag" python3 "$script"
  fi
}

run 01_fetch_inputs   python3 "$ROOT/code/fetch_data.py"
run 02_build_graph    python3 "$ROOT/code/build_graph.py"
run 03_verify_inputs  python3 verify_inputs.py
run 04_prep_neuropil  python3 prep_neuropil.py
hyp 05_h1_gap         h1_gap.py         h1  h1_sensitivity
hyp 06_h2_floor       h2_floor.py       h2  h2_sensitivity
hyp 07_h3_depth       h3_depth.py       h3  h3_full
hyp 08_h4_confidence  h4_confidence.py  h4  h4_full
hyp 09_h5_sign        h5_sign.py        h5  -
sha256sum "$MALECNS_RESULTS"/*.json > "$LOG/results.sha256"
python3 -c "
import json, glob
for f in sorted(glob.glob('$MALECNS_RESULTS/h[1-5].json')):
    r = json.load(open(f)); print(r['hypothesis'], r['verdict'])
" > "$LOG/verdicts.txt"
cat "$LOG/verdicts.txt" >> "$LOG/progress.log"
echo "ALL DONE" >> "$LOG/progress.log"
