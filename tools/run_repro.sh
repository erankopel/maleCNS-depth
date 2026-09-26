#!/usr/bin/env bash
# Phase 0 reproduction of the exploratory pipeline (README run order, steps 1 to 9).
# Touches none of H1 to H5. Logs go to logs/NN_<step>.log; stops at the first failure.
set -u
W=/home/claude/maleCNS-depth
export MALECNS_DATA=$W/data
export PYTHONUNBUFFERED=1
cd "$W"
mkdir -p logs

{
  echo "date_utc   $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "host       $(uname -srm)"
  echo "cpu        $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2 | xargs)  x $(nproc)"
  echo "mem        $(free -h | awk '/Mem:/{print $2}')"
  python3 - <<'EOF'
import sys, numpy, scipy, pandas, pyarrow
print('python    ', sys.version.split()[0])
for m in (numpy, scipy, pandas, pyarrow):
    print(f'{m.__name__:<10s}', m.__version__)
EOF
} > logs/00_environment.log 2>&1

run() {  # run <tag> <cmd...>
  local tag=$1; shift
  local log="logs/${tag}.log"
  echo "=== $tag: $* ===" > "$log"
  local t0=$(date +%s)
  python3 tools/peakrun.py "$@" >> "$log" 2>&1
  local rc=$?
  local t1=$(date +%s)
  echo "rc=$rc  wall=$((t1 - t0))s" >> "$log"
  echo "$(date -u +%H:%M:%S)  $tag  rc=$rc  wall=$((t1 - t0))s" >> logs/progress.log
  if [ $rc -ne 0 ]; then echo "STOPPED at $tag" >> logs/progress.log; exit $rc; fi
}

run 01_fetch_data        python3 code/fetch_data.py
run 02_build_graph       python3 code/build_graph.py
run 03_analyse           python3 code/analyse.py
run 04_spectra_scc       python3 code/spectra_scc.py
run 05_spectra_unsigned  python3 code/spectra.py unsigned
run 06_spectra_signed    python3 code/spectra.py signed
run 07_modes             python3 code/modes.py
run 08_slowmode_check    python3 code/slowmode_check.py
run 09_cert              python3 code/cert.py
echo "ALL DONE" >> logs/progress.log
