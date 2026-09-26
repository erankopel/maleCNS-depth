#!/usr/bin/env bash
# Download syn-partners (6.8 GB) and verify against the md5 that GCS reports for the object.
set -u
U=https://storage.googleapis.com/flyem-male-cns/v1.0/connectome-data/flat-connectome/syn-partners-male-cns-v1.0-minconf-0.5.feather
OUT=/home/claude/maleCNS-depth/data/syn-partners.feather
LOG=/home/claude/maleCNS-depth/logs/10_fetch_synpartners.log
{
  echo "=== syn-partners download $(date -u +%Y-%m-%dT%H:%M:%SZ) ==="
  curl -sI "$U" | grep -iE "content-length|x-goog-hash|last-modified|x-goog-generation"
  WANT=$(curl -sI "$U" | tr -d '\r' | awk -F'md5=' '/x-goog-hash: md5=/{print $2}' | python3 -c "import sys,base64; print(base64.b64decode(sys.stdin.read().strip()).hex())")
  echo "gcs md5 (hex): $WANT"
  t0=$(date +%s)
  curl -sS --retry 5 --retry-delay 5 -C - -o "$OUT.part" "$U"; rc=$?
  echo "curl rc=$rc  wall=$(( $(date +%s) - t0 ))s  bytes=$(stat -c %s "$OUT.part")"
  GOT=$(md5sum "$OUT.part" | cut -d' ' -f1)
  echo "local md5    : $GOT"
  if [ "$GOT" = "$WANT" ]; then mv "$OUT.part" "$OUT"; echo "OK md5 match -> $OUT"; else echo "MD5 MISMATCH; kept $OUT.part"; fi
} > "$LOG" 2>&1
