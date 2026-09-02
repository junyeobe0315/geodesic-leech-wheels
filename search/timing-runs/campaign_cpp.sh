#!/usr/bin/env bash
# Stage C': W_14 with the C++ engine (search/wheel_search_cpp), eight seeds,
# 18 h (64800 s) each, run to the limit. Replaces the CPython stage C of
# campaign.sh, which was stopped ~25 min in once the C++ port validated at
# ~87x CPython step throughput (see the notes in summary.tsv).
set -u
ROOT="/home/junyeop/projects/geodesic-leech-wheels"
OUT="$ROOT/search/timing-runs"
BIN="$ROOT/search/wheel_search_cpp"
SUMMARY="$OUT/summary.tsv"

echo "# stage C' (C++ engine) started: $(date -u +%FT%TZ)" >> "$SUMMARY"
t0=$(date +%s)
for s in 1 2 3 4 5 6 7 8; do
  nice -n 10 "$BIN" -n 14 --seed "$s" --time-limit 64800 --verbose > "$OUT/w14_cpp_seed${s}.log" 2>&1 &
done
wait
t1=$(date +%s)
for s in 1 2 3 4 5 6 7 8; do
  log="$OUT/w14_cpp_seed${s}.log"
  if grep -qE "^W_14 +PASS" "$log"; then st=PASS; else st=not_found; fi
  printf 'Ccpp\tW_14\t%s\t%s\t%s\t%s\n' "$s" "$st" "$((t1 - t0))" \
    "$(grep -E '^W_14 ' "$log" | head -1)" >> "$SUMMARY"
done
echo "# stage C' finished: $(date -u +%FT%TZ)" >> "$SUMMARY"
echo "CPP CAMPAIGN COMPLETE"
