#!/usr/bin/env bash
# Timing campaign backing the search-cost figures reported in Section 4 of the
# paper. Three stages:
#   A: W_5..W_12, default seed 20260831, sequential single-process baselines.
#   B: W_13, eight seeds (1..8) in parallel, stop at the first verified labeling.
#   C: W_14, eight seeds (1..8) x 18 h (64800 s) each, run to the limit.
# Every run's full stdout/stderr is kept next to this script; one summary row
# per run is appended to summary.tsv as soon as that run is decided.
set -u

ROOT="/home/junyeop/projects/geodesic-leech-wheels"
OUT="$ROOT/search/timing-runs"
SCRIPT="$ROOT/search/search_wheel_labelings.py"
SUMMARY="$OUT/summary.tsv"
mkdir -p "$OUT"

{
  echo "# campaign started: $(date -u +%FT%TZ)"
  echo "# host: $(uname -sr); cpu: $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2- | xargs); nproc=$(nproc)"
  echo "# python: $(python3 --version 2>&1)"
  printf 'stage\twheel\tseed\tstatus\twall_s\tdetail\n'
} >> "$SUMMARY"

row() { printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "$4" "$5" "$6" >> "$SUMMARY"; }
final_line() { grep -E "^W_$1 " "$2" 2>/dev/null | head -1; }

# ---------------- Stage A: baselines, default seed, one process at a time ----
for n in 5 6 7 8 9 10 11 12; do
  log="$OUT/stageA_W${n}.log"
  t0=$(date +%s)
  nice -n 10 python3 -u "$SCRIPT" -n "$n" --seed 20260831 --time-limit 3600 --verbose >"$log" 2>&1
  t1=$(date +%s)
  if grep -qE "^W_${n} +PASS" "$log"; then st=PASS; else st=not_found; fi
  row A "W_$n" 20260831 "$st" "$((t1 - t0))" "$(final_line "$n" "$log")"
done

# ---------------- Stage B: W_13 race, 8 seeds, first PASS wins --------------
declare -a PID SEED
t0=$(date +%s)
i=0
for s in 1 2 3 4 5 6 7 8; do
  log="$OUT/w13_seed${s}.log"
  nice -n 10 python3 -u "$SCRIPT" -n 13 --seed "$s" --time-limit 21600 --verbose >"$log" 2>&1 &
  PID[$i]=$!; SEED[$i]=$s; i=$((i + 1))
done

winner=""
live=8
while [ "$live" -gt 0 ] && [ -z "$winner" ]; do
  sleep 5
  for j in 0 1 2 3 4 5 6 7; do
    [ -z "${PID[$j]}" ] && continue
    s=${SEED[$j]}
    if grep -qE "^W_13 +PASS" "$OUT/w13_seed${s}.log" 2>/dev/null; then
      winner=$s
      break
    fi
    if ! kill -0 "${PID[$j]}" 2>/dev/null; then
      t_end=$(date +%s)
      row B W_13 "$s" not_found "$((t_end - t0))" "$(final_line 13 "$OUT/w13_seed${s}.log")"
      PID[$j]=""; live=$((live - 1))
    fi
  done
done
t1=$(date +%s)
if [ -n "$winner" ]; then
  row B W_13 "$winner" PASS "$((t1 - t0))" "$(final_line 13 "$OUT/w13_seed${winner}.log")"
  for j in 0 1 2 3 4 5 6 7; do
    [ -z "${PID[$j]}" ] && continue
    s=${SEED[$j]}
    if [ "$s" != "$winner" ]; then
      kill "${PID[$j]}" 2>/dev/null
      row B W_13 "$s" stopped_after_winner "$((t1 - t0))" ""
    fi
  done
else
  row B W_13 - none_within_limit "$((t1 - t0))" ""
fi
wait 2>/dev/null

# ---------------- Stage C: W_14, 8 seeds x 18 h ------------------------------
t0=$(date +%s)
for s in 1 2 3 4 5 6 7 8; do
  log="$OUT/w14_seed${s}.log"
  nice -n 10 python3 -u "$SCRIPT" -n 14 --seed "$s" --time-limit 64800 --verbose >"$log" 2>&1 &
done
wait
t1=$(date +%s)
for s in 1 2 3 4 5 6 7 8; do
  log="$OUT/w14_seed${s}.log"
  if grep -qE "^W_14 +PASS" "$log"; then st=PASS; else st=not_found; fi
  row C W_14 "$s" "$st" "$((t1 - t0))" "$(final_line 14 "$log")"
done

echo "# campaign finished: $(date -u +%FT%TZ)" >> "$SUMMARY"
echo "CAMPAIGN COMPLETE"
