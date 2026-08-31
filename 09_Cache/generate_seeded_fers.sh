#!/usr/bin/env bash
# generate_seeded_fers.sh
#
# Generate FERS baseband H5s for N_MC seeds x 3 scenarios (landing,
# takeoff, orbit360), injecting <randomseed>SEED</randomseed> into
# each fersxml so every run has an independent noise realisation.
#
# Output layout:
#   seeds/seed_NNN/direct_<scenario>.h5
#   seeds/seed_NNN/echo_<scenario>.h5
#
# Usage:
#   ./generate_seeded_fers.sh [N_MC] [JOBS]
# Defaults: N_MC=3, JOBS=12
set -euo pipefail

N_MC="${1:-3}"
JOBS="${2:-12}"

FERS_BIN="/home/itumeleng/Documents/Academia/MscEng/FERS/build/src/fers"
ROOT="/home/itumeleng/Documents/Academia/MscEng/FM-PR-TrackingFilters"

if [[ ! -x "$FERS_BIN" ]]; then
  echo "FERS binary not found or not executable: $FERS_BIN" >&2
  exit 1
fi

mkdir -p "$ROOT/seeds"

run_one() {
  local seed="$1" scen_name="$2" scen_file="$3"
  local seed_pad; seed_pad=$(printf "%03d" "$seed")
  local seed_dir="$ROOT/seeds/seed_$seed_pad"
  local work="$seed_dir/$scen_name"

  # Idempotent: skip if both H5s already exist.
  if [[ -s "$seed_dir/direct_${scen_name}.h5" && -s "$seed_dir/echo_${scen_name}.h5" ]]; then
    echo "  [skip] seed=$seed_pad scen=$scen_name (cached)"
    return 0
  fi

  mkdir -p "$work"
  ln -sfn "$ROOT/Waveform" "$work/Waveform"

  # Resolve scenario path (flightScenarios or BackupScenarios).
  local scen_path="$ROOT/FERS/flightScenarios/$scen_file.fersxml"
  if [[ ! -f "$scen_path" ]]; then
    scen_path="$ROOT/FERS/BackupScenarios/$scen_file.fersxml"
  fi

  # Inject <randomseed> just before </parameters> (deterministic per seed).
  awk -v seed="$seed" '
    /<\/parameters>/ && !inj { print "    <randomseed>" seed "</randomseed>"; inj = 1 }
    { print }
  ' "$scen_path" > "$work/run.fersxml"

  ( cd "$work" && "$FERS_BIN" run.fersxml > fers.log 2>&1 )

  # Receiver names vary by scenario (direct.h5 for single-target,
  # direct_N3.h5 for the 3-target); glob for whatever FERS produced.
  local direct_h5 echo_h5
  direct_h5=$(find "$work" -maxdepth 1 -name '*direct*.h5' 2>/dev/null | head -1)
  echo_h5=$(find "$work" -maxdepth 1 -name '*echo*.h5'   2>/dev/null | head -1)
  if [[ -z "$direct_h5" || -z "$echo_h5" ]]; then
    echo "  [FAIL] seed=$seed_pad scen=$scen_name: h5 output missing (see $work/fers.log)"
    return 1
  fi

  mv "$direct_h5" "$seed_dir/direct_${scen_name}.h5"
  mv "$echo_h5"   "$seed_dir/echo_${scen_name}.h5"
  rm -rf "$work"

  echo "  [ok] seed=$seed_pad scen=$scen_name"
}
export -f run_one
export FERS_BIN ROOT

tasks() {
  for s in $(seq 1 "$N_MC"); do
    echo "$s landing  scenario_2_landingManeuver"
    echo "$s takeoff  scenario_3_takeoffManeuver"
    echo "$s orbit360 scenario_4_360"
    echo "$s 3targets scenario_N3_targets"
  done
}

total=$((N_MC * 4))
echo "Launching $total FERS runs, $JOBS in parallel..."
t0=$(date +%s)

tasks | xargs -n 3 -P "$JOBS" bash -c 'run_one "$@"' _

t1=$(date +%s)
echo "Done: $total runs in $((t1 - t0))s"
