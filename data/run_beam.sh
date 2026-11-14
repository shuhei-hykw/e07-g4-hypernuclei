#!/bin/bash
# Perpendicular MIP beam tracks with the v20 signal model in real views
set -e
cd "$(dirname "$0")/.."
./run.py singletrack --name beam_eval --events 20 --tracks 25 --beam \
  --seed 11
./run.py images data/sim/beam_eval.root --events 0-19 --centered \
  --background random --bg-split eval --seed 5100 \
  --model output/signal_model_v20.json --name beam_eval_v20_dz4 -j 2
\cp -f data/sim/beam_eval_tracks.csv \
  data/images/beam_eval_v20_dz4/single_tracks.csv
# line comparison: data/run_beamlines.sh
echo ALLDONE
