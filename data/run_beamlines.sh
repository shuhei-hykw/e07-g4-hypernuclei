#!/bin/bash
# Steep-line (beam) comparison: real backgrounds, v20 perpendicular MIP
# embedded into real views, whole-synthetic views
cd "$(dirname "$0")/.."
PY=~/work/.venv/bin/python
$PY -m module.beamlines data/images/single_eval_v20_dz4 --events 20 \
  -o output/beamlines_real.md 2>&1 | grep -v Warning
$PY -m module.beamlines data/images/beam_eval_v20_dz4 --events 20 \
  -o output/beamlines_beam_v20.md 2>&1 | grep -v Warning
$PY -m module.beamlines data/images/v19_synth_eval_dz4 --events 10 \
  --scramble 4 -o output/beamlines_v19.md 2>&1 | grep -v Warning
$PY -m module.beamlines data/images/v24_synth_eval_dz4 --events 10 \
  --scramble 4 -o output/beamlines_v24.md 2>&1 | grep -v Warning
echo ALLDONE
