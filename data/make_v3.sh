#!/bin/bash
# Datasets with image model v3 (output/image_model.json)
set -e
cd /home/hayakawa/work/e07-g4-hypernuclei
J=20
./run.py images data/sim/eval_nagara.root --events 0-4 -j 5 --background random --name v3_embed_nagara5_dz3.0 &
./run.py images data/sim/eval_nagara.root --events 0-4 -j 5 --background random --dz 1.5 --name v3_embed_nagara5_dz1.5 &
./run.py images data/sim/eval_nagara.root --events 0-4 -j 5 --name v3_synth_nagara5_dz3.0 &
./run.py images data/sim/eval_nagara.root --events 0-4 -j 5 --dz 1.5 --name v3_synth_nagara5_dz1.5 &
wait
echo CHECKSETS_DONE
./run.py images data/sim/eval_nagara.root --events 0-99 -j $J --background random --name v3_embed_nagara_dz3.0
./run.py images data/sim/eval_nagara.root --events 0-99 -j $J --background random --dz 1.5 --name v3_embed_nagara_dz1.5
./run.py images data/sim/eval_ll11be.root --events 0-49 -j $J --background random --name v3_embed_ll11be_dz3.0
./run.py images data/sim/eval_ibuki.root --events 0-49 -j $J --background random --name v3_embed_ibuki_dz3.0
./run.py images data/sim/nagara_replay.root --events 0-19 -j $J --background random --name v3_embed_nagara_replay_dz3.0
./run.py images data/sim/ibuki_replay.root --events 0-19 -j $J --background random --name v3_embed_ibuki_replay_dz3.0
./run.py images data/sim/nagara_replay.root --events 0-19 -j $J --name v3_synth_nagara_replay_dz3.0
./run.py images data/sim/ibuki_replay.root --events 0-19 -j $J --name v3_synth_ibuki_replay_dz3.0
echo ALL_DONE
