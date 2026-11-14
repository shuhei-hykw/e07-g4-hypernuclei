set -e
M=output/image_model.json
./run.py images data/sim/eval_nagara.root --events 0-9 --name v19_synth_tune_dz4 --seed 1000 --model $M -j 5
./run.py images data/sim/eval_nagara.root --events 40-49 --name v19_synth_eval2_dz4 --seed 7000 --model $M -j 5
./run.py images data/sim/eval_nagara.root --events 0-99 --name v19_embed_nagara_dz4 --seed 3000 --model $M --background random --bg-split eval -j 6
./run.py images data/sim/eval_nagara.root --events 0-99 --name v19_embed_nagara_dz3.6 --dz 3.6 --seed 3000 --model $M --background random --bg-split eval -j 6
echo ALLDONE
