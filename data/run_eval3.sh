set -e
# Fresh held-out sets for the fullscan point-cloud check (events 60-69
# are unused by every other synthetic set; model v19)
M=output/image_model.json
./run.py images data/sim/eval_nagara.root --events 60-69 --name v19_synth_eval3_dz4 --seed 9000 --model $M -j 4
./run.py images data/sim/eval_nagara.root --events 60-69 --name v19_embed_eval3_dz4 --seed 9100 --model $M --background random --bg-split eval -j 4
echo ALLDONE
