set -e
# Held-out pair eval4: the 10 events of 100-299 (unused by any image
# set) with the most capture / decay charged tracks shorter than 15 um
# (4.7 per event vs 3.2 on average); model v19
M=output/image_model.json
E=162,106,125,135,202,238,283,103,105,107
./run.py images data/sim/eval_nagara.root --events $E --name v19_synth_eval4_dz4 --seed 9200 --model $M -j 4
./run.py images data/sim/eval_nagara.root --events $E --name v19_embed_eval4_dz4 --seed 9300 --model $M --background random --bg-split eval -j 4
echo ALLDONE
