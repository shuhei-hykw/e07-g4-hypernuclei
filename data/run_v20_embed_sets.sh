set -e
# v20-signal embed sets for the fullscan point-placement study; fresh
# events of 100-299 not in eval4, eval-split real backgrounds
M=output/signal_model_v20.json
./run.py images data/sim/eval_nagara.root --events 100,101,102,104,108-113 --name v20_embed_tune_dz4 --seed 9400 --model $M --background random --bg-split eval -j 5
./run.py images data/sim/eval_nagara.root --events 114-123 --name v20_embed_evalA_dz4 --seed 9500 --model $M --background random --bg-split eval -j 5
echo ALLDONE
