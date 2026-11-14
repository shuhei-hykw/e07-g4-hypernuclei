set -e
M=output/image_model_v15_dz4.json
./run.py images data/sim/eval_nagara.root --events 0-9 --name v15_synth_tune_dz4 --dz 4.0 --seed 1000 --model $M -j 5
./run.py images data/sim/eval_nagara.root --events 10-19 --name v15_synth_eval_dz4 --dz 4.0 --seed 2000 --model $M -j 5
./run.py images data/sim/eval_nagara.root --events 0-9 --name v15_embed_nagara_dz4 --dz 4.0 --seed 3000 --model $M --background random -j 5
./run.py validate data/images/v15_synth_eval_dz4 -o output/realism_v15.md > data/validate_v15.log 2>&1
echo ALLDONE
