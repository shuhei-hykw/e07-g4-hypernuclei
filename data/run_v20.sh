set -e
M=output/signal_model_v20.json
./run.py images data/sim/single_eval.root --events 0-79 --centered --background random --bg-split eval --seed 5000 --model $M --name single_eval_v20_dz4 -j 8
\cp -f data/sim/single_eval_tracks.csv data/images/single_eval_v20_dz4/single_tracks.csv
~/work/.venv/bin/python -m module.trackprofile data/images/single_eval_v20_dz4 --events 10 --slice-step 6 -o output/trackprofile_v20.md
./run.py images data/sim/eval_nagara.root --events 0-99 --name v20_embed_nagara_dz4 --seed 3000 --model $M --background random --bg-split eval -j 8
echo ALLDONE
