set -e
M=output/image_model_v15_dz4.json
./run.py singletrack --name single_eval --events 80 --tracks 25 --seed 11
./run.py images data/sim/single_eval.root --events 0-79 --centered --background random --bg-split eval --dz 4.0 --seed 5000 --model $M --name single_eval_dz4 -j 6
\cp -f data/sim/single_eval_tracks.csv data/images/single_eval_dz4/single_tracks.csv
./run.py images data/sim/single_eval.root --events 0-79 --centered --background random --bg-split eval --dz 3.6 --seed 5000 --model $M --name single_eval_dz3.6 -j 6
\cp -f data/sim/single_eval_tracks.csv data/images/single_eval_dz3.6/single_tracks.csv
echo ALLDONE
