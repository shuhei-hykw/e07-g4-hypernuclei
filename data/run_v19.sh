set -e
./run.py images data/sim/eval_nagara.root --events 30-39 --name v19_synth_eval_dz4 --dz 4.0 --seed 6000 --model output/image_model_v19_dz4.json -j 5
./run.py validate data/images/v19_synth_eval_dz4 -o output/realism_v19_dz4.md
~/work/.venv/bin/python -m module.linecheck --real 20 --synth data/images/v19_synth_eval_dz4 data/images/v15_synth_eval_dz4 -o output/linecheck_v19.md
echo ALLDONE
