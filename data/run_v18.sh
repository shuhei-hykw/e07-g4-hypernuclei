set -e
./run.py images data/sim/eval_nagara.root --events 20-29 --name v18_synth_eval_dz4 --dz 4.0 --seed 4000 --model output/image_model_v18_dz4.json -j 5
./run.py validate data/images/v18_synth_eval_dz4 -o output/realism_v18_dz4.md
~/work/.venv/bin/python -m module.linecheck --real 20 --synth data/images/v18_synth_eval_dz4 data/images/v15_synth_eval_dz4 -o output/linecheck_v18.md
echo ALLDONE
