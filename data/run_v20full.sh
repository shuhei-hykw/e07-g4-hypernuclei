set -e
# v20 grain settings applied to the whole synthetic view (background
# tracks too), same events / seed as v19_synth_eval_dz4 for comparison
M=output/signal_model_v20.json
./run.py images data/sim/eval_nagara.root --events 30-39 --name v20_synth_eval_dz4 --seed 6000 --model $M -j 5
./run.py validate data/images/v20_synth_eval_dz4 -o output/realism_v20full_dz4.md
~/work/.venv/bin/python -m module.linecheck --real 20 --synth data/images/v20_synth_eval_dz4 data/images/v19_synth_eval_dz4 -o output/linecheck_v20full.md
echo ALLDONE
