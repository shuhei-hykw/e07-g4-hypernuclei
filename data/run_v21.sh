until grep -q "^-> \|Traceback" data/calib_v21_dz4.log; do sleep 60; done
grep -q Traceback data/calib_v21_dz4.log && exit 1
M=output/image_model_v21_dz4.json
./run.py images data/sim/eval_nagara.root --events 30-39 --name v21_synth_eval_dz4 --dz 4.0 --seed 6000 --model $M -j 5
./run.py validate data/images/v21_synth_eval_dz4 -o output/realism_v21_dz4.md
~/work/.venv/bin/python -m module.linecheck --real 20 --synth data/images/v21_synth_eval_dz4 data/images/v19_synth_eval_dz4 -o output/linecheck_v21.md
~/work/.venv/bin/python -m module.depthsplit data/images/v21_synth_eval_dz4/ev003{0,1} > data/depth_v21.log 2>&1
echo ALLDONE
