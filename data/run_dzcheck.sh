M=output/image_model.json
for dz in 3.4 3.6; do
  ./run.py images data/sim/eval_nagara.root --events 50-59 --name v19_synth_eval_dz$dz --dz $dz --seed 8000 --model $M -j 5
  ./run.py validate data/images/v19_synth_eval_dz$dz -o output/realism_v19_dz$dz.md
  ~/work/.venv/bin/python -m module.depthsplit data/images/v19_synth_eval_dz$dz/ev005{0,1} > data/depth_v19_dz$dz.log 2>&1
done
echo ALLDONE
