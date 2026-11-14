for s in 0.034 0.028 0.023; do
  ./run.py images data/sim/eval_nagara.root --events 0-1 --name defocus${s}_tune_dz4 --seed 1000 --model data/trials/model_v19_defocus$s.json -j 2 &
done
wait
for s in 0.034 0.028 0.023; do
  echo "== $s"
  ~/work/.venv/bin/python -m module.depthsplit data/images/defocus${s}_tune_dz4/ev000{0,1} 2>&1 | grep "dz4/ev\|  compact\|fog, compact"
done
echo ALLDONE
