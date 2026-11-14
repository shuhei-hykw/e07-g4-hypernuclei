for p in 1.0 1.4; do
  ./run.py images data/sim/eval_nagara.root --events 0-1 --name power${p}_tune_dz4 --seed 1000 --model data/trials/model_v19_power$p.json -j 2 &
done
wait
for p in 1.0 1.4; do
  echo "== $p"
  ~/work/.venv/bin/python -m module.depthsplit data/images/power${p}_tune_dz4/ev000{0,1} 2>&1 | grep "dz4/ev\|  compact\|fog, compact"
done
echo ALLDONE
