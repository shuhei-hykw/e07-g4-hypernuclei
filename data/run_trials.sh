set -e
for m in data/trials/model_k8_s0.35.json data/trials/model_k12_s0.35.json data/trials/model_k8_s0.35_g0.35.json; do
  t=$(basename $m .json)
  ./run.py images data/sim/single_tune.root --events 0-7 --centered --background random --bg-split train --seed 8000 --model $m --name trial_$t -j 8 > /dev/null 2>&1
  \cp -f data/sim/single_tune_tracks.csv data/images/trial_$t/single_tracks.csv
  ~/work/.venv/bin/python -m module.trackprofile data/images/trial_$t --events 8 --slice-step 6 -o data/trials/profile_$t.md > /dev/null 2>&1
  echo "done $t"
done
echo ALLDONE
