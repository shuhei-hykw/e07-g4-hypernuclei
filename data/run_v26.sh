set -e
# v26 = v25 with the signal track grains at weight 0.57 (= 0.289 /
# 0.507: track darkness per pixel of length as calibrated, matches the
# vertex darkness of 70 eye-judged real full-scan vertices,
# output/vertex_darkness_scan_v26w.txt). Same events, seeds and
# backgrounds as v25 (paired); the real-view inputs are reused.
M=output/signal_model_v20.json
R=data/sim/eval_nagara.root
O="--model $M --background random --bg-split eval --pixel-um 0.507"
O="$O --track-weight 0.57 -j 4"
run() {  # name events seed
  ./run.py images $R --events $2 --name v26_embed_$1 --seed $3 \
    --bg-same-as data/images/v25_embed_$1 ${=O}
}
run nagara_dz4 0-99 3000
run tune_dz4 100,101,102,104,108-113 9400
run evalA_dz4 114-123 9500
./run.py images data/sim/evalV_nagara.root --events 0-99 \
  --name v26_embed_evalV_dz4 --seed 9700 \
  --bg-same-as data/images/v25_embed_evalV_dz4 ${=O}
echo ALLDONE
