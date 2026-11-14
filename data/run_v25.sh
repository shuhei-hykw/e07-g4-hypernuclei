set -e
# v25 = the v20 signal model drawn at the full-scan pixel, 0.507 um/px
# (grain image kept in pixels, Geant4 geometry at the true scale).
# Same events, seeds and background views as the v20 embed sets
# (--bg-same-as), so old and new sets are paired event by event.
M=output/signal_model_v20.json
R=data/sim/eval_nagara.root
O="--model $M --background random --bg-split eval --pixel-um 0.507 -j 4"
run() {  # name events seed
  ./run.py images $R --events $2 --name v25_embed_$1 --seed $3 \
    --bg-same-as data/images/v20_embed_$1 ${=O}
}
run tune_dz4 100,101,102,104,108-113 9400
run evalA_dz4 114-123 9500
run nagara_dz4 0-99 3000
E=$(command ls data/images/v20_embed_evalU_dz4 | sed 's/^ev0*//; s/\/$//' \
  | paste -sd,)
run evalU_dz4 $E 9600
echo ALLDONE
