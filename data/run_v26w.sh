set -e
# Track darkness scan: v25 (v20 model at 0.507 um/px) with the signal
# grain darkness scaled by --track-weight, first 30 nagara events, same
# seeds and backgrounds as v25 (paired).
M=output/signal_model_v20.json
R=data/sim/eval_nagara.root
O="--model $M --background random --bg-split eval --pixel-um 0.507 -j 4"
for w in 0.57 0.75; do
  ./run.py images $R --events 0-29 --name v26w${w#0.}_embed_nagara30_dz4 \
    --seed 3000 --bg-same-as data/images/v20_embed_nagara_dz4 \
    --track-weight $w ${=O}
done
echo ALLDONE
