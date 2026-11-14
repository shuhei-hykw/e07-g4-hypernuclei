#!/bin/bash
# B' border margin on the dev set (v20_embed_nagara_dz4); evalU frozen
cd "$(dirname "$0")/.."
D=~/work/e07-fullscan/results/sim_eval
for B in 10 25; do
  ./run.py unionstudy data/images/v20_embed_nagara_dz4 \
    --a $D/v20_embed_nagara_dz4_methodA_production.csv \
    --b $D/v20_embed_nagara_dz4_bprime.csv \
    --a-bg $D/v19_bgonly_nagara_dz4_methodA_production.csv \
    --b-bg $D/v19_bgonly_nagara_dz4_bprime.csv --border-px $B 2>&1 \
    | grep -v Warning > output/unionstudy_v20_dz4_border$B.txt &
done
wait
echo ALLDONE
