#!/bin/sh
# One-shot gate of the alen arm on v25_embed_evalV_dz4, fixed before the
# evalV candidates were opened. Decision: ./run.py armcompare rule
# (union capture 12 um inclusive diff >= 0, lower 95 % > -5 %, B'-only
# rescues and LL event capture not falling), plus G of both arms.
D=$HOME/work/e07-fullscan/results/sim_eval
S=data/images/v25_embed_evalV_dz4
O=output/arms
for f in evalV_dz4.done evalV_dz4_alen.done; do
  [ -f "$D/v25_embed_$f" ] || { echo "missing $f"; exit 1; }
done
[ -f "$D/v25_embed_evalV_dz4_real_alen_methodA_production.csv" ] \
  || { echo "real background runs missing"; exit 1; }
./run.py unionstudy $S \
  --a $D/v25_embed_evalV_dz4_methodA_production.csv \
  --b $D/v25_embed_evalV_dz4_bprime.csv \
  --a-bg $D/v25_embed_evalV_dz4_real_methodA_production.csv \
  --b-bg $D/v25_embed_evalV_dz4_real_bprime.csv \
  --per-event $O/v25_evalV_base.csv 2>&1 | grep -v Warning \
  > $O/v25_evalV_base.txt
./run.py unionstudy $S \
  --a $D/v25_embed_evalV_dz4_alen_methodA_production.csv \
  --b $D/v25_embed_evalV_dz4_bprime.csv \
  --a-bg $D/v25_embed_evalV_dz4_real_alen_methodA_production.csv \
  --b-bg $D/v25_embed_evalV_dz4_real_bprime.csv \
  --per-event $O/v25_evalV_alen.csv 2>&1 | grep -v Warning \
  > $O/v25_evalV_alen.txt
./run.py armcompare $O/v25_evalV_base.csv $O/v25_evalV_alen.csv \
  > $O/armcompare_evalV.txt
