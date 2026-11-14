set -e
# Scoring of the v26 dev set (track weight 0.57) once the fullscan
# candidates and features have landed; fixed before they did.
D=$HOME/work/e07-fullscan/results/sim_eval
S=data/images/v26_embed_nagara_dz4
P=$HOME/work/.venv/bin/python
V=$D/v26_embed_nagara_dz4
BG=$D/v19_bgonly_nagara_dz4
[ -f $V.done ] && [ -f ${V}_features.parquet ] || { echo "not ready"; exit 1; }
./run.py unionstudy $S --a ${V}_methodA_production.csv \
  --b ${V}_bprime.csv --a-bg ${BG}_methodA_production.csv \
  --b-bg ${BG}_bprime.csv --per-event output/arms/v26_nagara_base.csv \
  2>&1 | grep -v Warning > output/unionstudy_v26_dz4.txt
F="--a ${V}_methodA_production.csv --b ${V}_bprime.csv"
F="$F --features ${V}_features.parquet"
R="--a-bg ${BG}_methodA_production.csv --b-bg ${BG}_bprime.csv"
R="$R --features-bg ${BG}_features.parquet"
$P -m module.candchance $S ${=F} ${=R} -o output/candchance_v26_nagara.md
$P -m module.candchance $S ${=F} ${=R} --drop mean_intens \
  -o output/candchance_v26_nagara_no_mean_intens.md
$P -m module.candleak $S ${=F} ${=R} -o output/candleak_v26_nagara.md
echo ALLDONE
