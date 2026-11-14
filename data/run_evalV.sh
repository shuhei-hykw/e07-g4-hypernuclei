set -e
# Fresh frozen gate set at the full-scan pixel (v25 model): 100 new
# Geant4 events (data/sim/evalV_nagara.root, seed 31) and eval-split
# real views used by no existing image dataset, one per event.
EX=$(command ls -d data/images/*/ | paste -sd,)
./run.py images data/sim/evalV_nagara.root --events 0-99 \
  --name v25_embed_evalV_dz4 --seed 9700 \
  --model output/signal_model_v20.json --pixel-um 0.507 \
  --background random --bg-split eval --bg-exclude $EX --bg-unique -j 4
echo ALLDONE
