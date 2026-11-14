set -e
# v23: background recalibrated with the v20 track appearance fixed
# (grain appearance, density, defocus and image blur not moved)
~/work/.venv/bin/python -m module.synth.calibrate --ridge --views 12 --iterations 200 --dz 4.0 --z-loss --only brightness,mottle_large,mottle_small,mottle_small_um,mottle_z_um,mottle_long_rel,mottle_long_z_um,min_transmission,noise_sigma,fixed_noise_sigma,beam_slope_max,beam_spread,fog_per_um3,fog_weight,beam_per_mm2,random_per_mm2,random_ion_max,blobs_per_mm2,flat_per_mm2,flat_dip_max,flat_ion_max,flat_ion_min,speckle_rel,speckle_scale_px --output output/image_model_v23_dz4.json > data/calib_v23_dz4.log 2>&1
M=output/image_model_v23_dz4.json
./run.py images data/sim/eval_nagara.root --events 30-39 --name v23_synth_eval_dz4 --seed 6000 --model $M -j 5
./run.py validate data/images/v23_synth_eval_dz4 -o output/realism_v23_dz4.md
~/work/.venv/bin/python -m module.linecheck --real 20 --synth data/images/v23_synth_eval_dz4 data/images/v19_synth_eval_dz4 -o output/linecheck_v23.md
echo ALLDONE
