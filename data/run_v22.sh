~/work/.venv/bin/python -m module.synth.calibrate --ridge --views 12 --iterations 200 --dz 4.0 --z-loss --only brightness,mottle_large,mottle_small,mottle_small_um,mottle_z_um,mottle_long_rel,mottle_long_z_um,min_transmission,noise_sigma,fixed_noise_sigma,beam_slope_max,beam_spread,grain_sigma_um,image_blur_px,grain_density,fog_per_um3,fog_weight,grain_spread,beam_per_mm2,random_per_mm2,random_ion_max,blobs_per_mm2,flat_per_mm2,flat_dip_max,flat_ion_max,flat_ion_min,speckle_rel,speckle_scale_px --output output/image_model_v22_dz4.json > data/calib_v22_dz4.log 2>&1 || exit 1
M=output/image_model_v22_dz4.json
./run.py images data/sim/eval_nagara.root --events 50-59 --name v22_synth_eval_dz4 --dz 4.0 --seed 8000 --model $M -j 5
./run.py validate data/images/v22_synth_eval_dz4 -o output/realism_v22_dz4.md
~/work/.venv/bin/python -m module.linecheck --real 20 --synth data/images/v22_synth_eval_dz4 data/images/v19_synth_eval_dz4 -o output/linecheck_v22.md
~/work/.venv/bin/python -m module.depthsplit data/images/v22_synth_eval_dz4/ev005{0,1} > data/depth_v22.log 2>&1
echo ALLDONE
