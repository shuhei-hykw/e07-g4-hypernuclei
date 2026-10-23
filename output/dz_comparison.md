# Physical focal step: calibration under dz 3.0 and 1.5 um

Same target statistics (6 real tohoku-v1 views, central 1024 px, slices at offsets -2..2, 5, 13), same start (model v3) and 250 iterations each (`./run.py calibrate --dz {3.0,1.5} --exclude n_elong_per_mpx`). Residuals: relative (log for counts, absolute for correlations / fractions).

| statistic | real | dz 3.0 | dz 1.5 | residual 3.0 | residual 1.5 |
|---|---|---|---|---|---|
| median | 203.3 | 231 | 229 | +0.136 | +0.126 |
| std | 23.36 | 25.1 | 24.55 | +0.075 | +0.051 |
| p1 | 145.5 | 161 | 154 | +0.107 | +0.058 |
| p5 | 163.3 | 180 | 180 | +0.102 | +0.102 |
| n_cc_per_mpx | 1.391e+04 | 1.505e+04 | 1.538e+04 | +0.079 | +0.101 |
| cc_area_p90 | 28 | 29 | 31 | +0.035 | +0.102 |
| dark_blobs_per_mpx | 3524 | 3416 | 3409 | -0.031 | -0.033 |
| dark_blob_area_median | 3 | 3 | 3 | +0.000 | +0.000 |
| dark_contrast_median | 0.1941 | 0.1472 | 0.18 | -0.047 | -0.014 |
| dark_contrast_p90 | 0.3245 | 0.3524 | 0.3124 | +0.028 | -0.012 |
| corr1 | 0.9268 | 0.9279 | 0.9145 | +0.001 | -0.012 |
| corr2 | 0.8337 | 0.8262 | 0.8358 | -0.008 | +0.002 |
| corr5 | 0.5637 | 0.516 | 0.5346 | -0.048 | -0.029 |
| corr13 | 0.1875 | 0.1872 | 0.1793 | -0.000 | -0.008 |
| fg_fraction | 0.336 | 0.1725 | 0.1819 | -0.163 | -0.154 |

| parameter | dz 3.0 | dz 1.5 | ratio 1.5/3.0 |
|---|---|---|---|
| brightness | 242.7 | 234.1 | 0.96 |
| mottle_large | 0.1118 | 0.1118 | 1.00 |
| mottle_small | 0.002514 | 0 | 0.00 |
| min_transmission | 0.07029 | 0.02685 | 0.38 |
| noise_sigma | 1.661 | 2.549 | 1.53 |
| fixed_noise_sigma | 1.924 | 2 | 1.04 |
| beam_slope_max | 0.01195 | 0.006999 | 0.59 |
| beam_spread | 0.007251 | 0.001545 | 0.21 |
| grain_sigma_um | 0.2449 | 0.2616 | 1.07 |
| defocus_slope | 0.06095 | 0.09721 | 1.59 |
| grain_density | 0.2 | 0.2319 | 1.16 |
| fog_per_um3 | 0.02776 | 0.02776 | 1.00 |
| grain_spread | 0.3486 | 0.423 | 1.21 |
| beam_per_mm2 | 1e+04 | 3396 | 0.34 |
| random_per_mm2 | 647 | 1952 | 3.02 |
| random_ion_max | 3.273 | 4 | 1.22 |
| blobs_per_mm2 | 3.963 | 3.835 | 0.97 |

Result: the current search could not distinguish the two assumptions.
Both reach similar losses (0.163 / 0.143) and match the slice
correlation at offsets 1-13, but both leave large residuals in the same
places, notably the foreground fraction after binarization (real 0.336,
both ~0.18) and the intensity percentiles; the model is therefore not
complete enough for a decision. In addition the slice statistics
constrain mainly products of the step with the beam angular spread and
the defocus blur. Independent inputs (beam density per mm2, objective NA,
beam divergence, measured layer thickness and shrinkage) are needed.
