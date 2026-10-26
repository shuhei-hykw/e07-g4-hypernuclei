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

## Grains along perpendicular beam tracks (`./run.py beamgrains`)

A handle that is not degenerate with the beam divergence and the defocus
blur: a MIP track has a fixed grain density along its length (14.8 /
100 um assumed; this value is a measurement of e07-fullscan's own image
processing, STATUS 14.77 with px 0.29 and grains merged within 15 px, not
an external calibration), i.e. 0.148 x Sh = 0.28 grains
per developed um along z for a track perpendicular to the sheets, so the
number of grains seen per slice scales with the physical step. Persistent
dark dots (perpendicular beam tracks, 150 per view) are found in the
background-subtracted stack and peaks are counted along their depth
profile. The same procedure is applied to 4 real views and to synthetic
views rendered with image model v11 at dz 3.0 and 1.5 um (3 each):

| peak prominence (gray) | real | synthetic dz 3.0 | synthetic dz 1.5 |
|---|---|---|---|
| 4 | 0.063 | 0.089 | 0.074 |
| 6 | 0.032 | 0.064 | 0.033 |
| 9 | 0.010 | 0.041 | 0.008 |
| 12 | 0.003 | 0.029 | 0.002 |

Peaks per slice (mean over all dots); the four real views agree within
0.003. For every threshold the real views match the 1.5 um rendering and
not the 3.0 um one. The detection is non-linear: the ratio of the two
synthetic sets is 1.2 / 1.9 / 5.1 / 14.5 for thresholds 4 / 6 / 9 / 12.

Checks (threshold 6 / 9): taking every second slice of the real stacks,
i.e. doubling their step, gives 0.067 / 0.041, as the 3.0 um rendering
(0.064 / 0.041) and the 1.5 um rendering thinned the same way (0.069 /
0.047): the real stacks respond like the 1.5 um model at their native
step and like the 3.0 um model at twice the step. With 50 or 300 dots
per view instead of 150 the rates stay 0.030-0.033 (real), 0.032-0.034
(1.5 um) and 0.062-0.066 (3.0 um). With the image model v11 fixed and under the assumptions of the
MIP grain density (14.8 / 100 um, self-measured) and Sh = 1.9, the data
favour a physical step of about 1.5 um per slice in the developed
emulsion. The grain density and the shrinkage enter as a product; if it
were off by a factor of two the conclusion would change, so an
independent value (nominal emulsion sensitivity or a measurement) is
needed to settle the step.

## Attempt at an in-plane MIP grain density (`./run.py mipdensity`)

Sh = 1.9 is recorded by the experiment in the image JSON of the special
events. For the MIP grain density, peaks were counted along long straight
components of single slices (flat tracks in focus, exact pixel scale,
independent of the step). The estimator does not work at this
resolution: doubling the MIP grain density of the synthetic images
leaves the counted density unchanged (median 0.84 vs 0.78 per um, the
counts saturate on the fine texture along the lines), so the real /
synthetic ratio (1.09) carries no information. The MIP grain density, and
with it the final decision on the step, still needs an independent value
(nominal sensitivity of the emulsion batch or a dedicated measurement,
e.g. grain counting on flat MIP tracks at higher magnification).

## Emulsion layer extent in deep stacks (`./run.py layerdepth`)

A geometric measurement that needs neither the grain density nor the
image model: the 200-slice stacks of the special events contain whole
emulsion layers between the base and the surface. The fraction of dark
background-subtracted pixels per slice drops to zero outside the
emulsion; the half-height crossings bound each layer
(`output/layerdepth.txt`):

| stack | layer (slices) | stage (um) | physical step (um) | base gap (slices) |
|---|---|---|---|---|
| D005 | 54.2 | 163 | 4.37 | 7.0 |
| D013 | 62.3 | 187 | 3.80 | 6.9 |
| T004 | 59.0 | 177 | 4.01 | 6.8 |
| T011_200 upper | 64.8 | 195 | 3.65 | 7.0 |
| T011_200 lower | 59.5 | 178 | 3.98 | |

With 450 um per side before development and Sh = 1.9 (237 um developed)
the physical step is 4.0 +- 0.2 um per slice for the 3 um stage step
(developed / stage 1.2-1.5, as the focal shift of a dry objective in a
medium of index ~1.5); the 40 um base takes ~7 slices (21 um stage).
The fullscan stacks of ~58 slices thus cover about one whole layer, as
their profiles show (surface at slice ~5, emulsion to the end). A 1.5 um
step would need ~160 slices per layer and is excluded; the beam-grain
result above therefore means that the visible grain density along the
perpendicular beam tracks is lower than the assumed 14.8 / 100 um, not a
smaller step. Remaining assumption: nominal thickness 450 um and Sh 1.9
(a layer 20 % thinner would give 3.2 um).

Image model v15 (`output/image_model_v15_dz4.json`) is calibrated at
dz 4.0, starting from v14 with the depth scales rescaled by 4 / 1.5
(defocus slope, beam slopes, fog and track densities per volume divided,
mottle depth lengths and the flat-track dip multiplied): start loss 5.6
against 277 for v14 used unchanged at dz 4.0.
