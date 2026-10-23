# Realism of image model v4 on held-out real views

20 random tohoku-v1 views not used in the calibration (and not the CNN-labelled ones) against 5 full-size v4 synthetic views (central 1024 px, slices at offsets -2..2, 5, 13). Mean +- std over views; z = (synthetic mean - real mean) / real std.

| statistic | real (20 views) | synthetic v4 (5) | z |
|---|---|---|---|
| median | 202.8 +- 3.4 | 189.2 +- 0.98 | -4.0 |
| std | 26.41 +- 8.1 | 26.08 +- 0.42 | -0.0 |
| p1 | 133.2 +- 30 | 128 +- 0.63 | -0.2 |
| p5 | 154.1 +- 29 | 146.2 +- 0.4 | -0.3 |
| frac_bright | 0.02524 +- 0.0054 | 0.01067 +- 0.002 | -2.7 |
| fg_fraction | 0.3351 +- 0.015 | 0.3783 +- 0.0054 | +2.8 |
| n_cc_per_mpx | 1.394e+04 +- 3.6e+02 | 1.156e+04 +- 72 | -6.6 |
| cc_area_p90 | 26.44 +- 3.6 | 27.6 +- 2 | +0.3 |
| n_elong_per_mpx | 18.41 +- 4.9 | 18.88 +- 5.5 | +0.1 |
| dark_blobs_per_mpx | 3257 +- 7.5e+02 | 3680 +- 1.8e+02 | +0.6 |
| dark_contrast_median | 0.1933 +- 0.037 | 0.2233 +- 0.0082 | +0.8 |
| d_skew | 0.9715 +- 0.058 | 0.9228 +- 0.035 | -0.8 |
| d_kurt | 5.299 +- 0.45 | 5.503 +- 0.098 | +0.5 |
| d_iqr | 14.6 +- 0.52 | 15.69 +- 0.11 | +2.1 |
| pow_low | 6.789 +- 0.15 | 6.97 +- 0.017 | +1.2 |
| pow_mid | 5.414 +- 0.035 | 5.595 +- 0.011 | +5.2 |
| pow_high | 3.876 +- 0.027 | 4.049 +- 0.018 | +6.5 |
| mp_detail_std | 11.93 +- 0.32 | 15.27 +- 0.24 | +10.5 |
| mp_dark_frac | 0.06338 +- 0.0047 | 0.0962 +- 0.0032 | +7.0 |
| corr1 | 0.9241 +- 0.0085 | 0.9233 +- 0.0035 | -0.1 |
| corr2 | 0.8332 +- 0.0057 | 0.8101 +- 0.0059 | -4.0 |
| corr5 | 0.5614 +- 0.017 | 0.5152 +- 0.0099 | -2.8 |
| corr13 | 0.1865 +- 0.024 | 0.3299 +- 0.048 | +5.9 |
