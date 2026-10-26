# Signal vs real background lines in single_eval_v20_dz4

e07-fullscan 869804c; medians per line (darkness in gray levels, beading = spread / mean along the line, gap = fraction below 5, width = dark px within 4 px per px of line).

- real background lines: n=56716 dark 22.9 beading 0.449 gap 0.031 width 5.34
- signal lines (all): n=1289 dark 24.9 beading 0.486 gap 0.031 width 5.78
- signal dE/dx 0-0.0015 MeV/um: n=36 dark 12.7 beading 0.679 gap 0.242 width 3.93
- signal dE/dx 0.0015-0.005 MeV/um: n=95 dark 21.7 beading 0.556 gap 0.0577 width 5.27
- signal dE/dx 0.005-0.02 MeV/um: n=273 dark 24.2 beading 0.504 gap 0.0426 width 5.8
- signal dE/dx 0.02-1 MeV/um: n=879 dark 25.7 beading 0.476 gap 0.026 width 5.87
- signal flat: n=460 dark 21.9 beading 0.476 gap 0.0321 width 5.29
- signal steep: n=564 dark 31 beading 0.516 gap 0.0349 width 6.53

At equal apparent darkness (real lines have no truth):

| darkness | real: n, beading, gap, width | signal: n, beading, gap, width | signal flat (dip < 15 deg) |
|---|---|---|---|
| 5-10 | 5167, 0.7, 0.404, 1.6 | 52, 0.714, 0.368, 2.21 | 28, 0.649, 0.319, 2.26 |
| 10-15 | 2325, 0.726, 0.203, 3.42 | 68, 0.629, 0.154, 3.78 | 42, 0.574, 0.141, 3.86 |
| 15-20 | 7401, 0.503, 0.0671, 4.75 | 133, 0.523, 0.0633, 4.98 | 79, 0.493, 0.0563, 5.25 |
| 20-30 | 34141, 0.414, 0.0169, 5.75 | 674, 0.46, 0.0241, 5.68 | 295, 0.449, 0.0199, 5.58 |
| 30-50 | 4557, 0.365, 0, 6.37 | 323, 0.483, 0.0145, 7.01 | 11, 0.505, 0.0294, 6.37 |
| 50-255 | 14, 0.246, 0, 6.49 | 27, 0.554, 0.0286, 9.02 | 0, -, -, - |
