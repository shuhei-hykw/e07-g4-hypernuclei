# Signal vs real background lines in single_eval_v20_dz4

e07-fullscan d080818; medians per line (darkness in gray levels, beading = spread / mean along the line, gap = fraction below 5, width = dark px within 4 px per px of line).

- real background lines: n=14258 dark 22.8 beading 0.453 gap 0.0326 width 5.33
- signal lines (all): n=335 dark 25.6 beading 0.478 gap 0.0281 width 5.83
- signal dE/dx 0-0.0015 MeV/um: n=12 dark 13.7 beading 0.671 gap 0.229 width 3.83
- signal dE/dx 0.0015-0.005 MeV/um: n=29 dark 24.2 beading 0.473 gap 0.0435 width 5.75
- signal dE/dx 0.005-0.02 MeV/um: n=55 dark 24.8 beading 0.502 gap 0.0374 width 5.78
- signal dE/dx 0.02-1 MeV/um: n=237 dark 26.2 beading 0.466 gap 0.0241 width 5.9
- signal flat: n=101 dark 21.5 beading 0.473 gap 0.039 width 5.06
- signal steep: n=152 dark 32.3 beading 0.488 gap 0.0247 width 6.6

At equal apparent darkness (real lines have no truth):

| darkness | real: n, beading, gap, width | signal: n, beading, gap, width |
|---|---|---|
| 5-10 | 1312, 0.687, 0.403, 1.55 | 15, 0.724, 0.364, 2.14 |
| 10-15 | 557, 0.678, 0.189, 3.44 | 19, 0.594, 0.149, 3.98 |
| 15-20 | 1918, 0.505, 0.0682, 4.78 | 32, 0.482, 0.0576, 5.18 |
| 20-30 | 8606, 0.416, 0.0179, 5.73 | 167, 0.461, 0.0243, 5.63 |
| 30-50 | 1079, 0.364, 0, 6.31 | 92, 0.479, 0.016, 6.96 |
| 50-255 | 1, 0.512, 0.0656, 5.59 | 8, 0.537, 0, 9.99 |
