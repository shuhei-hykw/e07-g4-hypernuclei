# Signal vs real background lines in single_eval_dz4

e07-fullscan 2485ade; medians per line (darkness in gray levels, beading = spread / mean along the line, gap = fraction below 5, width = dark px within 4 px per px of line).

- real background lines: n=14327 dark 22.9 beading 0.452 gap 0.0323 width 5.33
- signal lines (all): n=251 dark 25.8 beading 0.651 gap 0.0976 width 5.37
- signal dE/dx 0-0.0015 MeV/um: n=3 dark 4.83 beading 1.45 gap 0.64 width 0.523
- signal dE/dx 0.0015-0.005 MeV/um: n=14 dark 24.6 beading 0.675 gap 0.107 width 5.53
- signal dE/dx 0.005-0.02 MeV/um: n=44 dark 24.3 beading 0.703 gap 0.126 width 5.19
- signal dE/dx 0.02-1 MeV/um: n=189 dark 26.4 beading 0.637 gap 0.0851 width 5.39
- signal flat: n=66 dark 19.3 beading 0.723 gap 0.147 width 4.35
- signal steep: n=129 dark 29.5 beading 0.621 gap 0.0619 width 5.83

At equal apparent darkness (real lines have no truth):

| darkness | real: n, beading, gap, width | signal: n, beading, gap, width |
|---|---|---|
| 5-10 | 1315, 0.685, 0.4, 1.55 | 7, 0.891, 0.38, 2.26 |
| 10-15 | 572, 0.689, 0.193, 3.46 | 16, 0.696, 0.204, 3.47 |
| 15-20 | 1916, 0.506, 0.069, 4.78 | 30, 0.754, 0.174, 4.03 |
| 20-30 | 8691, 0.417, 0.0179, 5.73 | 125, 0.631, 0.0864, 5.38 |
| 30-50 | 1084, 0.363, 0, 6.31 | 63, 0.623, 0.0462, 6.41 |
| 50-255 | 1, 0.512, 0.0656, 5.59 | 2, 0.527, 0.0481, 8.46 |
