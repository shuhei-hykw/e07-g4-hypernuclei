# Steep-line checks (module.beamlines)

Track vs persistent-texture checks on four raw real views (V00000885,
V00000260, V00001399, V00000560; finder of module.beamlines):

| check | result |
|---|---|
| unique slices | 58 / 58 in every view |
| global stack shift (phase correlation, 5 slices apart) | 0, 0 px |
| lines per crop | 2335-2563 |
| lines with the slice order permuted within the view | 100-254 |
| lines with selection slices from different views | 0-3 |
| line positions coinciding across views (within 2 px) | 62, 72, 70 (chance 68, 72, 73) |
| darkness outside the emulsion (line / shifted) | 8.1 / 5.6 (inside 46 / 11) |
| depth profile along the line: CV, lag-1 autocorrelation | 0.25, 0.84 |

The lines need the depth order, do not repeat at fixed pixels across
views and fade outside the emulsion. This argues against a fixed sensor
pattern and simple stack artefacts; the lines are consistent with steep
tracks through the layer such as the nearly perpendicular beam (beam
purity and finder recall unknown).

Perpendicular minimum-ionizing tracks rendered with signal model v20
and embedded into real views (`beam_eval_v20_dz4`, pi- 1.8 GeV/c, tilt
0-3 deg, `output/beamlines_beam_v20.md`), measurement slices only:

| sample | n | darkness 16 / 50 / 84 % | contrast | HWHM (px) | lag-1 autocorr |
|---|---|---|---|---|---|
| real lines, same views | 41389 | 36 / 45.5 / 56.2 | 35.8 | 2 | 0.69 |
| v20 beam, found by the finder | 109 | 33.4 / 41.9 / 51 | 41.3 | 3 | 0.67 |
| v20 beam, all at the truth line | 147 | 29.8 / 38.9 / 49.3 | 36.8 | 3 | 0.66 |

The finder recovers 109 distinct tracks of the 147 rendered inside the
crop. Among finder-selected tracks the v20 minimum-ionizing track is
about 8 % fainter than the real lines (15 % for all tracks at the truth
line) and 1 px wider (HWHM), with a similar depth autocorrelation. This
suggests a broadly similar appearance within this selection, not
agreement with an identified real beam sample.
