# Renderer revision: design note (not started)

Status 2026-11-02: deferred (user away until 11-09; agreed with Codex).
Embedded evaluations use real backgrounds and are not affected.

## What the beam-line study adds (output/beamlines_*.md)

- Real views hold ~2,370 finder lines per 1024 px crop (steep, 0.8-2.3
  deg from the normal, darkness 46). Whole-synthetic v19 views: no
  rendered beam track matched; v24: 581 lines per crop, 397 matched to
  beam tracks, at 3.9 deg.
- The v20 signal grain model gives perpendicular MIP tracks broadly
  similar to the real lines (8 % fainter, 1 px wider among
  finder-selected tracks).
- v24 compensated with beam density at 12x v19 and min_transmission at
  its bound; 9 / 23 held-out statistics stay beyond 3 spreads.

## Candidate targeted revision (instead of the 5-8 day split)

1. Render background beam tracks with the v20 signal grain model
   (the beam is minimum ionizing; fog and junk keep v19).
2. Beam axis and spread from the real line angle distribution, not free
   in the calibration.
3. Beam density as a forward-model parameter: choose it so that the
   finder line count and darkness on synthetic views match the real
   ones (never insert the real count directly; purity and recall are
   unknown).
4. Recalibrate the remaining parameters; adopt only if both the
   held-out statistics and the line observables improve (linecheck,
   beamlines), with the same frozen held-out views.

Cost estimate: 1-2 days of work plus calibration CPU (after the fullscan
B' run). Only worth it if whole-synthetic images are needed (e.g.
training).
