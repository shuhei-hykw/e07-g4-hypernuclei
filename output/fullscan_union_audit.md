# Read-only audit of the full-scan union (2026-11-02)

Input: `~/work/e07-fullscan/results/fullscan_union/union.parquet`
(2,025 views; A 1,554,446 rows, B' 1,559,250 rows). The four
CNN-labelled views are left out of every statistic. Scripts:
`data/union_audit/` (local).

## 1. View overlap and the stage coordinates

- Neighbouring views share content at a fixed pixel step, measured from
  coincident points of the union and confirmed by phase correlation of
  the middle slices (after masking the zero-shift static pattern):
  VX+1: (+987.5, -2.6) px, VY+1: (-1.1, -987.4) px, 16-84 % of 1,972
  pairs within +-1 px across the plate.
- The stage fields step by 0.5 per view (`view_x_mm`, `view_y_mm`,
  JSON `x`, `y`). At 0.289 um/px (simulation, specials, fullscan
  config) 987.5 px is 285 um, so these fields are not millimetres at
  that scale (factor 1.75; or the pixel would be 0.51 um). Global
  coordinates should use the measured pixel step, not `view_x_mm` x
  1000.
- Overlap: 1,060 px (52 %) per axis; every location of the interior is
  covered by about 4 views.

## 2. Reproducibility in the overlaps (cross-view duplicates)

Points shifted into the neighbour's field (60 px margin), fraction with
a point of the same source within r, against a control shift of
(150, 150) px:

| source | < 2 px | < 5 px | < 10 px | < 20 px | < 50 px |
|---|---|---|---|---|---|
| A | 0.7 % (0.2) | 4.2 % (1.4) | 14.7 % (5.4) | 41.9 % (21.1) | 91.2 % (80.6) |
| B' | 7.0 % (0.2) | 26.1 % (0.9) | 43.3 % (3.0) | 51.3 % (10.6) | 68.9 % (44.8) |

- Method A candidates are weakly reproducible between two images of
  the same place (14.7 % within 10 px against 5.4 % by chance); B'
  reproduces better (43 % against 3 %), limited by the per-view top 770.
- Exact duplicates (3 px) with an edge neighbour: A 5.0 %, B' 27.1 %.
- Consequence: the full scan looks at each place ~4 times with largely
  different candidates, so plate-level efficiency (found in any view)
  can exceed the single-view efficiencies of the simulation, and so can
  the chance level; a plate-level comparison needs a matched global
  candidate density and a cross-view merge with the measured step.

## 3. Image-border pile-up of B'

- 12.3 % of the B' rows (94.5 per view) lie within 10 px of the image
  border, 6.3x the mean density (10-200 px from the border: 0.46-0.75x).
  By side: left 6.4 %, top 3.8 %, right 1.6 %, bottom 0.5 %; dark_rank
  median 140-236 (uniform 385): they are ranked high.
- The asymmetry points to the image itself (dark border columns / rows)
  as well as the clipped darkness window. These points sit where the
  neighbouring view covers the same place with full context, so they
  take ~12 % of the B' budget without adding coverage. A border margin
  before the top-770 cut is a candidate change, to be tuned on a tune
  or dev set, not on evalU.

## 4. Darkness across views

- 770th darkness per view: 1 / 16 / 50 / 84 / 99 % = 15.4 / 15.8 / 16.2
  / 16.8 / 19.3; plate-edge views (VX or VY = 0 or 44) 16.4 against 16.2
  inside.
- 51 views beyond 4 robust spreads, all darker (threshold 18-19),
  clustered in a few plate regions (e.g. VX 43-44 / VY 41-42, VX 39-40 /
  VY 15-18, VX 20-22 / VY 18-19): local plate features, not errors.
- dark_rank is 0..769 with 770 unique values in every view.

## 5. Open question: the full-scan pixel scale (for the user)

The full-scan views (station NLAB-PC13) carry no pixel scale (identity
affine); the specials (NLAB-PC06, a different station) carry 0.28889
um/px. With the measured view step of 987.5 px, either the full-scan
pixel is 0.289 um and the stage x / y fields (0.5 per view) are not
millimetres, or the stage moved 0.5 mm and the pixel is ~0.506 um. The
z field steps by 0.003 per slice, natural as mm. Image evidence,
full scan vs specials (same 52-slice window, `data/union_audit/`):

| ruler | full scan | specials | ratio | 0.289 / 0.506 predict |
|---|---|---|---|---|
| darkness autocorrelation HWHM (px) | 3.12-3.13 | 3.09-3.62 | ~0.9 | 1.0 / 0.57 |
| steep-line radial HWHM, sub-pixel (px) | 1.78 | 2.47-2.60 | 0.70 | 1.0 / 0.57 |
| steep lines per 1024 px crop | ~2,400 | 840-1,000 | 2.4-2.9 | 1 / 3.06 (if equal beam flux) |
| line drift (px / slice) | 0.35 | 0.48-0.51 | 0.7 | 1.0 / 0.57 (if equal angle) |

The rulers disagree (the optics and the local beam differ between the
stations and plates), so the scale is not settled from the images.
Numbers of this repository that assume 0.289 um/px for full-scan views:
the Geant4 -> pixel geometry of every embedded signal (track lengths
and vertex separations in px), the um tolerances of the efficiency
tables (6 / 12 um), the 3-D distance with z in um, the steep-line
density per mm2 (~27,000 at 0.289, ~8,800 at 0.506) and the in-plane
MIP grain density. The image model itself was calibrated on pixel
statistics. Relative results (union against A on the same images) keep
their sign, but absolute efficiencies and tolerances would need a
re-render at the correct scale.

## 6. Border margin on the dev set (candidate only)

`./run.py unionstudy ... --border-px 10 | 25` on v20_embed_nagara_dz4
(`output/unionstudy_v20_dz4_border10.txt`, `_border25.txt`; evalU
untouched): B' alone at 12 um 26 -> 28 % (capture), LL event edge 22 ->
26 %, LL decay edge 38 -> 40 %; union capture inclusive 12 um unchanged
at 71 %, G unchanged at +19 % [+12, +27]; no cell worse, 10 and 25 px
identical. A consistent but small gain (1-4 of 100 events, within the
statistical uncertainty); adopting it would need a fresh unused set.
