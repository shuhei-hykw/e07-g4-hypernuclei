# Blind crop labelling, round 2 (protocol, fixed before any crop exists)

Goal: precision of the review rankers on real plate sites, where the
simulation gives recall only. Compared: darkness, logistic
(output/models/review_logistic_v25.json) and logistic without
mean_intens (review_logistic_v25_no_mean_intens.json). Labels are per
site, so a later ranker (v26 refit) can be rescored on the same sites.

## Sites

- Pool: the production union (all method A + B' top 770 by darkness)
  of every full-scan view of the plate, scored by each ranker; one
  site = one candidate, sites closer than 25 um merged (keep the best).
- Excluded views: output/labels/excluded_background_views.txt (every
  view used as an embedding background in any set, so neither a
  training nor a sealed evalV view, plus the four CNN-labelled views).
- Global rank per ranker over the plate; bins 0-100, 100-1000,
  1000-5000 (5000 = review budget target).

## Allocation (300 crops)

| group | definition (top 5000 = global rank < 5000) | n |
|---|---|---|
| L-only | logistic top 5000, not darkness top 5000 | 60 |
| D-only | darkness top 5000, not logistic top 5000 | 60 |
| N-only | no-mean_intens top 5000, not logistic top 5000 | 35 |
| M-only | logistic top 5000, not no-mean_intens top 5000 | 35 |
| both | logistic and darkness top 5000 | 40 |
| random | uniform over the pool (precision floor) | 40 |
| repeat | 30 of the above, shown again | 30 |

Each discordant group is split evenly over the three rank bins of the
ranker that selected it. Sites may fall in several groups; draw the
groups in the order above without replacement.

## Crop and display

- 50 um square (99 px at 0.507 um/px) centred on the site, full z
  stack, nothing else shown: no ranker, score, source or group.
- Order shuffled with a fixed seed; repeats at least 100 crops apart.
- Key (site id -> view, px, z, groups, ranks, scores) in a separate
  file the labeller does not open.

## Rubric (one grade per crop)

| grade | meaning |
|---|---|
| 0 | nothing: noise, dust, grain clusters, a single track |
| 1 | ambiguous: two tracks crossing or a kink, no clear common point |
| 2 | star: three or more prongs from one point, not Xi-like |
| 3 | Xi-like: a track that stops in a star (capture), or a star with a short heavy prong and a second vertex nearby |

Also per crop: click the vertex (x, y) and the slice of best focus
when grade >= 2; free comment optional. "good" = grade >= 2; Xi-like
= grade 3.

## Analysis (fixed now)

- Precision per group with Wilson 95 % intervals, also per rank bin.
- L-only vs D-only and N-only vs M-only: difference of good fractions
  (this is the McNemar test of the paired top-5000 sets restricted to
  their discordant sites), exact p-value; the same for Xi-like.
- Repeat agreement (fraction and Cohen kappa) as the label noise.
- Decision: logistic replaces darkness if L-only good > D-only good
  with p < 0.05; mean_intens is kept only if M-only good > N-only good
  with p < 0.05, otherwise the no-mean_intens model is the default.

## Implementation (agreed with e07-fullscan before any crop)

- Scores: candscore apply on results/fullscan_union/features.parquet;
  darkness from features.parquet (union.parquet has NaN for A rows).
- Sites per ranker: greedy merge in plate xy within 25 um (z ignored,
  the crop shows the full stack), best score kept; global rank = order
  of the merged sites; "in X's top 5000" = within 25 um of one of X's
  top 5000 sites. Groups skip sites within 25 um of a drawn site;
  random is uniform over pool candidates.
- Display view: the overlapping view (own or neighbour) placing the
  site farthest from the edge; the key keeps both views.
- Crop 99 x 99 px, every slice, 5x upscaled, grey stretch 1-99 % per
  site. Repeats: 30 random picks of the 270 under new ids.
- Context: a 150 um minimum projection with the 50 um box drawn. The
  grade refers to a vertex inside the box; the context only helps to
  read the topology (a track ending in the box, a nearby second vertex).

## Deviation at drawing (before any label)

Short bins pass their remainder to the later bins of the same group:
L-only 20/20/20, D-only 0/30/30 (every darkness top-100 site is in
the logistic top 5000), N-only 1/18/16, M-only 12/12/11; both 40,
random 40, repeats 30 (minimum gap 106). The group comparisons are
therefore also reported over bins 100-5000 only, where both groups
are populated. Three D-only B' sites lie 0-4 px from the edge of
every usable view (border pile-up); D-only precision is shown with
and without them. Key: e07-fullscan results/label_round2/key.csv,
code scripts/label_round2.py, seed 20261114.
