### Ranking of B' candidates (`./run.py rankstudy`, same 20 events)

Branch points scored by multiplicity, by the number of branch points
within 5 um (density), by the summed multiplicity within 5 um, or by the
image darkness around the point (background-subtracted minimum
projection over +-3 slices, mean within 6 px); the top K per view are
kept. Efficiency (chance) at <6 / <12 um, one-to-one.

| score | top K | capture | LL decay |
|---|---|---|---|
| multiplicity | 770 | 5% (0%) / 15% (4%) | 0% (1%) / 0% (6%) |
| density | 770 | 0% (2%) / 10% (4%) | 0% (1%) / 0% (3%) |
| mult_density | 770 | 5% (0%) / 5% (3%) | 0% (1%) / 0% (3%) |
| darkness | 100 | 10% (0%) / 20% (0%) | 0% (0%) / 11% (0%) |
| darkness | 300 | 10% (0%) / 40% (2%) | 0% (1%) / 16% (3%) |
| darkness | 770 | 20% (0%) / 55% (4%) | 11% (1%) / 37% (6%) |
| darkness | 3000 | 30% (3%) / 75% (15%) | 21% (2%) / 79% (14%) |

At the candidate rate of method A (~770 per view) the darkness ranking
keeps capture vertices with 20 % / 55 % efficiency at a chance level of
0 % / 4 %, against 15 % (2 %) / 15 % (8 %) for method A with z_best and
20 % (10 %) / 45 % (36 %) with the inclusive convention on the same
events. Positions and multiplicity alone do not rank. Caveat: the
signal darkness comes from the image model (checked against the real
Ibuki heavy tracks, section on calibration in README); a darker
simulated signal than real would favour this ranking.

Same study on the dz 1.5 um dataset (v3_embed_nagara_dz1.5, the same
20 events), next to dz 3.0, at 770 candidates per view:

| score | dz 1.5 capture | dz 1.5 LL decay | dz 3.0 capture | dz 3.0 LL decay |
|---|---|---|---|---|
| multiplicity | 0% (1%) / 0% (7%) | 0% (1%) / 11% (6%) | 5% (0%) / 15% (4%) | 0% (1%) / 0% (6%) |
| density | 10% (2%) / 35% (8%) | 5% (2%) / 16% (7%) | 0% (2%) / 10% (4%) | 0% (1%) / 0% (3%) |
| darkness | 50% (3%) / 70% (7%) | 21% (1%) / 47% (6%) | 20% (0%) / 55% (4%) | 11% (1%) / 37% (6%) |
| combined | 20% (2%) / 35% (8%) | 11% (1%) / 32% (5%) | 5% (1%) / 45% (6%) | 5% (0%) / 16% (5%) |

The multiplicity does not rank the simulated vertices at either step,
while it does on the real special events: the difference lies in the
simulated signal (its branch-point topology), not in the focal step.
The darkness ranking works at both steps and better at 1.5 um.
