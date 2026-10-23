# Vertex-finding efficiency on simulated events

Edge convention, one-to-one greedy matching in the developed emulsion; chance level (ctrl) from 10 random shifts of the truths (module/efficiency.py).

| dataset | method | cuts | vertex | n | <6 um (ctrl) | <12 um (ctrl) | <30 um (ctrl) | dup <30 um | cand/event |
|---|---|---|---|---|---|---|---|---|---|
| v3_embed_nagara5_dz1.5 | methodA | production | capture | 5 | 0% (6%) | 40% (22%) | 100% (90%) | 19 | 771 |
| v3_embed_nagara5_dz1.5 | methodA | production | LL decay | 4 | 25% (0%) | 25% (18%) | 100% (85%) | 19 | 771 |
| v3_embed_nagara5_dz1.5 | methodA | production | LL event | 4 | 0% (0%) | 25% (5%) | 100% (82%) | 19 | 771 |
| v3_embed_nagara5_dz1.5 | methodA | relaxed | capture | 5 | 0% (4%) | 60% (22%) | 100% (98%) | 22 | 906 |
| v3_embed_nagara5_dz1.5 | methodA | relaxed | LL decay | 4 | 0% (2%) | 25% (20%) | 100% (95%) | 22 | 906 |
| v3_embed_nagara5_dz1.5 | methodA | relaxed | LL event | 4 | 0% (0%) | 25% (2%) | 100% (95%) | 22 | 906 |
| v3_embed_nagara5_dz3.0 | methodA | production | capture | 5 | 0% (4%) | 0% (14%) | 80% (54%) | 10 | 777 |
| v3_embed_nagara5_dz3.0 | methodA | production | LL decay | 4 | 0% (0%) | 0% (5%) | 100% (52%) | 10 | 777 |
| v3_embed_nagara5_dz3.0 | methodA | production | LL event | 4 | 0% (0%) | 0% (0%) | 75% (40%) | 10 | 777 |
| v3_embed_nagara5_dz3.0 | methodA | relaxed | capture | 5 | 20% (0%) | 20% (10%) | 60% (48%) | 5 | 907 |
| v3_embed_nagara5_dz3.0 | methodA | relaxed | LL decay | 4 | 0% (2%) | 0% (8%) | 25% (32%) | 5 | 907 |
| v3_embed_nagara5_dz3.0 | methodA | relaxed | LL event | 4 | 0% (0%) | 0% (2%) | 0% (25%) | 5 | 907 |
| v3_synth_nagara5_dz1.5 | methodA | production | capture | 5 | 40% (0%) | 40% (0%) | 40% (4%) | 0 | 41 |
| v3_synth_nagara5_dz1.5 | methodA | production | LL decay | 4 | 0% (0%) | 0% (0%) | 0% (5%) | 0 | 41 |
| v3_synth_nagara5_dz1.5 | methodA | production | LL event | 4 | 0% (0%) | 0% (0%) | 0% (0%) | 0 | 41 |
| v3_synth_nagara5_dz1.5 | methodA | relaxed | capture | 5 | 60% (0%) | 80% (0%) | 80% (10%) | 3 | 99 |
| v3_synth_nagara5_dz1.5 | methodA | relaxed | LL decay | 4 | 0% (2%) | 25% (2%) | 50% (15%) | 3 | 99 |
| v3_synth_nagara5_dz1.5 | methodA | relaxed | LL event | 4 | 0% (0%) | 25% (0%) | 50% (5%) | 3 | 99 |
| v3_synth_nagara5_dz3.0 | methodA | production | capture | 5 | 0% (0%) | 20% (2%) | 20% (6%) | 0 | 42 |
| v3_synth_nagara5_dz3.0 | methodA | production | LL decay | 4 | 0% (0%) | 25% (0%) | 50% (10%) | 0 | 42 |
| v3_synth_nagara5_dz3.0 | methodA | production | LL event | 4 | 0% (0%) | 0% (0%) | 25% (5%) | 0 | 42 |
| v3_synth_nagara5_dz3.0 | methodA | relaxed | capture | 5 | 20% (0%) | 40% (6%) | 80% (12%) | 0 | 101 |
| v3_synth_nagara5_dz3.0 | methodA | relaxed | LL decay | 4 | 0% (0%) | 25% (2%) | 75% (12%) | 0 | 101 |
| v3_synth_nagara5_dz3.0 | methodA | relaxed | LL event | 4 | 0% (0%) | 0% (0%) | 75% (8%) | 0 | 101 |
