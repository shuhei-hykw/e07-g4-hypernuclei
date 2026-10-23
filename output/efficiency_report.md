# Vertex-finding efficiency on simulated events

Depth convention: edge; one-to-one greedy matching in the developed emulsion; chance level (ctrl) from 10 random shifts of the truths (module/efficiency.py).

| dataset | method | cuts | vertex | n | <6 um (ctrl) | <12 um (ctrl) | <30 um (ctrl) | dup <30 um | cand/event |
|---|---|---|---|---|---|---|---|---|---|
| v3_embed_nagara_dz3.0 | methodA | production | capture | 100 | 2% (2%) | 18% (10%) | 69% (58%) | 181 | 731 |
| v3_embed_nagara_dz3.0 | methodA | production | LL decay | 96 | 0% (1%) | 17% (7%) | 74% (57%) | 181 | 731 |
| v3_embed_nagara_dz3.0 | methodA | production | LL event | 96 | 0% (0%) | 6% (1%) | 59% (47%) | 181 | 731 |
| v3_embed_nagara_dz3.0 | methodA | relaxed | capture | 100 | 3% (2%) | 20% (12%) | 64% (60%) | 211 | 854 |
| v3_embed_nagara_dz3.0 | methodA | relaxed | LL decay | 96 | 2% (2%) | 17% (9%) | 67% (59%) | 211 | 854 |
| v3_embed_nagara_dz3.0 | methodA | relaxed | LL event | 96 | 0% (0%) | 6% (3%) | 55% (50%) | 211 | 854 |
