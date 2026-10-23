# Vertex-finding efficiency on simulated events

Depth convention: best; one-to-one greedy matching in the developed emulsion; chance level (ctrl) from 10 random shifts of the truths (module/efficiency.py).

| dataset | method | cuts | vertex | n | <6 um (ctrl) | <12 um (ctrl) | <30 um (ctrl) | dup <30 um | cand/event |
|---|---|---|---|---|---|---|---|---|---|
| v3_embed_nagara_dz1.5 | methodA | production | capture | 100 | 7% (2%) | 28% (10%) | 82% (68%) | 403 | 767 |
| v3_embed_nagara_dz1.5 | methodA | production | LL decay | 95 | 7% (2%) | 27% (9%) | 87% (64%) | 403 | 767 |
| v3_embed_nagara_dz1.5 | methodA | production | LL event | 95 | 0% (0%) | 9% (1%) | 74% (54%) | 403 | 767 |
| v3_embed_nagara_dz1.5 | methodA | relaxed | capture | 100 | 12% (2%) | 33% (11%) | 91% (73%) | 433 | 897 |
| v3_embed_nagara_dz1.5 | methodA | relaxed | LL decay | 95 | 7% (1%) | 25% (9%) | 86% (73%) | 433 | 897 |
| v3_embed_nagara_dz1.5 | methodA | relaxed | LL event | 95 | 0% (0%) | 7% (2%) | 80% (63%) | 433 | 897 |
| v3_embed_nagara_dz3.0 | methodA | production | capture | 100 | 12% (1%) | 24% (6%) | 65% (44%) | 191 | 768 |
| v3_embed_nagara_dz3.0 | methodA | production | LL decay | 96 | 3% (1%) | 14% (5%) | 61% (41%) | 191 | 768 |
| v3_embed_nagara_dz3.0 | methodA | production | LL event | 96 | 0% (0%) | 1% (0%) | 52% (31%) | 191 | 768 |
| v3_embed_nagara_dz3.0 | methodA | relaxed | capture | 100 | 6% (1%) | 21% (6%) | 63% (50%) | 221 | 899 |
| v3_embed_nagara_dz3.0 | methodA | relaxed | LL decay | 96 | 6% (1%) | 15% (7%) | 59% (49%) | 221 | 899 |
| v3_embed_nagara_dz3.0 | methodA | relaxed | LL event | 96 | 0% (0%) | 4% (1%) | 52% (37%) | 221 | 899 |
