# Vertex-finding efficiency on simulated events

Depth convention: edge; one-to-one greedy matching in the developed emulsion; chance level (ctrl) from 10 random shifts of the truths (module/efficiency.py).

| dataset | method | cuts | vertex | n | <6 um (ctrl) | <12 um (ctrl) | <30 um (ctrl) | dup <30 um | cand/event |
|---|---|---|---|---|---|---|---|---|---|
| v3_embed_nagara_dz1.5 | methodA | production | capture | 100 | 11% (3%) | 36% (18%) | 96% (88%) | 403 | 767 |
| v3_embed_nagara_dz1.5 | methodA | production | LL decay | 95 | 11% (3%) | 35% (15%) | 92% (86%) | 403 | 767 |
| v3_embed_nagara_dz1.5 | methodA | production | LL event | 95 | 2% (0%) | 16% (4%) | 89% (81%) | 403 | 767 |
| v3_embed_nagara_dz1.5 | methodA | relaxed | capture | 100 | 13% (4%) | 32% (21%) | 96% (92%) | 433 | 897 |
| v3_embed_nagara_dz1.5 | methodA | relaxed | LL decay | 95 | 12% (2%) | 24% (17%) | 89% (91%) | 433 | 897 |
| v3_embed_nagara_dz1.5 | methodA | relaxed | LL event | 95 | 1% (0%) | 13% (5%) | 87% (87%) | 433 | 897 |
| v3_embed_nagara_dz3.0 | methodA | production | capture | 100 | 3% (2%) | 20% (11%) | 73% (62%) | 191 | 768 |
| v3_embed_nagara_dz3.0 | methodA | production | LL decay | 96 | 0% (1%) | 18% (8%) | 77% (60%) | 191 | 768 |
| v3_embed_nagara_dz3.0 | methodA | production | LL event | 96 | 0% (0%) | 6% (1%) | 62% (50%) | 191 | 768 |
| v3_embed_nagara_dz3.0 | methodA | relaxed | capture | 100 | 4% (2%) | 21% (13%) | 67% (63%) | 221 | 899 |
| v3_embed_nagara_dz3.0 | methodA | relaxed | LL decay | 96 | 2% (2%) | 17% (9%) | 70% (61%) | 221 | 899 |
| v3_embed_nagara_dz3.0 | methodA | relaxed | LL event | 96 | 0% (0%) | 6% (3%) | 57% (52%) | 221 | 899 |
