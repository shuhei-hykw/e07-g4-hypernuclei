# Vertex-finding efficiency on simulated events

Depth convention: inclusive; one-to-one greedy matching in the developed emulsion; chance level (ctrl) from 10 random shifts of the truths (module/efficiency.py).

| dataset | method | cuts | vertex | n | <6 um (ctrl) | <12 um (ctrl) | <30 um (ctrl) | dup <30 um | cand/event |
|---|---|---|---|---|---|---|---|---|---|
| v3_embed_nagara_dz1.5 | methodA | production | capture | 100 | 42% (9%) | 71% (35%) | 99% (90%) | 403 | 767 |
| v3_embed_nagara_dz1.5 | methodA | production | LL decay | 95 | 35% (10%) | 60% (29%) | 97% (90%) | 403 | 767 |
| v3_embed_nagara_dz1.5 | methodA | production | LL event | 95 | 9% (1%) | 41% (14%) | 97% (87%) | 403 | 767 |
| v3_embed_nagara_dz1.5 | methodA | relaxed | capture | 100 | 42% (12%) | 78% (44%) | 100% (96%) | 433 | 897 |
| v3_embed_nagara_dz1.5 | methodA | relaxed | LL decay | 95 | 39% (12%) | 60% (39%) | 99% (96%) | 433 | 897 |
| v3_embed_nagara_dz1.5 | methodA | relaxed | LL event | 95 | 9% (1%) | 44% (20%) | 99% (95%) | 433 | 897 |
| v3_embed_nagara_dz3.0 | methodA | production | capture | 100 | 28% (10%) | 52% (31%) | 89% (79%) | 191 | 768 |
| v3_embed_nagara_dz3.0 | methodA | production | LL decay | 96 | 14% (8%) | 43% (26%) | 91% (79%) | 191 | 768 |
| v3_embed_nagara_dz3.0 | methodA | production | LL event | 96 | 2% (1%) | 19% (12%) | 84% (73%) | 191 | 768 |
| v3_embed_nagara_dz3.0 | methodA | relaxed | capture | 100 | 29% (14%) | 69% (42%) | 96% (87%) | 221 | 899 |
| v3_embed_nagara_dz3.0 | methodA | relaxed | LL decay | 96 | 23% (11%) | 54% (34%) | 98% (87%) | 221 | 899 |
| v3_embed_nagara_dz3.0 | methodA | relaxed | LL event | 96 | 5% (2%) | 38% (18%) | 95% (84%) | 221 | 899 |
