# Chance-subtracted ranking on v25_embed_nagara_dz4

100 events, learned features without mean_intens. Capture found among the top K of the view: embedded (real views alone), G = embedded - real with 95 % bootstrap interval; plate total = K x 2025 views (before cross-view dedup).

## capture within 6 um

| ranker | K=1 (2025) | K=2 (4050) | K=4 (8100) | K=8 (16200) | K=16 (32400) |
|---|---|---|---|---|---|
| darkness | 0% (0%) G +0% [+0%, +0%] | 1% (0%) G +1% [+0%, +3%] | 2% (0%) G +2% [+0%, +5%] | 3% (0%) G +3% [+0%, +7%] | 8% (0%) G +8% [+3%, +14%] |
| logistic | 3% (0%) G +3% [+0%, +7%] | 3% (0%) G +3% [+0%, +7%] | 4% (0%) G +4% [+1%, +8%] | 5% (0%) G +5% [+1%, +10%] | 9% (0%) G +9% [+4%, +15%] |
| logistic - darkness | +3% [+0%, +7%] | +2% [+0%, +5%] | +2% [+0%, +5%] | +2% [-2%, +6%] | +1% [-5%, +7%] |

## capture within 12 um

| ranker | K=1 (2025) | K=2 (4050) | K=4 (8100) | K=8 (16200) | K=16 (32400) |
|---|---|---|---|---|---|
| darkness | 0% (0%) G +0% [+0%, +0%] | 2% (0%) G +2% [+0%, +5%] | 6% (0%) G +6% [+2%, +11%] | 8% (0%) G +8% [+3%, +14%] | 19% (0%) G +19% [+12%, +27%] |
| logistic | 3% (0%) G +3% [+0%, +7%] | 3% (0%) G +3% [+0%, +7%] | 5% (0%) G +5% [+1%, +9%] | 9% (0%) G +9% [+4%, +15%] | 15% (0%) G +15% [+8%, +23%] |
| logistic - darkness | +3% [+0%, +7%] | +1% [+0%, +3%] | -1% [-5%, +3%] | +1% [-5%, +7%] | -4% [-12%, +4%] |
