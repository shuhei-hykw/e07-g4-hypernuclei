# Darkness at the true vertex: embedded v25 vs real specials (0.507 um/px)

Local darkness (B' darkness, module/rankstudy.py) at the true vertex
and its percentile among the B' points of the same view
(`data/union_audit/vtxdark.py`). Embedded: first 40 events of
v25_embed_nagara_dz4 (capture vertex, exact). Specials: the 9 clicked
vertices of results/specials_0507 (click good to 15-30 um, i.e. 30-60
px, which biases them low; another scanner, resampled).

| | n | darkness median [25 %, 75 %] | percentile median [25 %, 75 %] |
|---|---|---|---|
| embedded | 40 | 18.6 [14.4, 24.9] | 98.8 % [93.7, 99.8] |
| specials | 9 | 13.2 [7.6, 15.2] | 94.6 % [52.8, 97.1] |

Per special: D005 98.0 % (15.2), D013 40.2 % (2.8), IBUKI 46.5 % (7.6),
IRRAWADY 94.6 % (13.2), KISO 82.8 % (11.1), MINO 95.2 % (18.7), NAGARA
52.8 % (4.3), T004 99.9 % (27.0), T011 97.1 % (14.0).

Consistent with rendered vertices ~1.4 times darker than real ones, but
not conclusive: the click offset alone can explain low values (D013,
IBUKI, NAGARA).
