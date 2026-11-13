# Veto check on the special events at 0.507 um/px

Per arm: candidates A / B', A match, best darkness rank of a matching B' point ('-' = none), union found (A, or B' within the top K). Match: tolerance tol_xy_um of gt.json in XY and 30 um in stage z; K = 770 scaled by the image area to a 2048 px view. Control: union found at random positions 250-450 px away (450 in all, common to the arms).

| event | K | production | zfix | alen |
|---|---|---|---|---|
| D005 | 250 | 255/11394 A y B' 14 FOUND | 255/9178 A y B' 16 FOUND | 288/11394 A y B' 14 FOUND |
| D013 | 250 | 166/2680 A y B' 17 FOUND | 166/2218 A y B' 19 FOUND | 199/2680 A y B' 17 FOUND |
| IBUKI | 250 | 279/2587 A y B' 269 FOUND | 279/2247 A y B' 124 FOUND | 287/2587 A y B' 269 FOUND |
| IRRAWADY | 250 | 153/1906 A y B' 40 FOUND | 153/1681 A y B' 136 FOUND | 193/1906 A y B' 40 FOUND |
| KISO | 250 | 181/6571 A y B' 646 FOUND | 181/5450 A y B' 501 FOUND | 177/6571 A y B' 646 FOUND |
| MINO | 250 | 259/2497 A y B' 118 FOUND | 259/2144 A y B' 309 FOUND | 282/2497 A y B' 118 FOUND |
| NAGARA | 250 | 324/8223 A y B' 911 FOUND | 324/6815 A y B' 624 FOUND | 321/8223 A y B' 911 FOUND |
| T004 | 250 | 238/7608 A y B' 3 FOUND | 238/5997 A y B' 9 FOUND | 260/7608 A y B' 3 FOUND |
| T011 | 250 | 76/2856 A y B' 13 FOUND | 76/2378 A y B' 45 FOUND | 107/2856 A y B' 13 FOUND |

| arm | found | control | lost vs production |
|---|---|---|---|
| production | 9/9 | 78% | none |
| zfix | 9/9 | 80% | none |
| alen | 9/9 | 83% | none |
