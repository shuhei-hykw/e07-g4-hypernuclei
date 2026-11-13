# Leak check of the ranking features on v25_embed_nagara_dz4

Embedded vs real candidates of the same views, not within 6 um of a truth vertex; near = 6-50 um from the capture point, far = beyond 200 um (both images identical there, expected 0.5). AUC out of fold (folds grouped by view); per feature the AUC of the feature alone (folded to >= 0.5).

- near: 8446 embedded / 5735 real candidates, adversary AUC 0.630
- far: 200000 embedded / 200000 real candidates, adversary AUC 0.516

| source | feature | near AUC | far AUC |
|---|---|---|---|
| A | mean_intens | 0.710 | 0.505 |
| A | n_tracks_max | 0.640 | 0.500 |
| Bprime | n_bp_20px | 0.597 | 0.517 |
| Bprime | darkness | 0.580 | 0.500 |
| A | n_slices | 0.571 | 0.502 |
| A | z_span_um | 0.564 | 0.500 |
| Bprime | no_other | 0.549 | 0.501 |
| A | angle_spread_best | 0.547 | 0.505 |
| A | darkness | 0.540 | 0.503 |
| Bprime | dist_other_um | 0.528 | 0.500 |
| A | angle_spread_max | 0.515 | 0.503 |
| A | dist_other_um | 0.514 | 0.503 |
| A | no_other | 0.508 | 0.504 |
| Bprime | mult | 0.504 | 0.504 |
