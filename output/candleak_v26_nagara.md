# Leak check of the ranking features on v26_embed_nagara_dz4

Embedded vs real candidates of the same views, not within 6 um of a truth vertex; near = 6-50 um from the capture point, far = beyond 200 um (both images identical there, expected 0.5). AUC out of fold (folds grouped by view); per feature the AUC of the feature alone (folded to >= 0.5).

- near: 7239 embedded / 5735 real candidates, adversary AUC 0.550
- far: 200000 embedded / 200000 real candidates, adversary AUC 0.518

| source | feature | near AUC | far AUC |
|---|---|---|---|
| A | n_tracks_max | 0.577 | 0.508 |
| A | mean_intens | 0.559 | 0.510 |
| Bprime | n_bp_20px | 0.553 | 0.517 |
| A | n_slices | 0.537 | 0.500 |
| A | angle_spread_best | 0.535 | 0.502 |
| Bprime | darkness | 0.534 | 0.501 |
| A | angle_spread_max | 0.532 | 0.500 |
| A | darkness | 0.520 | 0.502 |
| A | z_span_um | 0.519 | 0.502 |
| Bprime | no_other | 0.510 | 0.501 |
| Bprime | dist_other_um | 0.508 | 0.500 |
| A | no_other | 0.505 | 0.501 |
| Bprime | mult | 0.504 | 0.504 |
| A | dist_other_um | 0.502 | 0.503 |
