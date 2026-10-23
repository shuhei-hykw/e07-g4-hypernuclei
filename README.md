# e07-g4-hypernuclei

Geant4 track simulation of double-strangeness events in the J-PARC E07
nuclear emulsion. Step 1: enumerate every energetically allowed final
state of Xi- capture at rest on the light emulsion nuclei 12C, 14N, 16O.

## Channel enumeration

```
python3 channels/enumerate_channels.py [--measured-only] \
  [--max-bodies N] [--outdir output]
```

Conservation laws: baryon number A+1, charge Z-1, two Lambdas
(S = -2). Every multiset of particle-stable species satisfying them with

    Q = M(A) + M(Xi-) - B_Xi - sum M_f > 0

is listed (up to `--max-bodies`, default 6).

### Inputs (`channels/species.py`)

| Item | Source / assumption |
|------|---------------------|
| Nuclear masses | AME2020 mass excess (unbound 5He, 5Li, 6Be, 7He, 8Be, 9B used only as hypernuclear cores) |
| B_Lambda | Emulsion compilation (Juric 1973 / Davis 2005); 7He_L and 16N_L from JLab |
| B_Xi | Atomic 3D state: 0.126 (C), 0.174 (N), 0.231 (O) MeV (`TARGETS`) |
| Unmeasured B_Lambda | Isobar mean of measured B_Lambda, interpolated in A (flagged `estimated`) |
| B_LambdaLambda | Measured: 6He_LL = 6.91 MeV (Nagara). Others: 2 B_Lambda(A-1) + 0.67 MeV (flagged) |

Hypernuclei (measured or estimated) are kept only if they are stable
against any strong breakup into lighter species of the database (e.g.
6Li_L is rejected because 6Li_L -> 5He_L + p is open).

### Outputs (`output/`)

- `species.csv`: particle-stable species with masses.
- `channels_<target>.csv`: one row per channel with columns
  `category` (`double`, `twin`, `single+L`, `free-LL`), `channel`,
  `q_mev`, `n_bodies`, `n_charged`, `n_neutral`, `p_2body_mev`
  (CM momentum for two-body states), `uses_estimate`.

A single-body final state (e.g. `LL13B + gamma`) means compound formation
followed by radiative de-excitation.

### Limitations

- Only ground states; excited hypernuclei / fragments (e.g. 12Be_LL*)
  shift Q by their excitation energy and are not listed separately.
- Estimated B_Lambda and B_LambdaLambda carry ~1 MeV uncertainties;
  channels with small Q and `uses_estimate = 1` may be closed.
- Xi hypernucleus formation and weak decays of the products are not
  included yet.
