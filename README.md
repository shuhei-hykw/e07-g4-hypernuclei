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

`output/targets.csv` holds the initial mass M(A) + M(Xi-) - B_Xi used by
the Geant4 generator.

### Limitations

- Only ground states; excited hypernuclei / fragments (e.g. 12Be_LL*)
  shift Q by their excitation energy and are not listed separately.
- Estimated B_Lambda and B_LambdaLambda carry ~1 MeV uncertainties;
  channels with small Q and `uses_estimate = 1` may be closed.
- Xi hypernucleus formation and weak decays of the products are not
  included yet.

## Geant4 track simulation (`g4/`)

`emulsim` generates one channel of Xi- capture at rest with
`TGenPhaseSpace` (isotropic, unweighted n-body phase space, masses from
`output/species.csv` / `output/targets.csv`) and tracks the products in
a uniform emulsion block with the capture point at its center.

### Build and run

Requires Geant4 (tested with 11.4.3) and ROOT (6.40).

```
source <geant4-install>/bin/geant4.sh
cd g4 && mkdir -p build && cd build
cmake .. && make -j
./emulsim nagara.mac        # Xi- + 12C -> LL6He + 4He + t, 1000 events
```

Macro commands:

| Command | Default | Meaning |
|---------|---------|---------|
| `/emul/density` | 3.53 | Emulsion density (g/cm3), before `/run/initialize` |
| `/emul/maxStep` | 0.1 um | Step limit for fine track sampling |
| `/xi/speciesFile`, `/xi/targetsFile` | `../../output/*.csv` | Mass tables |
| `/xi/target` | `12C` | `12C`, `14N` or `16O` |
| `/xi/channel` | `LL6He,4He,t` | Products, comma separated, names as in `species.csv` |
| `/xi/seed` | 12345 | Seed of `gRandom` (TGenPhaseSpace) |
| `/out/file` | `emulsim.root` | Output file |

The generator checks conservation of A, Z and S and that Q > 0.

### Physics

- FTFP_BERT with `G4EmStandardPhysics_option4`, production cut 1 um,
  `G4StepLimiterPhysics`.
- Hypernuclei are `G4IonTable::GetIon(Z, A, nL)` ions; their dynamic mass
  is overridden with the value of `species.csv`. They are not decayed
  yet (they stop at the end of their range).
- Emulsion: AgBr + gelatin with approximate mass fractions
  (`DetectorConstruction.cc`).

### Output (`tree`)

One entry per event; per-primary vectors: `name`, `z`, `a`, `nl`,
`ke` (MeV), `px/py/pz` (MeV/c), `range` (start-to-end distance, um),
`length` (path length, um), step points `x/y/zpos` (um) and `edep`
(MeV).

### Scripts

- `g4/scripts/draw_event.py file.root --event N -o ev.png`: XY / XZ
  projections of one event.
- `g4/scripts/check_range.py file.root --kinema <kinema dir>
  --density 3.53`: compares Geant4 ranges with the kinema range-energy
  relation (`rangeenergy.py`) and checks energy-momentum conservation.

With the default emulsion model the Geant4 ranges of t and alpha are
about 5 % shorter than the kinema relation at the same density; the
density / composition should be calibrated (e.g. with 8.78 MeV alpha
from ThC') before quantitative use.
