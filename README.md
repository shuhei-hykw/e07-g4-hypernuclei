# e07-g4-hypernuclei

Geant4 track simulation of double-strangeness events in the J-PARC E07
nuclear emulsion: Xi- capture at rest on the light emulsion nuclei 12C,
14N, 16O, followed by the weak decays of the produced hypernuclei.

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
- Xi hypernucleus formation is not included.

## Weak-decay enumeration

```
python3 channels/enumerate_decays.py [--measured-only] \
  [--max-bodies N] [--outdir output]
```

For every hypernucleus in `species.csv` one Lambda decays at a time:

| Mode | Process | Products (A, Z, nL) |
|------|---------|---------------------|
| `pi-` | Lambda -> p pi- | (A, Z+1, nL-1) + pi- |
| `pi0` | Lambda -> n pi0 | (A, Z, nL-1) + pi0 |
| `nm` | Lambda N -> N N, Lambda Lambda -> Lambda N | (A, Z, nL-1), >= 2 bodies |

All partitions into particle-stable species with Q > 0 are written to
`output/decays.csv` (`parent`, `mode`, `channel`, `q_mev`, `n_bodies`,
`n_charged`, `uses_estimate`), up to `--max-bodies` (default 4,
including the pion). Double-Lambda hypernuclei decay into single-Lambda
hypernuclei (or a free Lambda), which decay in turn.

## Geant4 track simulation (`g4/`)

`emulsim` generates one channel of Xi- capture at rest with
`TGenPhaseSpace` (isotropic, unweighted n-body phase space, masses from
`output/species.csv` / `output/targets.csv`) at a random point in the
emulsion layers of an E07 emulsion module and tracks the products.
Stopped
hypernuclei decay at rest into a channel of `output/decays.csv`, again
with phase-space kinematics, so that chains such as
LL6He -> L5He + p + pi-, L5He -> 4He + p + pi- are followed.

### Build and run

Requires Geant4 (tested with 11.4.3) and ROOT (6.40).

```
python3 channels/enumerate_channels.py && python3 channels/enumerate_decays.py
source <geant4-install>/bin/geant4.sh
cd g4 && mkdir -p build && cd build
cmake .. && make -j
./emulsim nagara.mac   # Xi- + 12C -> LL6He + 4He + t, simple decay policy
./emulsim ibuki.mac    # Xi- + 14N -> L10Be + L5He, uniform decay policy
```

Macro commands:

| Command | Default | Meaning |
|---------|---------|---------|
| `/emul/density` | 3.53 | Emulsion density (g/cm3) |
| `/emul/sheetX`, `/emul/sheetY` | 345 mm, 350 mm | Sheet area |
| `/emul/thickEmulsion`, `/emul/thickBase` | 450 um, 40 um | Thick sheet: emulsion (per side), base |
| `/emul/thinEmulsion`, `/emul/thinBase` | 100 um, 180 um | Thin sheet: emulsion (per side), base |
| `/emul/numThick` | 11 | Number of thick sheets |
| `/emul/numThinUp`, `/emul/numThinDown` | 1, 1 | Thin sheets upstream / downstream |
| `/emul/minStep`, `/emul/maxStep` | 0.1 um, 10 um | Bounds of the step limit |
| `/emul/rangeFraction` | 0.01 | Step limit / estimated residual range |
| `/xi/speciesFile`, `/xi/targetsFile`, `/xi/decaysFile` | `../../output/*.csv` | Tables |
| `/xi/target` | `12C` | `12C`, `14N` or `16O` |
| `/xi/channel` | `LL6He,4He,t` | Products, comma separated, names as in `species.csv` |
| `/xi/seed` | 12345 | Seed of `gRandom` (TGenPhaseSpace) |
| `/xi/vertexMode` | `random` | `random`: uniform in the thick emulsion layers; `fixed`: `/xi/vertex` |
| `/xi/vertex` | 0 0 0 um | Capture point for the fixed mode |
| `/decay/policy` | `simple` | Choice of the decay channel (see below) |
| `/decay/channel <parent> <products>` | | Fix the channel of one parent, e.g. `LL6He L5He,p,pi-` |
| `/out/file` | `emulsim.root` | Output file |

`/emul/*` must precede `/run/initialize`. The generator checks
conservation of A, Z and S and that Q > 0.

Decay policies (branching ratios are **not** physical):

- `simple`: the pi- mode with the fewest bodies, highest Q among them
  (falls back to all modes if no pi- mode is open).
- `uniform`: any channel of `decays.csv` with equal probability.
- `uniformPiMinus`: any pi- channel with equal probability.

### Geometry

One E07 emulsion module (nominal thicknesses before development): thin
sheet + 11 thick sheets + thin sheet stacked along the beam axis z
without gaps, each sheet being emulsion | polystyrene base | emulsion,
area 345 x 350 mm2, total thickness 11.1 mm, centered at the origin
inside an air world. The thin-sheet thicknesses and the emulsion
thickness at exposure (the measured sheet thickness of about 1.1 mm in
kinema suggests ~530 um per side) should be confirmed.

A track that leaves the module is flagged (`escaped`) and killed.

### Physics

- FTFP_BERT with `G4EmStandardPhysics_option4`, production cut 1 um,
  `G4StepLimiterPhysics`. The step limit of charged particles is
  `rangeFraction` times a Bragg-Kleeman estimate of the residual range,
  clamped to [`minStep`, `maxStep`], so short tracks and track ends are
  sampled finely.
- Hypernuclei are `G4IonTable::GetIon(Z, A, nL)` ions; their dynamic mass
  is that of `species.csv`. `HypernucleusDecayPhysics` adds the at-rest
  process `HypernucleusDecay` to GenericIon and to the light hypernuclei
  predefined by Geant4 (hypertriton, hyperHe5, ...), whose own `Decay`
  is removed and which are given ion EM processes. Decays in flight are
  neglected (c tau ~ 8 cm vs. ranges of um).
- pi- stop and are absorbed by Bertini (sigma stars); free Lambdas decay
  via `G4Decay`.
- Emulsion: AgBr + gelatin with approximate mass fractions
  (`DetectorConstruction.cc`).

### Output (`tree`)

One entry per event with the capture point `vx/vy/vz` (um) and
`contained` (no charged capture or hypernuclear decay product left the
module). All tracks except e+-, gamma and neutrinos are stored as
vectors: `id`, `parent` (0 for capture products), `process` (creator:
`primary`, `HypernucleusDecay`, `hBertiniCaptureAtRest`, ...), `name`,
`z`, `a`, `nl`, `escaped`, `ke` (MeV), `px/py/pz` (MeV/c), `range`
(start-to-end distance, um), `length` (path length, um), step points
`x/y/zpos` (um) and `edep` (MeV). About 0.11 MB per event for
nagara.mac.

The tree `layers` lists the module layers (`sheet`, `emulsion`, `thick`,
`zmin`, `zmax`, `halfx`, `halfy` in um).

### Scripts

- `g4/scripts/draw_event.py file.root --event N [--zoom 40] [--neutral]
  -o ev.png`: event display relative to the capture point; top row
  +-zoom um, bottom row the whole capture / decay chain; columns XY
  (sheet plane), ZX, ZY with base films shaded and the outside of the
  module hatched. Capture products thick, decay products medium, other
  secondaries thin; escaping tracks end with an open circle.
- `g4/scripts/summary.py file.root -o summary.png`: capture-point depth
  of all / contained events and escape fraction per particle.
- `g4/scripts/check_range.py file.root --kinema <kinema dir>
  --density 3.53`: compares Geant4 ranges of the capture products with
  the kinema range-energy relation (`rangeenergy.py`) and checks
  energy-momentum conservation.

With the default emulsion model the Geant4 ranges of t and alpha are
about 5 % shorter than the kinema relation at the same density; the
density / composition should be calibrated (e.g. with 8.78 MeV alpha
from ThC') before quantitative use.
