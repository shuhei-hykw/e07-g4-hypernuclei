# e07-g4-hypernuclei

Geant4 track simulation of double-strangeness events in the J-PARC E07
nuclear emulsion: Xi- capture at rest on the light emulsion nuclei 12C,
14N, 16O, followed by the weak decays of the produced hypernuclei.

## Quick start

Everything is driven by `./run.py` (`./run.py --help`, `./run.py <command>
--help`). Requires Python 3 with numpy and PyROOT, ROOT (tested with
6.40) and Geant4 (tested with 11.4.3; `GEANT4_DIR`, default
`~/software/geant4/11.4.3`).

```
./run.py setup                       # tables + Geant4 build
./run.py sim --preset nagara --events 100
./run.py draw data/sim/nagara.root --event 0 -o ev0.png
./run.py summary data/sim/nagara.root -o summary.png
./run.py check-range data/sim/nagara.root
```

| Command | Module | Output |
|---------|--------|--------|
| `channels` | `module/enumerate_channels.py` | `output/channels_*.csv`, `species.csv`, `targets.csv` |
| `decays` | `module/enumerate_decays.py` | `output/decays.csv` |
| `build` | `module/build.py` | `g4/build/emulsim` |
| `sim` | `module/simulate.py` | `data/sim/<name>.{mac,log,root}` |
| `draw` | `module/draw_event.py` | event display PNG |
| `summary` | `module/summary.py` | containment summary PNG |
| `check-range` | `module/check_range.py` | range comparison with kinema |
| `replay` | `module/replay.py` | `output/replay/{nagara,ibuki}.csv` fixed kinematics |
| `images` | `module/synth/` | `data/images/<name>/ev<NNNN>/` microscope stacks |
| `compare` | `module/imagestats.py` | image statistics of real vs synthetic stacks |
| `calibrate` | `module/synth/calibrate.py` | `output/image_model.json` tuned image model |
| `flatfield` | `module/synth/flatfield.py` | `output/flatfield.npy` illumination profile |
| `validate` | `module/synth/validate.py` | `output/realism.md` realism on held-out real views |
| `efficiency` | `module/efficiency.py` | vertex-finding efficiency vs truth |
| `setup` | | `channels` + `decays` + `build` |

Layout: `module/` Python package, `g4/` Geant4 application (C++),
`output/` tables (tracked), `data/` generated files (not tracked).

## Channel enumeration

```
./run.py channels [--measured-only] [--max-bodies N]
```

Conservation laws: baryon number A+1, charge Z-1, two Lambdas
(S = -2). Every multiset of particle-stable species satisfying them with

    Q = M(A) + M(Xi-) - B_Xi - sum M_f > 0

is listed (up to `--max-bodies`, default 6).

### Inputs (`module/species.py`)

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
./run.py decays [--measured-only] [--max-bodies N]
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

### Running

`./run.py sim` writes a macro to `data/sim/<name>.mac` and runs emulsim
with the Geant4 environment:

```
./run.py sim --preset nagara --events 1000      # 12C -> LL6He + 4He + t
./run.py sim --preset ibuki                      # 14N -> L10Be + L5He
./run.py sim --target 16O --channel LL13B,4He --policy uniform \
  --name o16_ll13b --events 100
./run.py sim --preset nagara --fix LL6He=L5He,p,pi- \
  --set "/emul/thickEmulsion 530 um"
```

Presets `nagara_replay` and `ibuki_replay` reproduce the observed events
with their measured kinematics (see below).

`--set` passes any macro command below (`/emul/*` are placed before
`/run/initialize` automatically).

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
| `/xi/replayFile` | | Fixed momenta of the capture and of hypernuclear decays (overrides `/xi/channel` and the decay policy for the listed parents) |
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

### Plots and checks

- `./run.py draw file.root --event N [--zoom 40] [--neutral]
  -o ev.png`: event display relative to the capture point; top row
  +-zoom um, bottom row the whole capture / decay chain; columns XY
  (sheet plane), ZX, ZY with base films shaded and the outside of the
  module hatched. Capture products thick, decay products medium, other
  secondaries thin; escaping tracks end with an open circle.
- `./run.py summary file.root -o summary.png`: capture-point depth
  of all / contained events and escape fraction per particle.
- `./run.py check-range file.root [--kinema <dir>] [--density 3.53]`
  (`KINEMA_DIR`, default `~/work/kinema`): compares Geant4 ranges of the capture products with
  the kinema range-energy relation (`rangeenergy.py`) and checks
  energy-momentum conservation.

With the default emulsion model the Geant4 ranges of t and alpha are
about 5 % shorter than the kinema relation at the same density; the
density / composition should be calibrated (e.g. with 8.78 MeV alpha
from ThC') before quantitative use.

## Observed events with measured kinematics (`./run.py replay`)

`module/replay.py` converts the measured ranges and angles of Nagara
(kinema `runcard.json`, NagaraA/B/C, density 3.629) and Ibuki (kinema
`input/twin_mod047pl10/e07twin_mod047pl10_{A,B,C}.json`, density 3.605)
into momentum vectors with the kinema range-energy relation, and writes
`output/replay/<event>.csv` (vertex, name, px, py, pz). The capture
products and the listed hypernuclear decays are then generated with
exactly these momenta (`./run.py sim --preset nagara_replay`); the
capture point and the background are random as usual.

| Vertex | Nagara | Ibuki |
|--------|--------|-------|
| A (capture) | LL6He + 4He + t | L10Be + L5He |
| B | LL6He -> L5He + p + pi- (pi- energy from Q) | L10Be -> d + t + t + p + n (kinema best candidate, not unique) |
| C | L5He -> d + p + n + n (only balancing 2-prong channel) | L5He -> 4He + pi- + p |

Missing neutrons take the missing momentum (one) or share the missing
four-momentum isotropically (two, `--seed`). The Geant4 ranges of the
short tracks reproduce the measured ones within a few percent (e.g.
Nagara t 88.3 / 88.3 um, Ibuki L10Be 9.2 / 9.17 um); long p and pi-
tracks may end earlier by interactions or by leaving the module.

## Synthetic microscope images (`./run.py images`)

```
./run.py images data/sim/nagara.root --events 0-99 -j 20 --name nagara_dz3.0
./run.py images data/sim/nagara.root --events 0-99 -j 20 --dz 1.5 \
  --name nagara_dz1.5
```

Each simulated event becomes one view of the E07 full scan
(`module/synth/render.py` holds all model parameters):

- Geometry: 2048 x 2048 px at 0.289 um/px, 60 slices. `--dz` is the
  physical focal step per slice in the developed emulsion (default 3.0
  um; 1.5 um is the alternative under discussion); emulsion layers
  shrink along z by `--shrink` (1.9) on development, base films do not.
  The capture point is placed at a random x, y in the view (100 px from
  the edges) and the stack at a random depth inside its emulsion layer
  (`--centered` puts it in the middle).
- Grains: along each charged track, Poisson with a linear density
  saturating in dE/dx (14.8 grains / 100 um for a MIP, 2.5 / um at
  saturation), only inside emulsion layers; fog 6e-3 grains / um3;
  lognormal grain darkness.
- Optics: a grain at defocus dz is a Gaussian of width
  sqrt(s0^2 + (k dz)^2) (s0 = 0.35 um, k = 0.35) with constant
  integrated density; I = B (t_min + (1 - t_min) exp(-sum D)) with a
  mottled brightness B (two scales), t_min = 0.4 and pixel noise.
- Background (densities per mm2 of view, tuned by `calibrate`): nearly
  parallel beam MIP tracks around an axis close to z, isotropic straight
  tracks, flat tracks nearly in the sheet plane (the many thin lines of
  real single slices), large junk clusters; mottle at two scales, a
  static fine texture, temporal and fixed-pattern pixel noise; the flat
  field of the real microscope.

Files per event:

| File | Content |
|------|---------|
| `0000.png` ... | 8-bit grayscale slices |
| `image.json` | full-scan JSON (`SpngReader`): `AffineP2S` in mm, stage z per slice (`--stage-dz`, default -0.003 mm as tohoku-v1), `Sh` |
| `truth.json` | col / row / slice and `polyline` are image (developed) coordinates, `um` the Geant4 (undeveloped) position; vertices (`capture`, `HypernucleusDecay`, `hBertiniCaptureAtRest`, `Decay`) with col / row / slice and `in_view`; charged tracks with origin, start / end, `vertex_start` / `vertex_end`, 3-D `polyline` (2 um steps); dataset parameters |
| `labels.npz` | per slice: `labels` class of in-focus grains (1 capture, 2 decay, 3 other secondary, 4 background track, 5 junk) and `instance` (track id; -1-i background track i, unique over beam / isotropic / flat tracks; -100001-i junk cluster i); grains within half a slice of the focal plane are painted as 2 px disks, lower classes win where disks overlap, fog is not labelled |
| `grains.npz` | per-grain truth: `col`, `row`, `slice` (view coordinates), `x_um`, `y_um`, `z_um` (developed um from the view corner and slice 0), `weight` (darkness), `cls` / `instance` as in the labels (fog 0 / 0), `fog` flag, `focus_slice` |
| `background_tracks.json` | straight background tracks: instance id, kind (beam, isotropic, flat), grains per um and straight segments [[col, row, slice] start, end], one per emulsion layer crossed, inside the emulsion near the view |
| `preview.png` | minimum projection: crop around the capture point and the whole view |
| `g4_tracks.png` | Geant4 event display before imaging (as `./run.py draw`) |
| `g4_view.png` | Geant4 tracks in the image coordinates of the view (crop as `preview.png`, whole view, column vs slice) |

Rendering takes about 70 s and 170 MB per event (`-j` renders events in
parallel). `--tracks-only` writes only the two Geant4 track pictures for
an existing dataset (the view placement is reproduced from the seed).

### Embedding into real views (`--background`)

```
./run.py images data/sim/eval_nagara.root --events 0-4 -j 5 \
  --background random --name embed_nagara_dz3.0
```

The Geant4 signal is rendered alone (same grain and optics model) and
its transmission multiplies the slices of a real tohoku-v1 view
(`/home/hayakawa/microscope-image/.../AREA00/IMAGE00_AREA00`, read
only), so fog, beam tracks, junk and noise are real and the candidate
background of the full scan is preserved. The stack size and stage z
follow the real view, `truth.json` records `background_view`, and the
labels cover the signal only. Views used for CNN labels (V00000011,
V00000991, V00001012, V00001339) are never picked.

### Calibration against real views

`./run.py compare --real <view.json> --synth <event dir> -o cmp.png`
computes, for central slices and a central crop: intensity percentiles,
blobs of the darkest 5 % of pixels, connected components after a
full-scan-like preprocessing (Gaussian background subtraction, Otsu),
the radial power spectrum and the correlation of background-subtracted
slices. In real views this correlation stays high over many slices
(about 0.7 at 1-5 slices, 0.33 at 13, 0.05 at 40): the nearly
perpendicular, dense beam tracks appear at almost the same place in
consecutive slices. The calibration therefore matches the correlation
at offsets 1, 2, 5 and 13, which separates beam density / slope /
spread from the defocus blur of the grains.

`./run.py calibrate --views 6 --iterations 150` averages these
statistics over random real views (`output/real_image_stats.json`) and
searches the model parameters (brightness, mottle, residual
transmission, temporal and fixed-pattern noise, grain image width,
defocus blur, grain density, fog density, grain darkness spread) that
reproduce them; the result `output/image_model.json` is used by
`./run.py images` unless `--model none`. Several parameters are moved
at once in the search (e.g. wider grains need weaker defocus blur and
more fog together).

Grain darkness is constrained by the depth of the darkest blobs
relative to their local background (real: median 0.19, 90 % 0.32).
As a check of the signal, `module/synth/profile.py` measures the
transverse profile of tracks on minimum projections: a heavy track of
the real Ibuki image (specials_x20, another scanner) has depth 0.22 and
FWHM 4 px, the simulated Ibuki replay gives 0.17-0.19 / 3-4 px for t
and L5He and 0.36 / 4 px for alpha, while p and pi- (0.03-0.04) are as
faint as real MIP tracks. The tuned model (v3) puts the beam density at
about 1e4 tracks / mm2 nearly perpendicular to the sheets (slope ~0.01),
consistent with the long-range slice correlation.

Image model v14 (`output/image_model.json`, the default, calibrated for
a physical step of 1.5 um per slice, now the default `--dz`) adds to the
model flat tracks nearly in the sheet plane, fine mottle of two depth
correlation lengths (`mottle_z_um`, `mottle_long`), a static fine
texture and an optics / camera blur of the transmitted image (0.6 px);
straight background tracks have a fixed length (1 mm) and flat /
isotropic ones a density per volume, so that the grain density does not
depend on the view size or stack depth. The calibration matches,
besides the statistics above, elongated components, the raw power
spectrum in three bands, the shape of the background-subtracted
intensity distribution, the near-saturated fraction and the minimum
projection; its loss is the mean squared z score in units of the real
view-to-view spread, with patch statistics averaged over three seeds.

Since v16 the calibration also matches the response of point-cloud
methods, which the background statistics do not see: line-like pixels
of the Hessian of the background-subtracted slice (`line_thin`, scale
1.5 px, strong; `line_wide`, 3 px, weak), DoG maxima per Mpx
(`dog_per_mpx`) and their fraction on lines (`dog_on_line`). Real flat
tracks look like smooth continuous lines of varying width with many
isolated grains around, v15 rendered them as beaded, sharp lines with
fewer isolated points (DoG 1.0e4 vs 1.2e4 per Mpx, 0.21 vs 0.17 on
lines). The lower bound of the flat-track grain density (`flat_ion_min`)
is tuned as well. Model v15 (`output/image_model_v15_dz4.json`) is v14
recalibrated at dz 4.0 (held-out: 20 / 23 within 3 spreads,
`output/realism_v15_dz4.md`).

`./run.py validate data/images/<dataset>` compares a synthetic dataset
with 20 real views not used in the calibration (`output/realism.md`).
For v14 all 23 statistics are within 3 spreads of the real views
(largest deviations: slice correlation at 13 slices z = -2.8,
minimum-projection detail z = +2.6, high-frequency power z = +2.5). v14
is v13 refined with `--only` on the grain, blur, flat-track, junk,
texture and fog parameters.
With synthetic images method A of e07-fullscan finds 800-900
candidates per view, as on real views (~770; model v3: ~40). Heavy
signal tracks have depth 0.11-0.22 and width 3-4 px on minimum
projections (real Ibuki heavy track 0.22 / 4 px). The v3_* evaluation
datasets were made with model v3 (`data/image_model_v3.json`).

`./run.py flatfield` stores the median middle slice of 20 real views
(smoothed, normalized to 1 at the center) in `output/flatfield.npy`:
the circular field stop and the vignetting, applied to the brightness
of fully synthetic images.

## Vertex-finding efficiency (`./run.py efficiency`)

```
./run.py efficiency data/images/nagara_dz3.0 candidates.csv \
  --tolerances 6,12,30 --details rows.csv -o eff.png
```

`candidates.csv` holds the merged vertex candidates of a full-scan
method (columns `event` or `view_id`, `vx_px`, `vy_px`, `z_min`/`z_max`
or `z_mean` in stage mm, optional `n_tracks_max`). The matching
specification, agreed with the e07-fullscan analysis:

- Truth: in-view vertices of `truth.json`.
- Distance in the developed emulsion (`--space undeveloped`: z x Sh):
  x, y = pixels x 0.289 um, z = slice difference x `physical_dz_um`.
- Candidate z ranges: `edge` (primary) uses the distance to the nearer
  end of the range, `inclusive` counts the inside of the range as zero;
  point candidates use the plain distance.
- One-to-one greedy matching per tolerance; the optimal assignment is
  computed for comparison; unmatched candidates within the tolerance of
  a truth are reported as duplicates.
- Control: mean over 10 random shifts of 400-1000 px of the truths.
- Full-scan criterion (200 px XY, depth within 30 um of the z range) as
  a loose recovery measure; "LL event" = capture and double-Lambda decay
  vertex both found; rank by `n_tracks_max`.
- `./run.py beamgrains <view.json> ...`: grains per slice along
  perpendicular beam tracks, a step handle independent of the beam
  divergence and defocus (output/dz_comparison.md).
- `./run.py layerdepth <image.json> ...`: physical step from the number
  of slices spanned by whole emulsion layers in deep stacks (450 um / Sh
  developed thickness); independent of the grain density. On the
  200-slice special events it gives 3.6-4.4 um per slice (mean 4.0 um
  for the 3 um stage step, `output/layerdepth.txt`), which excludes the
  1.5 um step.
- `./run.py rankstudy <dataset> <bprime.csv>`: ranking of point
  candidates by multiplicity, local density or image darkness at fixed
  candidates per view (`output/bprime_ranking.md`), with event-level
  binomial uncertainties; ties are averaged in the rank fractions and
  broken at random (fixed seed) at the top-K boundary.
- `./run.py specialsrank`: the same rankings on the B' branch points of
  the real special events (`output/specials_ranking.md`): integer rank,
  number of matching points and XY / z offset of the best one, control
  distribution with positions common to all scores (no match = rank
  fraction 1), at a loose (200 px) and a strict (20 px) XY tolerance.
  The darkness model was checked on IBUKI, so these events are not a
  fully independent test.
- `./run.py zbias <dataset> <csv> --sim <root>`: offset z_best - z_true
  of the XY-nearest candidate, split by layer / stack position, daughter
  z asymmetry and n_tracks_max.
- `--background-candidates <csv>`: background control with the method
  run on the original real views (lists in `output/background_views_*`).
- `--nearest`: XY distance to the nearest candidate in XY and the edge
  z gap of that candidate, true vs control (separates XY accuracy from
  the depth convention).

`./run.py evalreport` writes the user-facing summary
`output/efficiency_report.md` (dz 3.0 / 1.5 comparison for method A with
the three depth conventions, XY accuracy, z_best bias, B' vs A).

`./run.py report [--pattern 'v3_*']` evaluates every candidate file
`<dataset>_<method>_<cuts>.csv` of the e07-fullscan analysis
(`~/work/e07-fullscan/results/sim_eval`) and writes
`output/efficiency_report.md`: efficiency (chance level) of the capture
vertex, the double-Lambda decay vertex and the LL event per tolerance,
duplicates and candidates per event. With real backgrounds (about 800
method-A candidates per view, z ranges of ~50 slices) the chance level
at 30 um reaches 50-90 %, so only the tight tolerances are meaningful.
