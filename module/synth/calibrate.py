"""Tune the image model to the statistics of real full-scan views.

  ./run.py calibrate --views 6 --iterations 120

Target statistics (module/imagestats.py) are averaged over randomly
chosen real views (central slices, central crop) and stored in
output/real_image_stats.json. Background-only synthetic patches (fog,
beam and random tracks, junk) are rendered with trial parameters and a
coordinate / random search minimizes the squared relative deviations.
The best parameters are written to output/image_model.json, which
./run.py images uses by default.
"""

import argparse
import json
import random
from dataclasses import replace
from pathlib import Path

import numpy as np

from .. import config
from .. import imagestats as st
from . import render

REAL_DIR = ("/home/hayakawa/microscope-image/tohoku/fullscan/E07/MOD108/"
            "PL12/tohoku-v1/AREA00/IMAGE00_AREA00")
CROP_PX = 1024
HALF_SLICES = 2
# Slice offsets from the central slice used for the statistics; the
# long lags probe the nearly perpendicular beam tracks
OFFSETS = [-2, -1, 0, 1, 2, 5, 13]
CORR_LAGS = (1, 2, 5, 13)
# Statistic -> (weight, compare in log)
TARGETS = {
  "median": (3., False), "std": (2., False), "p1": (1., False),
  "p99": (3., False), "frac_bright": (10., False),
  "p5": (1., False), "n_cc_per_mpx": (1., True),
  "cc_area_p90": (2., True), "dark_blobs_per_mpx": (1., True),
  "dark_blob_area_median": (1., True),
  "dark_contrast_median": (5., False), "dark_contrast_p90": (5., False),
  "corr1": (3., False),
  "corr2": (3., False), "corr5": (3., False), "corr13": (6., False),
  "fg_fraction": (5., False), "n_elong_per_mpx": (3., True),
  "pow_low": (3., False), "pow_mid": (3., False), "pow_high": (3., False),
  "d_skew": (3., False), "d_kurt": (2., True), "d_iqr": (3., False),
  "mp_detail_std": (5., False), "mp_dark_frac": (5., False),
  "line_thin": (3., False), "line_wide": (3., False),
  "dog_per_mpx": (3., True), "dog_on_line": (3., False),
}
# Parameter -> (object, field, low, high, log scale)
PARAMS = {
  "brightness": ("background", "brightness", 150., 235., False),
  "mottle_large": ("background", "mottle", 0.0, 0.25, False),
  "mottle_small": ("background", "mottle", 0.0, 0.3, False),
  "mottle_small_um": ("background", "mottle", 1.5, 15., True),
  "mottle_z_um": ("background", "mottle_z_um", 2., 200., True),
  "mottle_long_rel": ("background", "mottle_long", 0.0, 0.2, False),
  "mottle_long_z_um": ("background", "mottle_long", 10., 600., True),
  "min_transmission": ("background", "min_transmission", 0.0, 0.7, False),
  "noise_sigma": ("background", "noise_sigma", 0.3, 6., True),
  "fixed_noise_sigma": ("background", "fixed_noise_sigma", 0.0, 2., False),
  "beam_slope_max": ("background", "beam_slope_max", 0.005, 0.3, True),
  "beam_spread": ("background", "beam_spread", 0.001, 0.1, True),
  "grain_sigma_um": ("optics", "grain_sigma_um", 0.15, 1.0, True),
  # Optics / camera MTF: a blur of at least ~0.5 px is physical
  "image_blur_px": ("optics", "image_blur_px", 0.5, 0.8, False),
  # Lower bound: defocus blur of an objective with NA >= ~0.1
  "defocus_slope": ("optics", "defocus_slope", 0.03, 0.5, True),
  "grain_density": ("optics", "grain_density", 0.03, 4., True),
  "fog_per_um3": ("emulsion", "fog_per_um3", 5e-4, 4e-2, True),
  "fog_weight": ("emulsion", "fog_weight", 0.5, 6., True),
  "grain_spread": ("emulsion", "grain_spread", 0.0, 1.0, False),
  "beam_per_mm2": ("background", "beam_per_mm2", 100., 5e4, True),
  "random_per_mm2": ("background", "random_per_mm2", 10., 1e4, True),
  "random_ion_max": ("background", "random_track_ionization", 1.5, 20.,
                     True),
  "blobs_per_mm2": ("background", "blobs_per_mm2", 1., 3e3, True),
  "flat_per_mm2": ("background", "flat_per_mm2", 1., 2e4, True),
  "flat_dip_max": ("background", "flat_dip_max", 0.01, 0.5, True),
  "flat_ion_max": ("background", "flat_ionization", 1.5, 30., True),
  "flat_ion_min": ("background", "flat_ionization", 0.3, 10., True),
  "speckle_rel": ("background", "speckle_rel", 0.0, 0.3, False),
  "speckle_scale_px": ("background", "speckle_scale_px", 0.5, 4., True),
}
MOTTLE_SCALES_UM = (25., 6.)
MAX_MOVED = 3
PATCH_SEEDS = (1, 2, 3)
# Minimum spread used for z scores: fraction of the target value (or of
# Z_FLOOR_ABS for values near zero)
Z_FLOOR_REL = 0.001
Z_FLOOR_ABS = 0.01
# Statistics compared by absolute difference (fractions, correlations)
ABSOLUTE_STATS = ("corr1", "corr2", "corr5", "corr13", "fg_fraction",
                  "dark_contrast_median", "dark_contrast_p90",
                  "pow_low", "pow_mid", "pow_high", "d_skew",
                  "frac_bright", "mp_dark_frac", "line_thin", "line_wide",
                  "dog_on_line")
MODEL_FILE = "image_model.json"
# Real views labelled by hand for the full-scan CNN: never used here
EXCLUDED_VIEWS = ("V00000011", "V00000991", "V00001012", "V00001339")


def real_targets(n_views, seed):
  views = sorted(p for p in Path(REAL_DIR).glob("V*.json")
                 if not p.name.startswith(EXCLUDED_VIEWS))
  random.Random(seed).shuffle(views)
  stats = []
  for v in views[:n_views]:
    c = st.n_slices(v) // 2
    idx = [c + o for o in OFFSETS]
    stats.append(flatten(st.stack_stats(st.read_stack(v, idx), CROP_PX,
                                        OFFSETS)))
  means = {k: float(np.mean([s[k] for s in stats])) for k in TARGETS}
  stds = {k: float(np.std([s[k] for s in stats])) for k in TARGETS}
  return means, stds, [str(v) for v in views[:n_views]]


def flatten(s):
  out = {k: v for k, v in s.items() if not k.startswith("_")}
  corr = out.pop("slice_corr")
  for lag in CORR_LAGS:
    out[f"corr{lag}"] = corr.get(lag, 0.)
  return out


def apply_values(values, optics, emulsion, background):
  """Models with the tuned parameter values substituted."""
  objs = {"optics": optics, "emulsion": emulsion, "background": background}
  for name, v in values.items():
    obj, field, *_ = PARAMS[name]
    if name == "random_ion_max":
      objs[obj] = replace(objs[obj], random_track_ionization=(1., v))
    elif name == "flat_ion_max":
      lo = min(values.get("flat_ion_min", 1.), v)
      objs[obj] = replace(objs[obj], flat_ionization=(lo, v))
    elif name == "flat_ion_min":
      continue
    elif name in ("mottle_long_rel", "mottle_long_z_um"):
      continue
    elif name == "mottle_z_um" or not name.startswith("mottle"):
      objs[obj] = replace(objs[obj], **{field: v})
  objs["background"] = replace(objs["background"], mottle=(
    (values["mottle_large"], MOTTLE_SCALES_UM[0]),
    (values["mottle_small"],
     values.get("mottle_small_um", MOTTLE_SCALES_UM[1]))),
    mottle_long=(values.get("mottle_long_rel", 0.),
                 values.get("mottle_small_um", MOTTLE_SCALES_UM[1]),
                 values.get("mottle_long_z_um", 40.)))
  return objs["optics"], objs["emulsion"], objs["background"]


def load_model(path):
  """Tuned parameter values from output/image_model.json."""
  return json.loads(open(path).read())["params"]


DZ_UM = render.Optics.dz_um  # physical focal step (set by --dz)


def make_models(values):
  optics = replace(render.Optics(), size_px=CROP_PX,
                   n_slices=2 * max(OFFSETS) + 1, dz_um=DZ_UM)
  background = replace(render.Background(), n_dead_pixels=0)
  return apply_values(values, optics, render.Emulsion(), background)


class _Layer:
  def __init__(self, emulsion, zmin, zmax):
    self.emulsion, self.zmin, self.zmax = emulsion, zmin, zmax


def render_patch(values, seed):
  """Background-only stack with the trial parameters."""
  optics, emulsion, background = make_models(values)
  rng = np.random.default_rng(seed)
  # A thick emulsion layer around z = 0 (undeveloped um)
  devz = render.DevelopedZ([_Layer(True, -500., 500.)], optics.shrink)
  view = render.View(0., 0., float(devz.to_dev(np.array([0.]))[0]),
                     devz, optics)
  grains = [render.fog_grains(view, emulsion, rng)]
  n_fog = len(grains[0])
  for maker in (render.beam_grains, render.random_track_grains,
                render.flat_track_grains):
    g, _ = maker(view, background, emulsion, 0., rng)
    grains.append(view.to_view(g[devz.in_emulsion(g[:, 2])]))
  blobs, _ = render.blob_grains(view, background, rng)
  grains.append(blobs)
  allg = np.vstack(grains)
  weights = render.grain_weights(len(allg), emulsion, rng)
  weights[:n_fog] *= emulsion.fog_weight
  renderer = render.Renderer(optics)
  c = optics.n_slices // 2
  mottle = render.MottleVolume(optics, background, rng,
                               view.slice_zdev(c + min(OFFSETS)),
                               view.slice_zdev(c + max(OFFSETS)))
  fixed = render.fixed_pattern(optics, background, rng)
  return [render.compose(mottle.at(view.slice_zdev(c + o)),
                         renderer.density(allg, weights, c + o),
                         background, fixed, rng, optics.image_blur_px)
          for o in OFFSETS]


def patch_stats(values):
  """Statistics averaged over patches of several seeds: one patch
  fluctuates far more than the real view-to-view spread."""
  runs = [flatten(st.stack_stats(render_patch(values, seed), CROP_PX,
                                 OFFSETS)) for seed in PATCH_SEEDS]
  return {k: float(np.mean([r[k] for r in runs])) for k in runs[0]}


def loss(stats, targets, stds=None):
  """Sum of squared deviations in units of the view-to-view spread of the
  real statistic (z scores); without spreads, relative deviations."""
  if stds is not None:
    total = 0.
    for k in TARGETS:
      floor = Z_FLOOR_REL * max(abs(targets[k]), Z_FLOOR_ABS)
      z = (stats[k] - targets[k]) / max(stds[k], floor)
      total += z * z
    return total / len(TARGETS)
  total = 0.
  for k, (w, logc) in TARGETS.items():
    a, b = stats[k], targets[k]
    if logc:
      d = np.log(max(a, 1e-3)) - np.log(max(b, 1e-3))
    else:
      d = (a - b) if k in ABSOLUTE_STATS else (a - b) / max(abs(b), 1.)
    total += w * d * d
  return total


def start_values(path):
  """Previous tuning result if any, else the built-in defaults."""
  values = current_defaults()
  if path.exists():
    values.update({k: v for k, v in load_model(path).items()
                   if k in PARAMS})
  # Keep the start inside the (possibly changed) search ranges
  return {k: float(np.clip(v, PARAMS[k][2], PARAMS[k][3]))
          for k, v in values.items()}


def current_defaults():
  o, e, b = render.Optics(), render.Emulsion(), render.Background()
  return {"brightness": b.brightness, "mottle_large": b.mottle[0][0],
          "mottle_small": b.mottle[1][0],
          "mottle_small_um": b.mottle[1][1],
          "mottle_z_um": b.mottle_z_um if b.mottle_z_um > 0 else 30.,
          "mottle_long_rel": b.mottle_long[0],
          "mottle_long_z_um": b.mottle_long[2],
          "min_transmission": b.min_transmission,
          "noise_sigma": b.noise_sigma,
          "fixed_noise_sigma": b.fixed_noise_sigma,
          "grain_sigma_um": o.grain_sigma_um,
          "image_blur_px": o.image_blur_px,
          "defocus_slope": o.defocus_slope, "grain_density": o.grain_density,
          "fog_per_um3": e.fog_per_um3, "grain_spread": e.grain_spread,
          "fog_weight": e.fog_weight,
          "beam_per_mm2": b.beam_per_mm2,
          "beam_slope_max": b.beam_slope_max,
          "beam_spread": b.beam_spread,
          "random_per_mm2": b.random_per_mm2,
          "random_ion_max": b.random_track_ionization[1],
          "blobs_per_mm2": b.blobs_per_mm2,
          "flat_per_mm2": b.flat_per_mm2, "flat_dip_max": b.flat_dip_max,
          "flat_ion_max": b.flat_ionization[1],
          "flat_ion_min": b.flat_ionization[0],
          "speckle_rel": b.speckle_rel,
          "speckle_scale_px": b.speckle_scale_px}


def perturb(values, rng, scale, names=None):
  """Move one to MAX_MOVED parameters at once (joint moves are needed,
  e.g. wider grains with weaker defocus blur and more fog)."""
  out = dict(values)
  pool = list(names or PARAMS)
  n = int(rng.integers(1, min(MAX_MOVED, len(pool)) + 1))
  for name in rng.choice(pool, n, replace=False):
    _, _, lo, hi, logs = PARAMS[name]
    if logs:
      v = np.exp(np.log(out[name]) + rng.normal(0., scale)
                 * np.log(hi / lo))
    else:
      v = out[name] + rng.normal(0., scale) * (hi - lo)
    out[name] = float(np.clip(v, lo, hi))
  return out


def main():
  parser = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument("--views", type=int, default=6)
  parser.add_argument("--iterations", type=int, default=120)
  parser.add_argument("--seed", type=int, default=0)
  parser.add_argument("--dz", type=float, default=render.Optics.dz_um,
                      help="physical focal step per slice (um)")
  parser.add_argument("--z-loss", action="store_true",
                      help="loss in units of the real view-to-view spread")
  parser.add_argument("--only", default="",
                      help="comma separated parameters to move (default all)")
  parser.add_argument("--exclude", default="",
                      help="comma separated target statistics to ignore")
  parser.add_argument("--output", default=str(config.OUTPUT_DIR
                                              / MODEL_FILE),
                      help="model file (also the start point if present)")
  args = parser.parse_args()
  global DZ_UM
  DZ_UM = args.dz
  for k in filter(None, args.exclude.split(",")):
    TARGETS.pop(k)

  target_file = config.OUTPUT_DIR / "real_image_stats.json"
  targets, stds, views = real_targets(args.views, args.seed)
  target_file.write_text(json.dumps({"views": views, "stats": targets,
                                     "spread": stds}, indent=1))
  if not args.z_loss:
    stds = None
  print("targets:", {k: round(v, 3) for k, v in targets.items()})

  rng = np.random.default_rng(args.seed)
  best = start_values(Path(args.output))
  best_stats = patch_stats(best)
  best_loss = loss(best_stats, targets, stds)
  print(f"start loss {best_loss:.3f}")
  for it in range(args.iterations):
    scale = 0.3 * (1. - it / args.iterations) + 0.03
    trial = perturb(best, rng, scale,
                    [p for p in args.only.split(",") if p] or None)
    s = patch_stats(trial)
    l = loss(s, targets, stds)
    if l < best_loss:
      best, best_stats, best_loss = trial, s, l
      print(f"[{it}] loss {l:.3f} "
            + " ".join(f"{k}={v:.3g}" for k, v in best.items()), flush=True)
  print("final stats vs targets:")
  for k in TARGETS:
    print(f"  {k:22s} {best_stats[k]:10.4g} {targets[k]:10.4g}")
  Path(args.output).write_text(json.dumps(
    {"params": best, "loss": best_loss, "stats": best_stats,
     "targets": targets, "dz_um": DZ_UM}, indent=1))
  print(f"-> {args.output}")


if __name__ == "__main__":
  main()
