"""Depth profile of isolated grains split by the synthetic truth.

  ~/work/.venv/bin/python -m module.depthsplit <event dir> [...]

Runs the truth-free isolated-grain depth profile of e07-fullscan
(scripts/grain_depth_profile.py: DoG grains with no other DoG maximum
within 10 px over +-3 slices, DoG response at +-8 slices normalised to
the peak, FWHM in slices, open = not below half on both sides) and
labels each grain of a synthetic view by the heaviest truth grain
(grains.npz) within MATCH_PX in xy and MATCH_SLICES in depth: fog,
background track (beam, isotropic and flat), signal or none (embedded
real background, junk). Real views pass through unlabelled ("all").
This separates grain defocus from the depth extent of nearly
perpendicular tracks, which also pass the isolation test. A truth-free
cut does the same in real views: "compact" grains have fallen below
COMPACT_LEVEL of their peak at both ends of the +-8 slice window.
"""

import argparse
from pathlib import Path

import numpy as np
from scipy.ndimage import maximum_filter
from scipy.spatial import cKDTree

from . import linecheck
from .synth import render

MATCH_PX = 2.
MATCH_SLICES = 2.
OFFSETS_SHOWN = (-8, -4, 4, 8)
COMPACT_LEVEL = 0.2


def grain_profiles(gp, stack_dir):
  """Per-grain (x, y, z) and normalised profile, as in gp.measure."""
  ep = gp.ep
  p = Path(stack_dir)
  reader = gp.load_spng(p if p.suffix == ".json" else p / "image.json")
  r = gp.responses(reader)
  n, h, w = r.shape
  hits, _ = ep.dog_hits(reader, ep.DOG_K)
  peaks, grains = hits["dog"], hits["dog3d"].astype(int)
  tree = cKDTree(peaks[:, :2])
  e, half = gp.EDGE_PX, gp.PROFILE_HALF
  ok = ((grains[:, 0] >= e) & (grains[:, 0] < w - e)
        & (grains[:, 1] >= e) & (grains[:, 1] < h - e)
        & (grains[:, 2] >= half) & (grains[:, 2] < n - half))
  iso = np.zeros(len(grains), bool)
  for i in np.flatnonzero(ok):
    nb = tree.query_ball_point(grains[i, :2], gp.ISO_R_PX)
    dz = np.abs(peaks[nb, 2] - grains[i, 2])
    d = np.hypot(*(peaks[nb, :2] - grains[i, :2]).T)
    iso[i] = not ((dz <= gp.ISO_DZ) & (d > gp.PEAK_HALF_PX + 1)).any()
  idx = np.flatnonzero(iso)
  rng = np.random.default_rng(gp.SEED)
  if idx.size > gp.MAX_GRAINS:
    idx = np.sort(rng.choice(idx, gp.MAX_GRAINS, replace=False))
  k = 2 * gp.PEAK_HALF_PX + 1
  box = maximum_filter(r, size=(1, k, k))
  offs = np.arange(-half, half + 1)
  x, y, z = grains[idx, 0], grains[idx, 1], grains[idx, 2]
  prof = box[z[:, None] + offs[None, :], y[:, None], x[:, None]]
  prof = prof / prof[:, half:half + 1]
  return np.c_[x, y, z], prof


def labels(ev_dir, xyz):
  path = Path(ev_dir) / "grains.npz"
  if not path.exists():
    return np.full(len(xyz), "all")
  g = np.load(path)
  pts = np.c_[g["col"], g["row"]]
  tree = cKDTree(pts)
  out = []
  for x, y, z in xyz:
    nb = np.array(tree.query_ball_point((x, y), MATCH_PX), int)
    nb = nb[np.abs(g["focus_slice"][nb] - z) <= MATCH_SLICES]
    if not len(nb):
      out.append("none")
      continue
    j = nb[np.argmax(g["weight"][nb])]
    c = g["cls"][j]
    out.append("fog" if g["fog"][j] else
               "background track" if c == render.LABEL_BEAM else
               "junk" if c == render.LABEL_JUNK else "signal")
  return np.array(out)


def summary(gp, prof):
  w = np.array([gp.fwhm(p) for p in prof])
  med = np.median(prof, axis=0)
  half = gp.PROFILE_HALF
  return (f"n={len(prof)} FWHM {np.nanmedian(w):.1f} "
          f"({np.nanquantile(w, .25):.1f}-{np.nanquantile(w, .75):.1f}) "
          f"open {np.isnan(w).mean():.0%} profile "
          + " ".join(f"{o:+d}:{med[half + o]:.2f}" for o in OFFSETS_SHOWN))


def main():
  parser = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument("events", nargs="+")
  args = parser.parse_args()
  gp = linecheck.import_fullscan("grain_depth_profile")
  for ev in args.events:
    xyz, prof = grain_profiles(gp, ev)
    lab = labels(ev, xyz)
    compact = np.maximum(prof[:, 0], prof[:, -1]) < COMPACT_LEVEL
    print(f"{ev}: all {summary(gp, prof)}")
    print(f"  compact: {summary(gp, prof[compact])}")
    for name in sorted(set(lab) - {"all"}):
      m = lab == name
      print(f"  {name}: {summary(gp, prof[m])}")
      if (m & compact).sum():
        print(f"  {name}, compact: {summary(gp, prof[m & compact])}",
              flush=True)


if __name__ == "__main__":
  main()
