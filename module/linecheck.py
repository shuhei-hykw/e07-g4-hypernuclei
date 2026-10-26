"""Long-line response of the e07-fullscan oriented line detector on real
and synthetic views (realism check not used in the calibration).

  ~/work/.venv/bin/python -m module.linecheck --real 20 \
    --synth data/images/<dataset> [-o output/linecheck.md]

ridge_response / line_mask of e07-fullscan scripts/eval_pointcloud.py
(oriented ridge filter integrated along the line, centre lines above
k robust sigmas and longer than a minimum length) are applied to the
central slice; the statistics are taken in the central crop: centre-line
pixels and lines per Mpx, line length distribution and the fraction of
DoG points (module/imagestats.py) within ON_LINE_PX of a centre line.
Needs the e07-fullscan environment (cv2, skimage), e.g. ~/work/.venv.
"""

import argparse
import importlib
import random
import subprocess
import sys
from pathlib import Path

import numpy as np

from . import imagestats as st
from .synth import calibrate as cal
from .synth import validate as val

FULLSCAN = Path.home() / "work/e07-fullscan"
CROP_PX = 1024
ON_LINE_PX = 3
LENGTH_PCT = (50, 90)


def fullscan_commit():
  return subprocess.run(["git", "-C", str(FULLSCAN), "rev-parse", "--short",
                         "HEAD"], capture_output=True, text=True).stdout.strip()


def import_fullscan(name="eval_pointcloud"):
  """A script of e07-fullscan/scripts (default eval_pointcloud.py). Both
  projects name their package `module`: ours is set aside while the
  full-scan one is imported, then restored (the imported functions keep
  their own)."""
  ours = {m: sys.modules.pop(m) for m in list(sys.modules)
          if m == "module" or m.startswith("module.")}
  path = sys.path[:]
  sys.path[:0] = [str(FULLSCAN / "scripts"), str(FULLSCAN)]
  try:
    mod = importlib.import_module(name)
  finally:
    for m in [m for m in sys.modules
              if m == "module" or m.startswith("module.")]:
      del sys.modules[m]
    sys.modules.update(ours)
    sys.path[:] = path
  return mod


def ridge_stats(img, ep):
  """Centre-line pixels per Mpx and fraction of DoG points within
  ON_LINE_PX of a centre line, in the central CROP_PX of a slice."""
  import cv2
  # Crop first: calibration patches are CROP_PX wide
  lo = max(0, (img.shape[0] - CROP_PX) // 2)
  img = np.ascontiguousarray(img[lo:lo + CROP_PX, lo:lo + CROP_PX])
  c = ep.line_mask(ep.ridge_response(img)) > 0
  d = st.gaussian_blur(img.astype(float), st.BG_SIGMA_PX) - img
  dog = st.dog_points(d)
  near = cv2.dilate(c.astype(np.uint8), np.ones((2 * ON_LINE_PX + 1,) * 2,
                                                np.uint8)) > 0
  return {"ridge_px_per_mpx": c.sum() / (c.size / 1e6),
          "ridge_dog_near": dog[near].sum() / max(dog.sum(), 1)}


def view_stats(json_path, ep):
  import cv2
  img = st.read_stack(json_path, [st.n_slices(json_path) // 2])[0]
  centre = ep.line_mask(ep.ridge_response(img)) > 0
  lo = (img.shape[0] - CROP_PX) // 2
  sl = np.s_[lo:lo + CROP_PX, lo:lo + CROP_PX]
  c = centre[sl]
  n, lab = cv2.connectedComponents(c.astype(np.uint8))
  length = np.bincount(lab.ravel())[1:]
  d = st.gaussian_blur(img.astype(float), st.BG_SIGMA_PX) - img
  dog = st.dog_points(d[sl])
  near = cv2.dilate(c.astype(np.uint8), np.ones((2 * ON_LINE_PX + 1,) * 2,
                                                np.uint8)) > 0
  mpx = c.size / 1e6
  return {"line_px_per_mpx": c.sum() / mpx, "lines_per_mpx": (n - 1) / mpx,
          **{f"length_p{p}": float(np.percentile(length, p))
             if len(length) else 0. for p in LENGTH_PCT},
          "dog_near_line": dog[near].sum() / max(dog.sum(), 1)}


def main():
  parser = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument("--real", type=int, default=20,
                      help="number of held-out real views")
  parser.add_argument("--synth", nargs="+", required=True)
  parser.add_argument("-o", "--output")
  args = parser.parse_args()
  ep = import_fullscan()
  real_views = val.holdout_views(args.real)
  real = [view_stats(v, ep) for v in real_views]
  lines = [f"Fullscan line detector at e07-fullscan {fullscan_commit()}; "
           f"{len(real)} held-out real views, central slice, central "
           f"{CROP_PX} px; z = (synthetic - real mean) / real spread.", "",
           "| statistic | real | " + " | ".join(
             Path(s).name for s in args.synth) + " |",
           "|---|---|" + "---|" * len(args.synth)]
  synth = [[view_stats(d / "image.json", ep)
            for d in sorted(Path(s).glob("ev*"))] for s in args.synth]
  for k in real[0]:
    r = np.array([x[k] for x in real])
    cells = []
    for ss in synth:
      v = np.array([x[k] for x in ss])
      z = (v.mean() - r.mean()) / (r.std() or 1.)
      cells.append(f"{v.mean():.4g} (z {z:+.1f})")
    lines.append(f"| {k} | {r.mean():.4g} +- {r.std():.2g} | "
                 + " | ".join(cells) + " |")
  text = "\n".join(lines) + "\n"
  if args.output:
    Path(args.output).write_text(text)
  print(text)


if __name__ == "__main__":
  main()
