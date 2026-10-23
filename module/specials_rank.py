"""Darkness ranking of B' branch points on the real special events.

  ./run.py specialsrank [--points DIR] [--gt specials_gt.json]
    [--images DIR] [--events IBUKI,MINO]

For every special event the branch points of B' (<event>_points.csv:
vx_px, vy_px, slice, z_mm, mult) are scored by multiplicity and by the
image darkness of module/rankstudy.py (background-subtracted minimum
projection over +-3 slices, mean within 6 px), and by the sum of the two
rank fractions (combined). A branch point matches
the clicked true vertex if it lies within --tol-xy px in XY and its
stage z within --tol-z um of the clicked slice. The figure of merit is
the best rank of a matching point divided by the number of points
(smaller is better), compared with the same quantity for the truth
shifted at random inside the view (control, mean over --n-shifts).
"""

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from . import imagestats as st
from .rankstudy import BG_SIGMA_PX, DARK_HALF_SLICES, DARK_RADIUS_PX

FULLSCAN = Path.home() / "work/e07-fullscan"
DEFAULT_POINTS = FULLSCAN / "results/graphdet_specials_e07srv"
DEFAULT_GT = FULLSCAN / "tests/specials_gt.json"
DEFAULT_IMAGES = Path.home() / "work/specials_x20"
DEFAULT_TOL_XY_PX = 200.
DEFAULT_TOL_Z_UM = 30.
DEFAULT_N_SHIFTS = 20
SHIFT_MIN_PX = 400.
SHIFT_MAX_PX = 1000.
UM_PER_MM = 1000.


def load_points(path):
  with open(path, newline="") as f:
    rows = list(csv.DictReader(f))
  return (np.array([[float(r["vx_px"]), float(r["vy_px"])] for r in rows]),
          np.array([int(float(r["slice"])) for r in rows]),
          np.array([float(r["z_mm"]) for r in rows]),
          np.array([float(r["mult"]) for r in rows]))


def darkness(stack, xy, slices):
  n = stack.shape[0]
  r = DARK_RADIUS_PX
  cache, out = {}, np.zeros(len(xy))
  for i, ((x, y), k) in enumerate(zip(xy, slices)):
    k = int(np.clip(k, 0, n - 1))
    if k not in cache:
      lo, hi = max(0, k - DARK_HALF_SLICES), min(n, k + DARK_HALF_SLICES + 1)
      mp = stack[lo:hi].min(axis=0).astype(float)
      cache[k] = st.gaussian_blur(mp, BG_SIGMA_PX) - mp
    det = cache[k]
    yi, xi = int(round(y)), int(round(x))
    win = det[max(0, yi - r):yi + r + 1, max(0, xi - r):xi + r + 1]
    out[i] = win.mean() if win.size else 0.
  return out


def best_rank_fraction(score, xy, z_mm, truth_xy, truth_z, tol_xy, tol_z):
  order = np.argsort(-score)
  rank = np.empty(len(score), dtype=int)
  rank[order] = np.arange(1, len(score) + 1)
  match = (np.hypot(*(xy - truth_xy).T) <= tol_xy) \
    & (np.abs(z_mm - truth_z) * UM_PER_MM <= tol_z)
  if not match.any():
    return None
  return rank[match].min() / len(score)


def main():
  parser = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument("--points", default=str(DEFAULT_POINTS))
  parser.add_argument("--gt", default=str(DEFAULT_GT))
  parser.add_argument("--images", default=str(DEFAULT_IMAGES))
  parser.add_argument("--events", help="comma separated (default: all)")
  parser.add_argument("--tol-xy", type=float, default=DEFAULT_TOL_XY_PX)
  parser.add_argument("--tol-z", type=float, default=DEFAULT_TOL_Z_UM)
  parser.add_argument("--n-shifts", type=int, default=DEFAULT_N_SHIFTS)
  parser.add_argument("-o", "--output")
  args = parser.parse_args()

  gt = json.loads(Path(args.gt).read_text())["events"]
  names = args.events.split(",") if args.events else sorted(gt)
  rng = np.random.default_rng(0)
  lines = ["| event | points | multiplicity: true (control) | "
           "darkness: true (control) | combined: true (control) |",
           "|---|---|---|---|---|"]
  for name in names:
    pts = Path(args.points) / f"{name}_points.csv"
    img = Path(args.images) / name / "image.json"
    g = gt.get(name, {})
    if not pts.exists() or not img.exists() or g.get("z_slice") is None:
      continue
    xy, slices, z_mm, mult = load_points(pts)
    meta = json.loads(img.read_text())
    truth_z = float(meta["Images"][int(g["z_slice"])]["z"])
    truth_xy = np.array([g["vx"], g["vy"]], dtype=float)
    # 8-bit stack: 200-slice specials would need ~7 GB as floats
    stack = np.array(st.read_stack(img))
    size = stack.shape[1]
    dark = darkness(stack, xy, slices)
    # Combined: sum of the two rank fractions (smaller is better)
    rk = lambda v: np.argsort(np.argsort(-v)) / len(v)
    combined = -(rk(mult) + rk(dark))
    cells = []
    for score in (mult, dark, combined):
      true = best_rank_fraction(score, xy, z_mm, truth_xy, truth_z,
                                args.tol_xy, args.tol_z)
      ctrl = []
      for _ in range(args.n_shifts):
        r = rng.uniform(SHIFT_MIN_PX, SHIFT_MAX_PX)
        a = rng.uniform(0., 2. * np.pi)
        sxy = (truth_xy + r * np.array([np.cos(a), np.sin(a)])) % size
        f = best_rank_fraction(score, xy, z_mm, sxy, truth_z, args.tol_xy,
                               args.tol_z)
        if f is not None:
          ctrl.append(f)
      cells.append(f"{true:.4f}" if true is not None else "no match")
      cells[-1] += f" ({np.median(ctrl):.4f})" if ctrl else " (-)"
    lines.append(f"| {name} | {len(xy)} | {cells[0]} | {cells[1]} | "
                 f"{cells[2]} |")
    print(lines[-1], flush=True)
  text = ("Best rank of a branch point matching the clicked vertex divided "
          "by the number of points (smaller is better); control: median "
          f"over {args.n_shifts} random shifts of the truth.\n\n"
          + "\n".join(lines) + "\n")
  if args.output:
    Path(args.output).write_text(text)
  print(text)


if __name__ == "__main__":
  main()
