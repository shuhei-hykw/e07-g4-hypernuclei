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
the best rank of a matching point (ties averaged) divided by the number
of points (smaller is better), compared with the same quantity for the
truth shifted at random inside the view (control). The shifted
positions are common to all scores; a control without any matching
point counts as rank fraction 1. Several (XY, z) tolerances are
evaluated in one run (--tols "200:30,20:30", px:um).

Caveat: the image model behind the darkness score was checked against
the heavy tracks of IBUKI, so these events are not a fully independent
test of the darkness ranking.
"""

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from . import imagestats as st
from .rankstudy import local_darkness, rank_fraction

FULLSCAN = Path.home() / "work/e07-fullscan"
DEFAULT_POINTS = FULLSCAN / "results/graphdet_specials_e07srv"
DEFAULT_GT = FULLSCAN / "tests/specials_gt.json"
DEFAULT_IMAGES = Path.home() / "work/specials_x20"
DEFAULT_TOLS = "200:30,20:30"
DEFAULT_N_SHIFTS = 20
SHIFT_MIN_PX = 400.
SHIFT_MAX_PX = 1000.
UM_PER_MM = 1000.
PIXEL_UM = 0.289
SCORES = ("multiplicity", "darkness", "combined")
CTRL_PCT = (16, 50, 84)
TOP_K = 770   # candidates per view of method A


def load_points(path):
  with open(path, newline="") as f:
    rows = list(csv.DictReader(f))
  return (np.array([[float(r["vx_px"]), float(r["vy_px"])] for r in rows]),
          np.array([float(r["slice"]) for r in rows]),
          np.array([float(r["z_mm"]) for r in rows]),
          np.array([float(r["mult"]) for r in rows]))


def best_match(frac, xy, z_mm, truth_xy, truth_z, tol_xy, tol_z):
  """Best rank fraction among matching points (1 if none), its index and
  the number of matching points."""
  match = (np.hypot(*(xy - truth_xy).T) <= tol_xy) \
    & (np.abs(z_mm - truth_z) * UM_PER_MM <= tol_z)
  if not match.any():
    return 1., None, 0
  idx = np.flatnonzero(match)
  j = idx[np.argmin(frac[idx])]
  return frac[j], j, len(idx)


def evaluate(fracs, xy, z_mm, truth_xy, truth_z, shifts, tol_xy, tol_z):
  rows = {}
  for name, frac in fracs.items():
    f, j, n = best_match(frac, xy, z_mm, truth_xy, truth_z, tol_xy, tol_z)
    ctrl = np.array([best_match(frac, xy, z_mm, sxy, truth_z, tol_xy,
                                tol_z)[0] for sxy in shifts])
    rows[name] = dict(
      frac=f, rank=f * len(frac) + 1 if j is not None else None,
      n_match=n,
      dxy_um=np.hypot(*(xy[j] - truth_xy)) * PIXEL_UM if j is not None
      else None,
      dz_um=(z_mm[j] - truth_z) * UM_PER_MM if j is not None else None,
      ctrl=np.percentile(ctrl, CTRL_PCT), ctrl_none=int((ctrl >= 1.).sum()),
      p=float((ctrl <= f).mean()))
  return rows


def fmt(r):
  if r["rank"] is None:
    true = "no match"
  else:
    true = (f"{r['rank']:.0f} ({r['frac']:.4f}), n={r['n_match']}, "
            f"dXY {r['dxy_um']:.0f} um, dz {r['dz_um']:+.0f} um")
  c = r["ctrl"]
  return (f"{true} | {c[1]:.4f} [{c[0]:.4f}, {c[2]:.4f}], "
          f"none {r['ctrl_none']} | {r['p']:.2f}")


def main():
  parser = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument("--points", default=str(DEFAULT_POINTS))
  parser.add_argument("--gt", default=str(DEFAULT_GT))
  parser.add_argument("--images", default=str(DEFAULT_IMAGES))
  parser.add_argument("--events", help="comma separated (default: all)")
  parser.add_argument("--tols", default=DEFAULT_TOLS)
  parser.add_argument("--n-shifts", type=int, default=DEFAULT_N_SHIFTS)
  parser.add_argument("-o", "--output")
  args = parser.parse_args()

  tols = [tuple(float(v) for v in t.split(":")) for t in args.tols.split(",")]
  gt = json.loads(Path(args.gt).read_text())["events"]
  names = args.events.split(",") if args.events else sorted(gt)
  rng = np.random.default_rng(0)
  results = {t: [] for t in tols}
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
    dark = local_darkness(stack, xy, slices)
    fracs = {"multiplicity": rank_fraction(mult),
             "darkness": rank_fraction(dark)}
    fracs["combined"] = rank_fraction(
      -(fracs["multiplicity"] + fracs["darkness"]))
    # Order independence of the tie-averaged ranks
    perm = rng.permutation(len(xy))
    comb_p = rank_fraction(-(rank_fraction(mult[perm])
                             + rank_fraction(dark[perm])))
    assert np.allclose(comb_p, fracs["combined"][perm])
    # Control positions common to all scores and tolerances
    r = rng.uniform(SHIFT_MIN_PX, SHIFT_MAX_PX, args.n_shifts)
    a = rng.uniform(0., 2. * np.pi, args.n_shifts)
    shifts = (truth_xy + np.c_[r * np.cos(a), r * np.sin(a)]) % size
    for t in tols:
      rows = evaluate(fracs, xy, z_mm, truth_xy, truth_z, shifts, *t)
      results[t].append((name, len(xy), rows))
      print(name, t, {k: (v["rank"], round(v["p"], 2))
                      for k, v in rows.items()}, flush=True)

  text = [
    "Best rank (ties averaged) of a B' branch point matching the clicked "
    "vertex, rank fraction in brackets (smaller is better), number of "
    "matching points and XY / z offset of the best one; control: median "
    f"[16 %, 84 %] of the rank fraction over {args.n_shifts} random "
    "shifts of the truth (common to all scores; no match counts as 1, "
    "'none' = such controls); p = fraction of controls at least as good "
    "as the truth. At the strict tolerance most controls have no "
    "matching point at all, so the comparison with the control median is "
    "weak there; the integer rank against the method-A rate "
    f"({TOP_K} per view) is the more useful figure.", "",
    "Caveats: the darkness image model was checked against the IBUKI "
    "heavy tracks, so this is not a fully independent test; a match "
    "only means a high-ranked point near the vertex, not that the tracks "
    "of the vertex are connected to it (not checked here).", ""]
  for t in tols:
    text += [f"### Tolerance {t[0]:g} px ({t[0] * PIXEL_UM:.0f} um) in XY, "
             f"{t[1]:g} um in z", "",
             "| event | points | score | true | control | p |",
             "|---|---|---|---|---|---|"]
    better = {s: 0 for s in SCORES}
    in_top = {s: 0 for s in SCORES}
    for name, n, rows in results[t]:
      for s in SCORES:
        text.append(f"| {name} | {n} | {s} | {fmt(rows[s])} |")
        better[s] += rows[s]["frac"] < rows[s]["ctrl"][1]
        in_top[s] += rows[s]["rank"] is not None \
          and rows[s]["rank"] <= TOP_K
    text += ["", "Better than the control median: " + ", ".join(
      f"{s} {better[s]}/{len(results[t])}" for s in SCORES)
      + f"; true rank within the top {TOP_K}: " + ", ".join(
      f"{s} {in_top[s]}/{len(results[t])}" for s in SCORES), ""]
  out = "\n".join(text) + "\n"
  if args.output:
    Path(args.output).write_text(out)
  print(out)


if __name__ == "__main__":
  main()
