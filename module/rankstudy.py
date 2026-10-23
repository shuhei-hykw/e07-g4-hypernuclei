"""Ranking study for point candidates (e.g. B' branch points).

  ./run.py rankstudy data/images/<dataset> <bprime.csv> [--top 100,300,770]

Each candidate is scored by
  multiplicity   the n_tracks_max column (branch multiplicity)
  density        the number of candidates within --radius um in 3-D
                 (developed emulsion), i.e. how strongly branch points
                 cluster: several tracks meeting at a vertex produce
                 several nearby branch points
  mult_density   the summed multiplicity of the candidates within the
                 radius
  darkness       image darkness around the candidate: background-
                 subtracted minimum projection over +-DARK_HALF_SLICES
                 slices, averaged within DARK_RADIUS_PX (heavy, dark
                 tracks start at hypernuclear vertices, while most
                 background branch points are crossings of faint MIPs)
and only the top K per event are kept; the one-to-one efficiency at 6
and 12 um and its chance level (shifted truths) are reported for each K,
so that rankings can be compared at the candidate rate of method A
(~770 per view).
"""

import argparse

from pathlib import Path

import numpy as np

from . import efficiency as eff
from . import imagestats as st

DEFAULT_TOP = "100,300,770"
DEFAULT_RADIUS_UM = 5.
TOLS = [6., 12.]
SCORES = ("multiplicity", "density", "mult_density", "darkness")
DARK_HALF_SLICES = 3
DARK_RADIUS_PX = 6
BG_SIGMA_PX = 8.
SHIFTS = eff.random_shifts(eff.DEFAULT_N_SHIFTS, eff.DEFAULT_SHIFT_MIN_PX,
                           eff.DEFAULT_SHIFT_MAX_PX, np.random.default_rng(0))


def coords_um(truth, cands):
  pix, dz = truth["pixel_um"], truth["physical_dz_um"]
  z0, sdz = truth.get("stage_z0_mm", 0.), truth["stage_dz_mm"]
  return np.array([[c["col"] * pix, c["row"] * pix,
                    (c["z"][0] - z0) / sdz * dz] for c in cands])


def neighbour_sums(xyz, weights, radius):
  """Sum of weights of the points within radius (including itself)."""
  cell = np.floor(xyz / radius).astype(np.int64)
  offset = cell.min(axis=0)
  cell -= offset
  dims = cell.max(axis=0) + 3
  key = ((cell[:, 0] + 1) * dims[1] + cell[:, 1] + 1) * dims[2] \
    + cell[:, 2] + 1
  order = np.argsort(key)
  skey = key[order]
  out = np.zeros(len(xyz))
  for dx in (-1, 0, 1):
    for dy in (-1, 0, 1):
      for dz in (-1, 0, 1):
        nk = key + (dx * dims[1] + dy) * dims[2] + dz
        lo = np.searchsorted(skey, nk, side="left")
        hi = np.searchsorted(skey, nk, side="right")
        for i in np.flatnonzero(hi > lo):
          j = order[lo[i]:hi[i]]
          d = np.linalg.norm(xyz[j] - xyz[i], axis=1)
          out[i] += weights[j][d <= radius].sum()
  return out


def darkness_scores(event_dir, truth, cands):
  """Mean background-subtracted darkness of the local minimum projection
  around every candidate."""
  path = Path(event_dir) / "image.json"
  stack = np.array(st.read_stack(path)).astype(float)
  n = stack.shape[0]
  z0, sdz = truth.get("stage_z0_mm", 0.), truth["stage_dz_mm"]
  r = DARK_RADIUS_PX
  out = np.zeros(len(cands))
  # Minimum projections cached per central slice
  cache = {}
  for i, c in enumerate(cands):
    k = int(np.clip(round((c["z"][0] - z0) / sdz), 0, n - 1))
    if k not in cache:
      lo, hi = max(0, k - DARK_HALF_SLICES), min(n, k + DARK_HALF_SLICES + 1)
      mp = stack[lo:hi].min(axis=0)
      cache[k] = st.gaussian_blur(mp, BG_SIGMA_PX) - mp
    det = cache[k]
    y, x = int(round(c["row"])), int(round(c["col"]))
    win = det[max(0, y - r):y + r + 1, max(0, x - r):x + r + 1]
    out[i] = win.mean() if win.size else 0.
  return out


def main():
  parser = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument("dataset")
  parser.add_argument("candidates")
  parser.add_argument("--top", default=DEFAULT_TOP)
  parser.add_argument("--radius", type=float, default=DEFAULT_RADIUS_UM)
  args = parser.parse_args()

  events = eff.load_truth(args.dataset)
  cands = eff.load_candidates(args.candidates)
  events = {k: v for k, v in events.items() if k in cands}
  scored = {s: {} for s in SCORES}
  for name, truth in events.items():
    cs = cands[name]
    xyz = coords_um(truth, cs)
    mult = np.array([c["score"] for c in cs])
    dens = neighbour_sums(xyz, np.ones(len(cs)), args.radius)
    mdens = neighbour_sums(xyz, mult, args.radius)
    dark = darkness_scores(Path(args.dataset) / name, truth, cs)
    for s, val in (("multiplicity", mult), ("density", dens),
                   ("mult_density", mdens), ("darkness", dark)):
      scored[s][name] = [dict(c, score=float(v)) for c, v in zip(cs, val)]
  print(f"{len(events)} events, radius {args.radius:g} um; cells: "
        "efficiency (chance) at <6 / <12 um, one-to-one")
  print(f"{'score':<14} {'top K':>6} {'capture':>22} {'LL decay':>22}")
  for s in SCORES:
    for k in [int(t) for t in args.top.split(",")]:
      top = {n: sorted(v, key=lambda c: -c["score"])[:k]
             for n, v in scored[s].items()}
      rows, _ = eff.evaluate(events, top, TOLS, SHIFTS, "developed")
      cells = []
      for vtype in ("capture", "HypernucleusDecay"):
        sel = [r for r in eff._select(rows, vtype, False)
               if vtype == "capture" or r["parent"].startswith("LL")]
        ctrl = [r for r in eff._select(rows, vtype, True)
                if vtype == "capture" or r["parent"].startswith("LL")]
        cells.append(" / ".join(
          f"{np.mean([r[f'found_edge_{t:g}'] for r in sel]):.0%}"
          f"({np.mean([r[f'found_edge_{t:g}'] for r in ctrl]):.0%})"
          for t in TOLS))
      print(f"{s:<14} {k:>6} {cells[0]:>22} {cells[1]:>22}", flush=True)


if __name__ == "__main__":
  main()
