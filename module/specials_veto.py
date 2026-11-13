"""Veto check of candidate arms on the real special events at 0.507 um/px.

  ./run.py specialsveto [--arm name:A.csv:B.csv ...] [-o FILE]

The special events resampled to the full-scan pixel by e07-fullscan
(results/specials_0507/<ev>: image.json, gt.json with the clicked vertex
vx / vy in px, z_mm the stage z of the clicked slice, tol_xy_um) are
scored per arm like the union of module/unionstudy.py: every method-A
candidate plus the B' points ranked by image darkness, the top K of the
view, with K = 770 scaled by the image area to the 2048 px full-scan
view (same candidate density per area). A candidate matches when it
lies within tol_xy_um in XY and its stage z range [z_min, z_max] within
TOL_Z_UM of the clicked z. Control: the same test at N_SHIFTS positions
shifted at random by SHIFT_MIN_PX..SHIFT_MAX_PX (common to all arms),
so the chance rate of a match at this loose tolerance is shown next to
the truth. The click is good to about 15-30 um only, so this is a veto
check: an arm should not lose a special event that production finds.
"""

import argparse
import json
from pathlib import Path

import numpy as np

from . import efficiency as eff
from . import imagestats as st
from .rankstudy import local_darkness, top_k

FULLSCAN = Path.home() / "work/e07-fullscan/results"
SPECIALS = FULLSCAN / "specials_0507"
SIM_EVAL = FULLSCAN / "sim_eval"
DEFAULT_ARMS = [
  "production:specials_0507_methodA_production.csv:"
  "specials_0507_bprime.csv",
  "zfix:specials_0507_methodA_production.csv:specials_0507_zfix_bprime.csv",
  "alen:specials_0507_alen_methodA_production.csv:"
  "specials_0507_bprime.csv",
]
FULL_VIEW_PX = 2048
FULL_TOP = 770
PIXEL_UM = 0.507
TOL_Z_UM = 30.
UM_PER_MM = 1000.
N_SHIFTS = 50
SHIFT_MIN_PX = 250.
SHIFT_MAX_PX = 450.
SEED = 0


def match(cs, xy, z_mm, tol_px):
  """Indices of candidates within tol_px in XY and TOL_Z_UM in z."""
  if not cs:
    return np.array([], dtype=int)
  p = np.array([(c["col"], c["row"]) for c in cs])
  z = np.array([c["z"] for c in cs])
  dz = np.maximum(z[:, 0] - z_mm, z_mm - z[:, 1]).clip(min=0.)
  ok = (np.hypot(*(p - xy).T) <= tol_px) & (dz * UM_PER_MM <= TOL_Z_UM)
  return np.flatnonzero(ok)


def score_arm(a, b, dark, k, xy, z_mm, tol_px):
  """(A found, B' best darkness rank or None, union found)."""
  fa = len(match(a, xy, z_mm, tol_px)) > 0
  mb = match(b, xy, z_mm, tol_px)
  rank = None
  if len(mb):
    order = np.empty(len(b), dtype=int)
    order[top_k(dark, len(b))] = np.arange(1, len(b) + 1)
    rank = int(order[mb].min())
  return fa, rank, fa or (rank is not None and rank <= k)


def main():
  parser = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument("--arm", action="append",
                      help="name:A.csv:B.csv (relative to sim_eval)")
  parser.add_argument("-o", "--output")
  args = parser.parse_args()

  arms = [s.split(":") for s in (args.arm or DEFAULT_ARMS)]
  cands = {n: (eff.load_candidates(SIM_EVAL / a),
               eff.load_candidates(SIM_EVAL / b)) for n, a, b in arms}
  rng = np.random.default_rng(SEED)
  rows, ctrl = [], {n: 0 for n, _, _ in arms}
  n_ctrl = 0
  for d in sorted(SPECIALS.iterdir()):
    gt_path, img = d / "gt.json", d / "image.json"
    if not gt_path.exists() or not img.exists():
      continue
    gt = json.loads(gt_path.read_text())
    meta = json.loads(img.read_text())
    zs = np.array([float(i["z"]) for i in meta["Images"]])
    size = int(meta["ImageType"]["Width"])
    k = round(FULL_TOP * (size / FULL_VIEW_PX) ** 2)
    tol_px = gt["tol_xy_um"] / PIXEL_UM
    xy = np.array([gt["vx"], gt["vy"]], dtype=float)
    r = rng.uniform(SHIFT_MIN_PX, SHIFT_MAX_PX, N_SHIFTS)
    phi = rng.uniform(0., 2. * np.pi, N_SHIFTS)
    shifts = xy + np.c_[r * np.cos(phi), r * np.sin(phi)]
    shifts = shifts[((shifts >= 0) & (shifts < size)).all(1)]
    n_ctrl += len(shifts)
    stack = np.array(st.read_stack(img))
    dz = np.median(np.diff(zs))
    out = {}
    for n, _, _ in arms:
      a, b = (c.get(d.name, []) for c in cands[n])
      dark = local_darkness(stack, [(c["col"], c["row"]) for c in b],
                            [(c["z"][0] - zs[0]) / dz for c in b]) \
        if b else np.array([])
      out[n] = (len(a), len(b)) + score_arm(a, b, dark, k, xy,
                                            gt["z_mm"], tol_px)
      ctrl[n] += sum(score_arm(a, b, dark, k, s, gt["z_mm"], tol_px)[2]
                     for s in shifts)
    rows.append((d.name, k, out))
    print(d.name, out, flush=True)

  names = [n for n, _, _ in arms]
  text = [
    "# Veto check on the special events at 0.507 um/px", "",
    "Per arm: candidates A / B', A match, best darkness rank of a "
    "matching B' point ('-' = none), union found (A, or B' within the "
    f"top K). Match: tolerance tol_xy_um of gt.json in XY and "
    f"{TOL_Z_UM:g} um in stage z; K = {FULL_TOP} scaled by the image area "
    f"to a {FULL_VIEW_PX} px view. Control: union found at random "
    f"positions {SHIFT_MIN_PX:g}-{SHIFT_MAX_PX:g} px away ({n_ctrl} in "
    "all, common to the arms).", "",
    "| event | K | " + " | ".join(names) + " |",
    "|---|---|" + "---|" * len(names)]
  for ev, k, out in rows:
    cells = []
    for n in names:
      na, nb, fa, rank, fu = out[n]
      cells.append(f"{na}/{nb} A {'y' if fa else 'n'} B' "
                   f"{rank if rank else '-'} {'FOUND' if fu else 'miss'}")
    text.append(f"| {ev} | {k} | " + " | ".join(cells) + " |")
  base = names[0]
  text += ["", "| arm | found | control | lost vs " + base + " |",
           "|---|---|---|---|"]
  for n in names:
    found = [out[n][4] for _, _, out in rows]
    lost = [ev for ev, _, out in rows if out[base][4] and not out[n][4]]
    text.append(f"| {n} | {sum(found)}/{len(rows)} | "
                f"{ctrl[n] / max(n_ctrl, 1):.0%} | "
                f"{', '.join(lost) or 'none'} |")
  out = "\n".join(text) + "\n"
  print(out)
  if args.output:
    Path(args.output).write_text(out)


if __name__ == "__main__":
  main()
