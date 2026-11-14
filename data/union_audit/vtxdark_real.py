"""Darkness near real full-scan vertices (judged by eye, v6 crops) vs
embedded v25 captures, same statistic: local darkness (B' darkness of
module/rankstudy.py) maximised within RADIUS_PX in XY and over
+-Z_SLACK slices around the vertex slice; percentile against the same
statistic at N_CTRL random positions of the view (same slices)."""
import csv
import json
import sys
from pathlib import Path

import numpy as np
from scipy.ndimage import maximum_filter, uniform_filter

from module import imagestats as st
from module.rankstudy import BG_SIGMA_PX, DARK_HALF_SLICES, DARK_RADIUS_PX

VIEWS = Path.home() / ("microscope-image/tohoku/fullscan/E07/MOD108/PL12/"
                       "tohoku-v1/AREA00/IMAGE00_AREA00")
REAL = Path.home() / "work/e07-fullscan/results/real_vertices_good_v6.csv"
EMBED = Path("data/images/v25_embed_nagara_dz4")
CNN_VIEWS = ("V00000011", "V00000991", "V00001012", "V00001339")
RADIUS_PX = 20      # ~10 um at 0.507, the position uncertainty
Z_SLACK = 3
N_CTRL = 200
EDGE_PX = 100
SEED = 1


def stat(json_path, x, y, k, rng):
  stack = np.array(st.read_stack(json_path))
  n, size = stack.shape[0], stack.shape[1]
  ctrl = rng.integers(EDGE_PX, size - EDGE_PX, (N_CTRL, 2))
  best, best_c = -np.inf, np.full(N_CTRL, -np.inf)
  for kk in range(max(0, k - Z_SLACK), min(n, k + Z_SLACK + 1)):
    lo, hi = max(0, kk - DARK_HALF_SLICES), min(n, kk + DARK_HALF_SLICES + 1)
    mp = stack[lo:hi].min(axis=0).astype(float)
    det = st.gaussian_blur(mp, BG_SIGMA_PX) - mp
    s = maximum_filter(uniform_filter(det, 2 * DARK_RADIUS_PX + 1),
                       2 * RADIUS_PX + 1)
    best = max(best, s[int(round(y)), int(round(x))])
    best_c = np.maximum(best_c, s[ctrl[:, 1], ctrl[:, 0]])
  return best, (best_c < best).mean()


def main():
  rng = np.random.default_rng(SEED)
  rows = []
  with open(REAL, newline="") as f:
    for r in csv.DictReader(f):
      if r["view"].startswith(CNN_VIEWS):
        continue
      path = VIEWS / f"{r['view']}.json"
      z = np.array([i["z"] for i in json.loads(path.read_text())["Images"]])
      k = int(np.argmin(np.abs(z - float(r["z_best_mm"]))))
      rows.append(("real",) + stat(path, float(r["x_px"]),
                                   float(r["y_px"]), k, rng))
  for p in sorted(EMBED.glob("ev*/truth.json")):
    t = json.loads(p.read_text())
    v = [x for x in t["vertices"] if x["type"] == "capture"][0]
    rows.append(("embedded",) + stat(p.parent / "image.json", v["col"],
                                     v["row"], int(round(v["slice"])), rng))
  for kind in ("real", "embedded"):
    d = np.array([r[1] for r in rows if r[0] == kind])
    p = np.array([r[2] for r in rows if r[0] == kind])
    q = lambda a: " ".join(f"{v:.3g}" for v in np.percentile(a, [25, 50, 75]))
    print(f"{kind}: n={len(d)} darkness [25 50 75 %] {q(d)}; percentile "
          f"{q(100 * p)}; >= 99th pct {np.mean(p >= .99):.0%}")


if __name__ == "__main__":
  main()
