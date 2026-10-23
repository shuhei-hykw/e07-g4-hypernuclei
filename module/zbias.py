"""Depth bias of candidate z_best relative to the true vertex depth.

  ./run.py zbias data/images/<dataset> candidates.csv [--sim data/sim/x.root]

For every true vertex of the chosen types, the candidate nearest in XY
(within --max-xy um) is taken and its offset z_best - z_true is given in
developed um and in slices. Positive offsets point to larger slice
indices, i.e. larger developed z (downstream along the beam); in the
tohoku-v1 convention the stage z decreases with the slice index, so a
positive offset is a smaller stage z.

The offsets are split by
  (a) the depth of the vertex inside its emulsion layer (thirds, from
      the upstream face; needs --sim for the layer table),
  (b) the z asymmetry of the charged daughters: (path length going to
      larger slices - to smaller) / total, using the polylines in view,
  (c) the n_tracks_max of the candidate,
  (d) the position of the vertex in the image stack (thirds of the
      slice range; separates a pull to the stack centre from one to the
      layer centre).
"""

import argparse
import json
from pathlib import Path

import numpy as np

from . import efficiency as eff

DEFAULT_MAX_XY_UM = 12.
TYPES = ("capture", "HypernucleusDecay")
LAYER_THIRDS = ("upstream third", "middle third", "downstream third")
NTRACK_BINS = ((0, 3), (4, 5), (6, 1000))
ASYM_CUT = 0.33


def daughter_asymmetry(truth, vertex_index):
  """(visible path to larger slices - to smaller) / total, in [-1, 1]."""
  plus = minus = 0.
  n = truth["n_slices"]
  for t in truth["tracks"]:
    if t.get("vertex_start") != vertex_index:
      continue
    p = np.asarray(t["polyline"])
    if len(p) < 2:
      continue
    inside = (p[:, 2] >= 0) & (p[:, 2] <= n - 1)
    seg = np.diff(p, axis=0)
    ok = inside[:-1] & inside[1:]
    length = np.linalg.norm(seg[ok, :2], axis=1)
    dz = seg[ok, 2]
    plus += length[dz > 0].sum()
    minus += length[dz < 0].sum()
  total = plus + minus
  return (plus - minus) / total if total > 0 else 0.


def layer_third(layers, z_um):
  for emul, lo, hi in layers:
    if emul and lo <= z_um <= hi:
      return min(int(3 * (z_um - lo) / (hi - lo)), 2)
  return None


def load_layers(sim):
  import ROOT
  f = ROOT.TFile.Open(str(sim))
  out = [(bool(l.emulsion), l.zmin, l.zmax) for l in f.Get("layers")]
  f.Close()
  return out


def summary(vals):
  v = np.asarray(vals)
  if len(v) == 0:
    return "n=0"
  p16, p50, p84 = np.percentile(v, [16, 50, 84])
  return (f"n={len(v):4d} mean {v.mean():6.1f} median {p50:6.1f} "
          f"p16/p84 {p16:6.1f}/{p84:6.1f}")


def main():
  parser = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument("dataset")
  parser.add_argument("candidates")
  parser.add_argument("--sim", help="emulsim ROOT file (layer table)")
  parser.add_argument("--max-xy", type=float, default=DEFAULT_MAX_XY_UM)
  args = parser.parse_args()

  events = eff.load_truth(args.dataset)
  cands = eff.load_candidates(args.candidates)
  layers = load_layers(args.sim) if args.sim else None
  rows = []
  for name, truth in events.items():
    cs = [c for c in cands.get(name, []) if "z_best" in c]
    if not cs:
      continue
    cc = np.array([[c["col"], c["row"]] for c in cs])
    z0, sdz = truth["stage_z0_mm"], truth["stage_dz_mm"]
    dz_um = truth["physical_dz_um"]
    for v in truth["vertices"]:
      if v["type"] not in TYPES or not v["in_view"]:
        continue
      dxy = np.hypot(cc[:, 0] - v["col"], cc[:, 1] - v["row"]) \
        * truth["pixel_um"]
      j = int(np.argmin(dxy))
      if dxy[j] > args.max_xy:
        continue
      slice_best = (cs[j]["z_best"] - z0) / sdz
      off_sl = slice_best - v["slice"]
      rows.append({
        "type": v["type"], "off_slices": off_sl, "off_um": off_sl * dz_um,
        "asym": daughter_asymmetry(truth, v["index"]),
        "stack_third": min(int(3 * v["slice"] / truth["n_slices"]), 2),
        "third": layer_third(layers, v["um"][2]) if layers else None,
        "ntracks": cs[j]["score"]})
  print(f"offset z_best - z_true of the XY-nearest candidate within "
        f"{args.max_xy:g} um; + = larger slice index (larger developed z, "
        "smaller tohoku stage z)")
  for vtype in TYPES:
    sel = [r for r in rows if r["type"] == vtype]
    print(f"\n{vtype}")
    print(f"  all          um: {summary([r['off_um'] for r in sel])}")
    print(f"               slices: "
          f"{summary([r['off_slices'] for r in sel])}")
    for lo, hi, label in ((-1.01, -ASYM_CUT, "daughters mostly -"),
                          (-ASYM_CUT, ASYM_CUT, "daughters balanced"),
                          (ASYM_CUT, 1.01, "daughters mostly +")):
      s = [r["off_um"] for r in sel if lo <= r["asym"] < hi]
      print(f"  {label:20s} um: {summary(s)}")
    for k, label in enumerate(("stack first third", "stack middle third",
                               "stack last third")):
      s = [r["off_um"] for r in sel if r["stack_third"] == k]
      print(f"  {label:20s} um: {summary(s)}")
    if layers:
      for k, label in enumerate(LAYER_THIRDS):
        s = [r["off_um"] for r in sel if r["third"] == k]
        print(f"  {label:16s} um: {summary(s)}")
    for lo, hi in NTRACK_BINS:
      s = [r["off_um"] for r in sel if lo <= r["ntracks"] <= hi]
      print(f"  n_tracks {lo}-{hi if hi < 1000 else 'max'}  um: "
            f"{summary(s)}")


if __name__ == "__main__":
  main()
