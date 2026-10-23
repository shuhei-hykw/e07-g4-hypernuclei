"""Vertex-finding efficiency of a full-scan method on synthetic images.

  ./run.py efficiency data/images/nagara_dz3.0 candidates.csv \
    [--tolerances 6,12,30] [-o eff.png]

`candidates.csv` holds the merged vertex candidates of one method, one
row per candidate, with the columns (extra columns are ignored)

  event         event directory name, e.g. ev0003
  vx_px, vy_px  position in pixels (column, row)
  z_min, z_max  stage z range of the contributing slices (mm), or
  z_mean        a single stage z (mm)
  n_tracks_max  (optional) ranking score, higher first

True vertices come from truth.json (in-view vertices only). Truth and
candidates are matched one-to-one, closest pairs first, in 3-D distance
in the developed emulsion: x, y from pixels and z from the stage z via
the slice index and the physical focal step per slice of the dataset.
A candidate z range counts as distance zero for truths inside it. The
full-scan criterion of tests/specials_gt.json (200 px in XY and the true
depth inside [z_min - 30 um, z_max + 30 um] of stage z) is also given.

Chance coincidences are estimated by shifting all true vertices by
--shift-px in x (wrapping inside the view) and matching again.
"""

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

DEFAULT_TOLERANCES = "6,12,30"
DEFAULT_SHIFT_PX = 600.
FULLSCAN_XY_PX = 200.
FULLSCAN_Z_UM = 30.
UM_PER_MM = 1000.
VERTEX_TYPES = ["capture", "HypernucleusDecay", "hBertiniCaptureAtRest",
                "Decay"]
LL_PREFIX = "LL"
N_LENGTH_BINS = 5


def load_truth(dataset):
  events = {}
  for d in sorted(Path(dataset).glob("ev*")):
    path = d / "truth.json"
    if path.exists():
      events[d.name] = json.loads(path.read_text())
  if not events:
    raise SystemExit(f"no truth.json under {dataset}")
  return events


def load_candidates(path):
  out = {}
  with open(path, newline="") as f:
    for row in csv.DictReader(f):
      c = {"col": float(row["vx_px"]), "row": float(row["vy_px"])}
      if row.get("z_min") not in (None, "") and row.get("z_max"):
        c["z"] = (float(row["z_min"]), float(row["z_max"]))
      else:
        z = float(row["z_mean"])
        c["z"] = (z, z)
      c["score"] = float(row.get("n_tracks_max") or 0.)
      out.setdefault(row["event"], []).append(c)
  return out


def stage_z(truth, slice_index):
  """Stage z (mm) of a fractional slice index."""
  return slice_index * truth["stage_dz_mm"]


def distances(truth, verts, cands):
  """3-D distance (um) and full-scan criterion for every pair."""
  pix, dz_um, sdz = truth["pixel_um"], truth["physical_dz_um"], \
    truth["stage_dz_mm"]
  n, m = len(verts), len(cands)
  dist = np.full((n, m), np.inf)
  fullscan = np.zeros((n, m), dtype=bool)
  for i, v in enumerate(verts):
    z_true = stage_z(truth, v["slice"])
    for j, c in enumerate(cands):
      dxy_px = np.hypot(c["col"] - v["col"], c["row"] - v["row"])
      lo, hi = sorted(c["z"])
      # Distance from the true depth to the candidate z range, in slices
      gap_mm = max(lo - z_true, 0., z_true - hi)
      gap_um = abs(gap_mm / sdz) * dz_um
      dist[i, j] = np.hypot(dxy_px * pix, gap_um)
      fullscan[i, j] = (dxy_px <= FULLSCAN_XY_PX
                        and gap_mm * UM_PER_MM <= FULLSCAN_Z_UM)
  return dist, fullscan


def match(dist, tolerance):
  """Greedy one-to-one matching of pairs closer than tolerance."""
  n, m = dist.shape
  pairs = np.argwhere(dist <= tolerance)
  order = np.argsort(dist[pairs[:, 0], pairs[:, 1]])
  used_t, used_c, result = set(), set(), {}
  for i, j in pairs[order]:
    if i in used_t or j in used_c:
      continue
    used_t.add(i)
    used_c.add(j)
    result[int(i)] = int(j)
  return result


def shifted(verts, shift_px, size_px):
  out = []
  for v in verts:
    w = dict(v)
    w["col"] = (v["col"] + shift_px) % size_px
    out.append(w)
  return out


def ll_track_length(truth):
  """Range (um) of the double-Lambda hypernucleus, if any."""
  for t in truth["tracks"]:
    if t["name"].startswith(LL_PREFIX) and t["parent"] == 0:
      return t["range_um"]
  return None


def evaluate(events, cands, tolerances, shift_px):
  rows = []
  for name, truth in events.items():
    verts = [v for v in truth["vertices"] if v["in_view"]]
    cs = sorted(cands.get(name, []), key=lambda c: -c["score"])
    size = truth.get("size_px", 2048)
    for label, vs in (("true", verts),
                      ("shifted", shifted(verts, shift_px, size))):
      if not vs:
        continue
      if cs:
        dist, fullscan = distances(truth, vs, cs)
      else:
        dist = np.full((len(vs), 0), np.inf)
        fullscan = np.zeros((len(vs), 0), dtype=bool)
      matched = {tol: match(dist, tol) for tol in tolerances}
      for i, v in enumerate(vs):
        row = {"event": name, "sample": label, "type": v["type"],
               "parent": v["parent"], "ll_range_um": ll_track_length(truth),
               "n_cand": len(cs),
               "nearest_um": float(dist[i].min()) if len(cs) else None,
               "fullscan": bool(fullscan[i].any()) if len(cs) else False}
        for tol in tolerances:
          j = matched[tol].get(i)
          row[f"found_{tol:g}"] = j is not None
          row[f"rank_{tol:g}"] = j + 1 if j is not None else None
        rows.append(row)
  return rows


def summarize(rows, tolerances):
  print(f"{'type':<24} {'sample':<8} {'n':>4} "
        + " ".join(f"{f'<{t:g}um':>8}" for t in tolerances)
        + f" {'fullscan':>9}")
  for vtype in VERTEX_TYPES + ["LL event"]:
    for sample in ("true", "shifted"):
      if vtype == "LL event":
        # Both the capture and the double-Lambda decay vertex found
        sel = _ll_events(rows, sample, tolerances)
      else:
        sel = [r for r in rows if r["type"] == vtype
               and r["sample"] == sample]
      if not sel:
        continue
      effs = [np.mean([r[f"found_{t:g}"] for r in sel]) for t in tolerances]
      fs = np.mean([r["fullscan"] for r in sel])
      print(f"{vtype:<24} {sample:<8} {len(sel):>4} "
            + " ".join(f"{e:8.1%}" for e in effs) + f" {fs:9.1%}")


def _ll_events(rows, sample, tolerances):
  out = []
  for name in sorted({r["event"] for r in rows}):
    rs = [r for r in rows if r["event"] == name and r["sample"] == sample]
    cap = [r for r in rs if r["type"] == "capture"]
    dec = [r for r in rs if r["type"] == "HypernucleusDecay"
           and r["parent"].startswith(LL_PREFIX)]
    if not cap or not dec:
      continue
    row = {"fullscan": cap[0]["fullscan"] and dec[0]["fullscan"]}
    for t in tolerances:
      row[f"found_{t:g}"] = cap[0][f"found_{t:g}"] and dec[0][f"found_{t:g}"]
    out.append(row)
  return out


def plot(rows, tolerances, output):
  import ROOT
  ROOT.gROOT.SetBatch(True)
  ROOT.gStyle.SetOptStat(0)
  canvas = ROOT.TCanvas("c", "", 1200, 500)
  canvas.Divide(2, 1)
  keep = []
  canvas.cd(1)
  frame = ROOT.gPad.DrawFrame(0., 0., max(tolerances) * 1.1, 1.05)
  frame.SetTitle("efficiency vs tolerance;3-D tolerance (#mum);efficiency")
  legend = ROOT.TLegend(0.55, 0.15, 0.9, 0.4)
  for k, vtype in enumerate(VERTEX_TYPES[:3]):
    for sample, style in (("true", 1), ("shifted", 2)):
      sel = [r for r in rows if r["type"] == vtype and r["sample"] == sample]
      if not sel:
        continue
      g = ROOT.TGraph(len(tolerances))
      for n, t in enumerate(tolerances):
        g.SetPoint(n, t, np.mean([r[f"found_{t:g}"] for r in sel]))
      g.SetLineColor(k + 2)
      g.SetMarkerColor(k + 2)
      g.SetLineStyle(style)
      g.SetMarkerStyle(20)
      g.Draw("LP")
      legend.AddEntry(g, f"{vtype} ({sample})", "lp")
      keep.append(g)
  legend.Draw()
  canvas.cd(2)
  t_mid = tolerances[len(tolerances) // 2]
  sel = [r for r in rows if r["type"] == "capture" and r["sample"] == "true"
         and r["ll_range_um"] is not None]
  lengths = np.array([r["ll_range_um"] for r in sel])
  edges = np.unique(np.quantile(lengths, np.linspace(0., 1., N_LENGTH_BINS
                                                     + 1))) if sel else []
  if len(edges) >= 2:
    found = np.array([r[f"found_{t_mid:g}"] for r in sel], dtype=float)
    h_all = ROOT.TH1D("all", "", len(edges) - 1, edges)
    h_ok = ROOT.TH1D("ok", "", len(edges) - 1, edges)
    for x, f in zip(lengths, found):
      h_all.Fill(x)
      if f:
        h_ok.Fill(x)
    eff = ROOT.TEfficiency(h_ok, h_all)
    eff.SetTitle(f"capture vertex, tolerance {t_mid:g} #mum;"
                 "LL hypernucleus range (#mum);efficiency")
    eff.Draw("AP")
    keep += [h_all, h_ok, eff]
  canvas.SaveAs(output)


def main():
  parser = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument("dataset", help="data/images/<name>")
  parser.add_argument("candidates", help="CSV of vertex candidates")
  parser.add_argument("--tolerances", default=DEFAULT_TOLERANCES,
                      help="3-D tolerances in um, comma separated")
  parser.add_argument("--shift-px", type=float, default=DEFAULT_SHIFT_PX)
  parser.add_argument("--details", help="write per-vertex rows to CSV")
  parser.add_argument("-o", "--output", help="efficiency plot (PNG)")
  args = parser.parse_args()

  tolerances = [float(t) for t in args.tolerances.split(",")]
  events = load_truth(args.dataset)
  cands = load_candidates(args.candidates)
  rows = evaluate(events, cands, tolerances, args.shift_px)
  summarize(rows, tolerances)
  if args.details:
    with open(args.details, "w", newline="") as f:
      w = csv.DictWriter(f, fieldnames=list(rows[0]))
      w.writeheader()
      w.writerows(rows)
  if args.output:
    plot(rows, tolerances, args.output)
  print(f"{len(events)} events, "
        f"{sum(len(c) for c in cands.values())} candidates", file=sys.stderr)


if __name__ == "__main__":
  main()
