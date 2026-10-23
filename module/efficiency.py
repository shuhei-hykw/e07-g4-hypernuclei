"""Vertex-finding efficiency of a full-scan method on synthetic images.

  ./run.py efficiency data/images/nagara_dz3.0 candidates.csv \
    [--tolerances 6,12,30] [--details rows.csv] [-o eff.png]

`candidates.csv` holds the merged vertex candidates of one method, one
row per candidate, with the columns (extra columns are ignored)

  event         event directory name, e.g. ev0003 (or view_id, the path
                of the event's image.json as written by the full scan)
  vx_px, vy_px  position in pixels (column, row)
  z_min, z_max  stage z range of the contributing slices (mm), or
  z_mean        a single stage z (mm; point candidates, e.g. method B')
  n_tracks_max  (optional) ranking score, higher first

Truth: the vertices of truth.json with in_view true.

Distance: 3-D in the developed emulsion (image space), x and y from
pixels (x pixel_um), z from the stage z via the slice index (x
physical_dz_um); --space undeveloped multiplies z by Sh (Geant4 space).
For a candidate with a z range two depth conventions are reported:
  edge       distance to the nearer end of the range (primary; a star
             extends to one side of its vertex, so the vertex lies at an
             end of the range)
  inclusive  zero inside the range, distance to the nearer end outside
  best       the candidate's z_best column as a point depth (stage z of
             the slice where most tracks cross), when present
Point candidates use the plain distance in all.

Matching: per event and tolerance, pairs within the tolerance are taken
closest first (greedy one-to-one). The optimal assignment (most matches,
then smallest summed distance) is computed too and the difference is
reported. Candidates within the tolerance of a truth but left unmatched
are counted as duplicates, separately from fakes.

Control: the truths are shifted --n-shifts times by a random vector of
length --shift-min..--shift-max px (wrapping inside the view) and
matched again; the mean efficiency estimates chance coincidences.

The full-scan criterion of tests/specials_gt.json (any candidate within
200 px in XY with the true stage depth within 30 um of its z range) is
also given as a loose recovery measure.
"""

import argparse
import csv
import json
from pathlib import Path

import numpy as np

DEFAULT_TOLERANCES = "6,12,30"
DEFAULT_N_SHIFTS = 10
DEFAULT_SHIFT_MIN_PX = 400.
DEFAULT_SHIFT_MAX_PX = 1000.
FULLSCAN_XY_PX = 200.
FULLSCAN_Z_UM = 30.
UM_PER_MM = 1000.
VERTEX_TYPES = ["capture", "HypernucleusDecay", "hBertiniCaptureAtRest",
                "Decay"]
LL_PREFIX = "LL"
N_LENGTH_BINS = 5
METRICS = ("edge", "inclusive", "best")
DEFAULT_SIZE_PX = 2048


def load_truth(dataset):
  events = {}
  for d in sorted(Path(dataset).glob("ev*")):
    path = d / "truth.json"
    if path.exists():
      events[d.name] = json.loads(path.read_text())
  if not events:
    raise SystemExit(f"no truth.json under {dataset}")
  return events


def load_candidates(path, min_score=None):
  out = {}
  with open(path, newline="") as f:
    for row in csv.DictReader(f):
      c = {"col": float(row["vx_px"]), "row": float(row["vy_px"])}
      if row.get("z_min") not in (None, "") and row.get("z_max"):
        c["z"] = tuple(sorted((float(row["z_min"]), float(row["z_max"]))))
      else:
        z = float(row["z_mean"])
        c["z"] = (z, z)
      c["score"] = float(row.get("n_tracks_max") or 0.)
      if row.get("z_best") not in (None, ""):
        c["z_best"] = float(row["z_best"])
      if min_score is not None and c["score"] < min_score:
        continue
      event = row.get("event") or Path(row["view_id"]).parent.name
      out.setdefault(event, []).append(c)
  return out


def stage_z(truth, slice_index):
  """Stage z (mm) of a fractional slice index."""
  return truth.get("stage_z0_mm", 0.) + slice_index * truth["stage_dz_mm"]


def distances(truth, verts, cands, space):
  """Distance matrices (um) per metric and the full-scan criterion."""
  pix, dz_um, sdz = (truth["pixel_um"], truth["physical_dz_um"],
                     abs(truth["stage_dz_mm"]))
  if space == "undeveloped":
    dz_um *= truth["shrink"]
  vc = np.array([[v["col"], v["row"]] for v in verts])
  vz = np.array([stage_z(truth, v["slice"]) for v in verts])
  cc = np.array([[c["col"], c["row"]] for c in cands])
  lo = np.array([c["z"][0] for c in cands])
  hi = np.array([c["z"][1] for c in cands])
  dxy_px = np.hypot(vc[:, None, 0] - cc[None, :, 0],
                    vc[:, None, 1] - cc[None, :, 1])
  edge_mm = np.minimum(np.abs(vz[:, None] - lo[None, :]),
                       np.abs(vz[:, None] - hi[None, :]))
  inside = (vz[:, None] >= lo[None, :]) & (vz[:, None] <= hi[None, :])
  incl_mm = np.where(inside, 0., edge_mm)
  gaps = [("edge", edge_mm), ("inclusive", incl_mm)]
  if all("z_best" in c for c in cands):
    zb = np.array([c["z_best"] for c in cands])
    gaps.append(("best", np.abs(vz[:, None] - zb[None, :])))
  out = {}
  for metric, gap_mm in gaps:
    out[metric] = np.hypot(dxy_px * pix, gap_mm / sdz * dz_um)
  out["fullscan"] = (dxy_px <= FULLSCAN_XY_PX) \
    & (incl_mm * UM_PER_MM <= FULLSCAN_Z_UM)
  return out


def match_greedy(dist, tolerance):
  pairs = np.argwhere(dist <= tolerance)
  order = np.argsort(dist[pairs[:, 0], pairs[:, 1]])
  used_t, used_c, result = set(), set(), {}
  for i, j in pairs[order]:
    if i in used_t or j in used_c:
      continue
    used_t.add(int(i))
    used_c.add(int(j))
    result[int(i)] = int(j)
  return result


def match_optimal(dist, tolerance):
  """Most matches, then smallest summed distance (exhaustive search over
  the few truths of an event, candidates restricted to the tolerance)."""
  options = [list(np.flatnonzero(dist[i] <= tolerance))
             for i in range(dist.shape[0])]
  best = {"key": (-1, 0.), "assigned": {}}

  def search(i, used, assigned, total):
    if i == len(options):
      key = (len(assigned), -total)
      if key > best["key"]:
        best["key"] = key
        best["assigned"] = dict(assigned)
      return
    if len(assigned) + (len(options) - i) < best["key"][0]:
      return
    for j in options[i]:
      if j not in used:
        used.add(j)
        assigned[i] = int(j)
        search(i + 1, used, assigned, total + dist[i, j])
        del assigned[i]
        used.discard(j)
    search(i + 1, used, assigned, total)

  search(0, set(), {}, 0.)
  return best["assigned"]


def shifted(verts, shift, size_px):
  out = []
  for v in verts:
    w = dict(v)
    w["col"] = (v["col"] + shift[0]) % size_px
    w["row"] = (v["row"] + shift[1]) % size_px
    out.append(w)
  return out


def ll_track_length(truth):
  """Range (um) of the double-Lambda hypernucleus, if any."""
  for t in truth["tracks"]:
    if t["name"].startswith(LL_PREFIX) and t["parent"] == 0:
      return t["range_um"]
  return None


def random_shifts(n, lo, hi, rng):
  r = rng.uniform(lo, hi, n)
  a = rng.uniform(0., 2. * np.pi, n)
  return np.column_stack([r * np.cos(a), r * np.sin(a)])


def nearest_components(events, cands, shifts, space):
  """For every truth (true and shifted) the XY distance (um) to the
  nearest candidate in XY and that candidate's edge z gap (um)."""
  out = {"true": {}, "control": {}}
  for name, truth in events.items():
    verts = [v for v in truth["vertices"] if v["in_view"]]
    cs = cands.get(name, [])
    if not verts or not cs:
      continue
    pix, dz_um = truth["pixel_um"], truth["physical_dz_um"]
    if space == "undeveloped":
      dz_um *= truth["shrink"]
    sdz = abs(truth["stage_dz_mm"])
    size = truth.get("size_px", DEFAULT_SIZE_PX)
    cc = np.array([[c["col"], c["row"]] for c in cs])
    lo = np.array([c["z"][0] for c in cs])
    hi = np.array([c["z"][1] for c in cs])
    samples = [("true", verts)] + [("control", shifted(verts, s, size))
                                   for s in shifts]
    for label, vs in samples:
      for v in vs:
        dxy = np.hypot(cc[:, 0] - v["col"], cc[:, 1] - v["row"]) * pix
        j = int(np.argmin(dxy))
        z = stage_z(truth, v["slice"])
        gap = min(abs(z - lo[j]), abs(z - hi[j])) / sdz * dz_um
        zb = cs[j].get("z_best")
        gap_best = abs(z - zb) / sdz * dz_um if zb is not None else np.nan
        out[label].setdefault(v["type"], []).append((dxy[j], gap, gap_best))
  return out


def print_nearest(near):
  print("nearest candidate in XY: median XY distance / median z gap to "
        "the range edge / to z_best (um), fraction within 6 / 12 um in XY")
  for vtype in VERTEX_TYPES:
    cells = []
    for label in ("true", "control"):
      vals = np.array(near[label].get(vtype, []))
      if len(vals) == 0:
        cells.append(f"{label}: -")
        continue
      cells.append(f"{label}: {np.median(vals[:, 0]):6.1f} / "
                   f"{np.median(vals[:, 1]):6.1f} / "
                   f"{np.nanmedian(vals[:, 2]) if np.isfinite(vals[:, 2]).any() else float('nan'):6.1f}  "
                   f"({np.mean(vals[:, 0] < 6):.0%}, "
                   f"{np.mean(vals[:, 0] < 12):.0%})  n={len(vals)}")
    print(f"  {vtype:<22} " + "   ".join(cells))


def evaluate(events, cands, tolerances, shifts, space):
  rows = []
  extra = {"duplicates": {t: 0 for t in tolerances},
           "optimal_minus_greedy": {t: 0 for t in tolerances},
           "n_candidates": 0}
  for name, truth in events.items():
    verts = [v for v in truth["vertices"] if v["in_view"]]
    cs = sorted(cands.get(name, []), key=lambda c: -c["score"])
    extra["n_candidates"] += len(cs)
    size = truth.get("size_px", DEFAULT_SIZE_PX)
    samples = [("true", verts)] + [
      (f"shift{k}", shifted(verts, s, size)) for k, s in enumerate(shifts)]
    for label, vs in samples:
      if not vs:
        continue
      if cs:
        d = distances(truth, vs, cs, space)
      else:
        empty = np.full((len(vs), 0), np.inf)
        d = {"edge": empty, "inclusive": empty, "best": empty,
             "fullscan": np.zeros((len(vs), 0), dtype=bool)}
      metrics = [m for m in METRICS if m in d]
      matched = {(m, t): match_greedy(d[m], t)
                 for m in metrics for t in tolerances}
      if label == "true":
        for t in tolerances:
          g = matched[("edge", t)]
          o = match_optimal(d["edge"], t)
          extra["optimal_minus_greedy"][t] += len(o) - len(g)
          near = np.flatnonzero((d["edge"] <= t).any(axis=0))
          extra["duplicates"][t] += len(set(near) - set(g.values()))
      for i, v in enumerate(vs):
        row = {"event": name, "sample": label, "type": v["type"],
               "parent": v["parent"], "ll_range_um": ll_track_length(truth),
               "n_cand": len(cs),
               "nearest_edge_um": float(d["edge"][i].min()) if cs else None,
               "fullscan": bool(d["fullscan"][i].any()) if cs else False}
        for m in metrics:
          for t in tolerances:
            j = matched[(m, t)].get(i)
            row[f"found_{m}_{t:g}"] = j is not None
            if m == "edge":
              row[f"rank_{t:g}"] = j + 1 if j is not None else None
        rows.append(row)
  return rows, extra


def _select(rows, vtype, control):
  if control:
    return [r for r in rows if r["type"] == vtype
            and r["sample"].startswith("shift")]
  return [r for r in rows if r["type"] == vtype and r["sample"] == "true"]


def _ll_events(rows, control, tolerances):
  out = []
  keys = sorted({(r["event"], r["sample"]) for r in rows
                 if r["sample"].startswith("shift") == control})
  for event, sample in keys:
    rs = [r for r in rows if r["event"] == event and r["sample"] == sample]
    cap = [r for r in rs if r["type"] == "capture"]
    dec = [r for r in rs if r["type"] == "HypernucleusDecay"
           and r["parent"].startswith(LL_PREFIX)]
    if not cap or not dec:
      continue
    row = {"fullscan": cap[0]["fullscan"] and dec[0]["fullscan"]}
    for k in cap[0]:
      if k.startswith("found_") and k in dec[0]:
        row[k] = cap[0][k] and dec[0][k]
    out.append(row)
  return out


def summarize(rows, extra, tolerances):
  has_best = any(k.startswith("found_best_") for k in rows[0])
  metrics = ["edge", "inclusive"] + (["best"] if has_best else [])
  cols = [f"{m[:4]}<{t:g}" for m in metrics for t in tolerances] \
    + [f"c.{m[:4]}<{t:g}" for m in metrics for t in tolerances] \
    + ["fullscan"]
  print("columns: efficiency per depth convention, then the control "
        "(c.) from shifted truths")
  print(f"{'type':<22} {'n':>4} " + " ".join(f"{c:>9}" for c in cols))
  for vtype in VERTEX_TYPES + ["LL event"]:
    if vtype == "LL event":
      sel = _ll_events(rows, False, tolerances)
      ctrl = _ll_events(rows, True, tolerances)
    else:
      sel, ctrl = _select(rows, vtype, False), _select(rows, vtype, True)
    if not sel:
      continue
    vals = [np.mean([r[f"found_{m}_{t:g}"] for r in sel])
            for m in metrics for t in tolerances]
    vals += [np.mean([r[f"found_{m}_{t:g}"] for r in ctrl]) if ctrl
             else float("nan") for m in metrics for t in tolerances]
    vals.append(np.mean([r["fullscan"] for r in sel]))
    print(f"{vtype:<22} {len(sel):>4} "
          + " ".join(f"{v:9.1%}" for v in vals))
  print("duplicates (unmatched candidates within tolerance, edge): "
        + ", ".join(f"<{t:g} um: {extra['duplicates'][t]}"
                    for t in tolerances))
  print("optimal - greedy matches (edge): "
        + ", ".join(f"<{t:g} um: {extra['optimal_minus_greedy'][t]}"
                    for t in tolerances))


def plot(rows, tolerances, output):
  import ROOT
  ROOT.gROOT.SetBatch(True)
  ROOT.gStyle.SetOptStat(0)
  canvas = ROOT.TCanvas("c", "", 1200, 500)
  canvas.Divide(2, 1)
  keep = []
  canvas.cd(1)
  frame = ROOT.gPad.DrawFrame(0., 0., max(tolerances) * 1.1, 1.05)
  frame.SetTitle("efficiency vs tolerance (edge);3-D tolerance (#mum);"
                 "efficiency")
  legend = ROOT.TLegend(0.5, 0.15, 0.9, 0.4)
  for k, vtype in enumerate(VERTEX_TYPES[:3]):
    for control, style in ((False, 1), (True, 2)):
      sel = _select(rows, vtype, control)
      if not sel:
        continue
      g = ROOT.TGraph(len(tolerances))
      for n, t in enumerate(tolerances):
        g.SetPoint(n, t, np.mean([r[f"found_edge_{t:g}"] for r in sel]))
      g.SetLineColor(k + 2)
      g.SetMarkerColor(k + 2)
      g.SetLineStyle(style)
      g.SetMarkerStyle(20)
      g.Draw("LP")
      legend.AddEntry(g, f"{vtype} ({'control' if control else 'true'})",
                      "lp")
      keep.append(g)
  legend.Draw()
  canvas.cd(2)
  t_mid = tolerances[len(tolerances) // 2]
  sel = [r for r in _select(rows, "capture", False)
         if r["ll_range_um"] is not None]
  lengths = np.array([r["ll_range_um"] for r in sel])
  edges = np.unique(np.quantile(lengths, np.linspace(
    0., 1., N_LENGTH_BINS + 1))) if sel else []
  if len(edges) >= 2:
    h_all = ROOT.TH1D("all", "", len(edges) - 1, edges)
    h_ok = ROOT.TH1D("ok", "", len(edges) - 1, edges)
    for r in sel:
      h_all.Fill(r["ll_range_um"])
      if r[f"found_edge_{t_mid:g}"]:
        h_ok.Fill(r["ll_range_um"])
    eff = ROOT.TEfficiency(h_ok, h_all)
    eff.SetTitle(f"capture vertex, edge < {t_mid:g} #mum;"
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
  parser.add_argument("--space", choices=["developed", "undeveloped"],
                      default="developed",
                      help="emulsion space of the 3-D matching")
  parser.add_argument("--n-shifts", type=int, default=DEFAULT_N_SHIFTS)
  parser.add_argument("--shift-min", type=float,
                      default=DEFAULT_SHIFT_MIN_PX)
  parser.add_argument("--shift-max", type=float,
                      default=DEFAULT_SHIFT_MAX_PX)
  parser.add_argument("--seed", type=int, default=0)
  parser.add_argument("--min-score", type=float,
                      help="keep candidates with n_tracks_max >= this")
  parser.add_argument("--details", help="write per-vertex rows to CSV")
  parser.add_argument("--nearest", action="store_true",
                      help="XY distance and z gap of the nearest "
                      "candidate, true vs control")
  parser.add_argument("-o", "--output", help="efficiency plot (PNG)")
  args = parser.parse_args()

  tolerances = [float(t) for t in args.tolerances.split(",")]
  events = load_truth(args.dataset)
  cands = load_candidates(args.candidates, args.min_score)
  shifts = random_shifts(args.n_shifts, args.shift_min, args.shift_max,
                         np.random.default_rng(args.seed))
  rows, extra = evaluate(events, cands, tolerances, shifts, args.space)
  print(f"3-D matching in the {args.space} emulsion "
        "(x, y unchanged; z = slice x physical_dz_um"
        + (" x Sh" if args.space == "undeveloped" else "") + "); "
        f"{len(events)} events, {extra['n_candidates']} candidates, "
        f"control = mean of {args.n_shifts} random shifts "
        f"({args.shift_min:g}-{args.shift_max:g} px)")
  summarize(rows, extra, tolerances)
  if args.nearest:
    print_nearest(nearest_components(events, cands, shifts, args.space))
  if args.details:
    with open(args.details, "w", newline="") as f:
      w = csv.DictWriter(f, fieldnames=list(rows[0]))
      w.writeheader()
      w.writerows(rows)
  if args.output:
    plot(rows, tolerances, args.output)


if __name__ == "__main__":
  main()
