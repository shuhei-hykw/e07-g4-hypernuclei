"""Steep straight lines (beam tracks) in real views vs minimum-ionizing
signal tracks, with independent selection and measurement slices.

  ~/work/.venv/bin/python -m module.beamlines data/images/<single-track
    dataset> [--events 20] [-o output/beamlines_v20.md]

The beam crosses the plates nearly perpendicularly and is minimum
ionizing. Following the darkest pixel (module.mipcheck) did not isolate
it from the texture, so here a beam track must be a straight 3-D line:
dark local peaks are found in the even slices inside the emulsion
(selection) and joined into lines through anchor pairs; a line needs
peaks within HIT_PX of its prediction in MIN_HIT_FRACTION of the
selection slices and a lateral drift below MAX_DRIFT_PX per slice. The
odd slices (measurement) are never used for the selection: there each
line gives the darkness at its position, the gap fraction and the
radial profile, against the same line shifted by SHIFT_PX (background
level). Coincidence lines are counted in scrambled stacks whose
selection slices come from different views (no real continuity).

Lines of an embedded event are split by truth into signal (matched to
a steep signal track) and real (background view, away from all signal
grains); in a whole-synthetic view (no background view) into synthetic
beam (rendered beam grains at the line) and other.

Result (output/beamlines_real.md, output/beamlines_checks.md,
output/beamlines_beam_v20.md): real views hold ~27,000 steep lines per
mm2 (none in the scrambled control); perpendicular v20 minimum-ionizing
tracks embedded into real views are about 8 % fainter and 1 px wider
than the real lines under the same selection.
"""

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy import ndimage
from scipy.spatial import cKDTree

from . import mipcheck as mc

PEAK_SIGMA_PX = 1.
PEAK_SIZE_PX = 5          # local maximum filter
PEAK_LEVEL = 8.           # gray levels below the local background
HIT_PX = 2.
MIN_HIT_FRACTION = 0.6
MAX_DRIFT_PX = 8.         # per slice, ~30 deg from the normal at dz 4
ANCHOR_QUANTILES = ((0.2, 0.8), (0.3, 0.7), (0.15, 0.6), (0.4, 0.85))
DEDUP_PX = 4.
EDGE_PX = 12
SHIFT_PX = 25
PROFILE_PX = 8
MATCH_PX = 3.
MATCH_FRACTION = 0.5
LABEL_BEAM = 4            # grain class of rendered beam tracks
STEEP_DEG = mc.STEEP_DEG
MIP_DEDX_MAX = mc.MIP_DEDX_MAX
PCT = mc.PCT
DARK_LEVEL = mc.DARK_LEVEL
MM_PER_UM = 1e-3
TINY = 1e-9


def peaks(d):
  s = ndimage.gaussian_filter(d, PEAK_SIGMA_PX)
  m = (s == ndimage.maximum_filter(s, PEAK_SIZE_PX)) & (s > PEAK_LEVEL)
  ys, xs = np.nonzero(m)
  return np.column_stack([ys, xs]).astype(float)


def find_lines(stack, sel):
  """Lines (y0, x0, dy, dx at slice k0 = sel middle) from selection
  slices; stack[k] are darkness images (stack may be scrambled)."""
  pk = {k: peaks(stack[k]) for k in sel}
  trees = {k: cKDTree(pk[k]) for k in sel}
  ks = np.asarray(sel, float)
  k0 = float(sel[len(sel) // 2])
  found = []
  for qa, qb in ANCHOR_QUANTILES:
    a, b = sel[int(qa * (len(sel) - 1))], sel[int(qb * (len(sel) - 1))]
    if not len(pk[a]) or not len(pk[b]):
      continue
    pairs = trees[a].query_ball_tree(trees[b], MAX_DRIFT_PX * (b - a))
    ia = np.repeat(np.arange(len(pairs)), [len(p) for p in pairs])
    ib = np.concatenate([np.asarray(p, int) for p in pairs])
    if not len(ia):
      continue
    pa, pb = pk[a][ia], pk[b][ib]
    slope = (pb - pa) / (b - a)
    hit = np.zeros((len(ia), len(sel)), bool)
    pos = np.zeros((len(ia), len(sel), 2))
    for j, k in enumerate(sel):
      pred = pa + slope * (k - a)
      dist, idx = trees[k].query(pred, distance_upper_bound=HIT_PX)
      hit[:, j] = np.isfinite(dist)
      pos[hit[:, j], j] = pk[k][idx[hit[:, j]]]
    ok = hit.mean(axis=1) >= MIN_HIT_FRACTION
    for h, p in zip(hit[ok], pos[ok]):
      kk = ks[h] - k0
      fy = np.polyfit(kk, p[h, 0], 1)
      fx = np.polyfit(kk, p[h, 1], 1)
      if np.hypot(fy[0], fx[0]) > MAX_DRIFT_PX:
        continue
      # Support recounted around the fitted line
      line = (fy[1], fx[1], fy[0], fx[0])
      n = sum(np.isfinite(trees[k].query(
        line_pos(line, k, k0), distance_upper_bound=HIT_PX)[0])
        for k in sel)
      if n >= MIN_HIT_FRACTION * len(sel):
        found.append((*line, n))
  found.sort(key=lambda f: -f[4])
  lines = []
  for f in found:
    if all(np.hypot(f[0] - g[0], f[1] - g[1]) > DEDUP_PX for g in lines):
      lines.append(f)
  return np.array(lines).reshape(-1, 5), k0


def line_pos(line, k, k0):
  return line[0] + line[2] * (k - k0), line[1] + line[3] * (k - k0)


def measure(det, line, meas, k0, shift=0):
  """Darkness (window max), gap fraction and radial profile on the
  measurement slices; None if the line leaves the crop."""
  vals, prof = [], np.zeros(PROFILE_PX + 1)
  rr = np.hypot(*np.mgrid[-PROFILE_PX:PROFILE_PX + 1,
                          -PROFILE_PX:PROFILE_PX + 1])
  ring = np.rint(rr).astype(int)
  for k in meas:
    y, x = line_pos(line, k, k0)
    y, x = int(round(y)), int(round(x)) + shift
    if not (EDGE_PX <= y < mc.CROP_PX - EDGE_PX
            and EDGE_PX <= x < mc.CROP_PX - EDGE_PX):
      return None
    vals.append(mc.window_max(det[k], y, x))
    win = det[k, y - PROFILE_PX:y + PROFILE_PX + 1,
              x - PROFILE_PX:x + PROFILE_PX + 1]
    prof += [win[ring == r].mean() for r in range(PROFILE_PX + 1)]
  v = np.array(vals)
  prof /= len(meas)
  half = (prof[0] + prof[-1]) / 2.
  below = np.flatnonzero(prof < half)
  hwhm = float(below[0]) if len(below) else float(PROFILE_PX)
  q = v - v.mean()
  autocorr = float((q[:-1] * q[1:]).sum() / max((q * q).sum(), TINY))
  return dict(dark=float(np.median(v)), autocorr=autocorr,
              gaps=float((v < DARK_LEVEL).mean()),
              contrast=float(prof[0] - prof[-1]), hwhm=hwhm)


def inside_slices(det):
  prof = (det > mc.LAYER_LEVEL).mean(axis=(1, 2))
  return np.flatnonzero(prof > mc.EMPTY_FRACTION * prof.max())


def signal_tracks(ev_dir, table, lo):
  """Steep MIP signal tracks: (track id, grain rows, cols, slices)."""
  g = np.load(ev_dir / "grains.npz")
  ev = int(ev_dir.name[2:])
  out = []
  for t in table:
    if int(t["event"]) != ev:
      continue
    m = (g["cls"] == 1) & (g["instance"] == int(t["track_id"]))
    out.append((t, g["row"][m] - lo, g["col"][m] - lo,
                g["focus_slice"][m]))
  allsig = g["cls"] == 1
  return out, g["row"][allsig] - lo, g["col"][allsig] - lo


def matched(line, k0, rows, cols, slices):
  if not len(rows):
    return False
  y, x = line_pos(line, slices, k0)
  return (np.hypot(rows - y, cols - x) < MATCH_PX).mean() >= MATCH_FRACTION


def near_signal(line, k0, ks, rows, cols):
  if not len(rows):
    return False
  tree = cKDTree(np.column_stack([rows, cols]))
  y, x = line_pos(line, np.asarray(ks, float), k0)
  d, _ = tree.query(np.column_stack([y, x]))
  return bool((d < mc.SIGNAL_EXCLUSION_PX).any())


def drift_deg(line, pixel_um, dz_um):
  return float(np.degrees(np.arctan(
    np.hypot(line[2], line[3]) * pixel_um / dz_um)))


def beam_match(line, k0, g, lo, used):
  """Rendered beam track whose grains follow the line (MATCH_FRACTION
  within MATCH_PX), each track matched at most once; None if none."""
  b = g["cls"] == LABEL_BEAM
  m = b & (np.abs(g["focus_slice"] - k0) <= 1)
  y, x = line_pos(line, k0, k0)
  near = np.hypot(g["row"][m] - lo - y, g["col"][m] - lo - x) < MATCH_PX
  for i in np.unique(g["instance"][m][near]):
    if i in used:
      continue
    t = b & (g["instance"] == i)
    if matched(line, k0, g["row"][t] - lo, g["col"][t] - lo,
               g["focus_slice"][t]):
      used.add(i)
      return i
  return None


def crossing_beams(g, lo, inside):
  """Rendered beam tracks with grains inside the crop in at least half
  of the inside slices (the ones a line can recover)."""
  b = g["cls"] == LABEL_BEAM
  r, c = g["row"][b] - lo, g["col"][b] - lo
  ink = (r >= 0) & (r < mc.CROP_PX) & (c >= 0) & (c < mc.CROP_PX)
  inst, sl = g["instance"][b][ink], g["focus_slice"][b][ink]
  return sum(len(np.unique(sl[inst == i])) >= len(inside) / 2.
             for i in np.unique(inst))


def event_rows(ev_dir, table, pixel_um, dz_um, synthetic):
  det, lo = mc.darkness_stack(ev_dir / "image.json")
  inside = inside_slices(det)
  sel, meas = list(inside[0::2]), list(inside[1::2])
  lines, k0 = find_lines(det, sel)
  tracks, srows, scols = signal_tracks(ev_dir, table, lo)
  g = np.load(ev_dir / "grains.npz") if synthetic else None
  used = set()
  out = []
  for line in lines:
    kind = None
    for t, r, c, s in tracks:
      if matched(line, k0, r, c, s):
        kind = ("signal_mip" if float(t["dedx_mev_per_um"])
                < MIP_DEDX_MAX else "signal_dense")
        break
    if kind is None and synthetic:
      kind = ("synth_beam" if beam_match(line, k0, g, lo, used)
              is not None else "synth_other")
    if kind is None:
      if near_signal(line, k0, inside, srows, scols):
        continue
      kind = "real"
    m = measure(det, line, meas, k0)
    c = measure(det, line, meas, k0, SHIFT_PX)
    if m is None or c is None:
      continue
    out.append(dict(kind=kind, hits=line[4] / len(sel),
                    angle=drift_deg(line, pixel_um, dz_um), **m,
                    **{f"{key}_shift": v for key, v in c.items()}))
  n_beams = crossing_beams(g, lo, inside) if synthetic else 0
  return out, det, inside, n_beams, len(lines)


def scrambled_lines(stacks):
  """Line count per view when each selection slice of a view comes from
  another view (coincidences only)."""
  n = len(stacks)
  out = []
  for i, (det, inside) in enumerate(stacks):
    sel = list(inside[0::2])
    mixed = {}
    for j, k in enumerate(sel):
      other, oin = stacks[(i + 1 + j) % n]
      kk = oin[min(len(oin) - 1, int(j * len(oin) / len(sel)))]
      mixed[k] = other[kk]
    lines, _ = find_lines(mixed, sel)
    out.append(len(lines))
  return out


def pct(rows, key):
  if not rows:
    return "-"
  return " / ".join(f"{np.percentile([r[key] for r in rows], p):.3g}"
                    for p in PCT)


def main():
  parser = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument("dataset")
  parser.add_argument("--events", type=int, default=20)
  parser.add_argument("--scramble", type=int, default=8,
                      help="views used for the scrambled control")
  parser.add_argument("-o", "--output")
  args = parser.parse_args()
  ds = Path(args.dataset)
  table = []
  if (ds / "single_tracks.csv").exists():
    with open(ds / "single_tracks.csv") as f:
      table = [t for t in csv.DictReader(f)
               if abs(float(t["dip_deg"])) > STEEP_DEG]
  rows, stacks, n_lines, n_beams = [], [], [], []
  for ev_dir in sorted(ds.glob("ev*"))[:args.events]:
    truth = json.loads((ev_dir / "truth.json").read_text())
    synthetic = truth.get("background_view") is None
    r, det, inside, nb, nf = event_rows(ev_dir, table, truth["pixel_um"],
                                truth["physical_dz_um"], synthetic)
    rows += r
    n_lines.append(nf)
    n_beams.append(nb)
    area_mm2 = (mc.CROP_PX * truth["pixel_um"] * MM_PER_UM) ** 2
    if len(stacks) < args.scramble:
      stacks.append((det, inside))
    print(ev_dir.name, len(r), "lines", flush=True)
  scr = scrambled_lines(stacks)
  n_synth = sum(r["kind"] == "synth_beam" for r in rows)
  recall = (f" Rendered beam tracks crossing half the layer in the crop: "
            f"{np.mean(n_beams):.1f} per view, uniquely matched lines "
            f"{n_synth / len(n_beams):.1f} per view (after the edge "
            f"cuts)." if sum(n_beams) else "")
  out = [f"# Steep lines in {ds.name}", "",
         f"{len(n_lines)} events; lines found per view (no cuts) "
         f"{np.mean(n_lines):.1f} ({np.mean(n_lines) / area_mm2:.0f} "
         f"per mm2 of crop), scrambled-control lines per view "
         f"{np.mean(scr):.1f} ({len(scr)} views, same finder).{recall} "
         f"Table: lines inside the crop edges and away from signal "
         f"grains; cells 16 / 50 / 84 %; measurement slices only; "
         f"'shifted' = the same line moved by {SHIFT_PX} px.", "",
         "| sample | n | angle (deg) | hit fraction | darkness | shifted "
         "| gaps | shifted | contrast | HWHM (px) | lag-1 autocorr |",
         "|---|---|---|---|---|---|---|---|---|---|---|"]
  for kind in ("real", "signal_mip", "signal_dense", "synth_beam",
               "synth_other"):
    rs = [r for r in rows if r["kind"] == kind]
    out.append(f"| {kind} | {len(rs)} | {pct(rs, 'angle')} | "
               f"{pct(rs, 'hits')} | {pct(rs, 'dark')} | "
               f"{pct(rs, 'dark_shift')} | {pct(rs, 'gaps')} | "
               f"{pct(rs, 'gaps_shift')} | {pct(rs, 'contrast')} | "
               f"{pct(rs, 'hwhm')} | {pct(rs, 'autocorr')} |")
  text = "\n".join(out) + "\n"
  print(text)
  if args.output:
    Path(args.output).write_text(text)


if __name__ == "__main__":
  main()
