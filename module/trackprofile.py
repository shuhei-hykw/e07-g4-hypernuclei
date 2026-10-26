"""Appearance of embedded signal tracks vs real background lines in the
same images (calibration of the signal rendering).

  ~/work/.venv/bin/python -m module.trackprofile data/images/<dataset> \
    [--events 20] [--slice-step 6] [-o output/trackprofile.md]

On embedded single-track datasets (./run.py singletrack, real eval
backgrounds), centre lines of the e07-fullscan line detector are found
in every --slice-step-th slice. A line is "signal" if most of its
pixels lie within NEAR_PX of signal grains within SIGNAL_SLICES of the
slice (grains.npz), "background" if none does (real tracks of the view).
Along each line: mean darkness of the background-subtracted slice,
its relative spread along the line (beading), the fraction of gaps
(darkness below GAP_LEVEL) and the width (dark pixels within WIDTH_PX
of the line per unit length). Signal lines are split by the dE/dx of
their track (single_tracks.csv) and by dip.
"""

import argparse
import csv
from pathlib import Path

import numpy as np

from . import imagestats as st
from . import linecheck

MIN_LEN_PX = 40
NEAR_PX = 4
SIGNAL_SLICES = 6.
SIGNAL_FRACTION = 0.6
GAP_LEVEL = 5.        # gray levels of darkness
WIDTH_PX = 4
WIDTH_LEVEL = 8.      # gray levels of darkness
DEDX_BINS = (0., 0.0015, 0.005, 0.02, 1.)   # MeV / um (MIP ~ 0.0006)
DIP_FLAT_DEG, DIP_STEEP_DEG = 15., 60.
# Comparison at equal apparent darkness (no truth for real lines)
DARK_BINS = (5., 10., 15., 20., 30., 50., 255.)


def line_stats(pix, d, wide):
  v = d[pix[:, 0], pix[:, 1]]
  m = max(v.mean(), 1e-6)
  return {"dark": float(v.mean()), "beading": float(v.std() / m),
          "gap": float(np.mean(v < GAP_LEVEL)),
          "width": float(wide / len(pix))}


def event_lines(ev_dir, ep, step, tracks):
  import cv2
  path = ev_dir / "image.json"
  n = st.n_slices(path)
  g = np.load(ev_dir / "grains.npz")
  sig = g["cls"] == 1
  out = []
  for k in range(step // 2, n, step):
    img = st.read_stack(path, [k])[0]
    d = st.gaussian_blur(img.astype(float), st.BG_SIGMA_PX) - img
    centre = (ep.line_mask(ep.ridge_response(img)) > 0).astype(np.uint8)
    near = np.zeros_like(centre)
    owner = np.zeros(centre.shape, np.int32)
    sel = sig & (np.abs(g["slice"] - k) < SIGNAL_SLICES)
    r = np.rint(g["row"][sel]).astype(int)
    c = np.rint(g["col"][sel]).astype(int)
    ok = (r >= 0) & (r < img.shape[0]) & (c >= 0) & (c < img.shape[1])
    near[r[ok], c[ok]] = 1
    owner[r[ok], c[ok]] = g["instance"][sel][ok]
    kern = np.ones((2 * NEAR_PX + 1,) * 2, np.uint8)
    near = cv2.dilate(near, kern)
    owner = cv2.dilate(owner.astype(np.float32), kern).astype(np.int32)
    darkmask = (d > WIDTH_LEVEL).astype(np.uint8)
    wkern = np.ones((2 * WIDTH_PX + 1,) * 2, np.uint8)
    nl, lab = cv2.connectedComponents(centre, connectivity=8)
    for i in range(1, nl):
      pix = np.argwhere(lab == i)
      if len(pix) < MIN_LEN_PX:
        continue
      frac = near[pix[:, 0], pix[:, 1]].mean()
      if 0. < frac < SIGNAL_FRACTION:
        continue
      y0, x0 = pix.min(axis=0) - WIDTH_PX
      y1, x1 = pix.max(axis=0) + WIDTH_PX + 1
      y0, x0 = max(y0, 0), max(x0, 0)
      sub = (lab[y0:y1, x0:x1] == i).astype(np.uint8)
      wide = int((cv2.dilate(sub, wkern) & darkmask[y0:y1, x0:x1]).sum())
      rec = line_stats(pix, d, wide)
      if frac == 0.:
        rec["kind"] = "background"
      else:
        ids = owner[pix[:, 0], pix[:, 1]]
        ids = ids[ids > 0]
        t = tracks.get(int(np.bincount(ids).argmax())) if len(ids) else None
        if t is None:
          continue
        rec.update(kind="signal", dedx=float(t["dedx_mev_per_um"]),
                   dip=abs(float(t["dip_deg"])), species=t["species"])
      out.append(rec)
  return out


def summarize(rows):
  keys = ("dark", "beading", "gap", "width")
  if not rows:
    return "n=0"
  return f"n={len(rows)} " + " ".join(
    f"{k} {np.median([r[k] for r in rows]):.3g}" for k in keys)


def main():
  parser = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument("dataset")
  parser.add_argument("--events", type=int, default=20)
  parser.add_argument("--slice-step", type=int, default=6)
  parser.add_argument("-o", "--output")
  args = parser.parse_args()
  ep = linecheck.import_fullscan()
  ds = Path(args.dataset)
  with open(ds / "single_tracks.csv") as f:
    table = list(csv.DictReader(f))
  rows = []
  for ev_dir in sorted(ds.glob("ev*"))[:args.events]:
    ev = int(ev_dir.name[2:])
    tracks = {int(t["track_id"]): t for t in table if int(t["event"]) == ev}
    rows += event_lines(ev_dir, ep, args.slice_step, tracks)
    print(ev_dir.name, len(rows), flush=True)
  bg = [r for r in rows if r["kind"] == "background"]
  sg = [r for r in rows if r["kind"] == "signal"]
  lines = [f"# Signal vs real background lines in {ds.name}", "",
           f"e07-fullscan {linecheck.fullscan_commit()}; medians per line "
           f"(darkness in gray levels, beading = spread / mean along the "
           f"line, gap = fraction below {GAP_LEVEL:g}, width = dark px "
           f"within {WIDTH_PX} px per px of line).", "",
           f"- real background lines: {summarize(bg)}",
           f"- signal lines (all): {summarize(sg)}"]
  for lo, hi in zip(DEDX_BINS[:-1], DEDX_BINS[1:]):
    sel = [r for r in sg if lo <= r["dedx"] < hi]
    lines.append(f"- signal dE/dx {lo:g}-{hi:g} MeV/um: {summarize(sel)}")
  for name, cond in (("flat", lambda a: a < DIP_FLAT_DEG),
                     ("steep", lambda a: a > DIP_STEEP_DEG)):
    sel = [r for r in sg if cond(r["dip"])]
    lines.append(f"- signal {name}: {summarize(sel)}")
  lines += ["", "At equal apparent darkness (real lines have no truth):",
            "", "| darkness | real: n, beading, gap, width | "
            "signal: n, beading, gap, width | signal flat (dip < "
            f"{DIP_FLAT_DEG:g} deg) |", "|---|---|---|---|"]
  flat = [r for r in sg if r["dip"] < DIP_FLAT_DEG]
  for lo, hi in zip(DARK_BINS[:-1], DARK_BINS[1:]):
    cells = []
    for group in (bg, sg, flat):
      sel = [r for r in group if lo <= r["dark"] < hi]
      cells.append(f"{len(sel)}, " + ", ".join(
        f"{np.median([r[k] for r in sel]):.3g}" if sel else "-"
        for k in ("beading", "gap", "width")))
    lines.append(f"| {lo:g}-{hi:g} | " + " | ".join(cells) + " |")
  text = "\n".join(lines) + "\n"
  if args.output:
    Path(args.output).write_text(text)
  print(text)


if __name__ == "__main__":
  main()
