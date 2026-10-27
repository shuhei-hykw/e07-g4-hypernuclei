"""From Geant4 tracks to a microscope image, step by step (one PNG).

  ./run.py pipelinefig data/images/<dataset>/ev0000
    --root data/sim/eval_nagara.root --model output/signal_model_v20.json
    [--half 46] [-o docs/pipeline_ev0000.png]

The event directory must be rendered into a real background view
(--background) without --bg-augment, from the given ROOT file and model.
All panels show +-half um around the capture point (image rows down).
Top row: (1) Geant4 tracks in the sheet plane and (2) seen from the side
(original, undeveloped z), (3) the developed silver grains placed along
the charged tracks, (4) their transmission in the slice in focus at the
capture point (defocused grains of the neighbouring depths included).
Bottom row: (5) that slice of the real background view, (6) the
synthetic slice = real x transmission, (7) the minimum projection over
all slices of the synthetic stack and (8) of the real view alone.
"""

import argparse
import json
from dataclasses import replace
from pathlib import Path

import numpy as np
import ROOT

from . import draw_event as de
from . import imagestats as st
from .synth import calibrate, render
from .synth.event import EventFile

CANVAS_W = 2000
CANVAS_H = 1050
N_COLS = 4
N_ROWS = 2
DEFAULT_HALF_UM = 46.
SIDE_HALF_Z_UM = 30.
GRAIN_MARKER = 20
GRAIN_MARKER_SIZE = 0.25
VERTEX_MARKER = 24
VERTEX_MARKER_SIZE = 1.6
LINE_WIDTH = {"capture": 3, "decay": 2, "other": 1}
GRAY_PCT = (0.5, 99.5)   # display range of the gray images
CLIP_EPS = 1e-3
TITLE_SIZE = 0.05
PAD_MARGIN = 0.13
LEGEND_BOX = (0.62, 0.14, 0.86, 0.40)


def colors_for(names):
  colors = dict(de.FIXED_COLORS, **{de.RECOIL: de.RECOIL_COLOR})
  free = 0
  for n in names:
    if n not in colors:
      colors[n] = de.COLORS[free % len(de.COLORS)]
      free += 1
  return colors


def species(t):
  return de.RECOIL if t.a > de.MAX_FRAGMENT_A else t.name


def gray_hist(name, img, half_px, pixel_um, lo=None, hi=None):
  """Image crop as a TH2 in um around the center, rows down."""
  n = img.shape[0]
  if lo is None:
    lo, hi = np.percentile(img, GRAY_PCT)
  # Values outside [lo, hi] would be left blank by COL
  img = np.clip(img, lo + CLIP_EPS * (hi - lo), hi - CLIP_EPS * (hi - lo))
  h = ROOT.TH2F(name, "", n, -half_px * pixel_um, half_px * pixel_um,
                n, -half_px * pixel_um, half_px * pixel_um)
  # Bin (i, j) = column i, row n - 1 - j (first row at the top)
  for j in range(n):
    for i in range(n):
      h.SetBinContent(i + 1, j + 1, float(img[n - 1 - j, i]))
  h.SetMinimum(lo)
  h.SetMaximum(hi)
  return h


def crop(img, c0, r0, half):
  return img[r0 - half:r0 + half, c0 - half:c0 + half]


def transmission(ev_dir, truth, model, c0, r0, half, k):
  """Transmission of the signal grains in slice k over the crop."""
  optics = replace(render.Optics(), size_px=2 * half,
                   n_slices=truth["n_slices"],
                   dz_um=truth["physical_dz_um"], shrink=truth["shrink"])
  optics, _, background = calibrate.apply_values(
    calibrate.load_model(model), optics, render.Emulsion(),
    render.Background())
  optics = replace(optics, size_px=2 * half)
  g = np.load(ev_dir / "grains.npz")
  grains = np.column_stack([g["col"] - (c0 - half), g["row"] - (r0 - half),
                            g["slice"]])
  d = render.Renderer(optics).density(grains, g["weight"], k)
  trans = np.exp(-d)
  if optics.image_blur_px > 0.:
    trans = render.gaussian_blur(trans, optics.image_blur_px)
  t_min = background.min_transmission
  return t_min + (1. - t_min) * trans


def main():
  parser = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument("event_dir")
  parser.add_argument("--root", required=True,
                      help="Geant4 file the event was rendered from")
  parser.add_argument("--model", required=True,
                      help="image model the event was rendered with")
  parser.add_argument("--half", type=float, default=DEFAULT_HALF_UM,
                      help="half width of the panels (um)")
  parser.add_argument("-o", "--output", default="pipeline.png")
  args = parser.parse_args()

  ev_dir = Path(args.event_dir)
  truth = json.loads((ev_dir / "truth.json").read_text())
  if not truth.get("background_view") or truth["background_transform"]:
    raise SystemExit("needs a real background without augmentation")
  pix = truth["pixel_um"]
  half = int(round(args.half / pix))
  cap = truth["vertices"][0]
  c0, r0 = int(round(cap["col"])), int(round(cap["row"]))
  k = int(round(cap["slice"]))
  ev = EventFile(args.root).event(truth["event"])
  v0 = ev.vertex

  ROOT.gROOT.SetBatch(True)
  ROOT.gStyle.SetOptStat(0)
  ROOT.gStyle.SetTitleSize(TITLE_SIZE, "t")
  ROOT.gStyle.SetPalette(ROOT.kGreyScale)
  canvas = ROOT.TCanvas("c", "", CANVAS_W, CANVAS_H)
  canvas.Divide(N_COLS, N_ROWS)
  keep = []

  def pad(i, title, xt, yt, y_half=args.half):
    canvas.cd(i)
    for m in ("Left", "Right", "Top", "Bottom"):
      getattr(ROOT.gPad, f"Set{m}Margin")(PAD_MARGIN)
    frame = ROOT.gPad.DrawFrame(-args.half, -y_half, args.half, y_half)
    frame.SetTitle(f"{title};{xt};{yt}")
    keep.append(frame)

  def marker(x, y):
    m = ROOT.TMarker(x, y, VERTEX_MARKER)
    m.SetMarkerSize(VERTEX_MARKER_SIZE)
    m.Draw()
    keep.append(m)

  charged = [t for t in ev.tracks if t.charge and len(t.points) > 1]
  colors = colors_for([species(t) for t in charged])
  legend = ROOT.TLegend(*LEGEND_BOX)
  shown = set()
  # (1), (2): Geant4 steps relative to the capture point
  for i, (iv, title, yt, y_half) in enumerate((
      (1, "(1) Geant4 tracks (sheet plane)", "#minus#Deltay (#mum)",
       args.half),
      (2, "(2) Geant4 tracks (side, original z)", "#Deltaz (#mum)",
       SIDE_HALF_Z_UM))):
    pad(i + 1, title, "#Deltax (#mum)", yt, y_half)
    sign = -1. if iv == 1 else 1.
    for t in charged:
      p = t.points - v0
      g = ROOT.TGraph(len(p), p[:, 0].copy(), (sign * p[:, iv]).copy())
      g.SetLineColor(colors[species(t)])
      g.SetLineWidth(LINE_WIDTH[t.origin])
      g.Draw("L")
      keep.append(g)
      if species(t) not in shown and t.origin != "other":
        legend.AddEntry(g, species(t), "l")
        shown.add(species(t))
    for v in truth["vertices"]:
      if v["in_view"]:
        d = np.array(v["um"]) - v0
        marker(d[0], sign * d[iv])
  canvas.cd(1)
  legend.Draw()

  # (3): developed silver grains, colored by the species of the track
  pad(3, f"(3) silver grains (z shrunk #times1/{truth['shrink']:g})",
      "#Deltax (#mum)", "#minus#Deltay (#mum)")
  g = np.load(ev_dir / "grains.npz")
  name_of = {t.id: species(t) for t in ev.tracks}
  for name in sorted(set(name_of.values())):
    sel = np.isin(g["instance"],
                  [i for i, n in name_of.items() if n == name])
    if not sel.any():
      continue
    x = (g["col"][sel] - cap["col"]) * pix
    y = -(g["row"][sel] - cap["row"]) * pix
    gr = ROOT.TGraph(int(sel.sum()), x.astype(float), y.astype(float))
    gr.SetMarkerStyle(GRAIN_MARKER)
    gr.SetMarkerSize(GRAIN_MARKER_SIZE)
    gr.SetMarkerColor(colors.get(name, de.RECOIL_COLOR))
    gr.Draw("P")
    keep.append(gr)

  # (4)-(8): images
  real_k = crop(st.read_stack(truth["background_view"], [k])[0],
                c0, r0, half).astype(float)
  synth_k = crop(st.read_stack(ev_dir / "image.json", [k])[0],
                 c0, r0, half).astype(float)
  synth_min = crop(np.min(st.read_stack(ev_dir / "image.json"), axis=0),
                   c0, r0, half).astype(float)
  real_min = crop(np.min(st.read_stack(truth["background_view"]), axis=0),
                  c0, r0, half).astype(float)
  trans = transmission(ev_dir, truth, args.model, c0, r0, half, k)
  lo_k, hi_k = np.percentile(np.r_[real_k.ravel(), synth_k.ravel()],
                             GRAY_PCT)
  lo_m, hi_m = np.percentile(np.r_[real_min.ravel(), synth_min.ravel()],
                             GRAY_PCT)
  for i, (title, img, lo, hi) in enumerate((
      (f"(4) transmission of the grains, slice {k}", trans,
       float(trans.min()), 1.),
      (f"(5) real background view, slice {k}", real_k, lo_k, hi_k),
      ("(6) synthetic slice = (5) #times (4)", synth_k, lo_k, hi_k),
      ("(7) synthetic: minimum over slices", synth_min, lo_m, hi_m),
      ("(8) real view alone: minimum", real_min, lo_m, hi_m))):
    pad(i + 4, title, "#Deltax (#mum)", "#minus#Deltay (#mum)")
    h = gray_hist(f"h{i}", img, half, pix, lo, hi)
    h.Draw("COL SAME")
    ROOT.gPad.RedrawAxis()
    keep.append(h)

  canvas.cd()
  Path(args.output).parent.mkdir(parents=True, exist_ok=True)
  canvas.SaveAs(str(args.output))


if __name__ == "__main__":
  main()
