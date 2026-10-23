"""Geant4 tracks drawn in the pixel / slice coordinates of a view.

Panels: crop around the capture point (same as preview.png), the whole
view, and the side view (column vs slice) of the whole view. Coordinates
are those of the synthetic images (developed emulsion, image rows down).
"""

import numpy as np
import ROOT

from .. import draw_event as de

CANVAS_W = 1800
CANVAS_H = 650
CROP_HALF_PX = 256
LEGEND_X = 0.90
TITLE_Y = 0.94
TITLE_SIZE = 0.03
VERTEX_MARKER = 24
VERTEX_COLOR = ROOT.kBlack
LINE_WIDTH = {"capture": 3, "decay": 2, "other": 1}


def draw_view(ev, view, output):
  ROOT.gROOT.SetBatch(True)
  ROOT.gStyle.SetOptStat(0)
  o = view.optics
  tracks = []
  for t in ev.tracks:
    if t.charge == 0 or len(t.points) < 2:
      continue
    name = de.RECOIL if t.a > de.MAX_FRAGMENT_A else t.name
    tracks.append((name, t.origin, view.to_view(t.points)))
  colors = dict(de.FIXED_COLORS, **{de.RECOIL: de.RECOIL_COLOR})
  n_free = 0
  for name, _, _ in tracks:
    if name not in colors:
      colors[name] = de.COLORS[n_free % len(de.COLORS)]
      n_free += 1
  vtx = view.to_view(ev.vertex[None, :])[0]

  canvas = ROOT.TCanvas("c_view", "", CANVAS_W, CANVAS_H)
  title = ROOT.TLatex(0.01, (1. + TITLE_Y) / 2.,
                      f"Geant4 tracks in image coordinates: Xi^{{-}} + "
                      f"{ev.channel} (event {ev.index}), capture at col "
                      f"{vtx[0]:.0f}, row {vtx[1]:.0f}, slice {vtx[2]:.1f}")
  title.SetNDC()
  title.SetTextSize(TITLE_SIZE)
  title.SetTextAlign(12)
  title.Draw()
  grid = ROOT.TPad("grid", "", 0., 0., LEGEND_X, TITLE_Y)
  grid.Draw()
  grid.Divide(3, 1)
  legend = ROOT.TLegend(LEGEND_X, 0.3, 1., TITLE_Y)
  keep = [grid, legend, title]
  # Crop center clipped like preview.png
  c0, r0 = np.clip(vtx[:2], CROP_HALF_PX, o.size_px - CROP_HALF_PX)
  # Image rows grow downwards: rows are plotted as -row so that the
  # panels look like preview.png
  panels = [
    ("crop (as preview.png)", 0, 1,
     (c0 - CROP_HALF_PX, -r0 - CROP_HALF_PX, c0 + CROP_HALF_PX,
      -r0 + CROP_HALF_PX), "column (px)", "-row (px)"),
    ("whole view", 0, 1, (0., -o.size_px, o.size_px, 0.),
     "column (px)", "-row (px)"),
    ("side view", 0, 2, (0., -1., o.size_px, o.n_slices),
     "column (px)", "slice"),
  ]
  shown = set()
  for k, (label, ih, iv, (x0, y0, x1, y1), xt, yt) in enumerate(panels):
    grid.cd(k + 1)
    ROOT.gPad.SetLeftMargin(0.14)
    frame = ROOT.gPad.DrawFrame(x0, y0, x1, y1)
    frame.SetTitle(f"{label};{xt};{yt}")
    sign = -1. if iv == 1 else 1.
    for name, origin, v in tracks:
      ys = sign * v[:, iv]
      g = ROOT.TGraph(len(v), v[:, ih].astype(float), ys.astype(float))
      g.SetLineColor(colors[name])
      g.SetLineWidth(LINE_WIDTH[origin])
      g.Draw("L")
      keep.append(g)
      if name not in shown:
        legend.AddEntry(g, name, "l")
        shown.add(name)
    m = ROOT.TMarker(vtx[ih], sign * vtx[iv], VERTEX_MARKER)
    m.SetMarkerColor(VERTEX_COLOR)
    m.Draw()
    keep.append(m)
  canvas.cd()
  legend.Draw()
  canvas.SaveAs(str(output))
