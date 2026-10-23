"""Event display with the emulsion module layers.

  ./run.py draw data/sim/nagara.root --event 0 --zoom 40 -o event0.png

Top row: +-zoom um around the capture point; bottom row: whole event.
Columns: XY (sheet plane, microscope view), ZX and ZY (z = beam axis,
base films shaded grey). Coordinates are relative to the capture point.
Capture products are drawn thick, hypernuclear decay products medium and
other secondaries thin; neutral tracks are dashed and hidden unless
--neutral. Tracks leaving the module end with an open circle.
"""

import argparse

import ROOT

CANVAS_W = 1500
CANVAS_H = 1000
FRAME_MARGIN = 1.1
WIDTH_PRIMARY = 3
WIDTH_DECAY = 2
WIDTH_OTHER = 1
NEUTRAL_STYLE = 2  # dashed
ESCAPE_MARKER = 24  # open circle
ESCAPE_MARKER_SIZE = 1.2
BASE_COLOR = ROOT.kGray
OUTSIDE_COLOR = ROOT.kGray + 2
OUTSIDE_STYLE = 3354  # hatched
# Canvas layout (NDC): title strip on top, legend strip on the right
TITLE_Y = 0.95
LEGEND_X = 0.88
TITLE_SIZE = 0.025
PAD_LEFT_MARGIN = 0.16
DECAY_PROCESS = "HypernucleusDecay"
FIXED_COLORS = {
  "pi-": ROOT.kMagenta + 1, "proton": ROOT.kOrange + 7,
  "deuteron": ROOT.kYellow + 3, "triton": ROOT.kRed + 1,
  "He3": ROOT.kPink + 6, "alpha": ROOT.kBlue + 1,
  "neutron": ROOT.kGray + 1, "lambda": ROOT.kSpring - 6,
}
# Remaining species (hypernuclei, heavier fragments) in order of appearance
COLORS = [ROOT.kGreen + 2, ROOT.kCyan + 2, ROOT.kAzure + 7,
          ROOT.kViolet + 1, ROOT.kTeal - 6, ROOT.kOrange + 3]
# Nuclei heavier than this are drawn as "recoil" in one color
MAX_FRAGMENT_A = 16
RECOIL = "recoil (A>16)"
RECOIL_COLOR = ROOT.kGray + 2
# (horizontal, vertical, label) of each projection
VIEWS = [("x", "y", "XY (sheet plane)"), ("z", "x", "ZX"),
         ("z", "y", "ZY")]


def load_tracks(tree, with_neutral):
  v = {"x": tree.vx, "y": tree.vy, "z": tree.vz}
  tracks = []
  for j in range(len(tree.name)):
    if tree.z[j] == 0 and not with_neutral:
      continue
    name = str(tree.name[j])
    if tree.a[j] > MAX_FRAGMENT_A:
      name = RECOIL
    tracks.append({
      "name": name,
      "charge": int(tree.z[j]),
      "parent": int(tree.parent[j]),
      "process": str(tree.process[j]),
      "escaped": bool(tree.escaped[j]),
      "x": [p - v["x"] for p in tree.x[j]],
      "y": [p - v["y"] for p in tree.y[j]],
      "z": [p - v["z"] for p in tree.zpos[j]],
    })
  return tracks


def load_layers(f, vz):
  layers = []
  for l in f.Get("layers"):
    layers.append((bool(l.emulsion), l.zmin - vz, l.zmax - vz))
  return layers


def line_width(track):
  if track["parent"] == 0:
    return WIDTH_PRIMARY
  if track["process"] == DECAY_PROCESS:
    return WIDTH_DECAY
  return WIDTH_OTHER


def draw_layers(layers, lo, hi, keep):
  """Shade base films and the outside of the module along z."""
  def box(z0, z1, color, style):
    z0, z1 = max(z0, lo), min(z1, hi)
    if z0 >= z1:
      return
    b = ROOT.TBox(z0, lo, z1, hi)
    b.SetFillColor(color)
    b.SetFillStyle(style)
    b.Draw()
    keep.append(b)

  solid = 1001
  for emulsion, z0, z1 in layers:
    if not emulsion:
      box(z0, z1, BASE_COLOR, solid)
  box(lo, layers[0][1], OUTSIDE_COLOR, OUTSIDE_STYLE)
  box(layers[-1][2], hi, OUTSIDE_COLOR, OUTSIDE_STYLE)


def draw_panel(tracks, layers, view, half, colors, keep, legend=None):
  h, v, label = view
  frame = ROOT.gPad.DrawFrame(-half, -half, half, half)
  frame.SetTitle(f"{label};{h} (#mum);{v} (#mum)")
  if h == "z":
    draw_layers(layers, -half, half, keep)
    ROOT.gPad.RedrawAxis()
  shown = set()
  for t in tracks:
    n = len(t[h])
    graph = ROOT.TGraph(n)
    for i in range(n):
      graph.SetPoint(i, t[h][i], t[v][i])
    graph.SetLineColor(colors[t["name"]])
    graph.SetLineWidth(line_width(t))
    if t["charge"] == 0:
      graph.SetLineStyle(NEUTRAL_STYLE)
    graph.Draw("L")
    keep.append(graph)
    if t["escaped"]:
      m = ROOT.TMarker(t[h][-1], t[v][-1], ESCAPE_MARKER)
      m.SetMarkerColor(colors[t["name"]])
      m.SetMarkerSize(ESCAPE_MARKER_SIZE)
      m.Draw()
      keep.append(m)
    if legend is not None and t["name"] not in shown:
      legend.AddEntry(graph, t["name"], "l")
      shown.add(t["name"])


def main():
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument("root_file")
  parser.add_argument("--event", type=int, default=0)
  parser.add_argument("--zoom", type=float, default=40.,
                      help="half width of the top row in um")
  parser.add_argument("--neutral", action="store_true",
                      help="also draw neutral tracks")
  parser.add_argument("-o", "--output", default="event.png")
  args = parser.parse_args()
  draw(args.root_file, args.event, args.output, args.zoom, args.neutral)


def draw(root_file, event, output, zoom=40., neutral=False):
  """Write the event display of one event to `output`."""
  ROOT.gROOT.SetBatch(True)
  ROOT.gStyle.SetOptStat(0)
  f = ROOT.TFile.Open(str(root_file))
  tree = f.Get("tree")
  tree.GetEntry(event)

  tracks = load_tracks(tree, neutral)
  layers = load_layers(f, tree.vz)
  # The wide view covers the capture and hypernuclear decay products
  chain = [t for t in tracks
           if t["parent"] == 0 or t["process"] == DECAY_PROCESS]
  wide = FRAME_MARGIN * max(
    max(abs(p) for p in t["x"] + t["y"] + t["z"]) for t in chain)
  colors = dict(FIXED_COLORS, **{RECOIL: RECOIL_COLOR})
  n_free = 0
  for t in tracks:
    if t["name"] not in colors:
      colors[t["name"]] = COLORS[n_free % len(COLORS)]
      n_free += 1

  canvas = ROOT.TCanvas("c", "", CANVAS_W, CANVAS_H)
  status = "all charged tracks contained" if tree.contained else (
    "escaping track(s), marked o")
  title = ROOT.TLatex(0.01, (1. + TITLE_Y) / 2.,
                      f"Xi^{{-}} + {tree.channel}   (event {event}, "
                      f"vertex z = {tree.vz:.0f} #mum, {status})")
  title.SetNDC()
  title.SetTextAlign(12)
  title.SetTextSize(TITLE_SIZE)
  title.Draw()
  grid = ROOT.TPad("grid", "", 0., 0., LEGEND_X, TITLE_Y)
  grid.Draw()
  grid.Divide(len(VIEWS), 2)
  legend = ROOT.TLegend(LEGEND_X, 0.3, 1., TITLE_Y)
  keep = [grid, legend]
  for row, half in enumerate([zoom, wide]):
    for col, view in enumerate(VIEWS):
      grid.cd(row * len(VIEWS) + col + 1)
      ROOT.gPad.SetLeftMargin(PAD_LEFT_MARGIN)
      draw_panel(tracks, layers, view, half, colors, keep,
                 legend if row == 1 and col == 0 else None)
  canvas.cd()
  legend.Draw()
  canvas.SaveAs(str(output))


if __name__ == "__main__":
  main()
