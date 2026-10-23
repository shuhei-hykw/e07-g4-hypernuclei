"""Draw XY / XZ projections of the tracks of one event.

  python3 draw_event.py emulsim.root --event 0 --zoom 30 -o event0.png

Primaries are drawn with thick lines, hypernuclear decay products with
medium lines and other secondaries (pi- capture, recoils) thin. Neutral
particles are dashed and hidden unless --neutral is given.
"""

import argparse

import ROOT

CANVAS_W = 1200
CANVAS_H = 600
FRAME_MARGIN = 1.1
WIDTH_PRIMARY = 3
WIDTH_DECAY = 2
WIDTH_OTHER = 1
LEGEND_POS = (0.60, 0.60, 0.90, 0.90)
NEUTRAL_STYLE = 2  # dashed
DECAY_PROCESS = "HypernucleusDecay"
COLORS = [ROOT.kRed + 1, ROOT.kBlue + 1, ROOT.kGreen + 2,
          ROOT.kMagenta + 1, ROOT.kOrange + 7, ROOT.kCyan + 2,
          ROOT.kViolet + 1, ROOT.kSpring - 6, ROOT.kAzure + 7,
          ROOT.kPink + 6, ROOT.kGray + 2, ROOT.kYellow + 3]


def load_tracks(tree, with_neutral):
  tracks = []
  for j in range(len(tree.name)):
    if tree.z[j] == 0 and not with_neutral:
      continue
    tracks.append({
      "name": str(tree.name[j]),
      "z": int(tree.z[j]),
      "parent": int(tree.parent[j]),
      "process": str(tree.process[j]),
      "x": list(tree.x[j]), "y": list(tree.y[j]), "zpos": list(tree.zpos[j]),
    })
  return tracks


def line_width(track):
  if track["parent"] == 0:
    return WIDTH_PRIMARY
  if track["process"] == DECAY_PROCESS:
    return WIDTH_DECAY
  return WIDTH_OTHER


def main():
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument("root_file")
  parser.add_argument("--event", type=int, default=0)
  parser.add_argument("--zoom", type=float, default=0.,
                      help="half width of frame in um (0: auto)")
  parser.add_argument("--neutral", action="store_true",
                      help="also draw neutral tracks")
  parser.add_argument("-o", "--output", default="event.png")
  args = parser.parse_args()

  ROOT.gROOT.SetBatch(True)
  ROOT.gStyle.SetOptStat(0)
  f = ROOT.TFile.Open(args.root_file)
  tree = f.Get("tree")
  tree.GetEntry(args.event)

  tracks = load_tracks(tree, args.neutral)
  half = args.zoom or FRAME_MARGIN * max(
    max(abs(v) for v in t["x"] + t["y"] + t["zpos"]) for t in tracks)

  color_of = {}
  for t in tracks:
    color_of.setdefault(t["name"], COLORS[len(color_of) % len(COLORS)])

  canvas = ROOT.TCanvas("c", "", CANVAS_W, CANVAS_H)
  canvas.Divide(2, 1)
  keep = []
  for pad, label in enumerate(["y", "zpos"], start=1):
    canvas.cd(pad)
    frame = ROOT.gPad.DrawFrame(-half, -half, half, half)
    axis = "z" if label == "zpos" else label
    frame.SetTitle(f"{tree.channel}  (event {args.event});x (#mum);"
                   f"{axis} (#mum)")
    legend = ROOT.TLegend(*LEGEND_POS)
    shown = set()
    for t in tracks:
      graph = ROOT.TGraph(len(t["x"]))
      for n, (xv, yv) in enumerate(zip(t["x"], t[label])):
        graph.SetPoint(n, xv, yv)
      graph.SetLineColor(color_of[t["name"]])
      graph.SetLineWidth(line_width(t))
      if t["z"] == 0:
        graph.SetLineStyle(NEUTRAL_STYLE)
      graph.Draw("L")
      if t["name"] not in shown:
        legend.AddEntry(graph, t["name"], "l")
        shown.add(t["name"])
      keep.append(graph)
    legend.Draw()
    keep.append(legend)
  canvas.SaveAs(args.output)


if __name__ == "__main__":
  main()
