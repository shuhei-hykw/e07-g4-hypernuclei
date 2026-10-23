"""Draw XY / XZ projections of the primary tracks of one event.

  python3 draw_event.py emulsim.root --event 0 -o event0.png
"""

import argparse

import ROOT

CANVAS_W = 1200
CANVAS_H = 600
FRAME_MARGIN = 1.1
LINE_WIDTH = 2
LEGEND_POS = (0.62, 0.70, 0.90, 0.90)
COLORS = [ROOT.kRed + 1, ROOT.kBlue + 1, ROOT.kGreen + 2,
          ROOT.kMagenta + 1, ROOT.kOrange + 7, ROOT.kCyan + 2]
NEUTRAL_STYLE = 2  # dashed


def main():
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument("root_file")
  parser.add_argument("--event", type=int, default=0)
  parser.add_argument("--max-range", type=float, default=0.,
                      help="half width of frame in um (0: auto)")
  parser.add_argument("-o", "--output", default="event.png")
  args = parser.parse_args()

  ROOT.gROOT.SetBatch(True)
  ROOT.gStyle.SetOptStat(0)
  f = ROOT.TFile.Open(args.root_file)
  tree = f.Get("tree")
  tree.GetEntry(args.event)

  tracks = []
  for j in range(len(tree.name)):
    tracks.append((str(tree.name[j]), int(tree.z[j]), list(tree.x[j]),
                   list(tree.y[j]), list(tree.zpos[j])))
  half = args.max_range or FRAME_MARGIN * max(
    max(abs(v) for v in xs + ys + zs) for _, _, xs, ys, zs in tracks)

  canvas = ROOT.TCanvas("c", "", CANVAS_W, CANVAS_H)
  canvas.Divide(2, 1)
  keep = []
  for pad, (label, iv) in enumerate([("y", 3), ("z", 4)], start=1):
    canvas.cd(pad)
    frame = ROOT.gPad.DrawFrame(-half, -half, half, half)
    frame.SetTitle(f"{tree.channel}  (event {args.event});x (#mum);"
                   f"{label} (#mum)")
    legend = ROOT.TLegend(*LEGEND_POS)
    for k, trk in enumerate(tracks):
      name, z, xs = trk[0], trk[1], trk[2]
      graph = ROOT.TGraph(len(xs))
      for n, (xv, yv) in enumerate(zip(xs, trk[iv])):
        graph.SetPoint(n, xv, yv)
      graph.SetLineColor(COLORS[k % len(COLORS)])
      graph.SetLineWidth(LINE_WIDTH)
      if z == 0:
        graph.SetLineStyle(NEUTRAL_STYLE)
      graph.Draw("L")
      legend.AddEntry(graph, name, "l")
      keep.append(graph)
    legend.Draw()
    keep.append(legend)
  canvas.SaveAs(args.output)


if __name__ == "__main__":
  main()
