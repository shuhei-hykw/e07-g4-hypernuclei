"""Containment summary of capture and hypernuclear decay products.

  python3 summary.py emulsim.root -o summary.png

Left: capture-point depth for all events and for events whose charged
capture / decay products all stop inside the module (base films and the
thin sheets marked). Right: escape fraction per particle and origin.
"""

import argparse
import collections

import ROOT

CANVAS_W = 1400
CANVAS_H = 600
Z_BIN_UM = 200.
DECAY_PROCESS = "HypernucleusDecay"
BASE_COLOR = ROOT.kGray
FILL_ALL = ROOT.kAzure - 9
FILL_CONTAINED = ROOT.kAzure + 2
BAR_COLOR = ROOT.kOrange + 7
BAR_WIDTH = 0.6
BAR_OFFSET = 0.2
LABEL_SIZE = 0.045
LEFT_MARGIN = 0.38
Y_MARGIN = 1.15


def main():
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument("root_file")
  parser.add_argument("-o", "--output", default="summary.png")
  args = parser.parse_args()

  ROOT.gROOT.SetBatch(True)
  ROOT.gStyle.SetOptStat(0)
  f = ROOT.TFile.Open(args.root_file)
  layers = [(bool(l.emulsion), bool(l.thick), l.zmin, l.zmax)
            for l in f.Get("layers")]
  z_lo, z_hi = layers[0][2], layers[-1][3]
  n_bins = int((z_hi - z_lo) / Z_BIN_UM)
  h_all = ROOT.TH1D("all", "", n_bins, z_lo, z_hi)
  h_cont = ROOT.TH1D("contained", "", n_bins, z_lo, z_hi)

  total = collections.Counter()
  escaped = collections.Counter()
  tree = f.Get("tree")
  channel = ""
  for ev in tree:
    channel = str(ev.channel)
    h_all.Fill(ev.vz)
    if ev.contained:
      h_cont.Fill(ev.vz)
    for j in range(len(ev.name)):
      if ev.z[j] == 0:
        continue
      if ev.parent[j] == 0:
        origin = "capture"
      elif str(ev.process[j]) == DECAY_PROCESS:
        origin = "decay"
      else:
        continue
      key = f"{ev.name[j]} ({origin})"
      total[key] += 1
      escaped[key] += ev.escaped[j]

  canvas = ROOT.TCanvas("c", "", CANVAS_W, CANVAS_H)
  canvas.Divide(2, 1)
  canvas.cd(1)
  n_ev = tree.GetEntries()
  n_cont = int(h_cont.GetEntries())
  h_all.SetTitle(f"{channel}: {n_cont}/{n_ev} events contained;"
                 "capture point z (#mum);events")
  h_all.SetFillColor(FILL_ALL)
  h_all.SetMaximum(h_all.GetMaximum() * Y_MARGIN)
  h_all.Draw("hist")
  y_max = h_all.GetMaximum()
  keep = []
  for emulsion, thick, z0, z1 in layers:
    if not emulsion or not thick:
      b = ROOT.TBox(z0, 0., z1, y_max)
      b.SetFillColor(BASE_COLOR)
      b.Draw()
      keep.append(b)
  h_all.Draw("hist same")
  h_cont.SetFillColor(FILL_CONTAINED)
  h_cont.Draw("hist same")
  leg = ROOT.TLegend(0.55, 0.80, 0.9, 0.9)
  leg.AddEntry(h_all, "all", "f")
  leg.AddEntry(h_cont, "all charged tracks contained", "f")
  leg.Draw()
  ROOT.gPad.RedrawAxis()

  canvas.cd(2)
  ROOT.gPad.SetLeftMargin(LEFT_MARGIN)
  keys = sorted(total)
  h_esc = ROOT.TH1D("esc", "escape fraction of charged tracks;;",
                    len(keys), 0, len(keys))
  for i, k in enumerate(keys, start=1):
    h_esc.SetBinContent(i, escaped[k] / total[k])
    h_esc.GetXaxis().SetBinLabel(i, f"{k} [{escaped[k]}/{total[k]}]")
  h_esc.SetMaximum(1.)
  h_esc.SetFillColor(BAR_COLOR)
  h_esc.SetBarWidth(BAR_WIDTH)
  h_esc.SetBarOffset(BAR_OFFSET)
  h_esc.GetXaxis().SetLabelSize(LABEL_SIZE)
  h_esc.Draw("hbar")
  canvas.SaveAs(args.output)


if __name__ == "__main__":
  main()
