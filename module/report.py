"""Efficiency report over datasets and full-scan methods.

  ./run.py report [--candidates DIR] [--pattern 'v3_*'] [-o report.md]

Every candidate file DIR/<dataset>_<method>_<cuts>.csv (as written by
the e07-fullscan analysis, e.g. v3_embed_nagara_dz3.0_methodA_relaxed)
whose dataset exists under data/images is evaluated with the matching
of module/efficiency.py; the report lists per dataset and method the
efficiencies of the capture vertex, the double-Lambda decay vertex and
the complete LL event at each tolerance (edge convention), the chance
level from shifted truths, duplicates and candidates per event.
"""

import argparse
import re
from pathlib import Path

import numpy as np

from . import config
from . import efficiency as eff

DEFAULT_CANDIDATES = Path.home() / "work/e07-fullscan/results/sim_eval"
NAME_RE = re.compile(r"^(?P<dataset>.+?)_(?P<method>method[^_]+|bprime)"
                     r"(?:_(?P<cuts>[^.]+))?\.csv$")


def rates(rows, extra, tolerances, metric):
  out = {}
  for label, vtype in (("capture", "capture"),
                       ("LL decay", "HypernucleusDecay")):
    sel = [r for r in eff._select(rows, vtype, False)
           if vtype == "capture" or r["parent"].startswith(eff.LL_PREFIX)]
    ctrl = [r for r in eff._select(rows, vtype, True)
            if vtype == "capture" or r["parent"].startswith(eff.LL_PREFIX)]
    out[label] = (len(sel),
                  [np.mean([r[f"found_{metric}_{t:g}"] for r in sel])
                   if sel else np.nan for t in tolerances],
                  [np.mean([r[f"found_{metric}_{t:g}"] for r in ctrl])
                   if ctrl else np.nan for t in tolerances])
  ll = eff._ll_events(rows, False, tolerances)
  llc = eff._ll_events(rows, True, tolerances)
  out["LL event"] = (len(ll),
                     [np.mean([r[f"found_{metric}_{t:g}"] for r in ll])
                      if ll else np.nan for t in tolerances],
                     [np.mean([r[f"found_{metric}_{t:g}"] for r in llc])
                      if llc else np.nan for t in tolerances])
  return out


def main():
  parser = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument("--candidates", default=str(DEFAULT_CANDIDATES))
  parser.add_argument("--pattern", default="*",
                      help="glob on candidate file names")
  parser.add_argument("--tolerances", default=eff.DEFAULT_TOLERANCES)
  parser.add_argument("--min-score", type=float,
                      help="keep candidates with n_tracks_max >= this")
  parser.add_argument("--metric", default="edge",
                      choices=["edge", "inclusive", "best"],
                      help="depth convention (best needs z_best)")
  parser.add_argument("-o", "--output",
                      default=str(config.OUTPUT_DIR / "efficiency_report.md"))
  args = parser.parse_args()

  tolerances = [float(t) for t in args.tolerances.split(",")]
  shifts = eff.random_shifts(eff.DEFAULT_N_SHIFTS, eff.DEFAULT_SHIFT_MIN_PX,
                             eff.DEFAULT_SHIFT_MAX_PX,
                             np.random.default_rng(0))
  head = " | ".join(f"<{t:g} um (ctrl)" for t in tolerances)
  lines = ["# Vertex-finding efficiency on simulated events", "",
           f"Depth convention: {args.metric}; one-to-one greedy matching "
           "in the developed "
           "emulsion; chance level (ctrl) from 10 random shifts of the "
           "truths (module/efficiency.py).", "",
           f"| dataset | method | cuts | vertex | n | {head} | "
           f"dup <{tolerances[-1]:g} um | cand/event |",
           "|" + "---|" * (7 + len(tolerances))]
  files = sorted(Path(args.candidates).glob(args.pattern + ".csv"))
  for path in files:
    m = NAME_RE.match(path.name)
    if not m:
      continue
    dataset = config.IMAGE_DIR / m["dataset"]
    if not dataset.exists():
      continue
    events = eff.load_truth(dataset)
    rows, extra = eff.evaluate(events,
                               eff.load_candidates(path, args.min_score),
                               tolerances, shifts, "developed")
    per_event = extra["n_candidates"] / max(len(events), 1)
    if not any(k.startswith(f"found_{args.metric}_") for k in rows[0]):
      continue
    for label, (n, vals, ctrl) in rates(rows, extra, tolerances,
                                        args.metric).items():
      if n == 0:
        continue
      lines.append(
        f"| {m['dataset']} | {m['method']} | {m['cuts'] or ''} | {label} "
        f"| {n} | " + " | ".join(f"{v:.0%} ({c:.0%})"
                                 for v, c in zip(vals, ctrl))
        + f" | {extra['duplicates'][tolerances[-1]]} | {per_event:.0f} |")
  Path(args.output).write_text("\n".join(lines) + "\n")
  print("\n".join(lines))
  print(f"-> {args.output}")


if __name__ == "__main__":
  main()
