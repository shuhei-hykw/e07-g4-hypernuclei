"""Gain of the union of method A and darkness-ranked B' over A alone.

  ./run.py unionstudy data/images/<dataset> --a A.csv --b B.csv
    --a-bg A_bg.csv --b-bg B_bg.csv [--top 770] [--boot 2000]

The union keeps every method-A candidate plus the top K B' points per
view by image darkness (module/rankstudy.py; ties broken at random with
a fixed seed). Recovery is scored on the ungrouped union; for the
inspection workload a B' point is grouped with the nearest method-A
candidate (3-D distance, then the first A in file order) within
GROUP_XY_UM in XY and GROUP_Z_UM of its z_best (developed space), and
the grouped candidates per view and the share of grouped B' points are
reported. Every event of the dataset is kept. Next to each candidate
file <name>_methodA_production.csv / <name>_bprime.csv the fullscan run
writes <name>_events.csv with one row per finished event: every event
of the dataset must be listed there (otherwise the study stops), and a
listed event without candidate rows is an empty view. Without that file
(older runs) events without rows are counted as empty views and their
number is printed. The
same selection is made on the real views before
embedding (--a-bg / --b-bg, darkness from the real view) and matched to
the same true positions. Reported per vertex (capture, LL decay) and
per LL event (both found), at 6 / 12 um, for the depth conventions
inclusive and edge: efficiency of A, of B' alone and of the union on the
embedded and on the real views, the number of A misses rescued by the
union and of A finds lost by it (greedy one-to-one matching is not
monotonic), and the gain
  G = (union - A)_embedded - (union - A)_real
with a 95 % interval from a paired bootstrap over events.
"""

import argparse
import csv
from pathlib import Path

import numpy as np

from . import efficiency as eff
from .rankstudy import darkness_scores, top_k

DEFAULT_TOP = 770   # candidates per view of method A
DEFAULT_BOOT = 2000
BOOT_SEED = 1
CI_PCT = (2.5, 97.5)
TOLS = [6., 12.]
METRICS = ("inclusive", "edge")
SETS = ("A", "B'", "union")
NO_SHIFTS = []
GROUP_XY_UM = 2.
GROUP_Z_UM = 2.
RUN_SUFFIXES = ("_methodA_production.csv", "_bprime.csv")


def check_finished(path, events):
  """Stop if the run of a candidate file has not finished every event;
  False if the run left no per-event record."""
  p = str(path)
  for suffix in RUN_SUFFIXES:
    if p.endswith(suffix):
      ev_file = Path(p[:-len(suffix)] + "_events.csv")
      break
  else:
    return False
  if not ev_file.exists():
    return False
  with open(ev_file, newline="") as f:
    done = {r["event"] for r in csv.DictReader(f)}
  missing = sorted(set(events) - done)
  if missing:
    raise SystemExit(f"{ev_file}: {len(missing)} events not finished, "
                     f"e.g. {missing[:5]}")
  return True


def stage_to_dev_um(truth, z_mm):
  return (z_mm - truth.get("stage_z0_mm", 0.)) / truth["stage_dz_mm"] \
    * truth["physical_dz_um"]


def n_grouped(truth, a, b):
  """Number of B' points that join a method-A candidate."""
  if not a or not b:
    return 0
  pix = truth["pixel_um"]
  axy = np.array([[c["col"], c["row"]] for c in a]) * pix
  az = stage_to_dev_um(truth, np.array(
    [c.get("z_best", 0.5 * sum(c["z"])) for c in a]))
  bxy = np.array([[c["col"], c["row"]] for c in b]) * pix
  bz = stage_to_dev_um(truth, np.array([c["z"][0] for c in b]))
  dxy = np.hypot(*(bxy[:, None, :] - axy[None, :, :]).transpose(2, 0, 1))
  dz = np.abs(bz[:, None] - az[None, :])
  return int(((dxy <= GROUP_XY_UM) & (dz <= GROUP_Z_UM)).any(axis=1).sum())


def bprime_top(events, cands, k, image_of):
  out = {}
  for name, truth in events.items():
    cs = cands.get(name, [])
    if not cs:
      out[name] = []
      continue
    dark = darkness_scores(image_of(name, truth), truth, cs)
    out[name] = [cs[i] for i in top_k(dark, k)]
  return out


def found(events, cands):
  """found[(vertex, metric, tol)][event] -> bool for the LL vertices."""
  rows, _ = eff.evaluate(events, cands, TOLS, NO_SHIFTS, "developed")
  out = {}
  for r in rows:
    if r["type"] == "capture":
      vtype = "capture"
    elif r["type"] == "HypernucleusDecay" \
        and r["parent"].startswith(eff.LL_PREFIX):
      vtype = "LL decay"
    else:
      continue
    for m in METRICS:
      for t in TOLS:
        out.setdefault((vtype, m, t), {})[r["event"]] = \
          r[f"found_{m}_{t:g}"]
  for m in METRICS:
    for t in TOLS:
      cap, dec = out[("capture", m, t)], out[("LL decay", m, t)]
      out[("LL event", m, t)] = {e: cap[e] and dec[e]
                                 for e in cap if e in dec}
  return out


def gain_ci(a_e, u_e, a_r, u_r, n_boot, rng):
  n = len(a_e)
  diff = (u_e - a_e) - (u_r - a_r)
  boot = [diff[rng.integers(0, n, n)].mean() for _ in range(n_boot)]
  return diff.mean(), np.percentile(boot, CI_PCT)


def main():
  parser = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument("dataset")
  parser.add_argument("--a", required=True, help="method A candidates")
  parser.add_argument("--b", required=True, help="B' branch points")
  parser.add_argument("--a-bg", required=True,
                      help="method A on the real views before embedding")
  parser.add_argument("--b-bg", required=True,
                      help="B' on the real views before embedding")
  parser.add_argument("--top", type=int, default=DEFAULT_TOP)
  parser.add_argument("--boot", type=int, default=DEFAULT_BOOT)
  args = parser.parse_args()

  events = eff.load_truth(args.dataset)
  inputs = [eff.load_candidates(p) for p in
            (args.a, args.b, args.a_bg, args.b_bg)]
  paths = (args.a, args.b, args.a_bg, args.b_bg)
  for label, c, p in zip(("A", "B'", "A real", "B' real"), inputs, paths):
    n_empty = sum(k not in c for k in events)
    record = "run record checked" if check_finished(p, events) \
      else "no run record"
    print(f"{label}: {n_empty} events without candidates ({record})")
  a, b, a_bg, b_bg = [{k: c.get(k, []) for k in events} for c in inputs]
  b = bprime_top(events, b, args.top,
                 lambda n, t: Path(args.dataset) / n / "image.json")
  b_bg = bprime_top(events, b_bg, args.top,
                    lambda n, t: t["background_view"])
  res = {}
  for label, ca, cb in (("embedded", a, b), ("real", a_bg, b_bg)):
    union = {n: ca[n] + cb[n] for n in events}
    for s, c in zip(SETS, (ca, cb, union)):
      res[(label, s)] = found(events, c)
  n_a = np.mean([len(a[n]) for n in events])
  n_b = np.mean([len(b[n]) for n in events])
  n_g = np.mean([n_grouped(t, a[n], b[n]) for n, t in events.items()])
  print(f"{len(events)} events; candidates per view: A {n_a:.0f}, "
        f"B' top {args.top} ({n_b:.0f}), union {n_a + n_b:.0f}, grouped "
        f"{n_a + n_b - n_g:.0f} ({n_g / max(n_b, 1.):.0%} of the B' "
        "points join an A candidate)")
  print("cells: embedded (real views alone); rescued / lost = union "
        "found and A not / the reverse (embedded); G = gain of the union over A minus the same on "
        "the real views, 95 % bootstrap interval")
  rng = np.random.default_rng(BOOT_SEED)
  head = f"{'vertex':<9} {'depth':<9} {'tol':>4}" + "".join(
    f"{s:>12}" for s in SETS) + f"{'resc/lost':>10}  G"
  print(head)
  for vtype in ("capture", "LL decay", "LL event"):
    for m in METRICS:
      for t in TOLS:
        key = (vtype, m, t)
        names = sorted(res[("embedded", "A")][key])
        arr = {(lab, s): np.array([res[(lab, s)][key][n] for n in names],
                                  dtype=float)
               for lab in ("embedded", "real") for s in SETS}
        cells = "".join(
          f"{arr[('embedded', s)].mean():>6.0%}"
          f"({arr[('real', s)].mean():>3.0%})" for s in SETS)
        ua, aa = arr[("embedded", "union")], arr[("embedded", "A")]
        resc = f"{int(np.sum((ua > 0) & (aa == 0)))}/" \
          f"{int(np.sum((ua == 0) & (aa > 0)))}"
        g, (lo, hi) = gain_ci(arr[("embedded", "A")],
                              arr[("embedded", "union")],
                              arr[("real", "A")], arr[("real", "union")],
                              args.boot, rng)
        print(f"{vtype:<9} {m:<9} {t:>4g}{cells}{resc:>10}  "
              f"{g:+.0%} [{lo:+.0%}, {hi:+.0%}]", flush=True)


if __name__ == "__main__":
  main()
