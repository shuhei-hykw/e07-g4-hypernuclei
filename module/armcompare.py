"""Paired comparison of candidate variants (arms) on one image dataset.

  ./run.py armcompare base.csv arm1.csv [arm2.csv ...] [--boot 2000]

Inputs are the per-event tables of ./run.py unionstudy --per-event,
made on the same dataset with the production candidates (base) and with
a variant (arm). For every cell (vertex, depth convention, tolerance)
the union efficiency of the arm minus the base is reported with a 95 %
interval from a bootstrap over events (paired: the same events are
resampled for both), next to the differences of A alone and B' alone
and the number of B'-only rescues (union found, A not; embedded). The
chance-subtracted difference (embedded minus real views, arm minus
base, per event) is reported too: an arm with more candidates also
finds more chance matches on the real views.

Decision rule (fixed with the fullscan session before any arm was
scored), on the embedded views: adopt an arm only if the union capture
at 12 um (inclusive) has a paired difference >= 0 with a 95 % lower
bound > -5 %, the B'-only rescues of that cell do not fall and the LL
event capture (inclusive, 12 um) does not fall. Added on the dev set
before evalV was opened, after alen was seen to raise the chance rate
of A on the real views: the chance-subtracted difference of that cell
must also be >= 0 with a 95 % lower bound > -5 %.
"""

import argparse
import csv
from pathlib import Path

import numpy as np

DEFAULT_BOOT = 2000
BOOT_SEED = 1
CI_PCT = (2.5, 97.5)
PRIMARY = ("capture", "inclusive", "12")
LL_EVENT = ("LL event", "inclusive", "12")
MIN_LOWER = -0.05


def load(path):
  """(vertex, depth, tol) -> {event: row} with integer flags."""
  out = {}
  with open(path, newline="") as f:
    for r in csv.DictReader(f):
      key = (r["vertex"], r["depth"], r["tol_um"])
      out.setdefault(key, {})[r["event"]] = {
        k: int(v) for k, v in r.items()
        if k.startswith(("embedded_", "real_"))}
  return out


def flags(cell, events, col):
  return np.array([cell[e][col] for e in events], dtype=float)


def rescues(cell, events):
  return int(np.sum(flags(cell, events, "embedded_union")
                    > flags(cell, events, "embedded_A")))


def paired(base, arm, events, col, n_boot, rng):
  d = flags(arm, events, col) - flags(base, events, col)
  boot = [d[rng.integers(0, len(d), len(d))].mean() for _ in range(n_boot)]
  return d.mean(), np.percentile(boot, CI_PCT)


def adjusted(base, arm, events, n_boot, rng):
  """Paired difference of (embedded - real) union, arm minus base."""
  def net(cell):
    return (flags(cell, events, "embedded_union")
            - flags(cell, events, "real_union"))
  d = net(arm) - net(base)
  boot = [d[rng.integers(0, len(d), len(d))].mean() for _ in range(n_boot)]
  return d.mean(), np.percentile(boot, CI_PCT)


def main():
  parser = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument("base")
  parser.add_argument("arms", nargs="+")
  parser.add_argument("--boot", type=int, default=DEFAULT_BOOT)
  args = parser.parse_args()

  base = load(args.base)
  for path in args.arms:
    arm = load(path)
    rng = np.random.default_rng(BOOT_SEED)
    print(f"## {Path(path).name} - {Path(args.base).name}")
    print("cell: union base -> arm, paired difference [95 %]; A and B' "
          "alone differences; B'-only rescues base -> arm (embedded)")
    verdict = {}
    for key in sorted(base):
      events = sorted(base[key])
      if sorted(arm.get(key, {})) != events:
        raise SystemExit(f"{path}: events of {key} differ from the base")
      d, (lo, hi) = paired(base[key], arm[key], events, "embedded_union",
                           args.boot, rng)
      da = (flags(arm[key], events, "embedded_A")
            - flags(base[key], events, "embedded_A")).mean()
      db = (flags(arm[key], events, "embedded_B'")
            - flags(base[key], events, "embedded_B'")).mean()
      rb, ra = rescues(base[key], events), rescues(arm[key], events)
      ub = flags(base[key], events, "embedded_union").mean()
      ua = flags(arm[key], events, "embedded_union").mean()
      print(f"{' '.join(key):<26} {ub:>4.0%} -> {ua:>4.0%} {d:+.0%} "
            f"[{lo:+.0%}, {hi:+.0%}]  A {da:+.0%}  B' {db:+.0%}  "
            f"rescues {rb} -> {ra}", flush=True)
      dn, (lon, hin) = adjusted(base[key], arm[key], events, args.boot,
                                rng)
      print(f"{'':<26} chance-subtracted {dn:+.0%} [{lon:+.0%}, "
            f"{hin:+.0%}]", flush=True)
      verdict[key] = (d, lo, rb, ra, dn, lon)
    d, lo, rb, ra, dn, lon = verdict[PRIMARY]
    ll = verdict[LL_EVENT][0]
    ok = d >= 0. and lo > MIN_LOWER and ra >= rb and ll >= 0. \
      and dn >= 0. and lon > MIN_LOWER
    print(f"decision: {'ADOPT' if ok else 'reject'} (capture 12 um "
          f"{d:+.0%}, lower {lo:+.0%}; chance-subtracted {dn:+.0%}, "
          f"lower {lon:+.0%}; rescues {rb} -> {ra}; LL event "
          f"{ll:+.0%})\n")


if __name__ == "__main__":
  main()
