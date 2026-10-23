"""Enumerate weak-decay channels of single- and double-Lambda hypernuclei.

One Lambda decays at a time:
  pi-  : Lambda -> p pi-          (A, Z+1, nL-1) + pi-
  pi0  : Lambda -> n pi0          (A, Z,   nL-1) + pi0
  nm   : Lambda N -> N N (or Lambda Lambda -> Lambda N), no pion
"""

import argparse
import csv
import os
import sys

import enumerate_channels as ec
import species as sp

MAX_DECAY_BODIES = 4   # including the pion
MODES = [("pi-", 1, sp.M_PI_MINUS), ("pi0", 0, sp.M_PI_ZERO),
         ("nm", 0, 0.)]


def enumerate_decays(parent, pool, max_bodies):
  rows = []
  for mode, dz, m_pion in MODES:
    has_pion = m_pion > 0.
    n_max = max_bodies - (1 if has_pion else 0)
    mass_max = parent.mass - m_pion
    for final in ec.partitions(pool, parent.a, parent.z + dz,
                               parent.n_lambda - 1, mass_max, n_max):
      # A decay into a single product is forbidden for no-pion modes
      if not has_pion and len(final) < 2:
        continue
      q = mass_max - sum(s.mass for s in final)
      if q < ec.EPS_MEV:
        continue
      names = [s.name for s in sorted(
        final, key=lambda s: (-s.n_lambda, -s.a, -s.z))]
      if has_pion:
        names.append(mode)
      rows.append({
        "parent": parent.name,
        "mode": mode,
        "channel": ",".join(names),
        "q_mev": f"{q:.3f}",
        "n_bodies": len(names),
        "n_charged": sum(1 for s in final if s.z > 0) + (mode == "pi-"),
        "uses_estimate": int(parent.estimated
                             or any(s.estimated for s in final)),
      })
  return rows


def main():
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument("--measured-only", action="store_true")
  parser.add_argument("--max-bodies", type=int, default=MAX_DECAY_BODIES)
  parser.add_argument("--outdir", default="output")
  args = parser.parse_args()

  pool = ec.stable_pool(not args.measured_only)
  parents = sorted((s for s in pool if s.n_lambda > 0 and s.a > 1),
                   key=lambda s: (s.n_lambda, s.z, s.a))
  rows = []
  for parent in parents:
    rows += enumerate_decays(parent, pool, args.max_bodies)
  os.makedirs(args.outdir, exist_ok=True)
  path = os.path.join(args.outdir, "decays.csv")
  with open(path, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)
  print(f"{len(parents)} hypernuclei, {len(rows)} decay channels -> "
        f"{path}", file=sys.stderr)


if __name__ == "__main__":
  main()
