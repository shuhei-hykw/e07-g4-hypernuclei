"""Enumerate energetically allowed final states of Xi- capture at rest.

  Xi- + A(Z) -> sum of nuclei / Lambda hypernuclei / free Lambdas
Conserved: baryon number A+1, charge Z-1, strangeness -2 (two Lambdas).
"""

import argparse
import csv
import math
import os
import sys

import species as sp

# Xi- binding energy in the atomic 3D state (MeV), as used in E07 analyses
TARGETS = {
  "12C": (12, 6, 0.126),
  "14N": (14, 7, 0.174),
  "16O": (16, 8, 0.231),
}

MAX_BODIES = 6
EPS_MEV = 1e-9


def partitions(pool, a, z, n_lam, mass_max, max_bodies, start=0, acc=()):
  """Yield multisets of species summing to (a, z, n_lam) below mass_max."""
  if a == 0:
    if z == 0 and n_lam == 0:
      yield acc
    return
  if len(acc) >= max_bodies:
    return
  used = sum(s.mass for s in acc)
  for i in range(start, len(pool)):
    s = pool[i]
    if s.a > a or s.z > z or s.n_lambda > n_lam:
      continue
    # Lower bound: remaining baryons at the lightest mass per baryon
    if used + s.mass + (a - s.a) * MIN_MASS_PER_A > mass_max:
      continue
    yield from partitions(pool, a - s.a, z - s.z, n_lam - s.n_lambda,
                          mass_max, max_bodies, i, acc + (s,))


def stable_pool(use_estimates):
  """Keep hypernuclei that are stable against strong particle emission."""
  global MIN_MASS_PER_A
  nuclei, lam, singles, doubles = sp.build_candidates(use_estimates)
  pool = nuclei + [lam]
  MIN_MASS_PER_A = min(s.mass / s.a for s in pool + singles + doubles)
  for group in (singles, doubles):
    for s in sorted(group, key=lambda x: x.a):
      breakups = partitions(pool, s.a, s.z, s.n_lambda, s.mass,
                            s.a, 0, ())
      # Any breakup with >= 2 bodies lighter than s means s is unbound
      if not any(len(b) >= 2 for b in breakups):
        pool.append(s)
  pool.sort(key=lambda x: (x.n_lambda, x.a, x.z), reverse=True)
  return pool


def classify(final):
  n_ll = sum(1 for s in final if s.n_lambda == 2)
  n_single = sum(1 for s in final if s.n_lambda == 1 and s.a > 1)
  if n_ll:
    return "double"
  if n_single == 2:
    return "twin"
  if n_single == 1:
    return "single+L"
  return "free-LL"


def two_body_momentum(m0, m1, m2):
  s = m0 * m0
  arg = (s - (m1 + m2) ** 2) * (s - (m1 - m2) ** 2)
  return math.sqrt(max(arg, 0.0)) / (2.0 * m0)


def enumerate_target(key, pool, max_bodies):
  a, z, b_xi = TARGETS[key]
  m_init = sp.nuclear_mass(a, z) + sp.M_XI_MINUS - b_xi
  rows = []
  for final in partitions(pool, a + 1, z - 1, 2, m_init, max_bodies):
    q = m_init - sum(s.mass for s in final)
    if q < EPS_MEV:
      continue
    names = [s.name for s in sorted(
      final, key=lambda s: (-s.n_lambda, -s.a, -s.z))]
    # A single-body final state de-excites only by gamma emission
    if len(final) == 1:
      names.append("gamma")
    p2 = ""
    if len(final) == 2:
      p2 = f"{two_body_momentum(m_init, final[0].mass, final[1].mass):.2f}"
    rows.append({
      "target": key,
      "category": classify(final),
      "channel": " + ".join(names),
      "q_mev": f"{q:.3f}",
      "n_bodies": len(final),
      "n_charged": sum(1 for s in final if s.z > 0),
      "n_neutral": sum(1 for s in final if s.z == 0),
      "p_2body_mev": p2,
      "uses_estimate": int(any(s.estimated for s in final)),
    })
  order = {"double": 0, "twin": 1, "single+L": 2, "free-LL": 3}
  rows.sort(key=lambda r: (order[r["category"]], -float(r["q_mev"])))
  return rows


def main():
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument("--measured-only", action="store_true",
                      help="drop single-L hypernuclei with estimated B_L")
  parser.add_argument("--max-bodies", type=int, default=MAX_BODIES)
  parser.add_argument("--outdir", default="output")
  args = parser.parse_args()

  pool = stable_pool(not args.measured_only)
  os.makedirs(args.outdir, exist_ok=True)

  with open(os.path.join(args.outdir, "species.csv"), "w",
            newline="") as f:
    w = csv.writer(f)
    w.writerow(["name", "A", "Z", "n_lambda", "mass_mev", "estimated"])
    for s in sorted(pool, key=lambda x: (x.n_lambda, x.z, x.a)):
      w.writerow([s.name, s.a, s.z, s.n_lambda, f"{s.mass:.3f}",
                  int(s.estimated)])

  for key in TARGETS:
    rows = enumerate_target(key, pool, args.max_bodies)
    path = os.path.join(args.outdir, f"channels_{key}.csv")
    with open(path, "w", newline="") as f:
      w = csv.DictWriter(f, fieldnames=list(rows[0]))
      w.writeheader()
      w.writerows(rows)
    counts = {}
    for r in rows:
      counts[r["category"]] = counts.get(r["category"], 0) + 1
    print(f"{key}: {len(rows)} channels {counts} -> {path}",
          file=sys.stderr)


if __name__ == "__main__":
  main()
