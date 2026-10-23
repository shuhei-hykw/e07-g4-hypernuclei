"""Compare Geant4 ranges with the kinema range-energy relation.

  ./run.py check-range data/sim/nagara.root [--density 3.53]
"""

import argparse
import statistics
import sys
import types

import ROOT

from . import config

MAX_PRINT_EVENTS = 3
DEFAULT_DENSITY = 3.53
SECANT_STEP = 1e-3
SECANT_TOL = 1e-10
SECANT_MAX_ITER = 100


def fsolve(func, x0):
  """Scalar secant solver standing in for scipy.optimize.fsolve."""
  x_prev, x = x0, x0 + SECANT_STEP
  f_prev = func(x_prev)
  for _ in range(SECANT_MAX_ITER):
    f = func(x)
    if abs(f) < SECANT_TOL or f == f_prev:
      break
    x_prev, x, f_prev = x, x - f * (x - x_prev) / (f - f_prev), f
  return [x]


def load_range_energy(kinema_dir):
  # rangeenergy.py needs matplotlib only for drawing and scipy only for
  # a scalar root finder
  for mod in ("matplotlib", "matplotlib.pyplot", "scipy",
              "scipy.optimize"):
    sys.modules.setdefault(mod, types.ModuleType(mod))
  sys.modules["scipy.optimize"].fsolve = fsolve
  sys.path.insert(0, kinema_dir)
  import rangeenergy
  return rangeenergy.RangeEnergy()


def main():
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument("root_file")
  parser.add_argument("--kinema", default=str(config.KINEMA_DIR))
  parser.add_argument("--density", type=float, default=DEFAULT_DENSITY)
  args = parser.parse_args()

  re = load_range_energy(args.kinema)
  f = ROOT.TFile.Open(args.root_file)
  tree = f.Get("tree")

  ratios = {}
  for i, ev in enumerate(tree):
    # Only the products of the Xi- capture
    prim = [j for j in range(len(ev.name)) if ev.parent[j] == 0]
    names = {j: str(ev.name[j]) for j in prim}
    e_sum = sum(ev.ke[j] for j in prim)
    p_sum = [sum(v[j] for j in prim) for v in (ev.px, ev.py, ev.pz)]
    if i < MAX_PRINT_EVENTS:
      print(f"event {ev.event}: {ev.channel}  Q={ev.q:.3f}"
            f"  sum KE={e_sum:.3f}  sum p=({p_sum[0]:.2e},"
            f" {p_sum[1]:.2e}, {p_sum[2]:.2e})")
    for j, name in names.items():
      ke, z = ev.ke[j], ev.z[j]
      if z <= 0 or ke <= 0:
        continue
      p2 = ev.px[j] ** 2 + ev.py[j] ** 2 + ev.pz[j] ** 2
      mass = (p2 - ke * ke) / (2. * ke)
      r_ref = re.RangeFromKE(mass, ke, z, args.density)
      ratios.setdefault(name, []).append((ev.range[j], r_ref, ke))
      if i < MAX_PRINT_EVENTS:
        print(f"  {name:12s} KE={ke:7.3f} MeV  R_G4={ev.range[j]:8.2f} um"
              f"  R_kinema={r_ref:8.2f} um  L_G4={ev.length[j]:8.2f} um")

  print("\nmean over events: R_G4 / R_kinema")
  for name, vals in ratios.items():
    r = [g / k for g, k, _ in vals if k > 0]
    print(f"  {name:12s} n={len(r):5d}  mean={statistics.mean(r):.3f}"
          f"  rms={statistics.pstdev(r):.3f}")


if __name__ == "__main__":
  main()
