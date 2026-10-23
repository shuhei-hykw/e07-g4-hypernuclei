"""Fixed kinematics of observed events (Nagara, Ibuki) for emulsim.

  ./run.py replay            # writes output/replay/{nagara,ibuki}.csv

Each charged track is given by its measured range, theta and phi
(theta from the z axis perpendicular to the sheets, phi in the sheet
plane; kinema convention) and its identification. The kinetic energy is
obtained from the kinema range-energy relation at the measured emulsion
density. A track flagged `from_q` takes its energy from energy
conservation instead (used for long pi- tracks, for which kinema's
relation is unreliable). Missing neutrons take the missing momentum
(one neutron) or share the missing four-momentum isotropically (two).

Sources (kinema): Nagara from runcard.json (NagaraA/B/C); Ibuki from
input/twin_mod047pl10/e07twin_mod047pl10_{A,B,C}.json. Identifications:
Nagara A, B and Ibuki A, C as published / kinema; Nagara C
(L5He -> d + p + 2n) is the only channel with two charged tracks that
balances with missing neutrons; Ibuki B (L10Be -> d + t + t + p + n) is
the kinema candidate with the smallest mass gap and is not unique.

Output CSV rows: vertex (capture or the decaying hypernucleus), name,
px, py, pz (MeV/c). /xi/replayFile in emulsim reproduces them.
"""

import argparse
import csv
import math
import sys

import numpy as np

from . import config
from .check_range import load_range_energy

EVENTS = {
  "nagara": {
    "target": "12C", "density": 3.629,
    "vertices": {
      "capture": [("LL6He", 8.1, 44.9, 337.5), ("4He", 3.2, 57.7, 174.9),
                  ("t", 88.3, 156.2, 143.0)],
      "LL6He": [("L5He", 9.1, 77.7, 115.9), ("p", 82.1, 122.8, 284.2),
                ("pi-", 13697., 81.1, 305.5, "from_q")],
      "L5He": [("d", 742.6, 138.5, 322.1), ("p", 5868., 52.2, 123.7),
               ("n", "missing"), ("n", "missing")],
    },
  },
  "ibuki": {
    "target": "14N", "density": 3.605,
    "vertices": {
      "capture": [("L10Be", 9.1696, 84.913, 298.521),
                  ("L5He", 89.44, 87.332, 121.576)],
      # kinema assignment (made on an earlier measurement of the same
      # tracks): d = short track, t = 25 um, t = 0.94 mm, p = 20 mm
      "L10Be": [("d", 4.148, 29.798, -50.906), ("t", 24.65, 109.6, 172.21),
                ("t", 935.9, 46.74, 189.14), ("p", 20195.9, 110.39, 325.86),
                ("n", "missing")],
      "L5He": [("4He", 18.67, 84.63, 307.53),
               ("pi-", 2007.9, 29.41, 234.61), ("p", 2134.1, 108.01, 120.45)],
    },
  },
}
MAX_MISSING = 2


def load_masses():
  table = {}
  with open(config.OUTPUT_DIR / "species.csv") as f:
    for r in csv.DictReader(f):
      table[r["name"]] = (float(r["mass_mev"]), int(r["Z"]))
  table["pi-"] = (139.57039, -1)
  with open(config.OUTPUT_DIR / "targets.csv") as f:
    initial = {r["name"]: float(r["m_init_mev"]) for r in csv.DictReader(f)}
  return table, initial


def direction(theta, phi):
  t, p = math.radians(theta), math.radians(phi)
  return np.array([math.sin(t) * math.cos(p), math.sin(t) * math.sin(p),
                   math.cos(t)])


def momentum_from_ke(ke, mass):
  return math.sqrt(ke * ke + 2. * ke * mass)


def two_body(m_sys, p_sys, e_sys, m1, m2, rng):
  """Isotropic decay of a system (mass, momentum, energy) into two."""
  pstar = math.sqrt(max((m_sys ** 2 - (m1 + m2) ** 2)
                        * (m_sys ** 2 - (m1 - m2) ** 2), 0.)) / (2. * m_sys)
  d = rng.normal(size=3)
  d /= np.linalg.norm(d)
  out = []
  beta = p_sys / e_sys
  gamma = e_sys / m_sys
  for sign, m in ((1., m1), (-1., m2)):
    p = sign * pstar * d
    e = math.sqrt(pstar ** 2 + m * m)
    bp = beta @ p
    b2 = beta @ beta
    p_lab = p + ((gamma - 1.) * bp / b2 if b2 > 0. else 0.) * beta \
      + gamma * e * beta
    out.append(p_lab)
  return out


def vertex_kinematics(parent_mass, products, masses, re, density, rng):
  """Momentum vectors of all products of one vertex."""
  p_sum = np.zeros(3)
  e_sum = 0.
  rows, missing, from_q = [], [], None
  for prod in products:
    name = prod[0]
    m, z = masses[name]
    if prod[1] == "missing":
      missing.append((name, m))
      continue
    rng_um, theta, phi = prod[1:4]
    d = direction(theta, phi)
    if len(prod) > 4 and prod[4] == "from_q":
      from_q = (name, m, d)
      continue
    ke = float(re.KEfromRange(m, rng_um, abs(z), density))
    p = momentum_from_ke(ke, m) * d
    rows.append((name, p, ke))
    p_sum += p
    e_sum += ke + m
  if from_q:
    name, m, d = from_q
    ke = parent_mass - e_sum - m - sum(mm for _, mm in missing)
    p = momentum_from_ke(ke, m) * d
    rows.append((name, p, ke))
    p_sum += p
    e_sum += ke + m
  if len(missing) == 1:
    name, m = missing[0]
    p = -p_sum
    rows.append((name, p, math.hypot(np.linalg.norm(p), m) - m))
  elif len(missing) == MAX_MISSING:
    e_miss = parent_mass - e_sum
    p_miss = -p_sum
    m_miss = math.sqrt(max(e_miss ** 2 - p_miss @ p_miss, 0.))
    (n1, m1), (n2, m2) = missing
    if m_miss < m1 + m2:
      raise ValueError("missing mass below the neutral masses")
    for (name, m), p in zip(missing, two_body(m_miss, p_miss, e_miss,
                                              m1, m2, rng)):
      rows.append((name, p, math.hypot(np.linalg.norm(p), m) - m))
  elif missing:
    raise ValueError("at most two missing particles are supported")
  return rows


def main():
  parser = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument("--events", default=",".join(EVENTS))
  parser.add_argument("--seed", type=int, default=0,
                      help="for sharing missing momentum between neutrons")
  args = parser.parse_args()

  masses, initial = load_masses()
  re = load_range_energy(str(config.KINEMA_DIR))
  rng = np.random.default_rng(args.seed)
  out_dir = config.OUTPUT_DIR / "replay"
  out_dir.mkdir(parents=True, exist_ok=True)
  for name in args.events.split(","):
    ev = EVENTS[name]
    path = out_dir / f"{name}.csv"
    with open(path, "w", newline="") as f:
      w = csv.writer(f)
      w.writerow(["vertex", "name", "px_mev", "py_mev", "pz_mev"])
      print(f"== {name} (Xi- + {ev['target']}, density {ev['density']})",
            file=sys.stderr)
      for vertex, products in ev["vertices"].items():
        m_parent = (initial[ev["target"]] if vertex == "capture"
                    else masses[vertex][0])
        rows = vertex_kinematics(m_parent, products, masses, re,
                                 ev["density"], rng)
        p_tot = sum(p for _, p, _ in rows)
        ke_tot = sum(ke for _, _, ke in rows)
        q = m_parent - sum(masses[n][0] for n, _, _ in rows)
        print(f"  {vertex:8s} -> {', '.join(n for n, _, _ in rows):28s}"
              f" Q {q:7.3f}  sum KE {ke_tot:7.3f}  |sum p| "
              f"{np.linalg.norm(p_tot):6.2f} MeV/c", file=sys.stderr)
        for n, p, ke in rows:
          w.writerow([vertex, n] + [f"{c:.4f}" for c in p])
    print(f"  -> {path}", file=sys.stderr)


if __name__ == "__main__":
  main()
