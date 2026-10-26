"""Single isolated tracks for the point-cloud evaluation.

  ./run.py singletrack --name single_eval --events 60 [--tracks 25]

For every event a reference point (the view and stack center) is put in
the middle of a random thick emulsion layer, and --tracks particles
(p, pi+, pi-, alpha, Li, C) start at random points of the view, at least
MIN_SEPARATION_UM apart. The momentum is log-uniform per species, from
nearly stopping to minimum ionizing; the dip angle (from the sheet
plane) is drawn from a mixture favouring flat and steep tracks. The gun
CSV (output/gun/<name>.csv) is simulated with emulsim (data/sim/<name>
.root) and, after the simulation, a per-track table data/sim/<name>
_tracks.csv lists event, Geant4 track id (= instance id of the labels),
species, dip, azimuth, kinetic energy, range in the module, mean
dE/dx, whether the track stops in the module, and its start point.

Render with the reference point centered, e.g.
  ./run.py images data/sim/single_eval.root --events 0-59 --centered \\
    --background random --bg-split eval --dz 4.0 --name single_eval_dz4
"""

import argparse
import csv
import subprocess
import sys
from pathlib import Path

import numpy as np

from . import config

SPECIES = {  # name: (p_min, p_max) MeV/c, log-uniform
  "p": (30., 1500.), "pi+": (10., 600.), "pi-": (10., 600.),
  "4He": (100., 3000.), "7Li": (200., 4000.), "12C": (400., 8000.),
}
DEFAULT_TRACKS = 25
DEFAULT_EVENTS = 60
VIEW_HALF_UM = 250.     # start points within the central part of a view
LAYER_MARGIN_UM = 10.
MIN_SEPARATION_UM = 60.
MAX_TRIES = 1000
# Dip mixture: flat (|dip| < FLAT_DEG), steep (|dip| > STEEP_DEG), any
FLAT_DEG, STEEP_DEG = 10., 60.
MIX = (1 / 3, 1 / 3, 1 / 3)
GEOMETRY_FILE = config.SIM_DIR / "eval_nagara.root"  # any emulsim output


def thick_layers():
  from .synth.event import EventFile
  layers = EventFile(str(GEOMETRY_FILE)).layers
  # Thick emulsion layers: the thickest emulsion layers of the module
  emu = [l for l in layers if l.emulsion]
  t = max(l.zmax - l.zmin for l in emu)
  return [l for l in emu if l.zmax - l.zmin > 0.9 * t]


def sample_dip(rng):
  u = rng.random()
  if u < MIX[0]:
    a = rng.uniform(0., FLAT_DEG)
  elif u < MIX[0] + MIX[1]:
    a = rng.uniform(STEEP_DEG, 90.)
  else:
    a = rng.uniform(0., 90.)
  return a * rng.choice([-1., 1.])


def species_mass(name):
  if name.startswith("pi"):
    return 139.570
  with open(config.OUTPUT_DIR / "species.csv") as f:
    for r in csv.DictReader(f):
      if r["name"] == name:
        return float(r["mass_mev"])
  raise KeyError(name)


def make_gun(n_events, n_tracks, rng):
  layers = thick_layers()
  rows, meta = [], []
  for ev in range(n_events):
    l = layers[int(rng.integers(len(layers)))]
    center = np.array([rng.uniform(-1e4, 1e4), rng.uniform(-1e4, 1e4),
                       0.5 * (l.zmin + l.zmax)])
    rows.append([ev, "center", *center, 0., 0., 0.])
    starts = []
    for _ in range(n_tracks):
      for _ in range(MAX_TRIES):
        xy = center[:2] + rng.uniform(-VIEW_HALF_UM, VIEW_HALF_UM, 2)
        if all(np.hypot(*(xy - s[:2])) >= MIN_SEPARATION_UM
               for s in starts):
          break
      z = rng.uniform(l.zmin + LAYER_MARGIN_UM, l.zmax - LAYER_MARGIN_UM)
      pos = np.array([xy[0], xy[1], z])
      starts.append(pos)
      name = rng.choice(list(SPECIES))
      p = np.exp(rng.uniform(*np.log(SPECIES[name])))
      dip, phi = np.radians(sample_dip(rng)), rng.uniform(0., 2. * np.pi)
      d = np.array([np.cos(dip) * np.cos(phi), np.cos(dip) * np.sin(phi),
                    np.sin(dip)])
      rows.append([ev, name, *pos, *(p * d)])
      m = species_mass(name)
      meta.append({"event": ev, "species": name,
                   "dip_deg": round(np.degrees(dip), 3),
                   "azimuth_deg": round(np.degrees(phi), 3),
                   "p_mev": round(p, 3),
                   "ke_mev": round(np.hypot(p, m) - m, 4),
                   "x_um": pos[0], "y_um": pos[1], "z_um": pos[2]})
  return rows, meta


def track_table(root, meta):
  """Per-track table: gun parameters joined with the Geant4 primaries
  (in the order of the gun rows of each event)."""
  from .synth.event import EventFile
  ef = EventFile(str(root))
  by_event = {}
  for m in meta:
    by_event.setdefault(m["event"], []).append(m)
  out = []
  for ev, ms in by_event.items():
    prim = sorted((t for t in ef.event(ev).tracks if t.parent == 0),
                  key=lambda t: t.id)
    for m, t in zip(ms, prim):
      length = float(np.linalg.norm(np.diff(t.points, axis=0),
                                    axis=1).sum()) if len(t.points) > 1 \
        else 0.
      out.append({**m, "track_id": t.id, "g4_name": t.name,
                  "range_um": round(t.range_um, 2),
                  "length_um": round(length, 2),
                  "dedx_mev_per_um": round(float(t.edep.sum())
                                           / max(length, 1e-9), 6),
                  "stops": not t.escaped})
  return out


def main():
  parser = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument("--name", required=True)
  parser.add_argument("--events", type=int, default=DEFAULT_EVENTS)
  parser.add_argument("--tracks", type=int, default=DEFAULT_TRACKS)
  parser.add_argument("--seed", type=int, default=1)
  args = parser.parse_args()
  rng = np.random.default_rng(args.seed)
  rows, meta = make_gun(args.events, args.tracks, rng)
  gun = config.OUTPUT_DIR / "gun" / f"{args.name}.csv"
  gun.parent.mkdir(parents=True, exist_ok=True)
  with open(gun, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["event", "name", "x_um", "y_um", "z_um", "px_mev",
                "py_mev", "pz_mev"])
    w.writerows(rows)
  r = subprocess.run([sys.executable, "-m", "module.simulate", "--gun",
                      str(gun), "--events", str(args.events), "--name",
                      args.name, "--seed", str(args.seed)],
                     cwd=Path(__file__).resolve().parent.parent)
  if r.returncode:
    return r.returncode
  table = track_table(config.SIM_DIR / f"{args.name}.root", meta)
  out = config.SIM_DIR / f"{args.name}_tracks.csv"
  with open(out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(table[0]))
    w.writeheader()
    w.writerows(table)
  print(f"{len(table)} tracks -> {out}")
  return 0


if __name__ == "__main__":
  raise SystemExit(main())
