"""Physical focal step from the emulsion layer extent in the stack.

  ./run.py layerdepth <image.json> [<image.json> ...] [--sh 1.9]

A geometric handle on the step, independent of the grain density: deep
stacks (e.g. the 200-slice special events) contain whole emulsion
layers, bounded by the base and by the surface. The fraction of dark
(background-subtracted) pixels per slice drops to zero outside the
emulsion; each layer is located by the half-height crossings of this
profile, and a layer whose two edges lie inside the stack gives

  physical step = (EMULSION_UM / Sh) / (number of slices of the layer)

for the developed layer (E07 thick plate, 450 um per side before
development). The base gap (BASE_UM, not shrunk) gives the ratio of the
geometric to the stage thickness in the base, as a cross-check of the
focal shift factor.
"""

import argparse
import json
from pathlib import Path

import numpy as np

from . import imagestats as st

CROP_PX = 1024
BG_SIGMA_PX = 8.
DARK_LEVEL = 15.      # gray levels below the local background
EMULSION_UM = 450.    # thick plate, one side, before development
BASE_UM = 40.         # thick plate base
DEFAULT_SH = 1.9
PLATEAU_SLICES = 5    # slices inside an edge that define its plateau
EMPTY_FRACTION = 0.1  # of the maximum: slice outside the emulsion
UM_PER_MM = 1000.


def dark_profile(json_path):
  n = st.n_slices(json_path)
  out = []
  for k in range(n):
    img = st.read_stack(json_path, [k])[0]
    lo = (img.shape[0] - CROP_PX) // 2
    c = img[lo:lo + CROP_PX, lo:lo + CROP_PX].astype(float)
    out.append(((st.gaussian_blur(c, BG_SIGMA_PX) - c) > DARK_LEVEL).mean())
  return np.array(out)


def edges(prof):
  """Half-height crossings (fractional slice) between emulsion and empty
  slices, with the direction +1 (entering emulsion) or -1."""
  inside = prof > EMPTY_FRACTION * prof.max()
  out = []
  for k in np.flatnonzero(inside[1:] != inside[:-1]):
    rising = inside[k + 1]
    sl = slice(k + 1, k + 1 + PLATEAU_SLICES) if rising \
      else slice(max(0, k + 1 - PLATEAU_SLICES), k + 1)
    half = 0.5 * prof[sl].mean()
    # Walk from the empty slice into the emulsion to the crossing
    step = 1 if rising else -1
    j = k if rising else k + 1
    while 0 <= j + step < len(prof) and prof[j + step] < half:
      j += step
    if 0 <= j + step < len(prof):
      a, b = prof[j], prof[j + step]
      out.append((j + step * (half - a) / (b - a), 1 if rising else -1))
  return out


def main():
  parser = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument("views", nargs="+")
  parser.add_argument("--sh", type=float, default=DEFAULT_SH)
  args = parser.parse_args()
  steps = []
  for v in args.views:
    meta = json.loads(Path(v).read_text())
    z = [im["z"] for im in meta["Images"]]
    stage_um = abs(z[1] - z[0]) * UM_PER_MM
    e = edges(dark_profile(Path(v)))
    for (a, da), (b, db) in zip(e[:-1], e[1:]):
      n = b - a
      if da > 0 and db < 0:
        step = EMULSION_UM / args.sh / n
        steps.append(step)
        print(f"{v}: layer {n:.1f} slices = {n * stage_um:.0f} um stage "
              f"-> physical step {step:.2f} um "
              f"(developed / stage {step / stage_um:.2f})")
      elif da < 0 and db > 0:
        print(f"{v}: base gap {n:.1f} slices = {n * stage_um:.0f} um stage "
              f"-> geometric / stage {BASE_UM / (n * stage_um):.2f}")
  if steps:
    print(f"physical step: mean {np.mean(steps):.2f} um, spread "
          f"{np.std(steps):.2f} um over {len(steps)} layers")


if __name__ == "__main__":
  main()
