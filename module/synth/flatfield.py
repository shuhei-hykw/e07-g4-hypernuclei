"""Illumination profile (flat field) of the full-scan microscope.

  ./run.py flatfield [--views 20]

The median of the middle slices of random real views, smoothed and
normalized to 1 in the central region, captures the circular field stop
and the vignetting. It is stored downsampled in output/flatfield.npy and
multiplies the brightness of fully synthetic images.
"""

import argparse
import random
from pathlib import Path

import numpy as np

from .. import config
from .. import imagestats as st
from . import calibrate

STORE_PX = 256
SMOOTH_PX = 2.          # on the stored grid
CENTER_FRACTION = 0.25  # half width of the normalization region
FLATFIELD_FILE = "flatfield.npy"


def compute(n_views, seed):
  views = [v for v in sorted(Path(calibrate.REAL_DIR).glob("V*.json"))
           if not v.name.startswith(calibrate.EXCLUDED_VIEWS)]
  random.Random(seed).shuffle(views)
  images = []
  for v in views[:n_views]:
    mid = st.n_slices(v) // 2
    img = st.read_stack(v, [mid])[0].astype(np.float32)
    f = img.shape[0] // STORE_PX
    images.append(img.reshape(STORE_PX, f, STORE_PX, f).mean(axis=(1, 3)))
  flat = st.gaussian_blur(np.median(images, axis=0), SMOOTH_PX)
  c, h = STORE_PX // 2, int(STORE_PX * CENTER_FRACTION)
  return flat / flat[c - h:c + h, c - h:c + h].mean()


def load(size_px):
  """Flat field at size_px resolution (ones if not computed)."""
  path = config.OUTPUT_DIR / FLATFIELD_FILE
  if not path.exists():
    return np.ones((size_px, size_px))
  flat = np.load(path)
  # Bilinear upsampling to the image grid
  x = (np.arange(size_px) + 0.5) * flat.shape[0] / size_px - 0.5
  x = np.clip(x, 0, flat.shape[0] - 1)
  i0 = np.floor(x).astype(int)
  i1 = np.minimum(i0 + 1, flat.shape[0] - 1)
  w = x - i0
  rows = flat[i0] * (1 - w)[:, None] + flat[i1] * w[:, None]
  return rows[:, i0] * (1 - w)[None, :] + rows[:, i1] * w[None, :]


def main():
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument("--views", type=int, default=20)
  parser.add_argument("--seed", type=int, default=1)
  args = parser.parse_args()
  flat = compute(args.views, args.seed)
  np.save(config.OUTPUT_DIR / FLATFIELD_FILE, flat.astype(np.float32))
  print(f"flat field min {flat.min():.3f} max {flat.max():.3f} -> "
        f"{config.OUTPUT_DIR / FLATFIELD_FILE}")


if __name__ == "__main__":
  main()
