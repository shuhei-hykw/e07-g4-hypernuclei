"""Grains along nearly perpendicular beam tracks: a handle on the
physical focal step per slice.

  ./run.py beamgrains <view.json> [<view.json> ...] [--n-tracks 150]

A MIP track has a known grain density along its length (14.8 / 100 um,
undeveloped), i.e. about 0.28 grains per developed um along z for a
track perpendicular to the sheets. The number of grains seen per slice
along such tracks therefore scales with the physical focal step per
slice, independently of the defocus blur and beam divergence that make
the slice statistics degenerate with the step.

Beam tracks are found as persistent dark dots: local minima of the mean
of background-subtracted slices that are dark in most slices. For each,
the depth profile of the background-subtracted darkness (maximum over a
small neighbourhood, following a slow drift) is searched for peaks; the
result is the number of peaks per slice. Real views and synthetic views
rendered with known steps are measured the same way.
"""

import argparse
from pathlib import Path

import numpy as np

from . import imagestats as st

CROP_PX = 1024
BG_SIGMA_PX = 8.
DEFAULT_N_TRACKS = 150
MIN_DOT_SEPARATION_PX = 15
NEIGHBOUR_PX = 2
PERSISTENCE = 0.6    # fraction of slices where the dot must be dark
DARK_LEVEL = 6.      # gray levels below the local background
PEAK_PROMINENCE = 6. # gray levels
SLICE_STEP = 1


def detail_stack(json_path):
  n = st.n_slices(json_path)
  stack = st.read_stack(json_path, list(range(0, n, SLICE_STEP)))
  lo = (stack[0].shape[0] - CROP_PX) // 2
  out = []
  for img in stack:
    c = img[lo:lo + CROP_PX, lo:lo + CROP_PX].astype(float)
    out.append(st.gaussian_blur(c, BG_SIGMA_PX) - c)  # darkness
  return np.array(out)


def find_dots(det, n_tracks):
  """Persistent dark dots (candidate perpendicular tracks)."""
  persistent = (det > DARK_LEVEL).mean(axis=0)
  mean = det.mean(axis=0)
  score = np.where(persistent >= PERSISTENCE, mean, -np.inf)
  order = np.argsort(score.ravel())[::-1]
  chosen = []
  for idx in order:
    if not np.isfinite(score.ravel()[idx]) or len(chosen) >= n_tracks:
      break
    y, x = divmod(int(idx), CROP_PX)
    if min(y, x) < NEIGHBOUR_PX + 1 or max(y, x) > CROP_PX - NEIGHBOUR_PX - 2:
      continue
    if all((y - cy) ** 2 + (x - cx) ** 2 >= MIN_DOT_SEPARATION_PX ** 2
           for cy, cx in chosen):
      chosen.append((y, x))
  return chosen


def depth_profile(det, y, x):
  """Darkness along z following the darkest pixel near (y, x)."""
  prof = []
  for k in range(det.shape[0]):
    win = det[k, y - NEIGHBOUR_PX:y + NEIGHBOUR_PX + 1,
              x - NEIGHBOUR_PX:x + NEIGHBOUR_PX + 1]
    prof.append(win.max())
  return np.array(prof)


def count_peaks(prof):
  """Local maxima with a prominence above PEAK_PROMINENCE."""
  n = 0
  for k in range(1, len(prof) - 1):
    if prof[k] >= prof[k - 1] and prof[k] > prof[k + 1]:
      # prominence against the minima within 3 slices on both sides
      lmin = prof[max(0, k - 3):k].min()
      rmin = prof[k + 1:k + 4].min()
      if prof[k] - max(lmin, rmin) >= PEAK_PROMINENCE:
        n += 1
  return n


def measure(json_path, n_tracks):
  det = detail_stack(json_path)
  dots = find_dots(det, n_tracks)
  rates = [count_peaks(depth_profile(det, y, x)) / det.shape[0]
           for y, x in dots]
  return np.array(rates), len(dots)


def main():
  parser = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument("views", nargs="+",
                      help="image.json of real or synthetic views")
  parser.add_argument("--n-tracks", type=int, default=DEFAULT_N_TRACKS)
  args = parser.parse_args()
  for v in args.views:
    rates, n = measure(Path(v), args.n_tracks)
    print(f"{v}: {n} dots, peaks per slice median {np.median(rates):.3f} "
          f"mean {rates.mean():.3f} +- {rates.std() / np.sqrt(len(rates)):.3f}")


if __name__ == "__main__":
  main()
