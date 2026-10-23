"""Transverse profile of tracks in image stacks.

The minimum projection over a slice range is sampled perpendicular to a
track segment; the background-normalized profile gives the track depth
(1 - I/I_bg at the center) and its full width at half depth.
"""

import numpy as np

from .. import imagestats as st

HALF_WIDTH_PX = 12
SAMPLE_STEP_PX = 1.
BG_EDGE_PX = 4  # outermost samples on each side used as background


def min_projection(json_path, slices):
  return np.min(np.array(st.read_stack(json_path, slices)), axis=0)


def profile(image, p0, p1):
  """Mean perpendicular profile along the segment p0 -> p1 (col, row)."""
  p0, p1 = np.asarray(p0, float), np.asarray(p1, float)
  d = p1 - p0
  length = np.linalg.norm(d)
  u = d / length
  n = np.array([-u[1], u[0]])
  offsets = np.arange(-HALF_WIDTH_PX, HALF_WIDTH_PX + 1, SAMPLE_STEP_PX)
  along = np.arange(0., length, SAMPLE_STEP_PX)
  pts = p0[None, None, :] + along[:, None, None] * u \
    + offsets[None, :, None] * n
  col = np.clip(np.rint(pts[..., 0]).astype(int), 0, image.shape[1] - 1)
  row = np.clip(np.rint(pts[..., 1]).astype(int), 0, image.shape[0] - 1)
  prof = image[row, col].astype(float).mean(axis=0)
  return offsets, prof


def depth_and_width(offsets, prof):
  bg = np.concatenate([prof[:BG_EDGE_PX], prof[-BG_EDGE_PX:]]).mean()
  depth_profile = 1. - prof / bg
  depth = depth_profile.max()
  above = offsets[depth_profile >= depth / 2.]
  return float(depth), float(above.max() - above.min() + SAMPLE_STEP_PX)


SEARCH_SHIFT_PX = 8
SEARCH_ANGLE_DEG = 4.
SEARCH_STEPS = 9


def refine(image, p0, p1):
  """Shift and rotate the segment around its middle to maximize depth."""
  p0, p1 = np.asarray(p0, float), np.asarray(p1, float)
  mid, half = (p0 + p1) / 2., (p1 - p0) / 2.
  u = half / np.linalg.norm(half)
  n = np.array([-u[1], u[0]])
  best = (-1., p0, p1)
  for shift in np.linspace(-SEARCH_SHIFT_PX, SEARCH_SHIFT_PX, 2 * SEARCH_SHIFT_PX + 1):
    for ang in np.radians(np.linspace(-SEARCH_ANGLE_DEG, SEARCH_ANGLE_DEG,
                                      SEARCH_STEPS)):
      c, s = np.cos(ang), np.sin(ang)
      h = np.array([c * half[0] - s * half[1], s * half[0] + c * half[1]])
      q0, q1 = mid + shift * n - h, mid + shift * n + h
      depth, _ = depth_and_width(*profile(image, q0, q1))
      if depth > best[0]:
        best = (depth, q0, q1)
  return best[1], best[2]
