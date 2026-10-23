"""Grain density of the faintest flat tracks: an in-plane, step-free
estimate of the MIP grain density.

  ./run.py mipdensity --real <view.json> ... --synth <event dir> ...

Long straight components of single binarized slices (flat tracks in
focus) are found, the background-subtracted darkness is sampled along
their axis and grain peaks are counted; the grain density (per um, from
the exact pixel scale) of each track is collected. The minimum ionizing
tracks form the low end of the distribution, summarized by its low
percentiles. Synthetic views (known MIP density, image model v14) are
measured identically to calibrate the estimator, which removes the bias
of merged or missed grains.

Result: at the 0.289 um/px, 20x resolution the counts saturate on the
fine texture along the lines and do not follow the grain density
(doubling it in the synthetic images leaves them unchanged), so this
estimator cannot provide the MIP density; see output/dz_comparison.md.
"""

import argparse
from pathlib import Path

import numpy as np

from . import imagestats as st

PIXEL_UM = 0.289
CROP_PX = 1024
BG_SIGMA_PX = 8.
MIN_LENGTH_PX = 60
MIN_ELONGATION = 8.
PEAK_PROMINENCE = 5.
MIN_PEAK_SEPARATION_PX = 2
PERCENTILES = (10, 25, 50)
SLICES = (-10, 0, 10)  # offsets from the central slice


def tracks_in_slice(img):
  c = img.astype(float)
  d = st.gaussian_blur(c, BG_SIGMA_PX) - c   # darkness
  fg = d > st.otsu(d)
  lab = st.label_components(fg)
  sel = lab[fg]
  ys, xs = np.nonzero(fg)
  densities = []
  for comp in np.unique(sel):
    m = sel == comp
    if m.sum() < MIN_LENGTH_PX:
      continue
    px, py = xs[m].astype(float), ys[m].astype(float)
    cov = np.cov(np.vstack([px, py]))
    w, v = np.linalg.eigh(cov)
    if w[0] <= 0 or np.sqrt(w[1] / w[0]) < MIN_ELONGATION:
      continue
    axis = v[:, 1]
    t = (px - px.mean()) * axis[0] + (py - py.mean()) * axis[1]
    length = t.max() - t.min()
    if length < MIN_LENGTH_PX:
      continue
    # Darkness profile along the axis (1 px steps, nearest pixel)
    s = np.arange(t.min(), t.max())
    qx = np.rint(px.mean() + s * axis[0]).astype(int)
    qy = np.rint(py.mean() + s * axis[1]).astype(int)
    ok = (qx >= 0) & (qx < d.shape[1]) & (qy >= 0) & (qy < d.shape[0])
    prof = d[qy[ok], qx[ok]]
    n, last = 0, -MIN_PEAK_SEPARATION_PX
    for k in range(1, len(prof) - 1):
      if prof[k] >= prof[k - 1] and prof[k] > prof[k + 1] \
         and k - last >= MIN_PEAK_SEPARATION_PX:
        lo = min(prof[max(0, k - 4):k].min(), prof[k + 1:k + 5].min())
        if prof[k] - lo >= PEAK_PROMINENCE:
          n += 1
          last = k
    densities.append(n / (len(prof) * PIXEL_UM))
  return densities


def view_densities(json_path):
  n = st.n_slices(json_path)
  idx = [n // 2 + o for o in SLICES]
  out = []
  for img in st.read_stack(json_path, idx):
    lo = (img.shape[0] - CROP_PX) // 2
    out += tracks_in_slice(img[lo:lo + CROP_PX, lo:lo + CROP_PX])
  return out


def summary(vals):
  v = np.asarray(vals)
  return (f"n={len(v)} " + " ".join(
    f"p{p}={np.percentile(v, p):.3f}" for p in PERCENTILES) + " /um")


def main():
  parser = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument("--real", nargs="*", default=[])
  parser.add_argument("--synth", nargs="*", default=[])
  args = parser.parse_args()
  real = sum((view_densities(Path(v)) for v in args.real), [])
  synth = sum((view_densities(Path(d) / "image.json") for d in args.synth),
              [])
  print("real ", summary(real))
  print("synth", summary(synth))
  for p in PERCENTILES:
    r, s = np.percentile(real, p), np.percentile(synth, p)
    print(f"p{p}: real / synthetic = {r / s:.2f}  -> MIP density estimate "
          f"{0.148 * r / s:.3f} /um ({14.8 * r / s:.1f} / 100 um)")


if __name__ == "__main__":
  main()
