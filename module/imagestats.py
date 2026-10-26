"""Image statistics to compare synthetic stacks with real full-scan views.

  ./run.py compare --real <view.json> --synth data/images/<name>/ev0000 \
    [--slices 3] -o compare.png

Both stacks are read the same way (the real one from its SPNG container)
and the same statistics are computed on the central slices:

- raw intensity: median, std, 1 / 5 percentiles
- darkest 5 % of pixels: number of connected blobs per Mpx and the blob
  area distribution (clustered grains vs pixel noise)
- after a full-scan-like preprocessing (Gaussian background subtraction
  of sigma ~8 px as cv2 ksize 51, Otsu threshold): foreground fraction
  and connected-component area / elongation distributions
- radially averaged power spectrum of the background-subtracted image
- correlation of background-subtracted slices k and k+1, k+2 (depth
  correlation of the defocus blur)
- response of point-cloud methods: line-like pixels of the Hessian of
  the background-subtracted slice (thin strong / wide weak lines), DoG
  maxima per Mpx and the fraction of them on lines (point-like vs
  line-like structure, e.g. smooth flat tracks vs beaded grains)
"""

import argparse
import json
import os
import tempfile
from pathlib import Path

import numpy as np

BG_SIGMA_PX = 8.0      # cv2.GaussianBlur ksize 51 -> sigma ~8 px
DARK_FRACTION = 0.05
N_SPECTRUM_BINS = 120
HIST_BINS = 64
MAX_LAG = 2
DEFAULT_CROP_PX = 1024
BRIGHT_LEVEL = 250  # near-saturated gray level
MP_DARK_LEVEL = 20.  # dark detail level in the minimum projection
# Bands of the raw radial power spectrum (cycles / px)
SPECTRUM_BANDS = {"pow_low": (0.004, 0.02), "pow_mid": (0.02, 0.08),
                  "pow_high": (0.08, 0.3)}
# Elongated components (track segments): elongation and minimum area
ELONG_MIN = 5.
ELONG_MIN_AREA = 30
# Line response: Hessian scale (px), scale-normalized level (gray) and
# maximum eigenvalue ratio |l_along| / |l_across| of a line
LINE_THIN = (1.5, 10.)
LINE_WIDE = (3., 3.)
LINE_ON = (3., 6.)       # line mask for the DoG points on lines
LINE_MAX_RATIO = 0.3
# DoG point detector: inner / outer sigma (px) and level (gray)
DOG_SIGMAS = (1., 3.)
DOG_LEVEL = 4.


def decode_png(data):
  """Decode PNG bytes to a 2-D uint8 array (OpenCV if available, else
  ROOT TASImage)."""
  try:
    import cv2
    return cv2.imdecode(np.frombuffer(data, np.uint8),
                        cv2.IMREAD_GRAYSCALE)
  except ImportError:
    pass
  import ROOT
  ROOT.gROOT.SetBatch(True)
  with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
    f.write(data)
    name = f.name
  try:
    img = ROOT.TASImage(name)
    w, h = img.GetWidth(), img.GetHeight()
    argb = np.frombuffer(img.GetArgbArray(), dtype=np.uint32, count=w * h)
    return (argb.reshape(h, w) & 0xff).astype(np.uint8)
  finally:
    os.unlink(name)


def read_stack(json_path, indices=None):
  """Slices of a full-scan view (SPNG or plain PNG entries)."""
  json_path = Path(json_path)
  meta = json.loads(json_path.read_text())
  images = meta["Images"]
  indices = range(len(images)) if indices is None else indices
  out = []
  for i in indices:
    path = images[i]["Path"]
    if "&" in path:
      name, offset, length = path.split("&")
      with open(json_path.parent / name, "rb") as f:
        f.seek(int(offset))
        data = f.read(int(length))
    else:
      data = (json_path.parent / path).read_bytes()
    out.append(decode_png(data))
  return out


def n_slices(json_path):
  return len(json.loads(Path(json_path).read_text())["Images"])


def gaussian_blur(img, sigma):
  n, m = img.shape
  fy = np.fft.fftfreq(n)[:, None]
  fx = np.fft.rfftfreq(m)[None, :]
  k = np.exp(-2. * np.pi ** 2 * sigma ** 2 * (fx ** 2 + fy ** 2))
  return np.fft.irfft2(np.fft.rfft2(img) * k, s=img.shape)


def otsu(values):
  hist, edges = np.histogram(values, bins=256)
  p = hist / hist.sum()
  centers = 0.5 * (edges[:-1] + edges[1:])
  w0 = np.cumsum(p)
  mu = np.cumsum(p * centers)
  between = (mu[-1] * w0 - mu) ** 2 / np.maximum(w0 * (1. - w0), 1e-12)
  return centers[np.argmax(between)]


def label_components(mask):
  """4-connected labels: minimum propagation with pointer jumping."""
  n, m = mask.shape
  big = n * m
  idx = np.arange(big).reshape(n, m)
  parent = np.where(mask, idx, big).ravel()
  parent = np.append(parent, big)  # sentinel for background
  flat_mask = mask.ravel()
  pairs = []
  for a, b in ((idx[1:, :], idx[:-1, :]), (idx[:, 1:], idx[:, :-1])):
    a, b = a.ravel(), b.ravel()
    ok = flat_mask[a] & flat_mask[b]
    pairs.append((a[ok], b[ok]))
  a = np.concatenate([p[0] for p in pairs])
  b = np.concatenate([p[1] for p in pairs])
  while True:
    ra, rb = parent[a], parent[b]
    lo = np.minimum(ra, rb)
    old = parent.copy()
    # Hook the larger root to the smaller one, then compress paths
    np.minimum.at(parent, ra, lo)
    np.minimum.at(parent, rb, lo)
    while True:
      jumped = parent[parent]
      if np.array_equal(jumped, parent):
        break
      parent = jumped
    if np.array_equal(parent, old):
      break
  return parent[:-1].reshape(n, m)


def blob_contrast(img, mask):
  """Per blob: (local background - minimum) / local background."""
  lab = label_components(mask)
  sel = lab[mask]
  _, inv = np.unique(sel, return_inverse=True)
  bg = gaussian_blur(img, BG_SIGMA_PX)[mask]
  vals = img[mask]
  n = inv.max() + 1 if len(inv) else 0
  vmin = np.full(n, np.inf)
  np.minimum.at(vmin, inv, vals)
  bmean = np.bincount(inv, bg, n) / np.maximum(np.bincount(inv, None, n), 1)
  return (bmean - vmin) / np.maximum(bmean, 1.)


def component_stats(mask):
  """Areas and elongations (sqrt of eigenvalue ratio) of components."""
  lab = label_components(mask)
  sel = lab[mask]
  ys, xs = np.nonzero(mask)
  uniq, inv, area = np.unique(sel, return_inverse=True, return_counts=True)
  sx = np.bincount(inv, xs)
  sy = np.bincount(inv, ys)
  sxx = np.bincount(inv, xs * xs)
  syy = np.bincount(inv, ys * ys)
  sxy = np.bincount(inv, xs * ys)
  mx, my = sx / area, sy / area
  cxx = sxx / area - mx ** 2
  cyy = syy / area - my ** 2
  cxy = sxy / area - mx * my
  tr = cxx + cyy
  det = cxx * cyy - cxy ** 2
  disc = np.sqrt(np.maximum(tr ** 2 / 4. - det, 0.))
  l1, l2 = tr / 2. + disc, np.maximum(tr / 2. - disc, 1e-3)
  return area, np.sqrt(l1 / l2)


def line_mask(d, sigma, level):
  """Pixels of dark lines (bright ridges of the darkness d) at scale
  sigma: scale-normalized curvature across above level, anisotropic."""
  g = gaussian_blur(d, sigma)
  gy, gx = np.gradient(g)
  gyy, gyx = np.gradient(gy)
  gxx = np.gradient(gx, axis=1)
  tr, det = gxx + gyy, gxx * gyy - gyx ** 2
  disc = np.sqrt(np.maximum(tr ** 2 / 4. - det, 0.))
  across, along = tr / 2. - disc, tr / 2. + disc
  return (sigma ** 2 * -across > level) \
    & (np.abs(along) < LINE_MAX_RATIO * np.abs(across))


def dog_points(d):
  """Local maxima (3x3) of the difference of Gaussians of d."""
  g = gaussian_blur(d, DOG_SIGMAS[0]) - gaussian_blur(d, DOG_SIGMAS[1])
  p = np.pad(g, 1, constant_values=-np.inf)
  n = g.shape
  mx = np.max([p[i:i + n[0], j:j + n[1]] for i in range(3)
               for j in range(3)], axis=0)
  return (g >= mx) & (g > DOG_LEVEL)


def radial_spectrum(img):
  f = np.abs(np.fft.fftshift(np.fft.fft2(img - img.mean()))) ** 2
  n = img.shape[0]
  y, x = np.indices(f.shape)
  r = np.hypot(x - n / 2, y - n / 2) / n
  edges = np.linspace(0., 0.5, N_SPECTRUM_BINS + 1)
  idx = np.digitize(r.ravel(), edges) - 1
  ok = (idx >= 0) & (idx < N_SPECTRUM_BINS)
  power = np.bincount(idx[ok], f.ravel()[ok], N_SPECTRUM_BINS)
  count = np.bincount(idx[ok], minlength=N_SPECTRUM_BINS)
  return 0.5 * (edges[:-1] + edges[1:]), power / np.maximum(count, 1)


def stack_stats(stack, crop, offsets=None):
  """Statistics of slices (central crop px). `offsets` gives the slice
  offset of each image from the central one (default: consecutive)."""
  n = stack[0].shape[0]
  lo = max(0, (n - crop) // 2)
  imgs = [s[lo:lo + crop, lo:lo + crop].astype(np.float64) for s in stack]
  if offsets is None:
    offsets = [i - len(imgs) // 2 for i in range(len(imgs))]
  c = offsets.index(0)
  mid = imgs[c]
  near = [i for i, off in zip(imgs, offsets) if abs(off) <= MAX_LAG]
  raw = np.concatenate([i.ravel() for i in near])
  diffs = [gaussian_blur(i, BG_SIGMA_PX) - i for i in imgs]
  d = diffs[c]
  thr = otsu(d)
  fg = d > thr
  zd = (d - d.mean()) / (d.std() or 1.)
  q25, q75 = np.percentile(d, [25, 75])
  area, elong = component_stats(fg)
  n_elong = int(((elong > ELONG_MIN) & (area >= ELONG_MIN_AREA)).sum())
  dark = mid <= np.quantile(mid, DARK_FRACTION)
  dark_area, _ = component_stats(dark)
  contrast = blob_contrast(mid, dark)
  k, power = radial_spectrum(d)
  # Minimum projection over all given slices (depth-integrated grains)
  mp = np.min(np.array(imgs), axis=0)
  dmp = gaussian_blur(mp, BG_SIGMA_PX) - mp
  # Spectrum of the raw image normalized by its mean (keeps the clouds)
  kr, praw = radial_spectrum(mid / mid.mean())
  bands = {}
  for name, (f0, f1) in SPECTRUM_BANDS.items():
    sel = (kr >= f0) & (kr < f1)
    bands[name] = float(np.log10(np.mean(praw[sel])))
  corr = {}
  for k, off in enumerate(offsets):
    if off > 0:
      corr[off] = float(np.corrcoef(diffs[c].ravel(),
                                    diffs[k].ravel())[0, 1])
  mpx = mid.size / 1e6
  dog = dog_points(d)
  on_line = line_mask(d, *LINE_ON)
  return {
    "median": float(np.median(raw)), "std": float(raw.std()),
    "p1": float(np.percentile(raw, 1)), "p5": float(np.percentile(raw, 5)),
    "p99": float(np.percentile(raw, 99)),
    "frac_bright": float(np.mean(raw >= BRIGHT_LEVEL)),
    "fg_fraction": float(fg.mean()), "n_cc_per_mpx": len(area) / mpx,
    "d_skew": float((zd ** 3).mean()), "d_kurt": float((zd ** 4).mean()),
    "d_iqr": float(q75 - q25),
    "cc_area_median": float(np.median(area)),
    "cc_area_p90": float(np.percentile(area, 90)),
    "cc_elong_p90": float(np.percentile(elong, 90)),
    "n_elong_per_mpx": n_elong / (mid.size / 1e6),
    "dark_blobs_per_mpx": len(dark_area) / mpx,
    "dark_blob_area_median": float(np.median(dark_area)),
    "dark_contrast_median": float(np.median(contrast)),
    "dark_contrast_p90": float(np.percentile(contrast, 90)),
    "slice_corr": corr,
    "mp_detail_std": float(dmp.std()),
    "mp_dark_frac": float(np.mean(dmp > MP_DARK_LEVEL)),
    **bands,
    "line_thin": float(line_mask(d, *LINE_THIN).mean()),
    "line_wide": float(line_mask(d, *LINE_WIDE).mean()),
    "dog_per_mpx": float(dog.sum() / mpx),
    "dog_on_line": float(dog[on_line].sum() / max(dog.sum(), 1)),
    "_raw": raw, "_area": area, "_elong": elong, "_k": k,
    "_power": power, "_dark_area": dark_area,
  }


def central(n, half):
  c = n // 2
  return list(range(max(0, c - half), min(n, c + half + 1)))


def plot(results, output):
  import ROOT
  ROOT.gROOT.SetBatch(True)
  ROOT.gStyle.SetOptStat(0)
  canvas = ROOT.TCanvas("c", "", 1600, 900)
  canvas.Divide(3, 2)
  colors = [ROOT.kBlack, ROOT.kRed + 1, ROOT.kBlue + 1, ROOT.kGreen + 2]
  keep = []

  def hist_panel(pad, key, title, bins, logx=False):
    canvas.cd(pad)
    if logx:
      ROOT.gPad.SetLogx()
    ROOT.gPad.SetLogy()
    leg = ROOT.TLegend(0.55, 0.75, 0.9, 0.9)
    for i, (label, s) in enumerate(results):
      h = ROOT.TH1D(f"h{pad}_{i}", title, len(bins) - 1,
                    np.asarray(bins, dtype=float))
      for v in s[key]:
        h.Fill(v)
      if h.Integral() > 0:
        h.Scale(1. / h.Integral())
      h.SetLineColor(colors[i % len(colors)])
      h.SetLineWidth(2)
      h.Draw("hist" if i == 0 else "hist same")
      leg.AddEntry(h, label, "l")
      keep.append(h)
    leg.Draw()
    keep.append(leg)

  hist_panel(1, "_raw", "raw intensity;gray level;fraction",
             np.linspace(0, 256, HIST_BINS + 1))
  area_bins = np.logspace(0, 4, 41)
  hist_panel(2, "_area", "components after preprocessing;area (px);"
             "fraction", area_bins, logx=True)
  hist_panel(3, "_dark_area", "darkest 5% blobs;area (px);fraction",
             area_bins, logx=True)
  hist_panel(4, "_elong", "component elongation;sqrt(l1/l2);fraction",
             np.logspace(0, 1.5, 31), logx=True)
  canvas.cd(5)
  ROOT.gPad.SetLogy()
  leg = ROOT.TLegend(0.55, 0.75, 0.9, 0.9)
  for i, (label, s) in enumerate(results):
    g = ROOT.TGraph(len(s["_k"]), s["_k"].astype(float),
                    s["_power"].astype(float))
    g.SetTitle("power spectrum (background subtracted);"
               "frequency (1/px);power")
    g.SetLineColor(colors[i % len(colors)])
    g.SetLineWidth(2)
    g.Draw("AL" if i == 0 else "L")
    leg.AddEntry(g, label, "l")
    keep.append(g)
  leg.Draw()
  canvas.SaveAs(output)


def main():
  parser = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument("--real", action="append", default=[],
                      help="real view JSON (repeatable)")
  parser.add_argument("--synth", action="append", default=[],
                      help="synthetic event directory (repeatable)")
  parser.add_argument("--slices", type=int, default=2,
                      help="half number of central slices used")
  parser.add_argument("--crop", type=int, default=DEFAULT_CROP_PX,
                      help="central square used for the statistics (px)")
  parser.add_argument("-o", "--output", help="comparison plot (PNG)")
  args = parser.parse_args()

  results = []
  for path in args.real:
    idx = central(n_slices(path), args.slices)
    results.append((f"real {Path(path).name[:9]}",
                    stack_stats(read_stack(path, idx), args.crop)))
  for d in args.synth:
    path = Path(d) / "image.json"
    idx = central(n_slices(path), args.slices)
    results.append((f"synth {Path(d).parent.name}/{Path(d).name}",
                    stack_stats(read_stack(path, idx), args.crop)))
  keys = [k for k in results[0][1] if not k.startswith("_")]
  print(f"{'':24s}" + "".join(f"{lab[:22]:>24s}" for lab, _ in results))
  for k in keys:
    vals = [s[k] for _, s in results]
    if isinstance(vals[0], dict):
      cells = [" ".join(f"{l}:{v:.2f}" for l, v in x.items()) for x in vals]
    else:
      cells = [f"{v:.4g}" for v in vals]
    print(f"{k:24s}" + "".join(f"{c:>24s}" for c in cells))
  if args.output:
    plot(results, args.output)


if __name__ == "__main__":
  main()
