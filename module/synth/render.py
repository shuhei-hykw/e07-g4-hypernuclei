"""Grain model and image formation of the emulsion microscope.

Coordinates: simulated tracks are in the undeveloped emulsion (um). The
emulsion layers shrink along z by the factor `shrink` on development
(base films do not), x and y are unchanged. Images are taken in the
developed emulsion: pixel (column, row) = (x, y) / pixel_um around the
view center, slices every dz_um in developed depth.

Image formation: each silver grain adds optical density
  D = d0 * (s0 / s)^2 * exp(-r^2 / 2 s^2),  s^2 = s0^2 + (k dz)^2
where dz is its distance from the focal plane, so a defocused grain is
wider and fainter. The intensity is I = B(x, y) exp(-sum D) + noise,
which saturates for dense (heavily ionizing) tracks; a residual
transmission keeps the darkest pixels at the level seen in real images:
  I = B (t_min + (1 - t_min) exp(-sum D)).
"""

from dataclasses import dataclass

import numpy as np

# Label classes of labels.npz (0: background / fog)
LABEL_CAPTURE = 1
LABEL_DECAY = 2
LABEL_OTHER = 3
LABEL_BEAM = 4    # beam and other straight background tracks
LABEL_JUNK = 5    # large grain clusters
LABEL_OF_ORIGIN = {"capture": LABEL_CAPTURE, "decay": LABEL_DECAY,
                   "other": LABEL_OTHER}


@dataclass(frozen=True)
class Optics:
  pixel_um: float = 0.289      # E07 full scan, 20x objective
  size_px: int = 2048
  n_slices: int = 60
  dz_um: float = 3.0           # focal step in the developed emulsion
  shrink: float = 1.9          # emulsion thickness ratio before / after
  grain_sigma_um: float = 0.35 # in-focus grain image (s0)
  defocus_slope: float = 0.35  # k: blur growth per um of defocus
  max_defocus_um: float = 6.0  # grains farther from focus are ignored
  n_defocus_bins: int = 12
  grain_density: float = 0.8   # d0: optical density of a focused grain
  pad_px: int = 64             # margin against FFT wrap-around
  label_radius_px: int = 2


@dataclass(frozen=True)
class Emulsion:
  grains_mip_per_um: float = 0.148  # MIP grain density (fullscan STATUS)
  grains_sat_per_um: float = 2.5    # saturation for heavy ionization
  dedx_mip_mev_per_um: float = 5.47e-4  # 1.55 MeV cm2/g x 3.53 g/cm3
  fog_per_um3: float = 6e-3         # fog grains per developed um3
  grain_jitter_um: float = 0.1
  grain_spread: float = 0.4         # lognormal sigma of grain darkness


@dataclass(frozen=True)
class Background:
  brightness: float = 167.
  # Relative rms and length scale of the slow brightness variations
  mottle: tuple = ((0.08, 25.), (0.06, 6.))
  min_transmission: float = 0.4     # darkest level relative to B
  noise_sigma: float = 4.
  n_beam_tracks: int = 100          # nearly parallel MIP tracks per view
  beam_slope_max: float = 0.1       # |dx/dz|, |dy/dz| of the beam axis
  beam_spread: float = 0.01         # rms slope around the beam axis
  n_random_tracks: int = 15         # isotropic straight background tracks
  random_track_ionization: tuple = (1., 4.)  # x MIP grain density
  n_blobs: int = 20                 # large dark grain clusters (junk)
  n_dead_pixels: int = 200          # pixels dark in every slice
  blob_radius_um: tuple = (1.5, 5.)
  blob_grains_per_um3: float = 0.5


class DevelopedZ:
  """Piecewise map from original z to developed z (emulsion shrinks)."""

  def __init__(self, layers, shrink):
    self.zmin = np.array([l.zmin for l in layers])
    self.zmax = np.array([l.zmax for l in layers])
    self.emulsion = np.array([l.emulsion for l in layers])
    scale = np.where(self.emulsion, 1. / shrink, 1.)
    self.scale = scale
    thick = (self.zmax - self.zmin) * scale
    self.dev_min = np.concatenate([[0.], np.cumsum(thick)[:-1]])
    self.dev_max = self.dev_min + thick

  def _layer(self, z, lo, hi):
    i = np.searchsorted(hi, z, side="left")
    i = np.clip(i, 0, len(hi) - 1)
    inside = (z >= lo[i]) & (z <= hi[i])
    return i, inside

  def to_dev(self, z):
    i, _ = self._layer(z, self.zmin, self.zmax)
    return self.dev_min[i] + (z - self.zmin[i]) * self.scale[i]

  def in_emulsion(self, z):
    i, inside = self._layer(z, self.zmin, self.zmax)
    return inside & self.emulsion[i]

  def dev_in_emulsion(self, zdev):
    i, inside = self._layer(zdev, self.dev_min, self.dev_max)
    return inside & self.emulsion[i]


class View:
  """Pixel / slice coordinates of one image stack.

  (cx, cy) is the view center in um and czdev the developed depth of the
  middle of the stack."""

  def __init__(self, cx, cy, czdev, devz, optics):
    self.optics = optics
    self.devz = devz
    self.cx, self.cy, self.czdev = cx, cy, czdev

  @classmethod
  def around(cls, vertex, devz, optics, rng, margin_px=0):
    """View containing `vertex` at a random position.

    The vertex is placed uniformly in x and y (keeping margin_px from the
    edges) and the stack is shifted in depth at random but kept inside
    the emulsion layer of the vertex when the layer is thick enough."""
    o = optics
    half_um = (o.size_px / 2. - margin_px) * o.pixel_um
    cx = vertex[0] + rng.uniform(-half_um, half_um)
    cy = vertex[1] + rng.uniform(-half_um, half_um)
    vdev = float(devz.to_dev(np.array([vertex[2]]))[0])
    half_stack = (o.n_slices - 1) / 2. * o.dz_um
    i, _ = devz._layer(np.array([vdev]), devz.dev_min, devz.dev_max)
    lay_lo, lay_hi = devz.dev_min[i[0]], devz.dev_max[i[0]]
    # Centers keeping the vertex in the stack ...
    lo, hi = vdev - half_stack, vdev + half_stack
    # ... and the stack inside the layer, if possible
    lo_in, hi_in = max(lo, lay_lo + half_stack), min(hi, lay_hi - half_stack)
    if lo_in <= hi_in:
      lo, hi = lo_in, hi_in
    return cls(cx, cy, rng.uniform(lo, hi), devz, optics)

  def to_view(self, points):
    """(n, 3) original um -> (n, 3) column, row, slice (float)."""
    o = self.optics
    col = (points[:, 0] - self.cx) / o.pixel_um + o.size_px / 2.
    row = (points[:, 1] - self.cy) / o.pixel_um + o.size_px / 2.
    zdev = self.devz.to_dev(points[:, 2])
    sl = (zdev - self.czdev) / o.dz_um + (o.n_slices - 1) / 2.
    return np.column_stack([col, row, sl])

  def slice_zdev(self, k):
    o = self.optics
    return self.czdev + (k - (o.n_slices - 1) / 2.) * o.dz_um

  def zdev_range(self):
    """Developed depth covered by the stack including defocus margin."""
    o = self.optics
    return (self.slice_zdev(0) - o.max_defocus_um,
            self.slice_zdev(o.n_slices - 1) + o.max_defocus_um)


def track_grains(points, edep, emulsion, rng):
  """Grain positions (original um) along one track."""
  a, b = points[:-1], points[1:]
  seg = b - a
  length = np.linalg.norm(seg, axis=1)
  de = edep[1:]
  ok = length > 0.
  a, seg, length, de = a[ok], seg[ok], length[ok], de[ok]
  if len(length) == 0:
    return np.zeros((0, 3))
  dedx = de / length
  e = emulsion
  density = e.grains_sat_per_um * (1. - np.exp(
    -e.grains_mip_per_um / e.grains_sat_per_um
    * dedx / e.dedx_mip_mev_per_um))
  n = rng.poisson(density * length)
  idx = np.repeat(np.arange(len(n)), n)
  t = rng.random(len(idx))
  grains = a[idx] + t[:, None] * seg[idx]
  return grains + rng.normal(0., e.grain_jitter_um, grains.shape)


def straight_grains(view, emulsion, vertex_z, n_tracks, direction_fn,
                    density_fn, rng):
  """Grains of straight tracks crossing the view (original um) and the
  index of the track each grain belongs to."""
  o = view.optics
  half = o.size_px * o.pixel_um / 2.
  lo, hi = view.zdev_range()
  span = (hi - lo) * o.shrink  # original-z extent, generous
  grains, owners = [], []
  for i in range(n_tracks):
    mid = np.array([view.cx + rng.uniform(-half, half),
                    view.cy + rng.uniform(-half, half),
                    vertex_z + rng.uniform(-span / 2., span / 2.)])
    d = direction_fn()
    # Long enough to cross the whole view in any direction
    length = 2. * np.hypot(2. * half, span)
    p0 = mid - d * length / 2.
    n = rng.poisson(density_fn() * length)
    g = p0 + rng.random(n)[:, None] * d * length
    grains.append(g + rng.normal(0., emulsion.grain_jitter_um, g.shape))
    owners.append(np.full(n, i))
  if not grains:
    return np.zeros((0, 3)), np.zeros(0, dtype=int)
  return np.vstack(grains), np.concatenate(owners)


def beam_grains(view, background, emulsion, vertex_z, rng):
  """Nearly parallel MIP tracks around a common beam axis close to z."""
  axis = rng.uniform(-background.beam_slope_max,
                     background.beam_slope_max, 2)

  def direction():
    s = axis + rng.normal(0., background.beam_spread, 2)
    d = np.array([s[0], s[1], 1.])
    return d / np.linalg.norm(d)
  return straight_grains(view, emulsion, vertex_z,
                         background.n_beam_tracks, direction,
                         lambda: emulsion.grains_mip_per_um, rng)


def random_track_grains(view, background, emulsion, vertex_z, rng):
  """Isotropic straight tracks of various ionization."""
  def direction():
    d = rng.normal(size=3)
    return d / np.linalg.norm(d)
  lo, hi = background.random_track_ionization
  return straight_grains(
    view, emulsion, vertex_z, background.n_random_tracks, direction,
    lambda: emulsion.grains_mip_per_um * rng.uniform(lo, hi), rng)


def blob_grains(view, background, rng):
  """Large dark clusters of grains (dust, crystals), in view coords, and
  the index of the cluster of each grain."""
  o = view.optics
  lo, hi = view.zdev_range()
  out, owners = [], []
  for i in range(background.n_blobs):
    r = rng.uniform(*background.blob_radius_um)
    n = rng.poisson(background.blob_grains_per_um3 * 4. / 3. * np.pi
                    * r ** 3)
    d = rng.normal(size=(n, 3))
    d *= (r * rng.random(n) ** (1. / 3.)
          / np.linalg.norm(d, axis=1))[:, None]
    center = np.array([rng.uniform(0., o.size_px * o.pixel_um),
                       rng.uniform(0., o.size_px * o.pixel_um),
                       rng.uniform(lo, hi)])
    p = center + d
    out.append(np.column_stack([
      p[:, 0] / o.pixel_um, p[:, 1] / o.pixel_um,
      (p[:, 2] - view.czdev) / o.dz_um + (o.n_slices - 1) / 2.]))
    owners.append(np.full(n, i))
  if not out:
    return np.zeros((0, 3)), np.zeros(0, dtype=int)
  return np.vstack(out), np.concatenate(owners)


def grain_weights(n, emulsion, rng):
  """Relative darkness of grains (lognormal, mean 1)."""
  s = emulsion.grain_spread
  return rng.lognormal(-s * s / 2., s, n)


def fog_grains(view, emulsion, rng):
  """Uniform fog in the developed emulsion of the view (view coords)."""
  o = view.optics
  lo, hi = view.zdev_range()
  volume = (o.size_px * o.pixel_um) ** 2 * (hi - lo)
  n = rng.poisson(emulsion.fog_per_um3 * volume)
  zdev = rng.uniform(lo, hi, n)
  keep = view.devz.dev_in_emulsion(zdev)
  col = rng.uniform(0., o.size_px, n)[keep]
  row = rng.uniform(0., o.size_px, n)[keep]
  sl = (zdev[keep] - view.czdev) / o.dz_um + (o.n_slices - 1) / 2.
  return np.column_stack([col, row, sl])


def mottle_field(optics, background, rng):
  """Slowly varying brightness B(x, y) as smoothed noise of two scales."""
  n = optics.size_px
  fy = np.fft.fftfreq(n)[:, None]
  fx = np.fft.rfftfreq(n)[None, :]
  total = np.ones((n, n))
  for rel, scale_um in background.mottle:
    sigma_px = scale_um / optics.pixel_um
    kernel = np.exp(-2. * np.pi ** 2 * sigma_px ** 2
                    * (fx ** 2 + fy ** 2))
    field = np.fft.irfft2(np.fft.rfft2(rng.normal(size=(n, n))) * kernel,
                          s=(n, n))
    total += rel * field / (field.std() or 1.)
  return background.brightness * total


class Renderer:
  def __init__(self, optics):
    self.o = optics
    n = optics.size_px + 2 * optics.pad_px
    self.shape = (n, n)
    fy = np.fft.fftfreq(n)[:, None]
    fx = np.fft.rfftfreq(n)[None, :]
    self.f2 = fx ** 2 + fy ** 2
    edges = np.linspace(0., optics.max_defocus_um,
                        optics.n_defocus_bins + 1)
    self.bin_edges = edges
    s0 = optics.grain_sigma_um / optics.pixel_um
    self.kernels = []
    for dz in 0.5 * (edges[:-1] + edges[1:]):
      s = np.hypot(optics.grain_sigma_um, optics.defocus_slope * dz)
      s /= optics.pixel_um
      # FT of d0 (s0/s)^2 exp(-r^2/2s^2): d0 2 pi s0^2 exp(-2 pi^2 s^2 f^2)
      self.kernels.append(optics.grain_density * 2. * np.pi * s0 ** 2
                          * np.exp(-2. * np.pi ** 2 * s ** 2 * self.f2))

  def _splat(self, col, row, weight):
    """Bilinear deposit of weights on the padded canvas."""
    n = self.shape[0]
    x = col + self.o.pad_px
    y = row + self.o.pad_px
    x0, y0 = np.floor(x).astype(int), np.floor(y).astype(int)
    fx, fy = x - x0, y - y0
    idx, val = [], []
    for dx, dy, w in ((0, 0, (1 - fx) * (1 - fy)), (1, 0, fx * (1 - fy)),
                      (0, 1, (1 - fx) * fy), (1, 1, fx * fy)):
      xi, yi = x0 + dx, y0 + dy
      ok = (xi >= 0) & (xi < n) & (yi >= 0) & (yi < n)
      idx.append(yi[ok] * n + xi[ok])
      val.append((w * weight)[ok])
    img = np.bincount(np.concatenate(idx), np.concatenate(val),
                      minlength=n * n)
    return img.reshape(self.shape)

  def density(self, grains, weights, k):
    """Optical density of slice k from grains (col, row, slice) with
    relative darkness weights."""
    o = self.o
    dz = np.abs(grains[:, 2] - k) * o.dz_um
    sel = dz < o.max_defocus_um
    g, w, dz = grains[sel], weights[sel], dz[sel]
    bins = np.digitize(dz, self.bin_edges) - 1
    spectrum = np.zeros((self.shape[0], self.shape[1] // 2 + 1),
                        dtype=complex)
    for b in np.unique(bins):
      m = bins == b
      spectrum += np.fft.rfft2(self._splat(g[m, 0], g[m, 1], w[m])) \
        * self.kernels[b]
    d = np.fft.irfft2(spectrum, s=self.shape)
    p = o.pad_px
    return d[p:p + o.size_px, p:p + o.size_px]

  def labels(self, grains, classes, instances, k):
    """Class and instance id of in-focus grains (within half a slice),
    painted as disks; lower classes win where grains overlap."""
    o = self.o
    cls = np.zeros((o.size_px, o.size_px), dtype=np.uint8)
    ins = np.zeros((o.size_px, o.size_px), dtype=np.int32)
    sel = np.abs(grains[:, 2] - k) < 0.5
    if not sel.any():
      return cls, ins
    col = np.rint(grains[sel, 0]).astype(int)
    row = np.rint(grains[sel, 1]).astype(int)
    c, d = classes[sel], instances[sel]
    r = o.label_radius_px
    xs, ys, cs, ds = [], [], [], []
    for dy in range(-r, r + 1):
      for dx in range(-r, r + 1):
        if dx * dx + dy * dy <= r * r:
          xs.append(col + dx)
          ys.append(row + dy)
          cs.append(c)
          ds.append(d)
    x, y = np.concatenate(xs), np.concatenate(ys)
    c, d = np.concatenate(cs), np.concatenate(ds)
    ok = (x >= 0) & (x < o.size_px) & (y >= 0) & (y < o.size_px)
    # Paint high classes first so that lower classes overwrite them
    order = np.argsort(-c[ok], kind="stable")
    cls[y[ok][order], x[ok][order]] = c[ok][order]
    ins[y[ok][order], x[ok][order]] = d[ok][order]
    return cls, ins
