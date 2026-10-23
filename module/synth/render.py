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
  grain_sigma_um: float = 0.5  # in-focus grain image (s0)
  defocus_slope: float = 0.08  # k: blur growth per um of defocus
  max_defocus_cap_um: float = 40.  # upper bound of the defocus range
  n_defocus_bins: int = 12
  image_blur_px: float = 0.    # optics / camera MTF on the intensity
  grain_density: float = 0.6   # d0: optical density of a focused grain
  pad_px: int = 64             # margin against FFT wrap-around
  label_radius_px: int = 2


# A grain is ignored where its peak density drops below this fraction
DEFOCUS_AMPLITUDE_CUT = 0.03


def max_defocus_um(optics):
  """Defocus beyond which grains are negligible: (s0/s)^2 < cut."""
  ratio = np.sqrt(1. / DEFOCUS_AMPLITUDE_CUT - 1.)
  return min(optics.max_defocus_cap_um,
             ratio * optics.grain_sigma_um / optics.defocus_slope)


@dataclass(frozen=True)
class Emulsion:
  grains_mip_per_um: float = 0.148  # MIP grain density (fullscan STATUS)
  grains_sat_per_um: float = 2.5    # saturation for heavy ionization
  dedx_mip_mev_per_um: float = 5.47e-4  # 1.55 MeV cm2/g x 3.53 g/cm3
  fog_per_um3: float = 2e-3         # fog grains per developed um3
  grain_jitter_um: float = 0.1
  grain_spread: float = 0.4         # lognormal sigma of grain darkness


@dataclass(frozen=True)
class Background:
  brightness: float = 210.
  # Relative rms and length scale of the slow brightness variations
  mottle: tuple = ((0.08, 25.), (0.06, 6.))
  # Depth correlation length of the mottle (developed um, 0: static)
  mottle_z_um: float = 0.
  min_transmission: float = 0.2     # darkest level relative to B
  noise_sigma: float = 2.           # temporal pixel noise (per slice)
  fixed_noise_sigma: float = 3.     # fixed-pattern noise (all slices)
  # Densities per mm2 of view area (counts scale with the view size)
  beam_per_mm2: float = 1100.       # nearly parallel MIP tracks
  beam_slope_max: float = 0.1       # |dx/dz|, |dy/dz| of the beam axis
  beam_spread: float = 0.01         # rms slope around the beam axis
  random_per_mm2: float = 170.      # isotropic straight tracks
  random_track_ionization: tuple = (1., 4.)  # x MIP grain density
  # Straight tracks nearly in the sheet plane (|dz/dl| < flat_dip_max)
  flat_per_mm2: float = 50.
  flat_dip_max: float = 0.15
  flat_ionization: tuple = (1., 3.)  # x MIP grain density
  # Static fine texture (gelatin), smoothed noise common to all slices
  speckle_rel: float = 0.0
  speckle_scale_px: float = 1.5
  blobs_per_mm2: float = 230.       # large dark grain clusters (junk)
  n_dead_pixels: int = 0            # pixels dark in every slice (the
                                    # tohoku-v1 dead-pixel map is empty)
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
    m = max_defocus_um(o)
    return (self.slice_zdev(0) - m, self.slice_zdev(o.n_slices - 1) + m)


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


# Straight background tracks have this fixed length and are placed in
# the view widened by half of it, so that the grain density inside the
# view does not depend on the view size
STRAIGHT_TRACK_UM = 1000.
DEPTH_REF_UM = 100.  # depth unit of the flat / isotropic track densities


def view_count(view, per_mm2, rng, margin_um=0.):
  """Poisson count of objects with a density per mm2 in the view widened
  by margin_um on each side."""
  o = view.optics
  width_mm = (o.size_px * o.pixel_um + 2. * margin_um) * 1e-3
  return int(rng.poisson(per_mm2 * width_mm ** 2))


def straight_grains(view, emulsion, vertex_z, per_mm2, direction_fn,
                    density_fn, rng, per_depth=False):
  """Grains of straight tracks of fixed length around the view (original
  um) and the index of the track each grain belongs to.

  per_mm2 counts tracks per mm2 of view area; with per_depth it is per
  mm2 and per DEPTH_REF_UM of (undeveloped) depth, for tracks spread in
  depth (flat, isotropic) rather than crossing the whole stack (beam)."""
  o = view.optics
  half = o.size_px * o.pixel_um / 2.
  lo, hi = view.zdev_range()
  span = (hi - lo) * o.shrink  # original-z extent of the view, generous
  length = max(STRAIGHT_TRACK_UM, 1.2 * span)
  margin = length / 2.
  scale = span / DEPTH_REF_UM if per_depth else 1.
  n_tracks = view_count(view, per_mm2 * scale, rng, margin)
  grains, owners = [], []
  for i in range(n_tracks):
    mid = np.array([view.cx + rng.uniform(-half - margin, half + margin),
                    view.cy + rng.uniform(-half - margin, half + margin),
                    vertex_z + rng.uniform(-span / 2., span / 2.)])
    d = direction_fn()
    p0 = mid - d * length / 2.
    n = rng.poisson(density_fn() * length)
    g = p0 + rng.random(n)[:, None] * d * length
    # Keep the grains near the view only (memory)
    near = (np.abs(g[:, 0] - view.cx) < half + o.pad_px * o.pixel_um) \
      & (np.abs(g[:, 1] - view.cy) < half + o.pad_px * o.pixel_um)
    g = g[near]
    grains.append(g + rng.normal(0., emulsion.grain_jitter_um, g.shape))
    owners.append(np.full(len(g), i))
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
                         background.beam_per_mm2, direction,
                         lambda: emulsion.grains_mip_per_um, rng)


def random_track_grains(view, background, emulsion, vertex_z, rng):
  """Isotropic straight tracks of various ionization."""
  def direction():
    d = rng.normal(size=3)
    return d / np.linalg.norm(d)
  lo, hi = background.random_track_ionization
  return straight_grains(
    view, emulsion, vertex_z, background.random_per_mm2, direction,
    lambda: emulsion.grains_mip_per_um * rng.uniform(lo, hi), rng,
    per_depth=True)


def flat_track_grains(view, background, emulsion, vertex_z, rng):
  """Straight tracks nearly parallel to the sheets, any azimuth."""
  def direction():
    a = rng.uniform(0., 2. * np.pi)
    dip = rng.uniform(-background.flat_dip_max, background.flat_dip_max)
    d = np.array([np.cos(a), np.sin(a), dip])
    return d / np.linalg.norm(d)
  lo, hi = background.flat_ionization
  return straight_grains(
    view, emulsion, vertex_z, background.flat_per_mm2, direction,
    lambda: emulsion.grains_mip_per_um * rng.uniform(lo, hi), rng,
    per_depth=True)


def blob_grains(view, background, rng):
  """Large dark clusters of grains (dust, crystals), in view coords, and
  the index of the cluster of each grain."""
  o = view.optics
  lo, hi = view.zdev_range()
  out, owners = [], []
  for i in range(view_count(view, background.blobs_per_mm2, rng)):
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


def _mottle_layer(optics, background, rng):
  """Relative brightness variation (mean 0) of the two mottle scales."""
  n = optics.size_px
  fy = np.fft.fftfreq(n)[:, None]
  fx = np.fft.rfftfreq(n)[None, :]
  total = np.zeros((n, n))
  for rel, scale_um in background.mottle:
    sigma_px = scale_um / optics.pixel_um
    kernel = np.exp(-2. * np.pi ** 2 * sigma_px ** 2
                    * (fx ** 2 + fy ** 2))
    field = np.fft.irfft2(np.fft.rfft2(rng.normal(size=(n, n))) * kernel,
                          s=(n, n))
    total += rel * field / (field.std() or 1.)
  return total


def mottle_field(optics, background, rng):
  """Static slowly varying brightness B(x, y)."""
  return background.brightness * (1. + _mottle_layer(optics, background,
                                                     rng))


class MottleVolume:
  """Brightness B(x, y, z) whose mottle decorrelates over mottle_z_um in
  developed depth: independent layers at nodes spaced mottle_z_um,
  interpolated with variance-preserving weights."""

  def __init__(self, optics, background, rng, zdev_min, zdev_max):
    self.b = background
    self.lz = background.mottle_z_um
    if self.lz <= 0.:
      self.static = mottle_field(optics, background, rng)
      return
    self.z0 = zdev_min
    n_nodes = int(np.ceil((zdev_max - zdev_min) / self.lz)) + 2
    self.layers = [_mottle_layer(optics, background, rng)
                   for _ in range(n_nodes)]

  def at(self, zdev):
    if self.lz <= 0.:
      return self.static
    t = (zdev - self.z0) / self.lz
    i = int(np.clip(np.floor(t), 0, len(self.layers) - 2))
    f = t - i
    w0, w1 = 1. - f, f
    norm = np.hypot(w0, w1)
    rel = (w0 * self.layers[i] + w1 * self.layers[i + 1]) / norm
    return self.b.brightness * (1. + rel)


def fixed_pattern(optics, background, rng):
  """Per-pixel noise common to all slices of a view (additive, gray
  levels) and the static fine texture (relative, smoothed)."""
  n = optics.size_px
  fixed = rng.normal(0., background.fixed_noise_sigma, (n, n))
  if background.speckle_rel > 0.:
    fy = np.fft.fftfreq(n)[:, None]
    fx = np.fft.rfftfreq(n)[None, :]
    k = np.exp(-2. * np.pi ** 2 * background.speckle_scale_px ** 2
               * (fx ** 2 + fy ** 2))
    sp = np.fft.irfft2(np.fft.rfft2(rng.normal(size=(n, n))) * k, s=(n, n))
    sp *= background.speckle_rel * background.brightness / (sp.std() or 1.)
    fixed = fixed + sp
  return fixed


def gaussian_blur(img, sigma_px):
  n, m = img.shape
  fy = np.fft.fftfreq(n)[:, None]
  fx = np.fft.rfftfreq(m)[None, :]
  k = np.exp(-2. * np.pi ** 2 * sigma_px ** 2 * (fx ** 2 + fy ** 2))
  return np.fft.irfft2(np.fft.rfft2(img) * k, s=img.shape)


def compose(bright, density, background, fixed, rng, blur_px=0.):
  """Slice intensity from brightness, optical density and noise; the
  transmitted image is blurred by the optics / camera (blur_px) before
  the sensor noise is added."""
  t_min = background.min_transmission
  trans = np.exp(-density)
  if blur_px > 0.:
    trans = gaussian_blur(trans, blur_px)
  img = bright * (t_min + (1. - t_min) * trans) + fixed \
    + rng.normal(0., background.noise_sigma, density.shape)
  return np.clip(np.rint(img), 0, 255).astype(np.uint8)


class Renderer:
  def __init__(self, optics):
    self.o = optics
    n = optics.size_px + 2 * optics.pad_px
    self.shape = (n, n)
    fy = np.fft.fftfreq(n)[:, None]
    fx = np.fft.rfftfreq(n)[None, :]
    self.f2 = fx ** 2 + fy ** 2
    self.max_defocus = max_defocus_um(optics)
    # Bins equally spaced in blur width
    s_edges = np.linspace(optics.grain_sigma_um,
                          np.hypot(optics.grain_sigma_um,
                                   optics.defocus_slope * self.max_defocus),
                          optics.n_defocus_bins + 1)
    edges = np.sqrt(np.maximum(s_edges ** 2 - optics.grain_sigma_um ** 2,
                               0.)) / optics.defocus_slope
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
    sel = dz < self.max_defocus
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
