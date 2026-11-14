"""Image shift vs encoder step for neighbouring full-scan views."""
import glob, json
import numpy as np
from module import imagestats as st
D = ("/home/hayakawa/microscope-image/tohoku/fullscan/E07/MOD108/PL12/"
     "tohoku-v1/AREA00/IMAGE00_AREA00/")
def view(x, y):
  return glob.glob(D + f"V*_L0_VX{x:04d}_VY{y:04d}_0_058.json")[0]
def img(p):
  k = st.n_slices(p) // 2
  a = st.read_stack(p, [k])[0].astype(float)
  return st.gaussian_blur(a, 8) - a, json.load(open(p))["Images"][k]
def shift(a, b):
  f = np.fft.fft2(a - a.mean()) * np.conj(np.fft.fft2(b - b.mean()))
  r = np.fft.ifft2(f / np.maximum(np.abs(f), 1e-9)).real
  for s in ((slice(None, 20), slice(None, 20)), (slice(None, 20), slice(-20, None)),
            (slice(-20, None), slice(None, 20)), (slice(-20, None), slice(-20, None))):
    r[s] = 0
  k = np.unravel_index(np.argmax(r), r.shape)
  n = r.shape
  sub = []
  for i in (0, 1):
    e = [0, 0]; e[i] = 1
    m = r[(k[0] - e[0]) % n[0], (k[1] - e[1]) % n[1]]
    p = r[(k[0] + e[0]) % n[0], (k[1] + e[1]) % n[1]]
    c = r[k]
    d = 0.5 * (m - p) / (m - 2 * c + p)
    sub.append(((k[i] + n[i] // 2) % n[i]) - n[i] // 2 + d)
  return sub
pairs = [((0, 0), (1, 0)), ((1, 0), (2, 0)), ((5, 5), (6, 5)), ((20, 20), (21, 20))]
pairs += [((x, 0), (x, 1)) for x in (3, 10, 17)]
pairs += [((x, 10), (x, 11)) for x in (3, 10, 17)]
for p, q in pairs:
  a, ja = img(view(*p))
  b, jb = img(view(*q))
  s = shift(a, b)
  ax = 1 if p[1] == q[1] else 0
  enc = (jb["x"] - ja["x"]) if ax else (jb["y"] - ja["y"])
  pix = abs(s[ax])
  print(p, q, "enc step", round(enc, 4), "pixel shift", np.round(s, 2),
        "um/px", round(abs(enc) * 1000 / pix, 4), flush=True)
