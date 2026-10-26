"""Minimum-ionizing signal vs the real beam tracks of the same images.

  ~/work/.venv/bin/python -m module.mipcheck data/images/<single-track
    dataset> [--events 40] [-o output/mipcheck_v20.md]

The beam (K- / pi- at ~1.8 GeV/c, minimum ionizing) crosses the plates
nearly perpendicularly, so every real background view carries MIP
tracks seen as dark dots persisting through the layer. Steep embedded
signal tracks of minimum-ionizing particles (dip above STEEP_DEG,
dE/dx below MIP_DEDX_MAX) should look the same. In each slice inside
the emulsion the background-subtracted darkness is taken as the maximum
within +-WINDOW_PX of the track position: for the signal the position
of its in-focus grains (truth), for a real dot the position followed
from slice to slice (start at the darkest persistent pixel, move to the
darkest pixel of the window when it is dark). Per track the median
darkness and the persistence (fraction of slices darker than
DARK_LEVEL) are compared. Real dots are selected by persistence, so the
signal is compared both before and after the same selection; real
dots near any signal grain are excluded. Following drifts towards dark
pixels, so the signal is also measured followed from its truth position
in the middle slice, and followed random positions give the level of
the background alone.

Result (single_eval_v20_dz4, output/mipcheck_v20.md): inconclusive. At
the truth position a minimum-ionizing signal track is only ~3 gray
levels darker than fixed random positions (median 14 vs 11), and the
real beam tracks cannot be isolated this way: 68 % of fixed and 85 % of
followed random positions pass the persistence cut, so the real
"dots" are the darkest texture rather than a beam sample.
"""

import argparse
import csv
from pathlib import Path

import numpy as np

from . import imagestats as st

CROP_PX = 1024
WINDOW_PX = 2
DARK_LEVEL = 6.          # gray levels below the local background
LAYER_LEVEL = 15.        # dark pixels defining slices inside the emulsion
EMPTY_FRACTION = 0.1     # of the maximum: slice outside the emulsion
PERSISTENCE = 0.6
MIN_SLICES = 10
MIN_DOT_SEPARATION_PX = 15
SIGNAL_EXCLUSION_PX = 15
FOCUS_SLICES = 1         # grains within this many slices give the position
STEEP_DEG = 60.
MIP_DEDX_MAX = 0.0005    # MeV/um, about 1.8 x the lowest value of the gun
DEDX_BINS = (0., MIP_DEDX_MAX, 0.0015, 1.)
PCT = (16, 50, 84)
N_RANDOM = 200           # followed random positions per event
RANDOM_SEED = 0


def darkness_stack(path):
  out = []
  for k in range(st.n_slices(path)):
    img = st.read_stack(path, [k])[0]
    lo = (img.shape[0] - CROP_PX) // 2
    c = img[lo:lo + CROP_PX, lo:lo + CROP_PX].astype(np.float32)
    out.append(st.gaussian_blur(c, st.BG_SIGMA_PX) - c)
  return np.array(out, np.float32), lo


def window_max(d, y, x):
  return d[max(0, y - WINDOW_PX):y + WINDOW_PX + 1,
           max(0, x - WINDOW_PX):x + WINDOW_PX + 1].max()


def track_profile(det, inside, cols, rows, slices):
  """Darkness per inside slice at the in-focus grains of a track."""
  vals = []
  for k in inside:
    m = np.abs(slices - k) <= FOCUS_SLICES
    if not m.any():
      continue
    x, y = int(round(cols[m].mean())), int(round(rows[m].mean()))
    if WINDOW_PX <= x < CROP_PX - WINDOW_PX \
       and WINDOW_PX <= y < CROP_PX - WINDOW_PX:
      vals.append(window_max(det[k], y, x))
  return np.array(vals)


def follow_profile(det, inside, y0, x0):
  """Darkness per inside slice following a dot from the middle slice."""
  mid = len(inside) // 2
  vals = {}
  for order in (inside[mid:], inside[:mid][::-1]):
    y, x = y0, x0
    for k in order:
      lo_y, lo_x = max(0, y - WINDOW_PX), max(0, x - WINDOW_PX)
      win = det[k, lo_y:y + WINDOW_PX + 1, lo_x:x + WINDOW_PX + 1]
      dy, dx = np.unravel_index(np.argmax(win), win.shape)
      vals[k] = win[dy, dx]
      if vals[k] > DARK_LEVEL:
        y = int(np.clip(lo_y + dy, WINDOW_PX, CROP_PX - WINDOW_PX - 1))
        x = int(np.clip(lo_x + dx, WINDOW_PX, CROP_PX - WINDOW_PX - 1))
  return np.array([vals[k] for k in inside])


def real_dots(det, inside, exclude):
  sub = det[inside]
  persistent = (sub > DARK_LEVEL).mean(axis=0)
  mean = sub.mean(axis=0)
  ok = (persistent >= PERSISTENCE) & ~exclude
  ok[:WINDOW_PX + 1], ok[-WINDOW_PX - 1:] = False, False
  ok[:, :WINDOW_PX + 1], ok[:, -WINDOW_PX - 1:] = False, False
  ys, xs = np.nonzero(ok)
  order = np.argsort(mean[ys, xs])[::-1]
  chosen = np.empty((0, 2), int)
  for y, x in zip(ys[order], xs[order]):
    if len(chosen) and (((chosen - (y, x)) ** 2).sum(axis=1)
                        < MIN_DOT_SEPARATION_PX ** 2).any():
      continue
    chosen = np.vstack([chosen, (y, x)])
  return chosen


def track_summary(v):
  return float(np.median(v)), float((v > DARK_LEVEL).mean())


def event_rows(ev_dir, tracks):
  det, lo = darkness_stack(ev_dir / "image.json")
  prof = (det > LAYER_LEVEL).mean(axis=(1, 2))
  inside = np.flatnonzero(prof > EMPTY_FRACTION * prof.max())
  g = np.load(ev_dir / "grains.npz")
  sig = g["cls"] == 1
  cols, rows = g["col"] - lo, g["row"] - lo
  out = []
  for t in tracks:
    m = sig & (g["instance"] == int(t["track_id"]))
    fs = g["focus_slice"][m]
    v = track_profile(det, inside, cols[m], rows[m], fs)
    if len(v) < MIN_SLICES:
      continue
    rec = dict(kind="signal", dedx=float(t["dedx_mev_per_um"]),
               dip=abs(float(t["dip_deg"])), n=len(v),
               **dict(zip(("dark", "persist"), track_summary(v))))
    # Followed from the in-focus position at the middle slice
    k = inside[len(inside) // 2]
    near = np.abs(fs - k) <= FOCUS_SLICES
    if near.any():
      y = int(round(rows[m][near].mean()))
      x = int(round(cols[m][near].mean()))
      if WINDOW_PX <= x < CROP_PX - WINDOW_PX \
         and WINDOW_PX <= y < CROP_PX - WINDOW_PX:
        f = follow_profile(det, inside, y, x)
        rec.update(zip(("dark_f", "persist_f"), track_summary(f)))
    out.append(rec)
  exclude = np.zeros((CROP_PX, CROP_PX), bool)
  r = np.rint(rows[sig]).astype(int)
  c = np.rint(cols[sig]).astype(int)
  ok = (r >= 0) & (r < CROP_PX) & (c >= 0) & (c < CROP_PX)
  exclude[r[ok], c[ok]] = True
  exclude = st.gaussian_blur(exclude.astype(float),
                             SIGNAL_EXCLUSION_PX / 2.) > 1e-3
  for y, x in real_dots(det, inside, exclude):
    v = follow_profile(det, inside, int(y), int(x))
    out.append(dict(kind="real", n=len(v),
                    **dict(zip(("dark", "persist"), track_summary(v)))))
  rng = np.random.default_rng(RANDOM_SEED)
  for y, x in rng.integers(WINDOW_PX, CROP_PX - WINDOW_PX, (N_RANDOM, 2)):
    v = follow_profile(det, inside, int(y), int(x))
    out.append(dict(kind="random", n=len(v),
                    **dict(zip(("dark", "persist"), track_summary(v)))))
    # Same position in every slice, as the truth position of the signal
    v = np.array([window_max(det[k], int(y), int(x)) for k in inside])
    out.append(dict(kind="fixed", n=len(v),
                    **dict(zip(("dark", "persist"), track_summary(v)))))
  return out


def pct(rows, key):
  if not rows:
    return "-"
  return " / ".join(f"{np.percentile([r[key] for r in rows], p):.3g}"
                    for p in PCT)


def main():
  parser = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument("dataset")
  parser.add_argument("--events", type=int, default=40)
  parser.add_argument("-o", "--output")
  args = parser.parse_args()
  ds = Path(args.dataset)
  with open(ds / "single_tracks.csv") as f:
    table = [t for t in csv.DictReader(f)
             if abs(float(t["dip_deg"])) > STEEP_DEG]
  rows = []
  for ev_dir in sorted(ds.glob("ev*"))[:args.events]:
    ev = int(ev_dir.name[2:])
    rows += event_rows(ev_dir, [t for t in table if int(t["event"]) == ev])
    print(ev_dir.name, len(rows), flush=True)
  real = [r for r in rows if r["kind"] == "real"]
  rand = [r for r in rows if r["kind"] == "random"]
  fixed = [r for r in rows if r["kind"] == "fixed"]
  lines = [f"# MIP signal vs real beam dots in {ds.name}", "",
           f"Steep tracks (dip > {STEEP_DEG:g} deg), {args.events} events; "
           f"per track the median darkness (max within +-{WINDOW_PX} px, "
           "gray levels) and the persistence (fraction of slices darker "
           f"than {DARK_LEVEL:g}); cells: {PCT[0]} / {PCT[1]} / {PCT[2]} %. "
           "Real dots are selected by persistence >= "
           f"{PERSISTENCE:g}; 'passing' applies the same cut to the "
           "signal; 'followed' measures the signal like the real dots "
           "(followed from its position in the middle slice).", "",
           "| sample | n | darkness | persistence | passing | darkness "
           "(passing) | followed: darkness, persistence |",
           "|---|---|---|---|---|---|---|",
           f"| real persistent dots | {len(real)} | {pct(real, 'dark')} | "
           f"{pct(real, 'persist')} | - | - | - |",
           f"| followed random positions | {len(rand)} | "
           f"{pct(rand, 'dark')} | {pct(rand, 'persist')} | "
           f"{sum(r['persist'] >= PERSISTENCE for r in rand) / len(rand):.0%}"
           " | - | - |",
           f"| fixed random positions | {len(fixed)} | "
           f"{pct(fixed, 'dark')} | {pct(fixed, 'persist')} | "
           f"{np.mean([r['persist'] >= PERSISTENCE for r in fixed]):.0%}"
           " | - | - |"]
  for lo, hi in zip(DEDX_BINS[:-1], DEDX_BINS[1:]):
    sel = [r for r in rows if r["kind"] == "signal"
           and lo <= r["dedx"] < hi]
    ok = [r for r in sel if r["persist"] >= PERSISTENCE]
    frac = f"{len(ok) / len(sel):.0%}" if sel else "-"
    fol = [r for r in sel if "dark_f" in r]
    lines.append(f"| signal dE/dx {lo:g}-{hi:g} MeV/um | {len(sel)} | "
                 f"{pct(sel, 'dark')} | {pct(sel, 'persist')} | {frac} | "
                 f"{pct(ok, 'dark')} | {pct(fol, 'dark_f')}, "
                 f"{pct(fol, 'persist_f')} |")
  text = "\n".join(lines) + "\n"
  if args.output:
    Path(args.output).write_text(text)
  print(text)


if __name__ == "__main__":
  main()
