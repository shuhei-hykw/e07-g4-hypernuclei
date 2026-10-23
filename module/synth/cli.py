"""Render microscope image stacks of simulated events.

  ./run.py images data/sim/nagara.root --events 0-4
  ./run.py images data/sim/nagara.root --events 3 --dz 1.5 --name nagara_dz1.5

For every event a directory data/images/<name>/ev<NNNN>/ is written:
  0000.png ...   8-bit slices, readable by e07-fullscan's SpngReader
  image.json     stack metadata in the full-scan JSON format (stage mm)
  truth.json     vertices, tracks and 3-D polylines in pixel / slice
                 coordinates (see README)
  labels.npz     per slice and pixel, for grains within half a slice of
                 the focal plane: `labels` class (1 capture products,
                 2 hypernuclear decay products, 3 other secondaries,
                 4 straight background tracks, 5 junk clusters, 0 none /
                 fog) and `instance` (track id; -1-i background track i,
                 -100001-i junk cluster i, 0 none)
  preview.png    min projection: crop around the vertex | whole view
  g4_tracks.png  Geant4 event display before imaging (./run.py draw)
  g4_view.png    Geant4 tracks in the image coordinates of the view

The capture point is placed at a random position in the view and the
stack at a random depth inside the emulsion layer (--centered disables).
"""

import argparse
import json
import multiprocessing
import sys
import time
from dataclasses import replace

import numpy as np

from .. import config
from . import render
from .. import draw_event
from .event import EventFile
from .pngio import write_gray_png
from .trackview import draw_view

CHAIN_PROCESSES = ("HypernucleusDecay", "hBertiniCaptureAtRest", "Decay")
PREVIEW_CROP_PX = 512
UM_PER_MM = 1000.
GRAY_MAX = 255
VERTEX_MARGIN_PX = 100
POLYLINE_STEP_UM = 2.
# Stage z of the first slice and its step (tohoku-v1: decreasing)
STAGE_Z0_MM = 0.
DEFAULT_STAGE_DZ_MM = -0.003
BACKGROUND_ID0 = -1
JUNK_ID0 = -100001


def parse_events(spec, n_events):
  """'3', '0-4' or '0,2,5' -> list of indices."""
  out = []
  for part in spec.split(","):
    if "-" in part:
      lo, hi = map(int, part.split("-"))
      out.extend(range(lo, hi + 1))
    else:
      out.append(int(part))
  return [i for i in out if 0 <= i < n_events]


def collect_grains(ev, view, devz, emulsion, background, rng):
  """Labelled grains (view coords) with class and instance, and fog."""
  o = view.optics
  lo = -o.pad_px
  hi = o.size_px + o.pad_px
  margin = o.max_defocus_um / o.dz_um + 1.

  def in_view(v):
    return ((v[:, 0] >= lo) & (v[:, 0] < hi) & (v[:, 1] >= lo)
            & (v[:, 1] < hi) & (v[:, 2] > -margin)
            & (v[:, 2] < o.n_slices - 1 + margin))

  grains, classes, instances, counts = [], [], [], {}

  def add(points_um, cls, ids):
    keep = devz.in_emulsion(points_um[:, 2])
    v = view.to_view(points_um[keep])
    ids = ids[keep]
    sel = in_view(v)
    grains.append(v[sel])
    classes.append(np.full(sel.sum(), cls))
    instances.append(ids[sel])
    return int(sel.sum())

  for t in ev.tracks:
    if t.charge == 0 or len(t.points) < 2:
      continue
    g = render.track_grains(t.points, t.edep, emulsion, rng)
    counts[t.id] = add(g, render.LABEL_OF_ORIGIN[t.origin],
                       np.full(len(g), t.id))
  beam, owner = render.beam_grains(view, background, emulsion,
                                   ev.vertex[2], rng)
  add(beam, render.LABEL_BEAM, BACKGROUND_ID0 - owner)
  rand, owner2 = render.random_track_grains(view, background, emulsion,
                                            ev.vertex[2], rng)
  offset = background.n_beam_tracks
  add(rand, render.LABEL_BEAM, BACKGROUND_ID0 - offset - owner2)
  blobs, owner3 = render.blob_grains(view, background, rng)
  grains.append(blobs)
  classes.append(np.full(len(blobs), render.LABEL_JUNK))
  instances.append(JUNK_ID0 - owner3)
  fog = render.fog_grains(view, emulsion, rng)
  return (np.vstack(grains), np.concatenate(classes),
          np.concatenate(instances), fog, counts)


def polyline(view, points):
  """Track points thinned to ~POLYLINE_STEP_UM, in view coordinates."""
  step = np.linalg.norm(np.diff(points, axis=0), axis=1)
  s = np.concatenate([[0.], np.cumsum(step)])
  keep = np.unique(np.concatenate([
    np.searchsorted(s, np.arange(0., s[-1], POLYLINE_STEP_UM)),
    [len(points) - 1]]))
  v = view.to_view(points[keep])
  return [[round(float(c), 2), round(float(r), 2), round(float(k), 3)]
          for c, r, k in v]


def truth_record(ev, view, counts):
  o = view.optics

  def point(p):
    v = view.to_view(p[None, :])[0]
    inside = (0 <= v[0] < o.size_px and 0 <= v[1] < o.size_px
              and 0 <= v[2] <= o.n_slices - 1)
    return {"col": round(float(v[0]), 2), "row": round(float(v[1]), 2),
            "slice": round(float(v[2]), 3), "in_view": bool(inside),
            "um": [round(float(c), 3) for c in p]}

  by_id = {t.id: t for t in ev.tracks}
  vertices = [dict(index=0, type="capture", parent="Xi-", parent_id=0,
                   children=[t.id for t in ev.tracks if t.parent == 0],
                   **point(ev.vertex))]
  start_vertex = {t.id: 0 for t in ev.tracks if t.parent == 0}
  end_vertex = {}
  groups = {}
  for t in ev.tracks:
    if t.parent > 0 and t.process in CHAIN_PROCESSES:
      groups.setdefault((t.parent, t.process), []).append(t)
  for (parent, process), children in sorted(groups.items()):
    index = len(vertices)
    vertices.append(dict(
      index=index, type=process,
      parent=by_id[parent].name if parent in by_id else "?",
      parent_id=parent, children=[c.id for c in children],
      **point(children[0].points[0])))
    end_vertex[parent] = index
    for c in children:
      start_vertex[c.id] = index

  tracks = []
  for t in ev.tracks:
    if t.charge == 0 or counts.get(t.id, 0) == 0:
      continue
    tracks.append({
      "id": t.id, "parent": t.parent, "name": t.name, "origin": t.origin,
      "process": t.process, "charge": t.charge, "a": t.a,
      "ke_mev": round(t.ke, 4), "range_um": round(t.range_um, 2),
      "escaped": t.escaped, "n_grains": counts[t.id],
      "vertex_start": start_vertex.get(t.id),
      "vertex_end": end_vertex.get(t.id),
      "start": point(t.points[0]), "end": point(t.points[-1]),
      "polyline": polyline(view, t.points)})
  return {"event": ev.index, "channel": ev.channel,
          "contained": ev.contained, "vertices": vertices,
          "tracks": tracks}


def view_json(view, stage_dz_mm, out_dir):
  o = view.optics
  pix_mm = o.pixel_um / UM_PER_MM
  x0 = view.cx / UM_PER_MM - o.size_px / 2. * pix_mm
  y0 = view.cy / UM_PER_MM - o.size_px / 2. * pix_mm
  images = [{"Path": f"{k:04d}.png", "x": round(view.cx / UM_PER_MM, 6),
             "y": round(view.cy / UM_PER_MM, 6),
             "z": round(STAGE_Z0_MM + k * stage_dz_mm, 6)}
            for k in range(o.n_slices)]
  meta = {
    "DocumentType": {"ver": "20180521", "unit": "mm",
                     "comment": "e07-g4-hypernuclei synthetic; physical "
                     f"focal step {o.dz_um} um per slice in developed "
                     "emulsion"},
    "Stage": {"Location": "synthetic", "ID": "emulsim", "Objective": "20"},
    "EmulsionType": {"Sh": o.shrink},
    "ImageType": {"Height": o.size_px, "Width": o.size_px, "Depth": 8},
    "AffineP2S": [pix_mm, 0., 0., pix_mm, x0, y0],
    "Images": images,
  }
  (out_dir / "image.json").write_text(json.dumps(meta, indent=1))


def preview(stack_min, vertex, out_dir):
  """Min projection: crop around the vertex next to the whole view."""
  n = stack_min.shape[0]
  half = PREVIEW_CROP_PX // 2
  c = np.clip(np.rint(vertex[:2]).astype(int), half, n - half)
  crop = stack_min[c[1] - half:c[1] + half, c[0] - half:c[0] + half]
  step = max(1, n // PREVIEW_CROP_PX)
  whole = stack_min[::step, ::step][:crop.shape[0], :crop.shape[1]]
  write_gray_png(out_dir / "preview.png", np.hstack([crop, whole]))


def render_event(ev, layers, args, optics, emulsion, background, out_dir):
  rng = np.random.default_rng(args.seed + ev.index)
  devz = render.DevelopedZ(layers, optics.shrink)
  if args.centered:
    vdev = float(devz.to_dev(np.array([ev.vertex[2]]))[0])
    view = render.View(ev.vertex[0], ev.vertex[1], vdev, devz, optics)
  else:
    view = render.View.around(ev.vertex, devz, optics, rng,
                              VERTEX_MARGIN_PX)
  out_dir.mkdir(parents=True, exist_ok=True)
  draw_event.draw(args.root_file, ev.index, out_dir / "g4_tracks.png")
  draw_view(ev, view, out_dir / "g4_view.png")
  if args.tracks_only:
    return
  labelled, classes, instances, fog, counts = collect_grains(
    ev, view, devz, emulsion, background, rng)
  all_grains = np.vstack([labelled, fog])
  weights = render.grain_weights(len(all_grains), emulsion, rng)
  renderer = render.Renderer(optics)
  bright = render.mottle_field(optics, background, rng)
  dead = rng.integers(0, optics.size_px, (background.n_dead_pixels, 2))

  stack_min = np.full((optics.size_px, optics.size_px), GRAY_MAX,
                      dtype=np.uint8)
  cls_stack, ins_stack = [], []
  t_min = background.min_transmission
  for k in range(optics.n_slices):
    d = renderer.density(all_grains, weights, k)
    img = bright * (t_min + (1. - t_min) * np.exp(-d)) \
      + rng.normal(0., background.noise_sigma, d.shape)
    img[dead[:, 1], dead[:, 0]] = 0.
    img = np.clip(np.rint(img), 0, GRAY_MAX).astype(np.uint8)
    write_gray_png(out_dir / f"{k:04d}.png", img)
    np.minimum(stack_min, img, out=stack_min)
    if not args.no_labels:
      cls, ins = renderer.labels(labelled, classes, instances, k)
      cls_stack.append(cls)
      ins_stack.append(ins)
  if cls_stack:
    np.savez_compressed(out_dir / "labels.npz", labels=np.array(cls_stack),
                        instance=np.array(ins_stack))
  view_json(view, args.stage_dz, out_dir)
  truth = truth_record(ev, view, counts)
  truth.update(
    n_track_grains=int(len(labelled)), n_fog_grains=int(len(fog)),
    physical_dz_um=optics.dz_um, stage_dz_mm=args.stage_dz,
    shrink=optics.shrink, pixel_um=optics.pixel_um,
    size_px=optics.size_px, n_slices=optics.n_slices,
    view_center_um=[round(view.cx, 3), round(view.cy, 3)],
    dead_pixels=dead.tolist())
  (out_dir / "truth.json").write_text(json.dumps(truth, indent=1))
  vtx = truth["vertices"][0]
  preview(stack_min, np.array([vtx["col"], vtx["row"]]), out_dir)


def main():
  parser = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument("root_file")
  parser.add_argument("--events", default="0", help="e.g. 0-4 or 0,3")
  parser.add_argument("--name", help="output name (default: file stem)")
  parser.add_argument("--size", type=int, default=render.Optics.size_px)
  parser.add_argument("--slices", type=int,
                      default=render.Optics.n_slices)
  parser.add_argument("--dz", type=float, default=render.Optics.dz_um,
                      help="physical focal step per slice in the "
                      "developed emulsion (um)")
  parser.add_argument("--stage-dz", type=float,
                      default=DEFAULT_STAGE_DZ_MM,
                      help="stage z step per slice in image.json (mm)")
  parser.add_argument("--shrink", type=float, default=render.Optics.shrink)
  parser.add_argument("--fog", type=float,
                      default=render.Emulsion.fog_per_um3,
                      help="fog grains per developed um3")
  parser.add_argument("--beam-tracks", type=int,
                      default=render.Background.n_beam_tracks)
  parser.add_argument("--random-tracks", type=int,
                      default=render.Background.n_random_tracks)
  parser.add_argument("--blobs", type=int,
                      default=render.Background.n_blobs)
  parser.add_argument("--centered", action="store_true",
                      help="vertex at the view center, stack centered")
  parser.add_argument("--seed", type=int, default=0)
  parser.add_argument("--no-labels", action="store_true")
  parser.add_argument("--tracks-only", action="store_true",
                      help="only write g4_tracks.png and g4_view.png "
                      "(same view placement as the full rendering)")
  parser.add_argument("-j", "--jobs", type=int, default=1,
                      help="events rendered in parallel")
  args = parser.parse_args()

  optics = replace(render.Optics(), size_px=args.size,
                   n_slices=args.slices, dz_um=args.dz,
                   shrink=args.shrink)
  emulsion = replace(render.Emulsion(), fog_per_um3=args.fog)
  background = replace(render.Background(),
                       n_beam_tracks=args.beam_tracks,
                       n_random_tracks=args.random_tracks,
                       n_blobs=args.blobs)

  name = args.name or args.root_file.rsplit("/", 1)[-1].rsplit(".", 1)[0]
  indices = parse_events(args.events, len(EventFile(args.root_file)))
  tasks = [(i, args, optics, emulsion, background, name) for i in indices]
  if args.jobs > 1:
    # spawn: ROOT state is not fork-safe
    ctx = multiprocessing.get_context("spawn")
    with ctx.Pool(args.jobs) as pool:
      for _ in pool.imap_unordered(_work, tasks):
        pass
  else:
    for task in tasks:
      _work(task)


_EVENT_FILE = {}


def _work(task):
  """Render one event; the ROOT file is opened once per process."""
  i, args, optics, emulsion, background, name = task
  if args.root_file not in _EVENT_FILE:
    _EVENT_FILE[args.root_file] = EventFile(args.root_file)
  ev_file = _EVENT_FILE[args.root_file]
  start = time.time()
  out_dir = config.IMAGE_DIR / name / f"ev{i:04d}"
  render_event(ev_file.event(i), ev_file.layers, args, optics, emulsion,
               background, out_dir)
  print(f"event {i}: {out_dir} ({time.time() - start:.1f} s)",
        file=sys.stderr, flush=True)


if __name__ == "__main__":
  main()
