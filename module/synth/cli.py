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
                 -100001-i junk cluster i, 0 none); grains within
                 label_radius_px (2 px) of a pixel paint it as a disk,
                 lower classes win where disks overlap, and within a
                 class the grain listed later wins; fog is not labelled
  grains.npz     per-grain truth: col, row, slice (view coordinates),
                 x_um, y_um, z_um (developed um from the view corner and
                 slice 0), weight (darkness), cls / instance as in the
                 labels (fog 0 / 0), fog flag and focus_slice
  background_tracks.json  straight background tracks (beam, isotropic,
                 flat): instance id, kind, grains per um and straight
                 segments [[col, row, slice] start, end], one per
                 emulsion layer crossed, inside the emulsion near the view
                 (same pad and defocus margin as the rendered grains)
  preview.png    min projection: crop around the vertex | whole view
  g4_tracks.png  Geant4 event display before imaging (./run.py draw)
  g4_view.png    Geant4 tracks in the image coordinates of the view

The capture point is placed at a random position in the view and the
stack at a random depth inside the emulsion layer (--centered disables).

Model parameters tuned by ./run.py calibrate (output/image_model.json)
are used when present (--model none for the built-in defaults).

--background random|<view.json> embeds the signal into a real full-scan
view instead: only the Geant4 tracks are rendered and their transmission
multiplies the real slices (fog, beam tracks, junk and noise are real);
the stack size and stage z follow the real view and the labels cover the
signal only.
"""

import argparse
import json
import multiprocessing
import sys
import time
from dataclasses import replace
from pathlib import Path

import numpy as np

from .. import config
from .. import draw_event
from .. import imagestats
from . import calibrate
from . import flatfield
from . import render
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
BG_SAMPLE_STEP_UM = 1.  # sampling to find the segment ends


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


def collect_grains(ev, view, devz, emulsion, background, rng,
                   signal_only=False, tracks=None):
  """Labelled grains (view coords) with class and instance, and fog.
  If tracks is a list, (instance id, kind, start, end, grains per um) of
  the straight background tracks (original um) is appended to it."""
  o = view.optics
  lo = -o.pad_px
  hi = o.size_px + o.pad_px
  margin = render.max_defocus_um(o) / o.dz_um + 1.

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
  if signal_only:
    return (np.vstack(grains), np.concatenate(classes),
            np.concatenate(instances), np.zeros((0, 3)), counts)
  # Instance ids of background tracks: BACKGROUND_ID0 - running index
  # over beam, isotropic and flat tracks (unique across the three kinds)
  offset = 0
  for kind, fn in (("beam", render.beam_grains),
                   ("isotropic", render.random_track_grains),
                   ("flat", render.flat_track_grains)):
    geometry = []
    g, owner = fn(view, background, emulsion, ev.vertex[2], rng,
                  geometry=geometry)
    add(g, render.LABEL_BEAM, BACKGROUND_ID0 - offset - owner)
    if tracks is not None:
      tracks += [(BACKGROUND_ID0 - offset - i, kind) + geo
                 for i, geo in enumerate(geometry)]
    offset += len(geometry)
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


def background_track_records(tracks, view, devz):
  """Straight background tracks as straight segments [[col, row, slice]
  start, end] (one per emulsion layer crossed) inside the emulsion and
  near the view; tracks not reaching it are dropped."""
  o = view.optics
  lo, hi = -o.pad_px, o.size_px + o.pad_px
  # Same depth margin as the grains kept by collect_grains
  margin = render.max_defocus_um(o) / o.dz_um + 1.
  out = []
  for tid, kind, start, end, density in tracks:
    length = np.linalg.norm(end - start)
    n = max(2, int(np.ceil(length / BG_SAMPLE_STEP_UM)) + 1)
    pts = start + np.linspace(0., 1., n)[:, None] * (end - start)
    v = view.to_view(pts)
    keep = devz.in_emulsion(pts[:, 2]) & (v[:, 0] >= lo) & (v[:, 0] < hi) \
      & (v[:, 1] >= lo) & (v[:, 1] < hi) & (v[:, 2] > -margin) \
      & (v[:, 2] < o.n_slices - 1 + margin)
    if not keep.any():
      continue
    # Contiguous runs of kept samples: straight segments (one per layer)
    edge = np.diff(np.r_[0, keep.astype(int), 0])
    first, last = np.flatnonzero(edge == 1), np.flatnonzero(edge == -1) - 1
    out.append({"id": int(tid), "kind": kind,
                "grains_per_um": round(float(density), 5),
                "segments": [np.round(v[[i, j]], 2).tolist()
                             for i, j in zip(first, last)]})
  return out


def grain_table(grains, weights, classes, instances, n_fog, optics):
  """Per-grain truth: view coordinates, developed um relative to the
  view corner and slice 0, weight, class / instance (fog: 0 / 0)."""
  n = len(grains)
  return dict(
    col=grains[:, 0].astype(np.float32), row=grains[:, 1].astype(np.float32),
    slice=grains[:, 2].astype(np.float32),
    x_um=(grains[:, 0] * optics.pixel_um).astype(np.float32),
    y_um=(grains[:, 1] * optics.pixel_um).astype(np.float32),
    z_um=(grains[:, 2] * optics.dz_um).astype(np.float32),
    weight=weights.astype(np.float32),
    cls=np.concatenate([classes, np.zeros(n_fog)]).astype(np.uint8),
    instance=np.concatenate([instances, np.zeros(n_fog)]).astype(np.int32),
    fog=np.arange(n) >= n - n_fog,
    focus_slice=np.rint(grains[:, 2]).astype(np.int16))


def view_json(view, stage_z_mm, out_dir):
  o = view.optics
  pix_mm = o.pixel_um / UM_PER_MM
  x0 = view.cx / UM_PER_MM - o.size_px / 2. * pix_mm
  y0 = view.cy / UM_PER_MM - o.size_px / 2. * pix_mm
  images = [{"Path": f"{k:04d}.png", "x": round(view.cx / UM_PER_MM, 6),
             "y": round(view.cy / UM_PER_MM, 6), "z": round(z, 8)}
            for k, z in enumerate(stage_z_mm)]
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


def pick_background(args, rng):
  """Real view JSON used as background, or None."""
  if not args.background:
    return None
  if args.background != "random":
    return Path(args.background)
  views = [v for v in sorted(Path(calibrate.REAL_DIR).glob("V*.json"))
           if not v.name.startswith(calibrate.EXCLUDED_VIEWS)]
  return views[int(rng.integers(len(views)))]


def render_event(ev, layers, args, optics, emulsion, background, out_dir):
  rng = np.random.default_rng(args.seed + ev.index)
  bg_json = pick_background(args, rng)
  if bg_json:
    real = imagestats.read_stack(bg_json)
    stage_z = [im["z"] for im in json.loads(bg_json.read_text())["Images"]]
    optics = replace(optics, n_slices=len(real), size_px=real[0].shape[0])
  else:
    stage_z = [STAGE_Z0_MM + k * args.stage_dz
               for k in range(optics.n_slices)]
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
  bg_tracks = []
  labelled, classes, instances, fog, counts = collect_grains(
    ev, view, devz, emulsion, background, rng, signal_only=bool(bg_json),
    tracks=bg_tracks)
  all_grains = np.vstack([labelled, fog])
  weights = render.grain_weights(len(all_grains), emulsion, rng)
  if len(fog):
    weights[-len(fog):] *= emulsion.fog_weight
  if not args.no_labels:
    np.savez_compressed(out_dir / "grains.npz", **grain_table(
      all_grains, weights, classes, instances, len(fog), optics))
    (out_dir / "background_tracks.json").write_text(json.dumps(
      background_track_records(bg_tracks, view, devz)))
  renderer = render.Renderer(optics)
  mottle = render.MottleVolume(optics, background, rng, view.slice_zdev(0),
                               view.slice_zdev(optics.n_slices - 1))
  flat = flatfield.load(optics.size_px)
  fixed = render.fixed_pattern(optics, background, rng)
  dead = rng.integers(0, optics.size_px, (background.n_dead_pixels, 2))

  stack_min = np.full((optics.size_px, optics.size_px), GRAY_MAX,
                      dtype=np.uint8)
  cls_stack, ins_stack = [], []
  t_min = background.min_transmission
  for k in range(optics.n_slices):
    d = renderer.density(all_grains, weights, k)
    if bg_json:
      trans = np.exp(-d)
      if optics.image_blur_px > 0.:
        trans = render.gaussian_blur(trans, optics.image_blur_px)
      img = real[k] * (t_min + (1. - t_min) * trans)
      img = np.clip(np.rint(img), 0, GRAY_MAX).astype(np.uint8)
    else:
      img = render.compose(mottle.at(view.slice_zdev(k)) * flat, d,
                           background, fixed, rng, optics.image_blur_px)
      img[dead[:, 1], dead[:, 0]] = 0
    write_gray_png(out_dir / f"{k:04d}.png", img)
    np.minimum(stack_min, img, out=stack_min)
    if not args.no_labels:
      cls, ins = renderer.labels(labelled, classes, instances, k)
      cls_stack.append(cls)
      ins_stack.append(ins)
  if cls_stack:
    np.savez_compressed(out_dir / "labels.npz", labels=np.array(cls_stack),
                        instance=np.array(ins_stack))
  view_json(view, stage_z, out_dir)
  truth = truth_record(ev, view, counts)
  truth.update(
    n_track_grains=int(len(labelled)), n_fog_grains=int(len(fog)),
    physical_dz_um=optics.dz_um,
    stage_z0_mm=stage_z[0],
    stage_dz_mm=(stage_z[-1] - stage_z[0]) / max(len(stage_z) - 1, 1),
    background_view=str(bg_json) if bg_json else None,
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
  parser.add_argument("--beam-density", type=float,
                      help="beam tracks per mm2 (default: model)")
  parser.add_argument("--random-density", type=float,
                      help="isotropic tracks per mm2 (default: model)")
  parser.add_argument("--blob-density", type=float,
                      help="junk clusters per mm2 (default: model)")
  parser.add_argument("--centered", action="store_true",
                      help="vertex at the view center, stack centered")
  parser.add_argument("--seed", type=int, default=0)
  parser.add_argument("--no-labels", action="store_true")
  parser.add_argument("--tracks-only", action="store_true",
                      help="only write g4_tracks.png and g4_view.png "
                      "(same view placement as the full rendering)")
  parser.add_argument("-j", "--jobs", type=int, default=1,
                      help="events rendered in parallel")
  parser.add_argument("--model", default=str(config.OUTPUT_DIR
                                             / calibrate.MODEL_FILE),
                      help="tuned model parameters ('none' for defaults)")
  parser.add_argument("--background",
                      help="'random' or a real view JSON to embed into")
  args = parser.parse_args()

  optics = replace(render.Optics(), size_px=args.size,
                   n_slices=args.slices, dz_um=args.dz,
                   shrink=args.shrink)
  emulsion = replace(render.Emulsion(), fog_per_um3=args.fog)
  background = render.Background()
  if args.model != "none" and Path(args.model).exists():
    optics, emulsion, background = calibrate.apply_values(
      calibrate.load_model(args.model), optics, emulsion, background)
    print(f"model parameters from {args.model}", file=sys.stderr)
  for field, value in (("beam_per_mm2", args.beam_density),
                       ("random_per_mm2", args.random_density),
                       ("blobs_per_mm2", args.blob_density)):
    if value is not None:
      background = replace(background, **{field: value})

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
