"""Load simulated events and the module layers from an emulsim file."""

from dataclasses import dataclass, field

import numpy as np
import ROOT

DECAY_PROCESS = "HypernucleusDecay"


@dataclass
class Track:
  id: int
  parent: int
  name: str
  process: str
  charge: int
  a: int
  ke: float          # MeV
  range_um: float
  escaped: bool
  points: np.ndarray  # (n, 3) um, original (undeveloped) coordinates
  edep: np.ndarray    # (n,) MeV deposited in the step ending at points[i]

  @property
  def origin(self):
    """capture, decay or other (pi- capture, recoils, ...)."""
    if self.parent == 0:
      return "capture"
    if self.process == DECAY_PROCESS:
      return "decay"
    return "other"


@dataclass
class Event:
  index: int
  channel: str
  vertex: np.ndarray  # (3,) um
  contained: bool
  tracks: list = field(default_factory=list)


@dataclass
class Layer:
  emulsion: bool
  zmin: float  # um
  zmax: float


class EventFile:
  def __init__(self, path):
    self._file = ROOT.TFile.Open(str(path))
    if not self._file or self._file.IsZombie():
      raise OSError(f"cannot open {path}")
    self._tree = self._file.Get("tree")
    self.layers = [Layer(bool(l.emulsion), l.zmin, l.zmax)
                   for l in self._file.Get("layers")]

  def __len__(self):
    return int(self._tree.GetEntries())

  def event(self, index):
    t = self._tree
    t.GetEntry(index)
    ev = Event(index, str(t.channel), np.array([t.vx, t.vy, t.vz]),
               bool(t.contained))
    for j in range(len(t.name)):
      points = np.column_stack([np.asarray(t.x[j], dtype=np.float64),
                                np.asarray(t.y[j], dtype=np.float64),
                                np.asarray(t.zpos[j], dtype=np.float64)])
      ev.tracks.append(Track(
        int(t.id[j]), int(t.parent[j]), str(t.name[j]), str(t.process[j]),
        int(t.z[j]), int(t.a[j]), float(t.ke[j]), float(t.range[j]),
        bool(t.escaped[j]), points,
        np.asarray(t.edep[j], dtype=np.float64)))
    return ev
