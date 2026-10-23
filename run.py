#!/usr/bin/env python3
"""E07 Geant4 hypernuclear event simulation entry point.

Usage:
  ./run.py <command> [args...]
  ./run.py <command> --help

Commands:
  channels     enumerate Xi- capture channels  -> output/channels_*.csv
  decays       enumerate hypernuclear weak decays -> output/decays.csv
  build        build the Geant4 application g4/build/emulsim
  sim          run emulsim for one channel     -> data/sim/<name>.root
  draw         event display of one event (PNG)
  summary      containment summary plot (PNG)
  check-range  compare Geant4 ranges with the kinema range-energy relation
  replay       fixed kinematics of Nagara / Ibuki -> output/replay/
  images       render microscope image stacks  -> data/images/<name>/
  compare      image statistics of real vs synthetic stacks
  calibrate    tune the image model to real views -> output/image_model.json
  flatfield    illumination profile of real views -> output/flatfield.npy
  validate     realism of synthetic stacks on held-out real views
  beamgrains   grains per slice on perpendicular beam tracks (step)
  rankstudy    ranking study for point candidates (B')
  specialsrank darkness ranking of B' points on real special events
  mipdensity   in-plane MIP grain density from flat tracks
  efficiency   vertex-finding efficiency of candidates vs truth
  report       efficiency table over datasets and methods
  zbias        depth bias of candidate z_best vs truth
  evalreport   markdown summary of the full-scan evaluation
  setup        channels + decays + build

Typical:
  ./run.py setup
  ./run.py sim --preset nagara --events 100
  ./run.py draw data/sim/nagara.root --event 0 -o ev0.png
  ./run.py images data/sim/nagara.root --events 0-4
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent

COMMANDS: dict[str, str] = {
  "channels":    "module.enumerate_channels",
  "decays":      "module.enumerate_decays",
  "build":       "module.build",
  "sim":         "module.simulate",
  "draw":        "module.draw_event",
  "summary":     "module.summary",
  "check-range": "module.check_range",
  "replay":      "module.replay",
  "images":      "module.synth.cli",
  "compare":     "module.imagestats",
  "calibrate":   "module.synth.calibrate",
  "flatfield":   "module.synth.flatfield",
  "validate":    "module.synth.validate",
  "beamgrains":  "module.beamgrains",
  "rankstudy":   "module.rankstudy",
  "mipdensity":  "module.mipdensity",
  "specialsrank": "module.specials_rank",
  "efficiency":  "module.efficiency",
  "report":      "module.report",
  "zbias":       "module.zbias",
  "evalreport":  "module.evalreport",
}
SEQUENCES: dict[str, list[str]] = {
  "setup": ["channels", "decays", "build"],
}


def _usage() -> None:
  print(__doc__)


def _run(cmd: str, rest: list[str]) -> int:
  argv = [sys.executable, "-m", COMMANDS[cmd], *rest]
  return subprocess.run(argv, cwd=str(ROOT)).returncode


def main() -> int:
  if len(sys.argv) < 2 or sys.argv[1] in ("-h", "--help"):
    _usage()
    return 0
  cmd = sys.argv[1]
  rest = sys.argv[2:]
  if cmd in SEQUENCES:
    for step in SEQUENCES[cmd]:
      code = _run(step, [])
      if code != 0:
        return code
    return 0
  if cmd not in COMMANDS:
    print(f"unknown command: {cmd}\n", file=sys.stderr)
    _usage()
    return 2
  return _run(cmd, rest)


if __name__ == "__main__":
  raise SystemExit(main())
