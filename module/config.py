"""Paths shared by the commands of run.py."""

import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUTPUT_DIR = ROOT / "output"          # mass / channel / decay tables
DATA_DIR = ROOT / "data"              # generated ROOT files and images
SIM_DIR = DATA_DIR / "sim"
IMAGE_DIR = DATA_DIR / "images"
G4_SOURCE_DIR = ROOT / "g4"
G4_BUILD_DIR = G4_SOURCE_DIR / "build"
EMULSIM = G4_BUILD_DIR / "emulsim"

# External installations, overridable by environment variables
GEANT4_DIR = Path(os.environ.get(
  "GEANT4_DIR", Path.home() / "software/geant4/11.4.3"))
KINEMA_DIR = Path(os.environ.get("KINEMA_DIR", Path.home() / "work/kinema"))


def geant4_env_command(command):
  """Shell command running `command` with the Geant4 environment."""
  return f"source {GEANT4_DIR / 'bin/geant4.sh'} && {command}"
