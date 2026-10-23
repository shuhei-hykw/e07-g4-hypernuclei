"""Configure and build the Geant4 application g4/build/emulsim."""

import argparse
import os
import subprocess
import sys

from . import config


def main():
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument("-j", "--jobs", type=int, default=os.cpu_count())
  parser.add_argument("--clean", action="store_true",
                      help="remove the CMake cache first")
  args = parser.parse_args()

  build = config.G4_BUILD_DIR
  build.mkdir(parents=True, exist_ok=True)
  cache = build / "CMakeCache.txt"
  if args.clean and cache.exists():
    cache.unlink()
  command = (f"cmake -S {config.G4_SOURCE_DIR} -B {build} "
             f"-DCMAKE_BUILD_TYPE=Release && "
             f"cmake --build {build} -j {args.jobs}")
  result = subprocess.run(["bash", "-c", config.geant4_env_command(command)])
  if result.returncode == 0:
    print(f"built {config.EMULSIM}", file=sys.stderr)
  return result.returncode


if __name__ == "__main__":
  raise SystemExit(main())
