"""Run emulsim for one Xi- capture channel.

  ./run.py sim --preset nagara --events 1000
  ./run.py sim --target 14N --channel L10Be,L5He --policy uniform \
    --name ibuki --events 100
  ./run.py sim --preset nagara --fix LL6He=L5He,p,pi- \
    --set "/emul/thickEmulsion 530 um"
  ./run.py sim --preset nagara_replay --events 100   # measured kinematics

Writes data/sim/<name>.{mac,log,root}.
"""

import argparse
import subprocess
import sys

from . import config

PRESETS = {
  # Xi- + 12C -> LL6He + 4He + t, LL6He -> L5He + p + pi-
  "nagara": {"target": "12C", "channel": "LL6He,4He,t",
             "policy": "simple", "fix": ["LL6He=L5He,p,pi-"]},
  # Xi- + 14N -> L10Be + L5He (twin)
  "ibuki": {"target": "14N", "channel": "L10Be,L5He",
            "policy": "simple", "fix": []},
  # Observed events with their measured kinematics (./run.py replay)
  "nagara_replay": {"target": "12C", "channel": "LL6He,4He,t",
                    "policy": "simple", "fix": [],
                    "replay": "nagara.csv"},
  "ibuki_replay": {"target": "14N", "channel": "L10Be,L5He",
                   "policy": "simple", "fix": [], "replay": "ibuki.csv"},
}
DEFAULT_EVENTS = 100
DEFAULT_SEED = 1


def make_macro(args):
  out = config.OUTPUT_DIR
  lines = ["/control/verbose 1", "/run/verbose 1"]
  # /emul/ settings must precede /run/initialize
  lines += [s for s in args.set if s.startswith("/emul/")]
  lines += [
    "/run/initialize",
    f"/xi/speciesFile {out / 'species.csv'}",
    f"/xi/targetsFile {out / 'targets.csv'}",
    f"/xi/decaysFile {out / 'decays.csv'}",
    f"/xi/target {args.target}",
    f"/xi/channel {args.channel}",
    *([f"/xi/replayFile {args.replay}"] if args.replay else []),
    f"/xi/seed {args.seed}",
    f"/decay/policy {args.policy}",
  ]
  for fix in args.fix:
    parent, products = fix.split("=")
    lines.append(f"/decay/channel {parent} {products}")
  lines += [s for s in args.set if not s.startswith("/emul/")]
  lines += [f"/out/file {config.SIM_DIR / (args.name + '.root')}",
            f"/run/beamOn {args.events}"]
  return "\n".join(lines) + "\n"


def main():
  parser = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument("--preset", choices=sorted(PRESETS))
  parser.add_argument("--target")
  parser.add_argument("--channel", help="e.g. LL6He,4He,t")
  parser.add_argument("--policy",
                      choices=["simple", "uniform", "uniformPiMinus"])
  parser.add_argument("--fix", action="append",
                      help="fix a decay, e.g. LL6He=L5He,p,pi-")
  parser.add_argument("--set", action="append", default=[],
                      help='extra macro command, e.g. "/emul/density 3.6"')
  parser.add_argument("--replay",
                      help="fixed kinematics CSV (output/replay/*.csv)")
  parser.add_argument("--events", type=int, default=DEFAULT_EVENTS)
  parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
  parser.add_argument("--name", help="output name (default: preset)")
  args = parser.parse_args()

  preset = PRESETS.get(args.preset, {})
  for key in ("target", "channel", "policy", "fix"):
    if getattr(args, key) is None:
      setattr(args, key, preset.get(key))
  if args.replay is None and "replay" in preset:
    args.replay = str(config.OUTPUT_DIR / "replay" / preset["replay"])
  args.policy = args.policy or "simple"
  args.fix = args.fix or []
  args.name = args.name or args.preset
  if not (args.target and args.channel and args.name):
    parser.error("give --preset or --target, --channel and --name")
  if not config.EMULSIM.exists():
    parser.error("emulsim not built; run ./run.py build")

  config.SIM_DIR.mkdir(parents=True, exist_ok=True)
  mac = config.SIM_DIR / f"{args.name}.mac"
  log = config.SIM_DIR / f"{args.name}.log"
  mac.write_text(make_macro(args))
  command = config.geant4_env_command(f"{config.EMULSIM} {mac}")
  with log.open("w") as f:
    result = subprocess.run(["bash", "-c", command], stdout=f,
                            stderr=subprocess.STDOUT)
  status = "done" if result.returncode == 0 else "FAILED"
  print(f"{status}: {config.SIM_DIR / (args.name + '.root')} (log {log})",
        file=sys.stderr)
  return result.returncode


if __name__ == "__main__":
  raise SystemExit(main())
