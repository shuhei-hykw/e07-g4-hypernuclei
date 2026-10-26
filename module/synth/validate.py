"""Realism of synthetic stacks on held-out real views.

  ./run.py validate data/images/<synthetic dataset> [--views 20]
    [-o output/realism.md]

Random tohoku-v1 views not used by the calibration (and not the
hand-labelled ones) give the mean and spread of every image statistic;
each synthetic event of the dataset is measured the same way, and the
deviation is given as z = (synthetic mean - real mean) / real spread.
"""

import argparse
import json
import random
from pathlib import Path

import numpy as np

from .. import config
from .. import imagestats as st
from . import calibrate as c

KEYS = ["median", "std", "p1", "p5", "frac_bright", "fg_fraction",
        "n_cc_per_mpx", "cc_area_p90", "n_elong_per_mpx",
        "dark_blobs_per_mpx", "dark_contrast_median", "d_skew", "d_kurt",
        "d_iqr", "pow_low", "pow_mid", "pow_high", "mp_detail_std",
        "mp_dark_frac", "corr1", "corr2", "corr5", "corr13"]
HOLDOUT_SEED = 7
Z_WARN = 3.


def stats_of(path):
  n = st.n_slices(path) // 2
  idx = [n + o for o in c.OFFSETS]
  return c.flatten(st.stack_stats(st.read_stack(path, idx), c.CROP_PX,
                                  c.OFFSETS))


def holdout_views(n):
  """Real views not used by the calibration nor labelled by hand, in a
  fixed random order."""
  used = set(json.loads((config.OUTPUT_DIR
                         / "real_image_stats.json").read_text())["views"])
  views = [v for v in sorted(Path(c.REAL_DIR).glob("V*.json"))
           if str(v) not in used
           and not v.name.startswith(c.EXCLUDED_VIEWS)]
  random.Random(HOLDOUT_SEED).shuffle(views)
  return views[:n]


def main():
  parser = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument("dataset")
  parser.add_argument("--views", type=int, default=20)
  parser.add_argument("-o", "--output",
                      default=str(config.OUTPUT_DIR / "realism.md"))
  args = parser.parse_args()

  real = [stats_of(v) for v in holdout_views(args.views)]
  synth = [stats_of(d / "image.json")
           for d in sorted(Path(args.dataset).glob("ev*"))]
  lines = [f"# Realism of {Path(args.dataset).name} on held-out real views",
           "", f"{len(real)} held-out tohoku-v1 views vs {len(synth)} "
           "synthetic views (central 1024 px, slices at offsets "
           f"{c.OFFSETS}); z = (synthetic - real mean) / real spread.", "",
           "| statistic | real | synthetic | z |", "|---|---|---|---|"]
  n_warn = 0
  for k in KEYS:
    r = np.array([x[k] for x in real])
    s = np.array([x[k] for x in synth])
    z = (s.mean() - r.mean()) / (r.std() or 1e-9)
    n_warn += abs(z) > Z_WARN
    lines.append(f"| {k} | {r.mean():.4g} +- {r.std():.2g} | "
                 f"{s.mean():.4g} +- {s.std():.2g} | {z:+.1f} |")
  lines += ["", f"{n_warn} of {len(KEYS)} statistics with |z| > {Z_WARN:g}."]
  Path(args.output).write_text("\n".join(lines) + "\n")
  print("\n".join(lines))


if __name__ == "__main__":
  main()
