"""Frozen logistic review ranker: fit once, apply to any feature table.

  ./run.py candscore fit model.json [--drop COL,...] DS A B F [...]
  ./run.py candscore apply model.json features.parquet out.parquet

fit: the logistic ranker of module/candrank.py (features of design(),
all positives and NEG_FRACTION of the negatives, balanced classes)
fitted on every candidate of the given datasets, written as JSON: the
design columns with their transforms, the standardisation (mean,
scale), coefficients and intercept. The score is
  intercept + sum(coef * (x - mean) / scale)
over the columns of design(), larger = review first.
apply: adds the column review_score to a feature table with the
fullscan schema (needs source, darkness, dist_other_um and the A / B'
features; rows of other sources are refused).
"""

import argparse
import datetime
import json
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

from . import candrank as cr

SCORE_COLUMN = "review_score"
SOURCES = ("A", "Bprime")


def git_commit():
  out = subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                       capture_output=True, text=True)
  return out.stdout.strip()


def fit(path, sets, drop):
  df = pd.concat([cr.load_set(*s, True) for s in sets], ignore_index=True)
  x, y = cr.design(df), (df["d_cap"] <= cr.LABEL_UM).to_numpy()
  rng = np.random.default_rng(cr.SEED)
  keep = y | (rng.random(len(y)) < cr.NEG_FRACTION)
  sc = StandardScaler().fit(x[keep])
  m = LogisticRegression(max_iter=cr.LR_MAX_ITER, class_weight="balanced")
  m.fit(sc.transform(x[keep]), y[keep])
  model = dict(
    columns=list(x.columns), mean=sc.mean_.tolist(),
    scale=sc.scale_.tolist(), coef=m.coef_[0].tolist(),
    intercept=float(m.intercept_[0]),
    transforms=dict(log1p=cr.LOG_FEATURES,
                    dist_cap_um=cr.DIST_CAP_UM,
                    no_other="1 if dist_other_um is not finite",
                    is_a="1 if source == 'A'",
                    missing="0 after the transforms"),
    dropped=sorted(drop), label_um=cr.LABEL_UM, neg_fraction=cr.NEG_FRACTION,
    trained_on=[Path(s[0]).name for s in sets],
    n_rows=int(len(df)), n_positive=int(y.sum()),
    commit=git_commit(),
    created=datetime.datetime.now().isoformat(timespec="seconds"))
  Path(path).write_text(json.dumps(model, indent=1) + "\n")
  print(f"{path}: {len(df)} rows, {y.sum()} positives")


def score(model, feats):
  bad = ~feats["source"].isin(SOURCES)
  if bad.any():
    raise SystemExit(f"{bad.sum()} rows with an unknown source")
  x = cr.design(feats)
  if list(x.columns) != model["columns"]:
    raise SystemExit("design columns differ from the model")
  z = (x.to_numpy() - np.array(model["mean"])) / np.array(model["scale"])
  return model["intercept"] + z @ np.array(model["coef"])


def main():
  parser = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
  sub = parser.add_subparsers(dest="cmd", required=True)
  f = sub.add_parser("fit")
  f.add_argument("model")
  f.add_argument("--drop", default="")
  f.add_argument("sets", nargs="+", help="DS A B F, repeated")
  a = sub.add_parser("apply")
  a.add_argument("model")
  a.add_argument("features")
  a.add_argument("output")
  args = parser.parse_args()

  if args.cmd == "fit":
    if len(args.sets) % 4:
      raise SystemExit("sets come in groups of DS A B F")
    drop = cr.drop_features(args.drop)
    fit(args.model, [args.sets[i:i + 4]
                     for i in range(0, len(args.sets), 4)], drop)
    return
  model = json.loads(Path(args.model).read_text())
  cr.drop_features(",".join(model.get("dropped", [])))
  feats = pd.read_parquet(args.features)
  feats[SCORE_COLUMN] = score(model, feats)
  feats.to_parquet(args.output)
  print(f"{args.output}: {len(feats)} rows scored")


if __name__ == "__main__":
  main()
