"""Leak check of the review-ranking features: embedded vs real views.

  ./run.py candleak data/images/<dataset> --a A.csv --b B.csv
    --features F.parquet --a-bg A_bg.csv --b-bg B_bg.csv
    --features-bg F_bg.parquet [-o output/candleak_<name>.md]

The same real views are scored before (--*-bg) and after embedding the
signal. An adversarial classifier (gradient-boosted trees, folds grouped
by view) tries to tell embedded from real candidates that are not
within VERTEX_UM of a truth vertex, in two regions around the capture
point: near (VERTEX_UM .. NEAR_UM, where the rendered signal tracks
change the image) and far (beyond FAR_UM, where both images are the
same real view, a sanity check expected at AUC 0.5). The out-of-fold
AUC of all learned features together and the AUC of each feature alone
are reported. Identical backgrounds make this a test of how the
rendered signal shifts the features of the surrounding candidates, not
of whether the embedded vertices look like real ones.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold

from . import candrank as cr
from . import efficiency as eff

VERTEX_UM = 6.
NEAR_UM = 50.
FAR_UM = 200.
MAX_FAR_ROWS = 200000
N_FOLDS = 5
SEED = 1
GBDT_ITER = 200
GBDT_LR = 0.05


def region(df, name):
  away = (df["d_cap"] > VERTEX_UM) & (df["d_ll"] > VERTEX_UM)
  if name == "near":
    return df[away & (df["d_cap"] <= NEAR_UM)]
  return df[df["d_cap"] > FAR_UM]


def adversary(x, y, groups):
  scores = np.zeros(len(y))
  for tr, te in GroupKFold(N_FOLDS).split(x, y, groups):
    m = HistGradientBoostingClassifier(max_iter=GBDT_ITER,
                                       learning_rate=GBDT_LR,
                                       random_state=SEED)
    m.fit(x.iloc[tr], y[tr])
    scores[te] = m.predict_proba(x.iloc[te])[:, 1]
  return roc_auc_score(y, scores)


def main():
  parser = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument("dataset")
  for k in ("a", "b", "features", "a-bg", "b-bg", "features-bg"):
    parser.add_argument(f"--{k}", required=True)
  parser.add_argument("-o", "--output")
  args = parser.parse_args()

  events = eff.load_truth(args.dataset)
  emb = cr.pool(events, eff.load_candidates(args.a),
                eff.load_candidates(args.b), pd.read_parquet(args.features))
  real = cr.pool(events, eff.load_candidates(args.a_bg),
                 eff.load_candidates(args.b_bg),
                 pd.read_parquet(args.features_bg))
  name = Path(args.dataset).name
  text = [f"# Leak check of the ranking features on {name}", "",
          "Embedded vs real candidates of the same views, not within "
          f"{VERTEX_UM:g} um of a truth vertex; near = {VERTEX_UM:g}-"
          f"{NEAR_UM:g} um from the capture point, far = beyond "
          f"{FAR_UM:g} um (both images identical there, expected 0.5). "
          "AUC out of fold (folds grouped by view); per feature the AUC "
          "of the feature alone (folded to >= 0.5).", ""]
  per_feature = {}
  for reg in ("near", "far"):
    e, r = region(emb, reg), region(real, reg)
    if reg == "far":
      e = e.sample(min(len(e), MAX_FAR_ROWS), random_state=SEED)
      r = r.sample(min(len(r), MAX_FAR_ROWS), random_state=SEED)
    both = pd.concat([e, r], ignore_index=True)
    y = np.r_[np.ones(len(e)), np.zeros(len(r))]
    x = cr.design(both)
    auc = adversary(x, y, both["view"])
    text.append(f"- {reg}: {len(e)} embedded / {len(r)} real "
                f"candidates, adversary AUC {auc:.3f}")
    for src in ("A", "Bprime"):
      m = (both["source"] == src).to_numpy()
      for c in x.columns:
        v = x.loc[m, c].to_numpy()
        if np.ptp(v) == 0.:
          continue
        a = roc_auc_score(y[m], v)
        per_feature.setdefault((src, c), {})[reg] = max(a, 1. - a)
  text += ["", "| source | feature | near AUC | far AUC |",
           "|---|---|---|---|"]
  for (src, c), v in sorted(per_feature.items(),
                            key=lambda t: -t[1].get("near", 0.)):
    text.append(f"| {src} | {c} | {v.get('near', np.nan):.3f} | "
                f"{v.get('far', np.nan):.3f} |")
  out = "\n".join(text) + "\n"
  print(out)
  if args.output:
    Path(args.output).write_text(out)


if __name__ == "__main__":
  main()
