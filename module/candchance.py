"""Chance-subtracted review ranking: learned ranker against darkness.

  ./run.py candchance data/images/<dataset> --a A.csv --b B.csv
    --features F.parquet --a-bg A_bg.csv --b-bg B_bg.csv
    --features-bg F_bg.parquet [--drop COL,...]
    [-o output/candchance_<name>.md]

The pool of module/candrank.py is built twice for the same events: on
the embedded views and on the real views alone (--*-bg, the background
views before embedding; distances to the same truth vertices, so a
found vertex there is a chance match). Per event and budget K the gain
G = found(embedded) - found(real) is reported per ranker, with the
logistic ranker scored out of fold (folds grouped by view; the real
view of an event is scored by the model that did not see that view).
The paired difference G(logistic) - G(darkness) has a 95 % bootstrap
interval over events. Per-view budgets map to plate totals by the
number of full-scan views (no cross-view merge: every synthetic event
lies in one view, so the total is an upper bound before dedup).
--drop removes learned features (robustness check against features
the rendering may shift, ./run.py candleak).
"""

import argparse
from pathlib import Path

import numpy as np
from sklearn.model_selection import GroupKFold

from . import candrank as cr

BUDGETS = [1, 2, 4, 8, 16]
N_VIEWS = 2025          # full-scan views per plate (45 x 45)
RANKERS = ("darkness", "logistic")


def logistic_both(emb, real):
  """Out-of-fold logistic scores of the embedded and real pools."""
  x, y = cr.design(emb).to_numpy(), (emb["d_cap"] <= cr.LABEL_UM).to_numpy()
  xr = cr.design(real).to_numpy()
  rng = np.random.default_rng(cr.SEED)
  se, sr = np.zeros(len(emb)), np.zeros(len(real))
  for tr, te in GroupKFold(cr.N_FOLDS).split(x, y, emb["view"]):
    rte = real["view"].isin(set(emb["view"].iloc[te])).to_numpy()
    s = cr.fit_predict(x[tr], y[tr], np.r_[x[te], xr[rte]], "logistic",
                       rng)
    se[te], sr[rte] = s[:len(te)], s[len(te):]
  return se, sr


def found(df, score, tol):
  """events x BUDGETS, events in sorted order."""
  f = cr.found_at(df, score, "d_cap", tol)
  return f[:, [cr.BUDGETS.index(k) for k in BUDGETS]]


def main():
  parser = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument("dataset")
  for k in ("a", "b", "features", "a-bg", "b-bg", "features-bg"):
    parser.add_argument(f"--{k}", required=True)
  parser.add_argument("--drop", default="")
  parser.add_argument("-o", "--output")
  args = parser.parse_args()

  drop = cr.drop_features(args.drop)

  emb = cr.load_set(args.dataset, args.a, args.b, args.features, False)
  real = cr.load_set(args.dataset, args.a_bg, args.b_bg, args.features_bg,
                     False)
  if sorted(emb["event"].unique()) != sorted(real["event"].unique()):
    raise SystemExit("embedded and real pools have different events")
  le, lr = logistic_both(emb, real)
  scores = {"darkness": (emb["darkness"].to_numpy(),
                         real["darkness"].to_numpy()),
            "logistic": (le, lr)}
  rng = np.random.default_rng(cr.SEED)
  name = Path(args.dataset).name
  text = [f"# Chance-subtracted ranking on {name}", "",
          f"{emb['event'].nunique()} events"
          + (f", learned features without {args.drop}" if drop else "")
          + ". Capture found among the top K "
          "of the view: embedded (real views alone), G = embedded - real "
          "with 95 % bootstrap interval; plate total = K x "
          f"{N_VIEWS} views (before cross-view dedup).", ""]
  for tol in cr.TOLS:
    text += [f"## capture within {tol:g} um", "",
             "| ranker | " + " | ".join(
               f"K={k} ({k * N_VIEWS})" for k in BUDGETS) + " |",
             "|---|" + "---|" * len(BUDGETS)]
    g = {}
    for r in RANKERS:
      fe, fr = found(emb, scores[r][0], tol), found(real, scores[r][1], tol)
      g[r] = fe - fr
      lo, hi = cr.boot_ci(g[r], rng)
      text.append(f"| {r} | " + " | ".join(
        f"{fe.mean(0)[i]:.0%} ({fr.mean(0)[i]:.0%}) G {g[r].mean(0)[i]:+.0%}"
        f" [{lo[i]:+.0%}, {hi[i]:+.0%}]" for i in range(len(BUDGETS)))
        + " |")
    d = g["logistic"] - g["darkness"]
    lo, hi = cr.boot_ci(d, rng)
    text += ["| logistic - darkness | " + " | ".join(
      f"{d.mean(0)[i]:+.0%} [{lo[i]:+.0%}, {hi[i]:+.0%}]"
      for i in range(len(BUDGETS))) + " |", ""]
  out = "\n".join(text)
  print(out)
  if args.output:
    Path(args.output).write_text(out)


if __name__ == "__main__":
  main()
