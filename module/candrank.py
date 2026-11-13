"""Ranking of the full-scan candidates (method A and B') for review.

  ./run.py candrank data/images/<dataset> --a A.csv --b B.csv
    --features F.parquet [-o output/candrank_<name>.md]

The review budget of the full scan is a few candidates per view, so the
question is which K candidates of a view to look at. Candidates of both
sources (every method-A candidate and every B' branch point) are put in
one pool per view and ranked by a score; the capture vertex (and the LL
decay vertex) counts as found at budget K when one of the top K lies
within the tolerance (3-D, developed space, inclusive depth convention
of module/efficiency.py). Rankers:
  random            uniform order (fixed seed)
  darkness          image darkness (same local darkness as the union)
  A tracks          method A by n_tracks_max only (B' never chosen)
  interleaved       A by n_tracks_max and B' by darkness alternately
  logistic, gbdt    learned on the features (logistic regression;
                    gradient-boosted trees), scored out of fold with
                    folds grouped by background view
Features come from the fullscan table (key event, source, cand_id = row
index within the event of the candidate file). Excluded from the
learned rankers: edge_px (synthetic capture points are placed at least
VERTEX_MARGIN_PX from the edge, a leak) and slice (absolute depth in the
stack). The training label is the capture vertex within LABEL_UM.
Efficiencies have a 95 % interval from a bootstrap over events; the
learned rankers are compared with darkness paired per event.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import StandardScaler

from . import efficiency as eff

BUDGETS = [1, 2, 4, 8, 16, 32, 64, 128, 256, 770, 1540]
TOLS = [6., 12.]
LABEL_UM = 6.
N_FOLDS = 5
NEG_FRACTION = 0.05     # negatives kept for training (all positives)
DIST_CAP_UM = 15.       # dist_other_um is inf beyond 15 um stage z
SEED = 1
N_BOOT = 2000
CI_PCT = (2.5, 97.5)
GBDT_ITER = 200
GBDT_LR = 0.05
LR_MAX_ITER = 1000
PAIRED_BUDGETS = (1, 2, 4, 8)
A_FEATURES = ["n_tracks_max", "n_slices", "mean_intens",
              "angle_spread_best", "angle_spread_max", "z_span_um"]
B_FEATURES = ["mult", "n_bp_20px"]
COMMON = ["darkness", "dark_rank_view", "dist_other_um"]
LOG_FEATURES = ["n_tracks_max", "dark_rank_view", "n_bp_20px", "mult"]


def ll_decays(verts):
  return [v for v in verts if v["type"] == "HypernucleusDecay"
          and v["parent"].startswith(eff.LL_PREFIX)]


def labels(events, cands):
  """Per event: distance (um) of each candidate to the capture vertex
  and to the nearest LL decay vertex."""
  out = {}
  for name, truth in events.items():
    cs = cands.get(name, [])
    if not cs:
      continue
    cap = [v for v in truth["vertices"] if v["type"] == "capture"]
    d_cap = eff.distances(truth, cap, cs, "developed")["inclusive"][0]
    ll = ll_decays(truth["vertices"])
    d_ll = eff.distances(truth, ll, cs, "developed")["inclusive"].min(0) \
      if ll else np.full(len(cs), np.inf)
    out[name] = (d_cap, d_ll)
  return out


def pool(events, a, b, feats):
  """One row per candidate of both sources with distances and features."""
  la, lb = labels(events, a), labels(events, b)
  rows = []
  for source, lab in (("A", la), ("Bprime", lb)):
    for name, (d_cap, d_ll) in lab.items():
      rows.append(pd.DataFrame(dict(
        event=name, source=source, cand_id=np.arange(len(d_cap)),
        d_cap=d_cap, d_ll=d_ll,
        view=events[name].get("background_view") or name)))
  df = pd.concat(rows, ignore_index=True)
  df = df.merge(feats, on=["event", "source", "cand_id"], how="left",
                validate="one_to_one")
  missing = df["darkness"].isna().sum()
  if missing:
    raise SystemExit(f"{missing} candidates without features")
  return df


def design(df):
  """Feature matrix of the learned rankers."""
  x = pd.DataFrame({"is_a": (df["source"] == "A").astype(float)})
  for c in COMMON + A_FEATURES + B_FEATURES:
    v = df[c].astype(float)
    if c == "dist_other_um":
      x["no_other"] = (~np.isfinite(v)).astype(float)
      v = v.clip(upper=DIST_CAP_UM).where(np.isfinite(v), DIST_CAP_UM)
    if c in LOG_FEATURES:
      v = np.log1p(v.clip(lower=0.))
    x[c] = v.fillna(0.)
  return x


def learned_scores(df, model):
  """Out-of-fold scores, folds grouped by background view."""
  x, y = design(df), (df["d_cap"] <= LABEL_UM).to_numpy()
  rng = np.random.default_rng(SEED)
  scores = np.zeros(len(df))
  for tr, te in GroupKFold(N_FOLDS).split(x, y, df["view"]):
    keep = tr[y[tr] | (rng.random(len(tr)) < NEG_FRACTION)]
    if model == "logistic":
      sc = StandardScaler().fit(x.iloc[keep])
      m = LogisticRegression(max_iter=LR_MAX_ITER, class_weight="balanced")
      m.fit(sc.transform(x.iloc[keep]), y[keep])
      scores[te] = m.decision_function(sc.transform(x.iloc[te]))
    else:
      m = HistGradientBoostingClassifier(
        max_iter=GBDT_ITER, learning_rate=GBDT_LR, random_state=SEED,
        class_weight="balanced")
      m.fit(x.iloc[keep], y[keep])
      scores[te] = m.predict_proba(x.iloc[te])[:, 1]
  return scores


def interleaved(df):
  """A by n_tracks_max and B' by darkness, taken alternately."""
  s = np.zeros(len(df))
  for src, key in (("A", "n_tracks_max"), ("Bprime", "darkness")):
    m = (df["source"] == src).to_numpy()
    r = df.loc[m].groupby("event")[key].rank(ascending=False,
                                              method="first")
    s[m] = -(2. * r.to_numpy() - (0. if src == "A" else 1.))
  return s


def found_at(df, score, dist_col, tol):
  """events x budgets: vertex within tol among the top K by score."""
  d = df.assign(_s=score)
  d = d[np.isfinite(d["_s"])]
  d["_r"] = d.groupby("event")["_s"].rank(ascending=False, method="first")
  hit = d[d[dist_col] <= tol].groupby("event")["_r"].min()
  ev = sorted(df["event"].unique())
  best = hit.reindex(ev).fillna(np.inf).to_numpy()
  return np.array([[b <= k for k in BUDGETS] for b in best], dtype=float)


def boot_ci(f, rng):
  n = len(f)
  b = [f[rng.integers(0, n, n)].mean(0) for _ in range(N_BOOT)]
  return np.percentile(b, CI_PCT, axis=0)


def main():
  parser = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument("dataset")
  parser.add_argument("--a", required=True)
  parser.add_argument("--b", required=True)
  parser.add_argument("--features", required=True)
  parser.add_argument("-o", "--output")
  args = parser.parse_args()

  events = eff.load_truth(args.dataset)
  a, b = eff.load_candidates(args.a), eff.load_candidates(args.b)
  feats = pd.read_parquet(args.features)
  df = pool(events, a, b, feats)
  rng = np.random.default_rng(SEED)
  rankers = {
    "random": rng.random(len(df)),
    "darkness": df["darkness"].to_numpy(),
    "A tracks": np.where(df["source"] == "A", df["n_tracks_max"], -np.inf),
    "interleaved": interleaved(df),
    "logistic": learned_scores(df, "logistic"),
    "gbdt": learned_scores(df, "gbdt"),
  }
  name = Path(args.dataset).name
  n_ev = df["event"].nunique()
  per_view = df.groupby("event")["source"].value_counts().unstack().mean()
  pos = df.groupby("event")["d_cap"].apply(lambda d: (d <= LABEL_UM).sum())
  text = [f"# Candidate ranking on {name}", "",
          f"{n_ev} events; candidates per view: A {per_view['A']:.0f}, "
          f"B' {per_view['Bprime']:.0f}; candidates within {LABEL_UM:g} um "
          f"of the capture vertex per event: median {pos.median():.0f}, "
          f"none in {(pos == 0).sum()} events (no ranker can find these)."
          " Cells: fraction of events with the vertex among the top K of "
          "the view, 95 % bootstrap interval in the paired table.", ""]
  flags = {}
  for vtx, col in (("capture", "d_cap"), ("LL decay", "d_ll")):
    for tol in TOLS:
      text += [f"## {vtx} within {tol:g} um", "",
               "| ranker | " + " | ".join(f"K={k}" for k in BUDGETS)
               + " |", "|---|" + "---|" * len(BUDGETS)]
      for r, s in rankers.items():
        f = found_at(df, s, col, tol)
        flags[(vtx, tol, r)] = f
        text.append(f"| {r} | " + " | ".join(
          f"{v:.0%}" for v in f.mean(0)) + " |")
      text.append("")
  text += ["## Learned minus darkness (paired per event, capture)", "",
           "| ranker | tol | " + " | ".join(
             f"K={k}" for k in PAIRED_BUDGETS) + " |",
           "|---|---|" + "---|" * len(PAIRED_BUDGETS)]
  cols = [BUDGETS.index(k) for k in PAIRED_BUDGETS]
  for r in ("logistic", "gbdt"):
    for tol in TOLS:
      d = flags[("capture", tol, r)] - flags[("capture", tol, "darkness")]
      lo, hi = boot_ci(d, rng)
      text.append(f"| {r} | {tol:g} | " + " | ".join(
        f"{d.mean(0)[i]:+.0%} [{lo[i]:+.0%}, {hi[i]:+.0%}]"
        for i in cols) + " |")
  text.append("")
  out = "\n".join(text)
  print(out)
  if args.output:
    Path(args.output).write_text(out)


if __name__ == "__main__":
  main()
