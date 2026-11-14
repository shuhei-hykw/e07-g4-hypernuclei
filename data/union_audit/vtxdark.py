"""Darkness at the true vertex, percentile among the B' points of the
view: embedded v25 nagara captures vs real specials at 0.507."""
import json
from pathlib import Path
import numpy as np, pandas as pd
from module import imagestats as st
from module.rankstudy import local_darkness
FS = Path.home() / "work/e07-fullscan/results"
def pct(img, xy, sl, feats_dark):
  stack = np.array(st.read_stack(img))
  d = local_darkness(stack, [xy], [sl])[0]
  return (feats_dark < d).mean(), d
out = []
f = pd.read_parquet(FS / "sim_eval/v25_embed_nagara_dz4_features.parquet")
f = f[f.source == "Bprime"]
for p in sorted(Path("data/images/v25_embed_nagara_dz4").glob("ev*/truth.json"))[:40]:
  t = json.loads(p.read_text())
  v = [x for x in t["vertices"] if x["type"] == "capture"][0]
  dk = f.loc[f.event == p.parent.name, "darkness"].to_numpy()
  out.append(("embedded",) + pct(p.parent / "image.json",
             (v["col"], v["row"]), v["slice"], dk))
f = pd.read_parquet(FS / "sim_eval/specials_0507_features.parquet")
f = f[f.source == "Bprime"]
for d in sorted((FS / "specials_0507").iterdir()):
  g = json.loads((d / "gt.json").read_text())
  dk = f.loc[f.event == d.name, "darkness"].to_numpy()
  r = pct(d / "image.json", (g["vx"], g["vy"]), g["z_slice"], dk)
  print(d.name, f"{r[0]:.1%} {r[1]:.1f}")
  out.append(("special",) + r)
df = pd.DataFrame(out, columns=["kind", "pct", "dark"])
print(df.groupby("kind").describe().T)
