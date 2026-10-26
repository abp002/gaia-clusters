"""More leaf runs of a scored region, from min_cluster_size 10 to 240, for the viewer's slider.

Exploration, not part of the score: the score is the four runs the protocols fix (score_region.py).
The minimum spanning tree is cached by score_region.py, so each size only re-cuts it.

uv run scripts/score_region.py <region>    # first: the four protocol runs and the cached tree
uv run scripts/sweep_region.py <region>    # labels in data/<region>_labels_sweep.parquet, summary in results/
"""
import sys
import time
from pathlib import Path

import pandas as pd

from gaia_clusters.cluster import features, hdbscan_runs
from gaia_clusters.region import REGIONS
from gaia_clusters.score import score, unmatched

SWEEP = [10, 15, 20, 30, 40, 60, 80, 120, 160, 240]

name = sys.argv[1]
box = REGIONS[name]
data, results = Path("data"), Path("results")
stars = pd.read_parquet(data / f"{name}_stars.parquet")
members = pd.read_parquet(data / f"{name}_hr23_members.parquet")
path = data / f"{name}_labels_sweep.parquet"

labels = pd.read_parquet(path) if path.exists() else pd.read_parquet(data / f"{name}_labels_leaf.parquet")
assert (labels.Source.to_numpy() == stars.Source.to_numpy()).all()
todo = [k for k in SWEEP if str(k) not in labels.columns]
if todo:
    t0 = time.time()
    X = features(stars, *box.centre, lon="GLON", lat="GLAT")
    runs = hdbscan_runs(X, todo, method="leaf", cache_dir=str(data / "hdbscan_cache"))
    labels = labels.assign(**{str(k): v for k, v in runs.items()})
    print(f"clustered {len(todo)} sizes of {len(stars):,} stars in {time.time() - t0:.0f} s", flush=True)
labels = labels[["Source"] + [str(k) for k in SWEEP]]
labels.to_parquet(path)

ref = pd.read_csv(results / f"{name}_scores_leaf.csv", index_col="Name")
ref = ref[ref.recoverable].index
rows = []
for k in SWEEP:
    lab = labels[str(k)]
    s = score(stars.Source, lab, members, ref)
    rows.append({"min_cluster_size": k, "recovered": int(s.recovered.sum()), "recoverable": len(ref),
                 "groups": int(lab.max() + 1), "unverified": len(unmatched(stars.Source, lab, members)),
                 "noise": round(float((lab < 0).mean()), 4)})
table = pd.DataFrame(rows)
table.to_csv(results / f"{name}_sweep.csv", index=False)
print(table.to_string(index=False))
