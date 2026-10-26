"""Write the files the 3D viewer loads for a scored region, into viewer/data/<region>/.

uv run scripts/fetch_region.py <region>    # stars, Bailer-Jones distances, Hunt & Reffert members
uv run scripts/score_region.py <region>    # leaf labels and scores
uv run scripts/sweep_region.py <region>    # optional: more leaf runs, for the min_cluster_size slider
uv run scripts/export_viewer.py <region>
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from gaia_clusters.region import REGIONS, SETUP
from gaia_clusters.score import DISTANCE_LABELS, MEMBERS_LABELS, score
from gaia_clusters.viewer import build, galactic_pm, rates, view, write

UNVERIFIED_RUN = 10  # the finest run; the README counts its groups that no Hunt & Reffert cluster owns

name = sys.argv[1]
box = REGIONS[name]
data = Path("data")
stars = pd.read_parquet(data / f"{name}_stars.parquet").merge(
    pd.read_parquet(data / f"{name}_distances.parquet"), on="Source", how="left")
stars["pm_l"], stars["pm_b"] = galactic_pm(stars.RA_ICRS, stars.DE_ICRS, stars.pmRA, stars.pmDE)
labels = pd.read_parquet(data / f"{name}_labels_leaf.parquet")
assert (labels.Source.to_numpy() == stars.Source.to_numpy()).all()
members = pd.read_parquet(data / f"{name}_hr23_members.parquet")
catalogue = pd.read_parquet(data / "hr23_clusters.parquet").set_index("Name")
table = pd.read_csv(Path("results") / f"{name}_scores_leaf.csv", index_col="Name").join(catalogue[["rtpc", "logAge16", "logAge50", "logAge84"]])
table["pm_l"], table["pm_b"] = galactic_pm(*catalogue.loc[table.index, ["RA_ICRS", "DE_ICRS", "pmRA", "pmDE"]].to_numpy().T)

built = build(stars, labels.drop(columns="Source"), members, table, SETUP[name]["named"], UNVERIFIED_RUN)
scored = table[table.recoverable]

# The slider's runs (sweep_region.py), each scored as score_region scores the protocol's
runs, sweep_path = {}, data / f"{name}_labels_sweep.parquet"
if sweep_path.exists():
    sweep = pd.read_parquet(sweep_path)
    assert (sweep.Source.to_numpy() == stars.Source.to_numpy()).all()
    sweep, by_run = sweep.drop(columns="Source"), scored.copy()
    for k in sweep.columns:
        s = score(stars.Source, sweep[k], members, scored.index).set_index("Name")
        by_run[f"recovered_{k}"], by_run[f"label_{k}"] = s.recovered, s.label
        by_run[f"recovery_{k}"], by_run[f"purity_{k}"] = s.recovery, s.purity
    runs = {int(k): view(stars, sweep, members, by_run, run=int(k)) for k in sweep.columns}
summary = {"region": name, "title": name.title(), "standing": SETUP[name]["methods"]["leaf"],
           "box": {"l": [box.l_min, box.l_max], "b": [box.b_min, box.b_max]},
           "recoverable": len(scored), "recovered": int(scored.recovered.sum()),
           "reliable": int(scored.reliable.sum()), "reliable_recovered": int((scored.reliable & scored.recovered).sum()),
           "unverified_run": UNVERIFIED_RUN,
           "rates": {"distance": rates(table, "distance", DISTANCE_LABELS), "members": rates(table, "members", MEMBERS_LABELS)}}
protocol = next(iter(SETUP[name]["methods"]))
if protocol != "leaf":  # the leaf runs are post hoc here: the page says what the protocol's own runs gave
    first = pd.read_csv(Path("results") / f"{name}_scores_{protocol}.csv")
    summary["protocol"] = {"selection": {"eom": "excess of mass"}[protocol],
                           "recovered": int(first[first.recoverable].recovered.sum())}
# Why members sit at their cluster's distance: how uncertain a member's own distance is, against how big
# a cluster is (Hunt & Reffert's radius holding half the members), both medians over recovered clusters
half_width = ((stars.r_hi_geo - stars.r_lo_geo) / 2).to_numpy()
summary["depth"] = {"sigma": round(float(np.median(half_width[built["member_index"]]))),
                    "r50": round(float(catalogue.r50pc.reindex(scored.index[scored.recovered]).median()), 1)}
summary["named"] = list(SETUP[name]["named"])
out = Path("viewer/data") / name
meta = write(out, built, summary, runs, protocol=[int(c) for c in labels.columns if c != "Source"])

print(f"{name}: {meta['stars']:,} stars, {stars.r_med_geo.isna().sum():,} without a distance; "
      f"field {meta['field']['count']:,} in {len(meta['field']['levels'])} levels, "
      f"{meta['members']['count']:,} members of {meta['recovered']} recovered clusters, "
      f"{meta['cmd']['count']:,} colour-magnitude points, {meta['unverified']['count']:,} unverified groups, "
      f"{len(runs)} single runs")
print(f"{sum(f.stat().st_size for f in out.iterdir()) / 1e6:.1f} MB in {out}")
