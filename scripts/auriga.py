"""Auriga region: four HDBSCAN runs scored against Hunt & Reffert, as fixed in docs/protocol-auriga.md.

uv run scripts/fetch_region.py   # once: downloads data/
uv run scripts/auriga.py         # clusters (cached in data/auriga_labels.parquet), scores, draws
"""
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from gaia_clusters.cluster import features, hdbscan_runs
from gaia_clusters.region import AURIGA
from gaia_clusters.score import score, unmatched

SIZES = [10, 20, 40, 80]
DISTANCE = pd.IntervalIndex.from_breaks([0, 1000, 2000, 3000, np.inf], closed="left")
DISTANCE_LABELS = ["< 1 kpc", "1–2 kpc", "2–3 kpc", "≥ 3 kpc"]
MEMBERS = pd.IntervalIndex.from_breaks([10, 30, 100, np.inf], closed="left")
MEMBERS_LABELS = ["10–29", "30–99", "≥ 100"]
MESSIER = {"NGC_1960": "M36", "NGC_2099": "M37", "NGC_1912": "M38"}
INK, RECOVERED, MISSED, FIELD = "#222222", "#1f3b73", "#c8553d", "#9a9a9a"

data, results, figures = Path("data"), Path("results"), Path("figures")
results.mkdir(exist_ok=True)
figures.mkdir(exist_ok=True)

stars = pd.read_parquet(data / "auriga_stars.parquet")
members = pd.read_parquet(data / "auriga_hr23_members.parquet")
clusters = pd.read_parquet(data / "hr23_clusters.parquet")
ref = clusters[AURIGA.contains(clusters.GLON, clusters.GLAT)].set_index("Name")

labels_path = data / "auriga_labels.parquet"
if labels_path.exists():
    labels = pd.read_parquet(labels_path)
else:
    t0 = time.time()
    runs = hdbscan_runs(features(stars, *AURIGA.centre, lon="GLON", lat="GLAT"), SIZES)
    labels = pd.DataFrame({"Source": stars.Source.to_numpy()} | {f"mcs{k}": v for k, v in runs.items()})
    labels.to_parquet(labels_path)
    print(f"clustered {len(stars):,} stars in {time.time() - t0:.0f} s")
assert (labels.Source.to_numpy() == stars.Source.to_numpy()).all()

# One row per reference cluster
table = ref[["GLON", "GLAT", "dist50", "N", "CST", "CMDCl50"]].copy()
table["reliable"] = (table.CST >= 5) & (table.CMDCl50 >= 0.5)
for k in SIZES:
    s = score(stars.Source, labels[f"mcs{k}"], members, ref.index).set_index("Name")
    table["n"], table["recoverable"] = s.n, s.recoverable
    table[f"recovered_{k}"], table[f"label_{k}"] = s.recovered, s.label
    table[f"recovery_{k}"], table[f"purity_{k}"] = s.recovery.round(3), s.purity.round(3)
table["recovered"] = table[[f"recovered_{k}" for k in SIZES]].any(axis=1)
table["distance"] = pd.cut(table.dist50, DISTANCE).cat.rename_categories(DISTANCE_LABELS)
table["members"] = pd.cut(table.n, MEMBERS).cat.rename_categories(MEMBERS_LABELS)
table.to_csv(results / "auriga_scores.csv")
rec = table[table.recoverable]


def share(df, col="recovered") -> str:
    k, n = int(df[col].sum()), len(df)
    return f"{k} of {n} ({k / n:.0%})" if n else "–"


lines = [
    f"Stars: {len(stars):,}. Reference clusters (centre in the box): {len(table)}. "
    f"Recoverable (≥ 10 members in our sample): {len(rec)}.",
    "",
    f"**Primary: {share(rec)} recovered by at least one run.**",
    "",
    "| Run | All recoverable | Highly reliable | Clusters found | Noise | No majority owner |",
    "|---|---|---|---|---|---|",
    f"| Any of the four | {share(rec)} | {share(rec[rec.reliable])} | | | |",
]
for k in SIZES:
    lab = labels[f"mcs{k}"]
    lines.append(f"| min_cluster_size {k} | {share(rec, f'recovered_{k}')} | "
                 f"{share(rec[rec.reliable], f'recovered_{k}')} | {lab.max() + 1:,} | {(lab < 0).mean():.0%} | "
                 f"{len(unmatched(stars.Source, lab, members)):,} |")
lines += ["", "| Distance | Recovered (any run) | Highly reliable | Under 10 members after our cuts |",
          "|---|---|---|---|"]
for c in DISTANCE_LABELS:
    lines.append(f"| {c} | {share(rec[rec.distance == c])} | {share(rec[(rec.distance == c) & rec.reliable])} | "
                 f"{int(((table.distance == c) & ~table.recoverable).sum())} |")
lines += ["", "| Members in our sample | Recovered (any run) | Highly reliable |", "|---|---|---|"]
for c in MEMBERS_LABELS:
    lines.append(f"| {c} | {share(rec[rec.members == c])} | {share(rec[(rec.members == c) & rec.reliable])} |")
summary = "\n".join(lines)
(results / "auriga.md").write_text(summary + "\n")
print(summary)

# Figure: the sky with every reference cluster, and recovery by distance and by members
plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK, "axes.labelcolor": INK, "xtick.color": INK,
                     "ytick.color": INK, "axes.spines.top": False, "axes.spines.right": False})
fig = plt.figure(figsize=(11, 8.2))
gs = fig.add_gridspec(2, 2, height_ratios=[1.35, 1], hspace=0.35, wspace=0.25)
ax = fig.add_subplot(gs[0, :])
ax.hexbin(stars.GLON, stars.GLAT, gridsize=(200, 100), bins="log", cmap="Greys", mincnt=1, linewidths=0, alpha=0.55)
miss = table[table.recoverable & ~table.recovered]
lost = table[~table.recoverable]
hit = table[table.recovered]
ax.scatter(hit.GLON, hit.GLAT, s=8 + 2.2 * np.sqrt(hit.n), c=RECOVERED, label=f"recovered ({len(hit)})", zorder=3)
ax.scatter(miss.GLON, miss.GLAT, s=8 + 2.2 * np.sqrt(miss.n), facecolors="white", edgecolors=MISSED, linewidths=1.3,
           label=f"missed ({len(miss)})", zorder=3)
ax.scatter(lost.GLON, lost.GLAT, s=14, marker="x", c=FIELD, linewidths=0.9,
           label=f"under 10 members after our cuts ({len(lost)})", zorder=3)
for name, messier in MESSIER.items():
    r = table.loc[name]
    ax.annotate(messier, (r.GLON, r.GLAT), xytext=(6, 6), textcoords="offset points", fontsize=9, color=INK)
ax.set(xlim=(AURIGA.l_max, AURIGA.l_min), ylim=(AURIGA.b_min, AURIGA.b_max), aspect="equal",
       xlabel="Galactic longitude l (°)", ylabel="Galactic latitude b (°)")
ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False, title="Hunt & Reffert clusters\n(size: members)",
          alignment="left")
ax.set_title(f"Auriga, l 170–180°: {share(rec)} of the recoverable Hunt & Reffert clusters recovered",
             loc="left", fontsize=11, color=INK)

for slot, (title, col, cats) in enumerate([("Distance", "distance", DISTANCE_LABELS),
                                           ("Members in our sample", "members", MEMBERS_LABELS)]):
    a = fig.add_subplot(gs[1, slot])
    x = np.arange(len(cats))
    groups = [rec[rec[col] == c] for c in cats]
    rate = [g.recovered.mean() if len(g) else 0 for g in groups]
    a.bar(x, rate, width=0.6, color=RECOVERED, alpha=0.85, label="any run")
    for i, g in enumerate(groups):
        a.text(i, rate[i] + 0.03, f"{int(g.recovered.sum())}/{len(g)}", ha="center", fontsize=9, color=INK)
        for j, k in enumerate(SIZES):
            a.plot(i - 0.21 + 0.14 * j, g[f"recovered_{k}"].mean() if len(g) else 0, "o", ms=3.5,
                   color=INK, alpha=0.35 + 0.15 * j, label=f"min_cluster_size {k}" if i == 0 else None)
    a.set(xticks=x, xticklabels=cats, ylim=(0, 1.12), yticks=[0, 0.25, 0.5, 0.75, 1],
          yticklabels=["0", "25 %", "50 %", "75 %", "100 %"], title=title)
    a.title.set(loc="left", fontsize=10, color=INK)
    if slot == 1:
        a.legend(frameon=False, fontsize=7.5, loc="upper left")
fig.savefig(figures / "auriga.png", dpi=130, bbox_inches="tight")

# Colour-magnitude diagrams of the Messier clusters and the farthest cluster recovered, from the run that did best
picks = list(MESSIER) + [hit.dist50.idxmax()] if len(hit) else list(MESSIER)
fig, axes = plt.subplots(1, len(picks), figsize=(3.1 * len(picks), 4), sharey=False)
ref_ids = members.groupby("Name").GaiaDR3.apply(set)
for a, name in zip(np.atleast_1d(axes), picks):
    r = table.loc[name]
    k = max(SIZES, key=lambda k: (r[f"recovered_{k}"], r[f"recovery_{k}"] * r[f"purity_{k}"]))
    ours = stars[labels[f"mcs{k}"].to_numpy() == r[f"label_{k}"]] if r[f"label_{k}"] >= 0 else stars.iloc[:0]
    theirs = ours.Source.isin(ref_ids.get(name, set()))
    a.scatter(ours["BP-RP"][theirs], ours.Gmag[theirs], s=4, c=RECOVERED, label="also in Hunt & Reffert")
    a.scatter(ours["BP-RP"][~theirs], ours.Gmag[~theirs], s=4, c=MISSED, label="ours only")
    a.invert_yaxis()
    status = "recovered" if r.recovered else "missed"
    a.set_title(f"{MESSIER.get(name, name.replace('_', ' '))} · {r.dist50 / 1000:.1f} kpc · {status}\n"
                f"min_cluster_size {k}: recovery {r[f'recovery_{k}']:.0%}, purity {r[f'purity_{k}']:.0%}",
                loc="left", fontsize=8.5, color=INK)
    a.set(xlabel="BP − RP", ylabel="G")
np.atleast_1d(axes)[0].legend(frameon=False, fontsize=7.5, loc="lower left")
fig.tight_layout()
fig.savefig(figures / "auriga_cmd.png", dpi=130, bbox_inches="tight")
