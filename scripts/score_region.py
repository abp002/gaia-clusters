"""Score a region: HDBSCAN runs against Hunt & Reffert, as fixed in docs/protocol-<region>.md.

Auriga's protocol used HDBSCAN's default excess-of-mass selection; its leaf runs are post hoc (see
the README). Carina's protocol fixes leaf selection in advance and states what Auriga's leaf rates
predict for it.

uv run scripts/fetch_region.py <region>              # once: downloads data/
uv run scripts/score_region.py <region> --predict    # the prediction alone, before any clustering
uv run scripts/score_region.py <region>              # clusters (labels cached in data/), scores, draws
"""
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from gaia_clusters.cluster import features, hdbscan_runs
from gaia_clusters.predict import predict_share
from gaia_clusters.region import REGIONS
from gaia_clusters.score import score, unmatched

# Per region: each selection method with its standing under that region's protocol (the protocol's
# own first), the clusters named on the sky map and drawn as colour-magnitude diagrams (the three
# richest in Hunt & Reffert for Carina), and the region whose leaf rates predict this one.
SETUP = {
    "auriga": {"methods": {"eom": "Protocol: excess of mass (library default)", "leaf": "Post hoc: leaf selection"},
               "named": {"NGC_1960": "M36", "NGC_2099": "M37", "NGC_1912": "M38"}},
    "carina": {"methods": {"leaf": "Protocol: leaf selection", "eom": "Secondary: excess of mass (library default)"},
               "named": {"NGC_3532": "NGC 3532", "IC_2714": "IC 2714", "Trumpler_19": "Trumpler 19"},
               "predicted_from": "auriga"},
}
SIZES = [10, 20, 40, 80]
DISTANCE = pd.IntervalIndex.from_breaks([0, 1000, 2000, 3000, np.inf], closed="left")
DISTANCE_LABELS = ["< 1 kpc", "1–2 kpc", "2–3 kpc", "≥ 3 kpc"]
MEMBERS = pd.IntervalIndex.from_breaks([10, 30, 100, np.inf], closed="left")
MEMBERS_LABELS = ["10–29", "30–99", "≥ 100"]
INK, RECOVERED, MISSED, FIELD = "#222222", "#1f3b73", "#c8553d", "#9a9a9a"
METHOD_COLOUR = {"eom": "#8a8a8a", "leaf": RECOVERED}

name = sys.argv[1]
BOX, METHODS, NAMED = REGIONS[name], SETUP[name]["methods"], SETUP[name]["named"]
data, results, figures = Path("data"), Path("results"), Path("figures")
results.mkdir(exist_ok=True)
figures.mkdir(exist_ok=True)

stars = pd.read_parquet(data / f"{name}_stars.parquet")
members = pd.read_parquet(data / f"{name}_hr23_members.parquet")
clusters = pd.read_parquet(data / "hr23_clusters.parquet")
ref = clusters[BOX.contains(clusters.GLON, clusters.GLAT)].set_index("Name")


def labels_for(method: str) -> pd.DataFrame:
    path = data / f"{name}_labels_{method}.parquet"
    if not path.exists():
        t0 = time.time()
        X = features(stars, *BOX.centre, lon="GLON", lat="GLAT")
        runs = hdbscan_runs(X, SIZES, method=method, cache_dir=str(data / "hdbscan_cache"))
        pd.DataFrame({"Source": stars.Source.to_numpy()} | {str(k): v for k, v in runs.items()}).to_parquet(path)
        print(f"{method}: clustered {len(stars):,} stars in {time.time() - t0:.0f} s", flush=True)
    labels = pd.read_parquet(path)
    assert (labels.Source.to_numpy() == stars.Source.to_numpy()).all()
    return labels


def scores(labels: pd.DataFrame) -> pd.DataFrame:
    """One row per reference cluster, with the result of every run and of any of them."""
    table = ref[["GLON", "GLAT", "dist50", "N", "CST", "CMDCl50"]].copy()
    table["reliable"] = (table.CST >= 5) & (table.CMDCl50 >= 0.5)
    for k in SIZES:
        s = score(stars.Source, labels[str(k)], members, ref.index).set_index("Name")
        table["n"], table["recoverable"] = s.n, s.recoverable
        table[f"recovered_{k}"], table[f"label_{k}"] = s.recovered, s.label
        table[f"recovery_{k}"], table[f"purity_{k}"] = s.recovery.round(3), s.purity.round(3)
    table["recovered"] = table[[f"recovered_{k}" for k in SIZES]].any(axis=1)
    table["distance"] = pd.cut(table.dist50, DISTANCE).cat.rename_categories(DISTANCE_LABELS)
    table["members"] = pd.cut(table.n, MEMBERS).cat.rename_categories(MEMBERS_LABELS)
    return table


def share(df, col="recovered") -> str:
    k, n = int(df[col].sum()), len(df)
    return f"{k} of {n} ({k / n:.0%})" if n else "–"


def prediction(table: pd.DataFrame) -> list[str]:
    """What the reference region's leaf rates by members in our sample predict for this region's mix."""
    source = SETUP[name].get("predicted_from")
    if not source:
        return []
    past = pd.read_csv(results / f"{source}_scores_leaf.csv")
    past = past[past.recoverable].groupby("members").recovered.agg(["sum", "count"]).reindex(MEMBERS_LABELS)
    now = table[table.recoverable].members.value_counts().reindex(MEMBERS_LABELS)
    mean, lo, hi = predict_share(now.to_numpy(), past["sum"].to_numpy(), past["count"].to_numpy())
    rows = " · ".join(f"{c}: {n} here, {s} of {k} in {source.title()}" for c, n, s, k in
                      zip(MEMBERS_LABELS, now, past["sum"], past["count"]))
    return [f"Prediction from {source.title()}'s leaf rates by members in our sample ({rows}): "
            f"**{mean:.1%}**, 95 % range {lo:.1%} to {hi:.1%}."]


def failure(r: pd.Series) -> str:
    """How a missed cluster fails, judged on the run that came closest (highest min(recovery, purity))."""
    k = max(SIZES, key=lambda k: min(r[f"recovery_{k}"], r[f"purity_{k}"]))
    if r[f"recovery_{k}"] > 0.5:
        return "inside a bigger group (purity ≤ 0.5)"
    if r[f"purity_{k}"] > 0.5:
        return "only partly found (recovery ≤ 0.5)"
    return "not found"


def summary(table: pd.DataFrame, labels: pd.DataFrame) -> list[str]:
    rec = table[table.recoverable]
    lines = [f"**{share(rec)} recovered by at least one run.**", "",
             "| Run | All recoverable | Highly reliable | Clusters found | Largest | Noise | No majority owner |",
             "|---|---|---|---|---|---|---|",
             f"| Any of the four | {share(rec)} | {share(rec[rec.reliable])} | | | | |"]
    for k in SIZES:
        lab = labels[str(k)]
        sizes = lab[lab >= 0].value_counts()
        lines.append(f"| min_cluster_size {k} | {share(rec, f'recovered_{k}')} | "
                     f"{share(rec[rec.reliable], f'recovered_{k}')} | {len(sizes):,} | {sizes.max():,} | "
                     f"{(lab < 0).mean():.0%} | {len(unmatched(stars.Source, lab, members)):,} |")
    lines += ["", "| Distance | Recovered (any run) | Highly reliable | Under 10 members after our cuts |",
              "|---|---|---|---|"]
    for c in DISTANCE_LABELS:
        lines.append(f"| {c} | {share(rec[rec.distance == c])} | {share(rec[(rec.distance == c) & rec.reliable])} | "
                     f"{int(((table.distance == c) & ~table.recoverable).sum())} |")
    lines += ["", "| Members in our sample | Recovered (any run) | Highly reliable |", "|---|---|---|"]
    for c in MEMBERS_LABELS:
        lines.append(f"| {c} | {share(rec[rec.members == c])} | {share(rec[(rec.members == c) & rec.reliable])} |")
    missed = rec[~rec.recovered]
    if len(missed):
        lines += ["", "| How the missed clusters fail (closest run) | Clusters | Median members in our sample |",
                  "|---|---|---|"]
        for how, g in missed.groupby(missed.apply(failure, axis=1)):
            lines.append(f"| {how} | {len(g)} | {int(g.n.median())} |")
    return lines


if "--predict" in sys.argv:  # n and the members bins do not depend on the labels: no clustering needed
    table = scores(pd.DataFrame(-1, index=stars.index, columns=[str(k) for k in SIZES]))
    print(f"Stars: {len(stars):,}. Reference clusters: {len(table)}. Recoverable: {int(table.recoverable.sum())}.")
    print("\n".join(prediction(table)))
    sys.exit()

labels = {m: labels_for(m) for m in METHODS}
tables = {m: scores(labels[m]) for m in METHODS}
for m, t in tables.items():
    t.to_csv(results / f"{name}_scores_{m}.csv")

first = tables[next(iter(METHODS))]
md = [f"Stars: {len(stars):,}. Reference clusters (centre in the box): {len(first)}. "
      f"Recoverable (≥ 10 members in our sample): {int(first.recoverable.sum())}."]
if expected := prediction(first):
    rec = first[first.recoverable]
    md += ["", expected[0] + f" Result: **{rec.recovered.mean():.1%}**."]
for m, title in METHODS.items():
    md += ["", f"### {title}", ""] + summary(tables[m], labels[m])
(results / f"{name}.md").write_text("\n".join(md) + "\n")
print("\n".join(md))

# Figure: the sky for each selection method, and recovery by distance and by members
plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK, "axes.labelcolor": INK, "xtick.color": INK,
                     "ytick.color": INK, "axes.spines.top": False, "axes.spines.right": False})
fig = plt.figure(figsize=(11, 10.5))
gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 0.8], hspace=0.5, wspace=0.18)
for col_, (m, title) in enumerate(METHODS.items()):
    t = tables[m]
    ax = fig.add_subplot(gs[0, col_])
    ax.hexbin(stars.GLON, stars.GLAT, gridsize=(200, 100), bins="log", cmap="Greys", mincnt=1, linewidths=0,
              alpha=0.55)
    hit, miss, lost = t[t.recovered], t[t.recoverable & ~t.recovered], t[~t.recoverable]
    ax.scatter(hit.GLON, hit.GLAT, s=8 + 2.2 * np.sqrt(hit.n), c=RECOVERED, label=f"recovered ({len(hit)})", zorder=3)
    ax.scatter(miss.GLON, miss.GLAT, s=8 + 2.2 * np.sqrt(miss.n), facecolors="white", edgecolors=MISSED,
               linewidths=1.3, label=f"missed ({len(miss)})", zorder=3)
    ax.scatter(lost.GLON, lost.GLAT, s=14, marker="x", c=FIELD, linewidths=0.9,
               label=f"< 10 members after cuts ({len(lost)})", zorder=3)
    for cl, label in NAMED.items():
        ax.annotate(label, (t.loc[cl, "GLON"], t.loc[cl, "GLAT"]), xytext=(6, 6), textcoords="offset points",
                    fontsize=9, color=INK)
    ax.set(xlim=(BOX.l_max, BOX.l_min), ylim=(BOX.b_min, BOX.b_max), aspect="equal",
           xlabel="Galactic longitude l (°)", ylabel="Galactic latitude b (°)")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.13), frameon=False, ncols=3, fontsize=8,
              handletextpad=0.3, columnspacing=1.0)
    rec = t[t.recoverable]
    ax.set_title(f"{title}\n{share(rec)} recovered", loc="left", fontsize=10.5, color=INK)

for slot, (title, col, cats) in enumerate([("By distance", "distance", DISTANCE_LABELS),
                                           ("By members in our sample", "members", MEMBERS_LABELS)]):
    a = fig.add_subplot(gs[1, slot])
    x = np.arange(len(cats))
    for j, m in enumerate(METHODS):
        rec = tables[m][tables[m].recoverable]
        groups = [rec[rec[col] == c] for c in cats]
        rate = np.array([g.recovered.mean() if len(g) else 0 for g in groups])
        xs = x + (j - 0.5) * 0.36
        a.bar(xs, rate, width=0.34, color=METHOD_COLOUR[m], label=METHODS[m].split(":")[0].lower() + f" ({m})")
        for xi, r, g in zip(xs, rate, groups):
            a.text(xi, r + 0.03, f"{int(g.recovered.sum())}/{len(g)}", ha="center", fontsize=7.5, color=INK)
    a.set(xticks=x, xticklabels=cats, ylim=(0, 1.12), yticks=[0, 0.25, 0.5, 0.75, 1],
          yticklabels=["0", "25 %", "50 %", "75 %", "100 %"])
    a.set_title(title, loc="left", fontsize=10, color=INK)
    if slot == 1:
        a.legend(frameon=False, fontsize=8, loc="upper left")
fig.savefig(figures / f"{name}.png", dpi=130, bbox_inches="tight")

# Colour-magnitude diagrams from the leaf runs: the named clusters, and the farthest recovered cluster
# with at least 30 members in our sample (fewer do not draw a sequence)
t, lab = tables["leaf"], labels["leaf"]
hit = t[t.recovered & (t.n >= 30)]
far = [hit.dist50.idxmax()] if len(hit) else []
picks = list(NAMED) + [cl for cl in far if cl not in NAMED]
ref_ids = members.groupby("Name").GaiaDR3.apply(set)
fig, axes = plt.subplots(1, len(picks), figsize=(3.6 * len(picks), 4.2))
for a, cl in zip(np.atleast_1d(axes), picks):
    r = t.loc[cl]
    k = max(SIZES, key=lambda k: (r[f"recovered_{k}"], r[f"recovery_{k}"] * r[f"purity_{k}"]))
    ours = stars[lab[str(k)].to_numpy() == r[f"label_{k}"]] if r[f"label_{k}"] >= 0 else stars.iloc[:0]
    theirs = ours.Source.isin(ref_ids.get(cl, set()))
    a.scatter(ours["BP-RP"][theirs], ours.Gmag[theirs], s=4, c=RECOVERED, label="also in Hunt & Reffert")
    a.scatter(ours["BP-RP"][~theirs], ours.Gmag[~theirs], s=4, c=MISSED, label="ours only")
    a.invert_yaxis()
    a.set_title(f"{NAMED.get(cl, cl.replace('_', ' '))} · {r.dist50 / 1000:.1f} kpc · "
                f"{'recovered' if r.recovered else 'missed'}\nleaf, min_cluster_size {k}: "
                f"recovery {r[f'recovery_{k}']:.0%}, purity {r[f'purity_{k}']:.0%}",
                loc="left", fontsize=8.5, color=INK)
    a.set(xlabel="BP − RP", ylabel="G")
np.atleast_1d(axes)[0].legend(frameon=False, fontsize=7.5, loc="upper right")
fig.tight_layout()
fig.savefig(figures / f"{name}_cmd.png", dpi=130, bbox_inches="tight")
