"""Our isochrone ages against Hunt & Reffert's, as docs/protocol-ages.md reports them.

uv run scripts/fit_ages.py auriga && uv run scripts/fit_ages.py carina
uv run scripts/compare_ages.py    # results/ages.md, figures/ages.png
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

INK, RULE = "#222222", "#c4c5bf"
REGION_HEX = {"carina": "#2a78d6", "auriga": "#eda100"}  # the viewer's region colours
BANDS = pd.IntervalIndex.from_breaks([-np.inf, 7, 8, 9, np.inf], closed="left")
BAND_LABELS = ["under 10 Myr", "10–100 Myr", "0.1–1 Gyr", "over 1 Gyr"]
WITHIN = 0.3

results = Path("results")
t = pd.concat([pd.read_csv(results / f"{r}_ages.csv").assign(region=r) for r in REGION_HEX], ignore_index=True)
given = t[t.given]
d = given.delta


def share(mask) -> str:
    k, n = int(mask.sum()), len(mask)
    return f"{k} of {n} ({k / n:.0%})" if n else "–"


rho = spearmanr(given.age, given.logAge50).statistic
agree = (d.abs() <= WITHIN).mean() >= 2 / 3 and rho >= 0.8
md = [f"Test clusters: {len(t)} ({', '.join(f'{(t.region == r).sum()} in {r.title()}' for r in REGION_HEX)}). "
      f"Given an age: {share(t.given)}.", "",
      f"**Primary: {share(d.abs() <= WITHIN)} within {WITHIN} dex of Hunt & Reffert's median age; Spearman "
      f"ρ = {rho:.2f}. The protocol's bar (two thirds, and ρ ≥ 0.8) is {'met' if agree else 'not met'}.**", "",
      f"Median Δ (ours − theirs): {d.median():+.2f} dex; median |Δ|: {d.abs().median():.2f} dex.", "",
      "| Region | Given an age | Within 0.3 dex | Median Δ | Spearman ρ |", "|---|---|---|---|---|"]
for r in REGION_HEX:
    g = given[given.region == r]
    md.append(f"| {r.title()} | {share(t[t.region == r].given)} | {share(g.delta.abs() <= WITHIN)} | "
              f"{g.delta.median():+.2f} | {spearmanr(g.age, g.logAge50).statistic:.2f} |")
md += ["", "| Hunt & Reffert age band | Clusters | Within 0.3 dex | Median Δ |", "|---|---|---|---|"]
band = pd.cut(given.logAge50, BANDS).cat.rename_categories(BAND_LABELS)
for b in BAND_LABELS:
    g = given[band == b]
    md.append(f"| {b} | {len(g)} | {share(g.delta.abs() <= WITHIN)} | {g.delta.median():+.2f} |" if len(g) else f"| {b} | 0 | – | – |")
md += ["", "| Stars in our group | Clusters | Within 0.3 dex | Median Δ |", "|---|---|---|---|"]
for lo, hi in [(0, 50), (50, 100), (100, 300), (300, np.inf)]:
    g = given[(given.n >= lo) & (given.n < hi)]
    md.append(f"| {lo}–{'' if hi == np.inf else hi - 1} | {len(g)} | {share(g.delta.abs() <= WITHIN)} | "
              f"{g.delta.median():+.2f} |" if len(g) else f"| {lo}–{hi} | 0 | – | – |")
md += ["", f"Extinction: median ours − theirs {(given.av - given.AV50).median():+.2f} mag "
       f"(Spearman ρ {spearmanr(given.av, given.AV50).statistic:.2f}). Distance modulus: median ours − theirs "
       f"{(given.dm - given.MOD50).median():+.2f} mag.", "",
       f"Δ against distance: Spearman ρ {spearmanr(given.dist50, given.delta).statistic:+.2f}; against Hunt & "
       f"Reffert's extinction: {spearmanr(given.AV50, given.delta).statistic:+.2f}; against stars in our group: "
       f"{spearmanr(given.n, given.delta).statistic:+.2f}.", ""]
no = t[~t.given]
if len(no):
    md += ["| Why no age | Clusters |", "|---|---|"] + [f"| {why} | {n} |" for why, n in no.reason.value_counts().items()]
(results / "ages.md").write_text("\n".join(md) + "\n")
print("\n".join(md))

plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK, "axes.labelcolor": INK, "xtick.color": INK,
                     "ytick.color": INK, "axes.spines.top": False, "axes.spines.right": False})
fig, (a, b) = plt.subplots(1, 2, figsize=(11, 5), gridspec_kw={"width_ratios": [1.15, 1]})
lim = (6.5, 10.2)
a.fill_between(lim, np.subtract(lim, WITHIN), np.add(lim, WITHIN), color=RULE, alpha=0.35, lw=0,
               label=f"within {WITHIN} dex")
a.plot(lim, lim, c=INK, lw=0.8)
for r, colour in REGION_HEX.items():
    g = given[given.region == r]
    a.errorbar(g.logAge50, g.age, xerr=[g.logAge50 - g.logAge16, g.logAge84 - g.logAge50],
               yerr=[g.age - g.age16, g.age84 - g.age], fmt="o", ms=3.5, c=colour, ecolor=colour, elinewidth=0.6,
               alpha=0.85, label=f"{r.title()} ({len(g)})")
a.set(xlim=lim, ylim=lim, aspect="equal", xlabel="Hunt & Reffert: log10(age / yr)", ylabel="Ours: log10(age / yr)")
a.legend(frameon=False, loc="upper left")
a.set_title(f"{share(d.abs() <= WITHIN)} within {WITHIN} dex, ρ = {rho:.2f}", loc="left", fontsize=10.5)
for r, colour in REGION_HEX.items():
    g = given[given.region == r]
    b.scatter(g.n, g.delta, s=12, c=colour, alpha=0.85)
b.axhspan(-WITHIN, WITHIN, color=RULE, alpha=0.35, lw=0)
b.axhline(0, c=INK, lw=0.8)
b.set(xscale="log", xlabel="Stars in our group (5 ≤ G ≤ 18)", ylabel="Δ = ours − Hunt & Reffert (dex)")
b.set_title("Δ by the number of stars", loc="left", fontsize=10.5)
fig.tight_layout()
fig.savefig(Path("figures") / "ages.png", dpi=130, bbox_inches="tight")
