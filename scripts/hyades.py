"""Select the Hyades' structure and score its tails against the published ones (docs/protocol-hyades.md).

uv run scripts/fetch_hyades.py    # first
uv run scripts/hyades.py          # results/hyades.md, results/hyades_members.csv, figures/hyades.png
"""
from pathlib import Path

import astropy.units as u
import hdbscan
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord

from gaia_clusters.hyades import HYADES, basis, heliocentric_icrs, residual, selected, space_velocity

INK, FIELD, OURS, THEIRS = "#222222", "#b4b4ae", "#1f3b73", "#c8553d"
data, results, figures = Path("data"), Path("results"), Path("figures")
V = space_velocity(*(HYADES[k] for k in ("ra", "dec", "plx", "pmra", "pmdec", "rv")))
centre = basis(HYADES["ra"], HYADES["dec"])[0] * HYADES["dist"]
tail_from = 2 * HYADES["rt"]
ahead = SkyCoord(l=90 * u.deg, b=0 * u.deg, frame="galactic").icrs.cartesian.xyz.value  # Galactic rotation, ICRS

sample = pd.read_parquet(data / "hyades_sample.parquet")
refs = pd.read_parquet(data / "hyades_references.parquet")  # every reference star, from DR3 by ID, any velocity
jer = set(pd.read_parquet(data / "hyades_jerabkova.parquet").GaiaEDR3)
roser = pd.read_parquet(data / "hyades_roser.parquet")
gcns = pd.read_parquet(data / "hyades_gcns.parquet")
gcns = set(gcns.GaiaEDR3[gcns.Cluster == "Hyades"])


def place(t):
    xyz = heliocentric_icrs(t.RA_ICRS, t.DE_ICRS, t.Plx)
    return xyz, np.linalg.norm(xyz - centre, axis=1), (xyz - centre) @ ahead


# 1. Convergent point
dv, sigma = residual(sample.RA_ICRS, sample.DE_ICRS, sample.Plx, sample.e_Plx, sample.pmRA, sample.pmDE, V)
kept = sample[selected(dv, sigma)].reset_index(drop=True)
xyz, r, along = place(kept)
# 2. HDBSCAN on position; the Hyades' group holds the most GCNS members
labels = hdbscan.HDBSCAN(min_cluster_size=20, min_samples=10, cluster_selection_method="eom").fit_predict(xyz)
named = pd.Series(labels[kept.Source.isin(gcns).to_numpy() & (labels >= 0)]).value_counts()
ours = labels == named.index[0]
tail = ours & (r > tail_from)
print(f"velocity cut keeps {len(kept):,} of {len(sample):,}; {labels.max() + 1} groups; the Hyades' holds {ours.sum()}, "
      f"{tail.sum()} beyond {tail_from:.1f} pc", flush=True)

# 3. Scores. A reference star counts if it passes the sample's cuts (refs carry them), and is a tail star past 17.2 pc
refs["xyz"], refs["r"], refs["along"] = list(place(refs)[0]), *place(refs)[1:]
base = refs[(refs.Plx > 2) & (refs.Plx > 10 * refs.e_Plx) & (refs.RUWE < 1.4)].set_index("Source")
found_tail, found = set(kept.Source[tail]), set(kept.Source[ours])


def recovery(ids, tail_only=True):
    ids = [i for i in ids if i in base.index and (not tail_only or base.r[i] > tail_from)]
    hit = sum(i in (found_tail if tail_only else found) for i in ids)
    return hit, len(ids)


def share(k, n):
    return f"{k} of {n} ({k / n:.0%})" if n else "–"


jk, jn = recovery(jer)
within = kept.Source[tail & (np.linalg.norm(xyz, axis=1) <= 500)]
pure = within.isin(jer).sum()
reproduced = jn and jk / jn >= 0.5 and pure / max(len(within), 1) >= 0.5
dr2 = set(roser.Source)
matched = len(dr2 & set(refs.Source))
lead, trail = (roser.Source[roser.Member == m] for m in ("Pr_tail", "Tr_tail"))


def reach(ids):
    """95th percentile of the distance from the centre, for stars ahead of it and behind it."""
    t = base.loc[[i for i in ids if i in base.index]]
    t = t[t.r > tail_from]
    return [np.percentile(t.r[side], 95) if side.any() else np.nan for side in (t.along > 0, t.along <= 0)]


ours_reach = [np.percentile(r[tail & side], 95) if (tail & side).any() else np.nan for side in (along > 0, along <= 0)]
jer_reach = reach(jer)
md = [f"Sample: {len(sample):,} stars within 10 km/s of the Hyades' tangential velocity (500 pc, all sky). "
      f"The velocity cut keeps {len(kept):,}; HDBSCAN finds {labels.max() + 1} groups; the Hyades' group holds "
      f"{ours.sum()} stars, {tail.sum()} of them more than {tail_from:.1f} pc from the centre.", "",
      f"**Primary: {share(jk, jn)} of Jerabkova et al.'s tail stars recovered; {share(pure, len(within))} of our "
      f"tail stars within 500 pc are among their 862. The protocol's bar (both at least half) is "
      f"{'met' if reproduced else 'not met'}.**", "",
      "| Reference | Recovered |", "|---|---|",
      f"| Jerabkova et al. (2021), tail | {share(jk, jn)} |",
      f"| Röser et al. (2019), leading tail | {share(*recovery(lead))} |",
      f"| Röser et al. (2019), trailing tail | {share(*recovery(trail))} |",
      f"| GCNS core (within {tail_from:.1f} pc) | {share(*recovery([i for i in gcns if i in base.index and base.r[i] <= tail_from], tail_only=False))} |",
      "", f"Röser et al.'s DR2 source IDs found in DR3 under the same ID: {share(matched, len(dr2))}.", "",
      "| 95th percentile of the distance from the centre | Ahead (leading) | Behind (trailing) |", "|---|---|---|",
      f"| Ours | {ours_reach[0]:.0f} pc | {ours_reach[1]:.0f} pc |",
      f"| Jerabkova et al. | {jer_reach[0]:.0f} pc | {jer_reach[1]:.0f} pc |"]
# Diagnostic, added after the run (not in the protocol): how close Jerabkova et al.'s tail stars are to the
# cluster's tangential velocity, and how many field stars are as close, by distance from the Sun
dist, is_jer = 1000 / sample.Plx, sample.Source.isin(jer)
md += ["", "Diagnostic, after the run: Δ is the difference from the Hyades' tangential velocity at the star's position.", "",
       "| Distance from the Sun | Jerabkova et al. tail stars | Their median Δ | Other stars with Δ < 3 km/s |",
       "|---|---|---|---|"]
for lo, hi in [(0, 100), (100, 200), (200, 300), (300, 500)]:
    m = (dist >= lo) & (dist < hi)
    md.append(f"| {lo}–{hi} pc | {int((m & is_jer).sum())} | {np.median(dv[m & is_jer]):.1f} km/s | {int((m & ~is_jer & (dv < 3)).sum()):,} |")
(results / "hyades.md").write_text("\n".join(md) + "\n")
print("\n".join(md))
kept.assign(label=labels, tail=tail, r=r.round(2), along=along.round(2)).loc[ours].to_csv(results / "hyades_members.csv", index=False)

# Figure: from above the Galactic plane, and edge-on, in Galactic Cartesian coordinates around the Sun
gal = lambda p: SkyCoord(x=p[:, 0], y=p[:, 1], z=p[:, 2], unit="pc", representation_type="cartesian",
                         frame="icrs").galactic.cartesian.xyz.value.T
g_kept, g_jer = gal(xyz), gal(np.stack(base.loc[[i for i in jer if i in base.index]].xyz.to_numpy()))
plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK, "axes.labelcolor": INK, "xtick.color": INK,
                     "ytick.color": INK, "axes.spines.top": False, "axes.spines.right": False})
fig, axes = plt.subplots(2, 1, figsize=(10, 9), gridspec_kw={"height_ratios": [2.2, 1]})
for a, (i, j, label) in zip(axes, [(0, 1, "y (pc), towards l = 90°"), (0, 2, "z (pc), north Galactic pole")]):
    a.scatter(g_kept[~ours, i], g_kept[~ours, j], s=1, c=FIELD, lw=0, label=f"velocity cut, other groups and noise ({(~ours).sum():,})")
    a.scatter(g_kept[ours, i], g_kept[ours, j], s=1, c=OURS, lw=0, alpha=0.25, label=f"our Hyades group ({ours.sum():,})")
    a.scatter(g_jer[:, i], g_jer[:, j], s=12, facecolors="none", edgecolors=THEIRS, lw=0.8, zorder=3,
              label=f"Jerabkova et al. (2021) tail members ({len(g_jer)})")
    a.scatter([0], [0], marker="+", c=INK, s=40)
    a.set(xlabel="x (pc), towards the Galactic centre", ylabel=label, aspect="equal")
axes[0].legend(frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1), markerscale=3)
axes[0].set_title(f"The velocity cut keeps {len(kept):,} stars; excess of mass puts {ours.sum():,} of them in one group",
                  loc="left", fontsize=10.5)
fig.tight_layout()
fig.savefig(figures / "hyades.png", dpi=130, bbox_inches="tight")
