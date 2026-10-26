"""Fit PARSEC isochrones to each recovered cluster's group, as fixed in docs/protocol-ages.md.

uv run scripts/fit_ages.py <region>                  # the test clusters: results/<region>_ages.csv
uv run scripts/fit_ages.py <region> --calibration    # only the calibration clusters, with a figure

Needs score_region.py's labels and scores. The isochrones are downloaded once into data/.
"""
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from joblib import Parallel, delayed

from gaia_clusters import fetch
from gaia_clusters.ages import A_BP, A_G, A_RP, Grid, fit
from gaia_clusters.score import shown_run

CALIBRATION = ["NGC_3532", "Trumpler_19", "NGC_2099"]
MIN_STARS = 30  # stars with a BP − RP colour in our group

name, calibration = sys.argv[1], "--calibration" in sys.argv
data, results, figures = Path("data"), Path("results"), Path("figures")
path = data / "parsec_solar.parquet"
if not path.exists():
    fetch.isochrones().to_parquet(path)
grid = Grid(pd.read_parquet(path))

stars = pd.read_parquet(data / f"{name}_stars.parquet")
labels = pd.read_parquet(data / f"{name}_labels_leaf.parquet")
assert (labels.Source.to_numpy() == stars.Source.to_numpy()).all()
sizes = [int(c) for c in labels.columns if c != "Source"]
table = pd.read_csv(results / f"{name}_scores_leaf.csv", index_col="Name")
catalogue = pd.read_parquet(data / "hr23_clusters.parquet").set_index("Name")


def group(r: pd.Series) -> np.ndarray:
    """Our group's stars with a colour: the run the cluster is shown with, as in the viewer."""
    k = shown_run(r, sizes)
    idx = np.flatnonzero(labels[str(k)].to_numpy() == r[f"label_{k}"])
    return idx[np.isfinite(stars["BP-RP"].to_numpy()[idx])]


rec = table[table.recovered & table.recoverable]
groups = {n: group(r) for n, r in rec.iterrows()}
picked = [n for n, g in groups.items() if len(g) >= MIN_STARS and (n in CALIBRATION) == calibration]
print(f"{name}: {len(rec)} recovered, {sum(len(g) >= MIN_STARS for g in groups.values())} with at least "
      f"{MIN_STARS} stars, {len(picked)} {'calibration' if calibration else 'test'} clusters", flush=True)


def one(n):
    s = stars.iloc[groups[n]]
    return n, fit(grid, s["BP-RP"].to_numpy(), s.Gmag.to_numpy(), s.Plx.to_numpy())


t0 = time.time()
fits = dict(Parallel(n_jobs=-1)(delayed(one)(n) for n in picked))
print(f"fitted in {time.time() - t0:.0f} s", flush=True)
out = pd.DataFrame([{"Name": n, "n": f.n, "age": f.age, "age16": f.lo, "age84": f.hi, "av": round(f.av, 3),
                     "dm": round(f.dm, 3), "given": f.given, "reason": f.reason} for n, f in fits.items()]).set_index("Name")
out = out.join(catalogue[["logAge16", "logAge50", "logAge84", "AV50", "MOD50", "dist50"]])
out["delta"] = (out.age - out.logAge50).round(3)
print(out.drop(columns=["dist50"]).to_string())

if not calibration:
    out.to_csv(results / f"{name}_ages.csv")
else:  # each calibration cluster's diagram with the best isochrone and Hunt & Reffert's
    fig, axes = plt.subplots(1, len(out), figsize=(4 * len(out), 5), squeeze=False)
    table_iso = pd.read_parquet(path)
    for a, (n, r) in zip(axes[0], out.iterrows()):
        s = stars.iloc[groups[n]]
        a.scatter(s["BP-RP"], s.Gmag, s=3, c="#1f3b73")
        for age, av, dm, colour, label in [(r.age, r.av, r.dm, "#c8553d", "ours"),
                                           (round(r.logAge50 * 20) / 20, r.AV50, r.MOD50, "#9a9a9a", "Hunt & Reffert")]:
            iso = table_iso[(table_iso.logAge.round(2) == round(age, 2)) & (table_iso.label <= grid.settings.phases)]
            a.plot(iso.G_BPmag - iso.G_RPmag + (A_BP - A_RP) * av, iso.Gmag + dm + A_G * av, c=colour, lw=1,
                   label=f"{label}: {age:.2f}, A_V {av:.2f}, dm {dm:.2f}")
        a.set(xlim=(s["BP-RP"].min() - 0.3, s["BP-RP"].max() + 0.3), ylim=(s.Gmag.max() + 0.5, s.Gmag.min() - 1),
              title=f"{n}: {r.n} stars", xlabel="BP − RP", ylabel="G")
        a.legend(fontsize=7, loc="lower left")
    fig.tight_layout()
    fig.savefig(figures / f"{name}_ages_calibration.png", dpi=110)
