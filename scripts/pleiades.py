"""Pleiades check: does HDBSCAN find the cluster on its own, and how close is it to Hunt & Reffert?

uv run scripts/pleiades.py [min_cluster_size] [min_samples]
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from gaia_clusters import fetch
from gaia_clusters.cluster import features, hdbscan_runs, tangent_plane
from gaia_clusters.score import score

RA, DEC, RADIUS = 56.75, 24.12, 6.0
DATA = Path("data")
CONE, MEMBERS = DATA / "pleiades_cone.parquet", DATA / "pleiades_hr23.parquet"

min_cluster_size = int(sys.argv[1]) if len(sys.argv) > 1 else 40
min_samples = int(sys.argv[2]) if len(sys.argv) > 2 else 10

DATA.mkdir(exist_ok=True)
if not CONE.exists():
    fetch.cone(RA, DEC, RADIUS).to_parquet(CONE)
if not MEMBERS.exists():
    fetch.hr23_members("Melotte_22").to_parquet(MEMBERS)
stars, ref = pd.read_parquet(CONE), pd.read_parquet(MEMBERS)

labels = hdbscan_runs(features(stars, RA, DEC), [min_cluster_size], min_samples)[min_cluster_size]
m = score(stars.Source, labels, ref, ["Melotte_22"]).iloc[0]
cl = stars[labels == m.label]

print(f"stars: {len(stars)}  clusters: {labels.max() + 1}  noise: {(labels < 0).mean():.1%}")
print(f"HR23 members: {len(ref)}, passing our cuts inside the cone: {m.n}")
print(f"Pleiades cluster: {m['size']} stars, {m.hits} in HR23")
print(f"  recovery {m.recovery:.1%}  purity {m.purity:.1%}")
print(f"  centre: plx {cl.Plx.median():.2f} mas, pm ({cl.pmRA.median():.1f}, {cl.pmDE.median():.1f}) mas/yr")

x, y = tangent_plane(stars.RA_ICRS, stars.DE_ICRS, RA, DEC)
inside = labels == m.label
in_ref = cl.Source.isin(ref.GaiaDR3)
fig, ax = plt.subplots(1, 3, figsize=(16, 5.2))
ax[0].scatter(x[~inside], y[~inside], s=0.2, c="0.8", rasterized=True)
ax[0].scatter(x[inside], y[inside], s=2, c="C3")
ax[0].set(xlabel="Δα (°)", ylabel="Δδ (°)", title="Sky")
ax[0].invert_xaxis()
ax[1].scatter(stars.pmRA[~inside], stars.pmDE[~inside], s=0.2, c="0.8", rasterized=True)
ax[1].scatter(cl.pmRA, cl.pmDE, s=2, c="C3")
ax[1].set(xlim=(-60, 60), ylim=(-80, 40), xlabel="μα* (mas/yr)", ylabel="μδ (mas/yr)", title="Proper motion")
abs_g = cl.Gmag + 5 * np.log10(cl.Plx / 100)
ax[2].scatter(cl["BP-RP"], abs_g, s=3, c=np.where(in_ref, "C3", "C0"))
ax[2].invert_yaxis()
ax[2].set(xlabel="BP−RP", ylabel="M_G", title="Colour-magnitude (red: in HR23, blue: ours only)")
fig.suptitle(f"Pleiades · HDBSCAN min_cluster_size={min_cluster_size}, min_samples={min_samples} · "
             f"recovers {m.recovery:.0%} of HR23, purity {m.purity:.0%}")
fig.tight_layout()
Path("figures").mkdir(exist_ok=True)
fig.savefig("figures/pleiades.png", dpi=110)
