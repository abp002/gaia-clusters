# gaia-clusters

Recover known open star clusters from Gaia DR3 with HDBSCAN, and measure how many come back.

Stars born in the same cloud share position, distance and motion. On the sky they blend into the
field; in five dimensions (position, parallax, proper motion) they form dense clumps. HDBSCAN
finds those clumps without being told how many there are and leaves the field as noise.

This is not a discovery project. [Hunt & Reffert (2023)](https://doi.org/10.1051/0004-6361/202346285)
already ran HDBSCAN over Gaia DR3 and catalogued about 7,000 clusters. Their catalogue is used here
as the reference: the score is the share of their members we recover, and the share of what we
find that they also list.

## First check: the Pleiades

A 6° cone around the Pleiades with parallax > 1 mas, RUWE < 1.4 and parallax/error > 5
(118,080 stars), clustered with `min_cluster_size=40, min_samples=10`, untuned:

| | |
|---|---|
| Clusters found | 2: the field and the Pleiades |
| Pleiades cluster | 1,497 stars |
| Recovery (Hunt & Reffert members passing our cuts that we find) | 96.9 % (1,424 of 1,470) |
| Purity (our members that Hunt & Reffert also list) | 95.1 % |

These numbers come from the `hdbscan` library, which the project uses because it scales to millions
of stars. The first run used scikit-learn's HDBSCAN and gave 96.3 % and 96.3 %.

![Pleiades](figures/pleiades.png)

The Pleiades are the easy case: close, rich and far from the field in proper motion. The number
that matters is the all-sky one, with distant and sparse clusters. Also, this is the same method
on the same data as the reference, so it measures agreement with it, not truth.

## Next: a 10° × 10° region

The Pleiades are one cluster. The next measurement is a whole region of the Galactic plane, Auriga
(l 170–180°, |b| < 5°), with 121 Hunt & Reffert clusters from 350 pc to 5.4 kpc. Its region,
cuts, parameters and matching rule are in [docs/protocol-auriga.md](docs/protocol-auriga.md),
committed before the region was clustered.

## Run

```bash
uv sync
uv run pytest
uv run scripts/pleiades.py       # downloads into data/ on first run (~15 min from VizieR)
uv run scripts/fetch_region.py   # Auriga stars from GAVO and Hunt & Reffert members from VizieR
uv run scripts/auriga.py         # four HDBSCAN runs, scores, figures
```
