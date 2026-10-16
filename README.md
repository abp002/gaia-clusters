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
| Pleiades cluster | 1,471 stars |
| Recovery (Hunt & Reffert members passing our cuts that we find) | 96.3 % (1,416 of 1,470) |
| Purity (our members that Hunt & Reffert also list) | 96.3 % |

![Pleiades](figures/pleiades.png)

The Pleiades are the easy case: close, rich and far from the field in proper motion. The number
that matters is the all-sky one, with distant and sparse clusters. Also, this is the same method
on the same data as the reference, so it measures agreement with it, not truth.

## Run

```bash
uv sync
uv run pytest
uv run scripts/pleiades.py   # downloads into data/ on first run (~15 min from VizieR)
```
