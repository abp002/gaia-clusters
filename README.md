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

## A whole region: Auriga, 10° × 10° of the Galactic plane

Galactic longitude 170–180°, latitude −5° to +5°, towards the Galactic anticentre:
576,442 Gaia DR3 stars with RUWE < 1.4 and parallax/error > 5, and 121 Hunt & Reffert clusters
centred in the box, from 350 pc to 5.4 kpc. 104 of them keep at least 10 members after our cuts;
those are the ones scored. Region, cuts, parameters and matching rule were fixed in
[docs/protocol-auriga.md](docs/protocol-auriga.md) and committed before the first run (`e2a63a1`).

A cluster counts as recovered when one of our groups holds more than half of its members and more
than half of that group are its members. HDBSCAN runs four times, with min_cluster_size 10, 20, 40
and 80 as Hunt & Reffert did, and a cluster recovered by any run counts.

**Under the protocol: 3 of 104 clusters recovered (3 %).** HDBSCAN's default cluster selection,
excess of mass, put most of the field into a single group in every run (545,214 of the 576,442
stars with min_cluster_size 10), and the clusters stayed inside it. The protocol missed that Hunt &
Reffert's first paper ([2021](https://doi.org/10.1051/0004-6361/202039341)), whose parameter
choices the catalogue paper reuses, found leaf selection "almost always superior" to excess of mass
and much better at recovering small clusters, which hold very few stars next to their field.

**Post hoc, switching to leaf selection and changing nothing else: 84 of 104 recovered (81 %)**,
and 62 of the 70 that Hunt & Reffert rate highly reliable (89 %). This number was obtained after
seeing the first one, so it needs a second region with leaf selection fixed in advance before it
counts as a result.

| | Protocol (excess of mass) | Post hoc (leaf) |
|---|---|---|
| Under 1 kpc | 3 of 18 | 17 of 18 |
| 1–2 kpc | 0 of 46 | 37 of 46 |
| 2–3 kpc | 0 of 22 | 18 of 22 |
| 3 kpc and beyond | 0 of 18 | 12 of 18 |
| 10–29 members in our sample | 1 of 45 | 27 of 45 |
| 30–99 members | 1 of 49 | 47 of 49 |
| 100 members or more | 1 of 10 | 10 of 10 |

![Auriga: recovered and missed clusters on the sky, and recovery by distance and by members](figures/auriga.png)

Both expectations written in the protocol hold in the leaf runs: recovery falls with distance
(94, 80, 82 and 67 %, not strictly), and highly reliable clusters come back more often than the rest
(89 % against 65 %). Beyond 3 kpc the cuts weigh as much as the clustering: 13 of the 31 clusters
there keep fewer than 10 members once parallax/error > 5 is applied.

Of the 20 misses, 11 sit inside a bigger group of ours (more than half of their members are there,
but they are less than half of the group), 7 are not found at all (12 members in our sample, median)
and 2 are only partly found. The leaf run with min_cluster_size 10 also finds 494 groups that no
Hunt & Reffert cluster owns. They are not discoveries: Hunt & Reffert discard groups like these with
a significance test and a colour-magnitude classifier, which this project does not have yet.

![Colour-magnitude diagrams of M36, M37, M38 and Berkeley 71 as recovered by the leaf runs](figures/auriga_cmd.png)

Full tables in [results/auriga.md](results/auriga.md), one row per cluster in `results/auriga_scores_*.csv`.

## Run

```bash
uv sync
uv run pytest
uv run scripts/pleiades.py              # downloads into data/ on first run (~15 min from VizieR)
uv run scripts/fetch_region.py auriga   # stars from GAVO, Hunt & Reffert members from VizieR
uv run scripts/score_region.py auriga   # HDBSCAN runs, scores, figures (same for carina)
```
