# Protocol: the Auriga region

Written and committed before the region was clustered. Anything done after the first run that
departs from it is reported as post hoc, next to the result it would replace, never instead of it.

## Region

Galactic longitude 170° to 180°, latitude −5° to +5°: the Galactic anticentre, in Auriga.

Chosen from the Hunt & Reffert (2023) catalogue alone, before any clustering. Every 10° window
along the plane holds between 61 and 224 of their clusters. This one looks outwards, away from the
Galactic centre, where the field is thinnest, and still has 121 clusters (119 open clusters, 2 moving
groups) from 350 pc to 5.4 kpc: 20 under 1 kpc, 46 at 1–2 kpc, 24 at 2–3 kpc and 31 beyond. Among
them are M36, M37 and M38 (NGC 1960, NGC 2099 and NGC 1912).

## Stars

Gaia DR3 from the GAVO Data Center (`gaia.dr3lite`), every source in the box with RUWE < 1.4 and
parallax_over_error > 5. No parallax floor (the Pleiades check had parallax > 1 mas, which would
drop every cluster beyond 1 kpc) and no magnitude cut.

Deviation from Hunt & Reffert: they kept sources with a Rybizki et al. (2022) fidelity_v1 of at
least 0.5, with no RUWE or parallax S/N cut, which keeps many more faint stars.

## Clustering

- Five axes, as Hunt & Reffert: position on the gnomonic tangent plane centred on (l, b) = (175°, 0°),
  parallax, pmRA and pmDE, each rescaled to the same median and interquartile range (`RobustScaler`).
- HDBSCAN (`hdbscan` 0.8.44, default excess-of-mass selection), min_samples = 10, four runs with
  min_cluster_size = 10, 20, 40 and 80: Hunt & Reffert's parameter grid. Their paper does not state
  the cluster selection method; we keep the library default.
- Deviation: one field for the whole region. Hunt & Reffert clustered fields of about 5.5° for
  clusters beyond 750 pc, of about 44° with parallax > 1 mas for nearer ones, and one all-sky field
  for the nearest 250 pc.

## Reference

The 121 Hunt & Reffert clusters whose centre lies in the box. Their members are every star the
catalogue lists for them (any membership probability, inside or outside the tidal radius) that is
in our sample, so a star removed by our cuts is never counted against us.

A cluster is **recoverable** when at least 10 of its members are in our sample: with fewer, no run
can find it, whatever the method.

## Matching

A run **recovers** a reference cluster when one of its clusters holds more than half of the
reference cluster's members in our sample (recovery > 0.5) and more than half of that cluster's
stars are members of the reference cluster (purity > 0.5). Both are strict, so the match is
one-to-one: a cluster split in two halves, or two clusters merged into one, does not count.
Noise never matches.

**Primary result:** the share of recoverable clusters recovered by at least one of the four runs.
Hunt & Reffert built their catalogue by merging the same four runs.

Secondary, reported in full whatever they show:
- the same share for each run on its own;
- the same shares restricted to their 73 "highly reliable" clusters in the box;
- broken down by distance (their `dist50`: under 1 kpc, 1–2, 2–3, 3 kpc and beyond) and by members
  in our sample (10–29, 30–99, 100 and more). "Highly reliable" is their own cut: CST ≥ 5 and median
  CMD class ≥ 0.5, 4,105 clusters over the whole sky;
- per run, the number of our clusters in which no Hunt & Reffert cluster, wherever its centre, holds
  more than half the stars. These are not presented as discoveries: they used the same method on the
  same data and rejected what failed their significance and CMD tests.

## Expectations

Written down so the result can contradict them: recovery falls with distance, because distant
clusters keep only their brightest members after the parallax S/N cut and sit where the field is
densest; and their "highly reliable" clusters are recovered more often than the rest.
