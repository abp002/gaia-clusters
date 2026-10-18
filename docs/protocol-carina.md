# Protocol: the Carina region

Written and committed before the region was clustered. It repeats the Auriga run
([protocol-auriga.md](protocol-auriga.md)) with the one change Auriga's result called for, leaf
cluster selection, now fixed in advance. Anything done after the first run that departs from it is
reported as post hoc, next to the result it would replace, never instead of it.

## Region

Galactic longitude 285° to 295°, latitude −5° to +5°: Carina, along the Carina arm.

Chosen from the Hunt & Reffert (2023) catalogue alone, before any clustering, by a rule: the
window with the most clusters among the 10° windows, every 5° along the plane, that the Auriga
protocol counted (61 to 224). This is the one with 224: 221 open clusters and 3 moving groups, 135
of them highly reliable, from 150 pc to 12 kpc (and one moving group with a negative parallax).
21 are under 1 kpc, 38 at 1–2 kpc, 82 at 2–3 kpc and 83 beyond, against Auriga's 20, 46, 24 and 31.
Among them are NGC 3532, IC 2602 (the Southern Pleiades), NGC 3766, NGC 3293, NGC 3603, and
Trumpler 14 and 16 in the Carina Nebula.

The field is denser than Auriga's, so the method is tested where it is harder: 2,521,944 stars pass
the cuts below, 4.4 times Auriga's 576,442. They were downloaded before this protocol was
committed, to quote that number and the prediction below; nothing had been clustered.

## Stars

As in Auriga: Gaia DR3 from the GAVO Data Center (`gaia.dr3lite`), every source in the box with
RUWE < 1.4 and parallax_over_error > 5, no parallax floor and no magnitude cut. The same deviation
from Hunt & Reffert: no Rybizki et al. (2022) fidelity cut.

## Clustering

As in Auriga, except for the cluster selection method:

- Five axes: position on the gnomonic tangent plane centred on (l, b) = (290°, 0°), parallax, pmRA
  and pmDE, each rescaled to the same median and interquartile range (`RobustScaler`).
- HDBSCAN (`hdbscan` 0.8.44) with **leaf** cluster selection, min_samples = 10, four runs with
  min_cluster_size = 10, 20, 40 and 80.
- Changed from Auriga: leaf instead of excess of mass. Hunt & Reffert's first paper (2021, A&A 646,
  A104) found leaf selection "almost always superior". In Auriga, excess of mass recovered 3 % of
  the clusters, and leaf, run after seeing that, 81 %.
- The same deviation: one field for the whole region, where Hunt & Reffert clustered fields of
  about 5.5° for clusters beyond 750 pc.
- If HDBSCAN cannot run on the whole region on this machine, the run stops there, and any split of
  the region is reported as a deviation.

## Reference and matching

As in Auriga. The reference is the 224 clusters centred in the box, with every member the catalogue
lists for them that is in our sample. A cluster is recoverable when at least 10 of its members are.
A run recovers a cluster when one of its groups holds more than half of the cluster's members in our
sample and more than half of that group's stars are its members.

**Primary result:** the share of recoverable clusters recovered by at least one of the four leaf runs.

## Prediction: when the replication holds

Auriga's leaf runs recovered 27 of the 45 clusters with 10–29 members in our sample, 47 of 49 with
30–99 and 10 of 10 with 100 or more. Carina has 202 recoverable clusters: 78, 94 and 30 in those
bins. If every bin keeps Auriga's rate, Carina should recover **81.5 %** of them, with a 95 % range
of 73.3 % to 89.1 % (148 to 180 clusters). The numbers come from
`uv run scripts/score_region.py carina --predict`: each Auriga rate is drawn from a Jeffreys Beta
distribution, then a binomial count over Carina's clusters (`src/gaia_clusters/predict.py`).

- **The Auriga result replicates if Carina recovers at least 148 of its 202 recoverable clusters
  (73.3 %).**
- Below that, it does not: something in Carina, most likely the denser field, costs more than its
  mix of cluster sizes explains, and the 81 % cannot be quoted beyond thin fields like Auriga's.

The README reports the result next to Auriga's, whatever it turns out to be.

## Secondary

Reported in full whatever they show, as in Auriga:
- the share for each leaf run on its own;
- the same shares restricted to the highly reliable clusters (CST ≥ 5 and median CMD class ≥ 0.5);
- broken down by distance (under 1 kpc, 1–2, 2–3, 3 kpc and beyond) and by members in our sample
  (10–29, 30–99, 100 and more);
- per run, the groups that no Hunt & Reffert cluster owns, which are not discoveries;
- the same four runs with excess of mass, to see whether Auriga's failure repeats.

## Expectations

Written down so the result can contradict them. They are not part of the test above.
- Excess of mass fails again: most of the field ends up in one group, and under 10 % of the
  clusters are recovered.
- With leaf, recovery falls with distance, and highly reliable clusters come back more often than
  the rest.
- In each members bin, recovery is lower than in Auriga, because the field is denser.
