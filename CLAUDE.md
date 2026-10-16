# gaia-clusters

Recover known open clusters from Gaia DR3 with HDBSCAN and score them against Hunt & Reffert
(2023). Written in English. The why and the decisions live in the vault note (Proyectos/gaia-clusters).

## QA
Nivel: activo

Issues: abp002/gaia-clusters

## How it works
- `fetch.py`: cones from VizieR (`I/355/gaiadr3`, ~15 min for 6°), regions from the GAVO Data Center
  (`gaia.dr3lite`, async TAP polled by hand: GAVO drops pyvo's long polls), Hunt & Reffert tables from
  VizieR (`J/A+A/673/A114`). ESA's Gaia TAP returned 500/408 on cone queries in Oct 2026.
- `region.py`: `GalacticBox`; the download is a 4-corner ICRS polygon, then filtered exactly in (l, b).
- Quality cuts live in `fetch.py`: parallax > 1 mas, RUWE < 1.4, parallax/error > 5.
- `cluster.py`: gnomonic tangent plane (degrees) + parallax + proper motion, `RobustScaler` (as
  Hunt & Reffert), HDBSCAN from the `hdbscan` library (Boruvka, scales to millions; scikit-learn's is
  quadratic). `hdbscan_runs` caches the tree so the four min_cluster_size runs share it.
- `score.py`: the matching rule of `docs/protocol-auriga.md` (recovery > 0.5 and purity > 0.5, over
  reference members that passed our cuts; recoverable = at least 10 of them). The best match is the
  cluster with most reference members, not the biggest: the field can come out as one huge cluster.
- Cluster selection: use `method="leaf"`. Excess of mass (the library default, used by the Auriga
  protocol) puts the field into one giant group and hides the clusters in it: 3 % recovered against
  81 % with leaf, post hoc. Hunt & Reffert's Paper I (A&A 646, A104) picked leaf for this reason.
- Protocols in `docs/` are fixed before a run. Never edit one after seeing results: report changes as post hoc.
- Never call it discovery: the reference used the same method on the same data.
- `data/` is gitignored and rebuilt by the scripts.

## Commands
- `uv run pytest` — projection and scoring tests.
- `uv run scripts/pleiades.py [min_cluster_size] [min_samples]` — Pleiades check, writes `figures/pleiades.png`.
- `uv run scripts/fetch_region.py` then `uv run scripts/auriga.py` — Auriga region; labels cached in
  `data/auriga_labels.parquet` (delete it to recluster), scores in `results/`.
