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
- `region.py`: `GalacticBox` and `REGIONS` (auriga, carina); the download is a 4-corner ICRS polygon,
  then filtered exactly in (l, b).
- Quality cuts live in `fetch.py`: parallax > 1 mas, RUWE < 1.4, parallax/error > 5.
- `cluster.py`: gnomonic tangent plane (degrees) + parallax + proper motion, `RobustScaler` (as
  Hunt & Reffert), HDBSCAN from the `hdbscan` library (Boruvka, scales to millions; scikit-learn's is
  quadratic). `hdbscan_runs` caches the tree so the four min_cluster_size runs share it.
- `score.py`: the matching rule of `docs/protocol-auriga.md` (recovery > 0.5 and purity > 0.5, over
  reference members that passed our cuts; recoverable = at least 10 of them). The best match is the
  cluster with most reference members, not the biggest: the field can come out as one huge cluster.
- `predict.py`: what one region's recovery rates by members bin predict for another (Jeffreys Beta
  per bin, then binomial). It is the replication test of `docs/protocol-carina.md`.
- `viewer.py` + `scripts/export_viewer.py`: the 3D viewer's data. Field stars at their Bailer-Jones
  geometric distance (`data/<region>_distances.parquet`, from GAVO `gedr3dist.litewithdist` with the
  star query's polygon and cuts); members of a recovered cluster at the cluster's `dist50`, because a
  member's own distance is ~25 times more uncertain than a cluster is big (decision in the README).
  The field is shuffled once (seed 0) and cut into nested levels of 65,536 × 2^k stars, quantized to
  uint16 between `quant.lo` and `quant.hi`; members are float32 (quantization steps show up close).
- `viewer/`: static page, three.js 0.170 from jsdelivr, no build step. The field is rendered as a
  density (stars per pixel counted in a half-float buffer, shown in grey on a log scale). `viewer/data/`
  is committed: Pages cannot rebuild it, since it needs `data/`. `.github/workflows/pages.yml` deploys
  `viewer/` and skips while the repo is private.
- Cluster selection: use `method="leaf"`. Excess of mass (the library default, used by the Auriga
  protocol) puts the field into one giant group and hides the clusters in it: 3 % recovered against
  81 % with leaf, post hoc. Carina, with leaf fixed in advance: 74 % (bar 73.3 %), excess of mass 1 %.
  Hunt & Reffert's Paper I (A&A 646, A104) picked leaf for this reason.
- Protocols in `docs/` are fixed before a run. Never edit one after seeing results: report changes as post hoc.
- Never call it discovery: the reference used the same method on the same data.
- `data/` is gitignored and rebuilt by the scripts.

## Commands
- `uv run pytest` — projection and scoring tests.
- `uv run scripts/pleiades.py [min_cluster_size] [min_samples]` — Pleiades check, writes `figures/pleiades.png`.
- `uv run scripts/export_viewer.py <region>` — the viewer's data into `viewer/data/<region>/`;
  `python3 -m http.server -d viewer` to look at it (any free port: another session may hold 8000).
- `uv run scripts/fetch_region.py <region>` then `uv run scripts/score_region.py <region>` — one region
  (`auriga`, `carina`); labels cached in `data/<region>_labels_<method>.parquet` (delete to recluster),
  scores in `results/`, figures in `figures/`. `--predict` prints the prediction alone, before clustering.
