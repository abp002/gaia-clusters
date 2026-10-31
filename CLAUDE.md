# gaia-clusters

Recover known open clusters from Gaia DR3 with HDBSCAN and score them against Hunt & Reffert
(2023). Written in English.

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
- `ages.py` + `scripts/fit_ages.py`, `scripts/compare_ages.py`: PARSEC isochrones (`fetch.isochrones`, CMD 3.9 form,
  `data/parsec_solar.parquet`) fitted per recovered group, IMF-weighted mixture likelihood, parallax prior on the
  distance modulus, bootstrap interval (`docs/protocol-ages.md`; settings frozen in `Settings`). The protocol's bar
  failed (56 % within 0.3 dex): distant groups truncated by the parallax cut fit as ~10 Myr. `--post-hoc`
  (`faint_quantile=0.9`) is reported beside it, never instead. The CMD service needs `truststore` (no intermediate cert).
- `hyades.py` + `scripts/fetch_hyades.py`, `scripts/hyades.py`: Hyades tails (`docs/protocol-hyades.md`). The download
  applies the convergent-point condition server-side in ADQL (`adql_residual_below`, tested against `residual`): 914k
  stars within 500 pc and 10 km/s. Reference stars are fetched by ID with no cuts (`fetch.by_source_id`). Result: eom
  swallows the volume; Jerabkova's tail stars drift 3–9 km/s off the cluster's velocity beyond 100 pc.
- `viewer.py` + `scripts/export_viewer.py`: the 3D viewer's data. Field stars at their Bailer-Jones
  geometric distance (`data/<region>_distances.parquet`, from GAVO `gedr3dist.litewithdist` with the
  star query's polygon and cuts); members of a recovered cluster at the cluster's `dist50`, because a
  member's own distance is ~25 times more uncertain than a cluster is big (decision in the README).
  The field is shuffled once (seed 0) and cut into nested levels of 65,536 × 2^k stars, quantized to
  uint16 between `quant.lo` and `quant.hi`; members are float32 (quantization steps show up close).
  Proper motions (galactic, int16 steps of 0.01 mas/yr) go in `field-<k>-pm.bin` and `members-pm.bin`,
  fetched only when the time slider moves; `meta.reference` is the field's median motion per 250 pc,
  which the page takes off by default. `drift()` in the shaders and `drifted()` in JS must stay twins.
  Each card's proper-motion diagram reads `cmd-pm.bin` (row for row with `cmd.bin`) and `near-pm.bin`
  (field stars around each cluster, `near()` in viewer.py; offset and count in `clusters[].near`).
- `scripts/sweep_region.py`: leaf runs from min_cluster_size 10 to 240 (`data/<region>_labels_sweep.parquet`,
  summary in `results/<region>_sweep.csv`), exploration outside the score. The export turns each into a
  single-run view (`view(run=k)` in viewer.py: `run-<k>.json` and `run-<k>-*.bin`); the page's slider swaps
  views (`makeView`/`setView` in main.js). `meta.json` stays the score's view: each cluster from its best run.
- `scripts/render.py`: Blender (5.2, Cycles on Metal) still from `viewer/data/<region>/`, run headless. Points through
  geometry nodes (Mesh to Points); an unlit transparent-emission material composites the faint field stars into a
  density on paper, as the page does. Default camera: perspective along the wedge at the recovered clusters'
  median distance. Views tried: from above (a 10° wedge is a thin sliver) and from behind the Sun (near stars blow up).
- `scripts/film.py` + `scripts/film_encode.py`: the one-minute film (night style) from `viewer/data/<region>/`. Stars drift
  through a `vel` attribute times a keyed clock in geometry nodes; shots are keyframes in `SHOTS`, titles (ffmpeg drawtext,
  Charter: New York's default instance has no hyphen) follow the same seconds. Members get a modelled depth, as deep as
  their cluster is wide on the sky; point sizes shrink for the close-up. Frames in `data/film/<region>/`, ~6 s each on an M4.
- `viewer/`: static page, three.js 0.170 from jsdelivr, no build step. Clusters are spheres (instanced,
  tidal radius `rtpc`, at least 4 CSS px, resized every frame in `drawSpheres`). The field is rendered as a
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
- `uv run scripts/sweep_region.py <region>` — the slider's runs (under a minute; reuses the cached tree), before the export.
- `blender -b -P scripts/render.py -- <region>` — the README's still, `figures/render-<region>.jpg` (after the export).
- `blender -b -P scripts/film.py -- <region> [--preview]` then `uv run scripts/film_encode.py <region> [--preview]` —
  the film, `figures/film-<region>.mp4` (~3 h for 1080p; the preview, 640x360 every third frame, ~25 min).
- `uv run scripts/export_viewer.py <region>` — the viewer's data into `viewer/data/<region>/`;
  `python3 -m http.server -d viewer` to look at it (any free port: another session may hold 8000).
- `uv run scripts/fetch_region.py <region>` then `uv run scripts/score_region.py <region>` — one region
  (`auriga`, `carina`); labels cached in `data/<region>_labels_<method>.parquet` (delete to recluster),
  scores in `results/`, figures in `figures/`. `--predict` prints the prediction alone, before clustering.
