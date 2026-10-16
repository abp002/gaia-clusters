# gaia-clusters

Recover known open clusters from Gaia DR3 with HDBSCAN and score them against Hunt & Reffert
(2023). Written in English. The why and the decisions live in the vault note (Proyectos/gaia-clusters).

## QA
Nivel: activo

Issues: abp002/gaia-clusters

## How it works
- `fetch.py` downloads from VizieR (`I/355/gaiadr3`, `J/A+A/673/A114/members`). ESA's Gaia TAP
  returned 500/408 on cone queries in Oct 2026; VizieR works but a 6° cone takes ~15 min.
- Quality cuts live in `fetch.py`: parallax > 1 mas, RUWE < 1.4, parallax/error > 5.
- `cluster.py`: gnomonic tangent plane (degrees) + parallax + proper motion, `RobustScaler`, HDBSCAN.
  `best_match` picks the cluster with the most reference members (not the biggest: the field can
  come out as one huge cluster). Recovery is over reference members that passed our cuts.
- Never call it discovery: the reference used the same method on the same data.
- `data/` is gitignored and rebuilt by the scripts.

## Commands
- `uv run pytest` — projection and scoring tests.
- `uv run scripts/pleiades.py [min_cluster_size] [min_samples]` — Pleiades check, writes `figures/pleiades.png`.
