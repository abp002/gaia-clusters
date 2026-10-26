"""Download the Hyades sample and the published tail members, as docs/protocol-hyades.md fixes them.

uv run scripts/fetch_hyades.py    # data/hyades_*.parquet; the sample takes the longest
"""
import time
from pathlib import Path

import pandas as pd
from astroquery.vizier import Vizier

from gaia_clusters import fetch
from gaia_clusters.hyades import HYADES, adql_residual_below, space_velocity

DATA = Path("data")
V = space_velocity(*(HYADES[k] for k in ("ra", "dec", "plx", "pmra", "pmdec", "rv")))
catalogue = lambda name: lambda: Vizier(row_limit=-1, columns=["**"], timeout=600).get_catalogs(name)[0].to_pandas()

for file, get in [("hyades_jerabkova", catalogue("J/A+A/647/A137/tableb2")),
                  ("hyades_roser", catalogue("J/A+A/621/L2/stars")),
                  ("hyades_gcns", catalogue("J/A+A/649/A6/hyacomb")),
                  ("hyades_sample", lambda: fetch.hyades_sample(adql_residual_below(V))),
                  ("hyades_references", lambda: fetch.by_source_id(  # every reference star, whatever its velocity
                      [*pd.read_parquet(DATA / "hyades_jerabkova.parquet").GaiaEDR3,
                       *pd.read_parquet(DATA / "hyades_roser.parquet").Source,  # DR2 IDs: most are the same in DR3
                       *pd.read_parquet(DATA / "hyades_gcns.parquet").GaiaEDR3]))]:
    path = DATA / f"{file}.parquet"
    if path.exists():
        continue
    t0 = time.time()
    t = get()
    t.to_parquet(path)
    print(f"{file}: {len(t)} rows in {time.time() - t0:.0f} s", flush=True)
