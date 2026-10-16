"""Download the Auriga region (l 170-180°, |b| < 5°): Gaia DR3 stars and Hunt & Reffert members.

uv run scripts/fetch_region.py
"""
import time
from pathlib import Path

from gaia_clusters import fetch
from gaia_clusters.region import AURIGA

DATA = Path("data")
DATA.mkdir(exist_ok=True)

for name, get in [("hr23_clusters", fetch.hr23_clusters),
                  ("auriga_hr23_members", lambda: fetch.hr23_members_in(AURIGA)),
                  ("auriga_stars", lambda: fetch.box(AURIGA))]:
    path = DATA / f"{name}.parquet"
    if path.exists():
        continue
    t0 = time.time()
    t = get()
    t.to_parquet(path)
    print(f"{name}: {len(t)} rows in {time.time() - t0:.0f} s", flush=True)
