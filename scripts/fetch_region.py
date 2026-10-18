"""Download a region: Gaia DR3 stars and their Bailer-Jones distances from GAVO, and Hunt & Reffert
members from VizieR.

uv run scripts/fetch_region.py auriga|carina
"""
import sys
import time
from pathlib import Path

from gaia_clusters import fetch
from gaia_clusters.region import REGIONS

name = sys.argv[1]
box = REGIONS[name]
DATA = Path("data")
DATA.mkdir(exist_ok=True)

for file, get in [("hr23_clusters", fetch.hr23_clusters),
                  (f"{name}_hr23_members", lambda: fetch.hr23_members_in(box)),
                  (f"{name}_stars", lambda: fetch.box(box)),
                  (f"{name}_distances", lambda: fetch.distances(box))]:
    path = DATA / f"{file}.parquet"
    if path.exists():
        continue
    t0 = time.time()
    t = get()
    t.to_parquet(path)
    print(f"{file}: {len(t)} rows in {time.time() - t0:.0f} s", flush=True)
