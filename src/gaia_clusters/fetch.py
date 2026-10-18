"""Download Gaia DR3 stars and Hunt & Reffert (2023) members.

ESA's Gaia TAP returned 500/408 on cone queries in October 2026 (archive release 3.15). Small
regions come from VizieR (a 6° cone takes about 15 minutes); big ones from the GAVO Data Center
in Heidelberg (`gaia.dr3lite`), through an async TAP job.
"""
import time

import astropy.units as u
import pandas as pd
import pyvo
import requests
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier

from gaia_clusters.region import GalacticBox, galactic

GAIA_DR3 = "I/355/gaiadr3"
HR23_CLUSTERS = "J/A+A/673/A114/clusters"
HR23_MEMBERS = "J/A+A/673/A114/members"
GAVO_TAP = "https://dc.zah.uni-heidelberg.de/tap"
COLUMNS = ["Source", "RA_ICRS", "DE_ICRS", "Plx", "e_Plx", "pmRA", "pmDE", "Gmag", "BP-RP", "RUWE"]

# Quality cuts: parallax > 1 mas (~1 kpc), RUWE < 1.4, parallax_over_error > 5
MIN_PARALLAX, MAX_RUWE, MIN_PARALLAX_SNR = 1, 1.4, 5


def cone(ra: float, dec: float, radius_deg: float) -> pd.DataFrame:
    v = Vizier(columns=COLUMNS, row_limit=-1, timeout=1200,
               column_filters={"Plx": f">{MIN_PARALLAX}", "RUWE": f"<{MAX_RUWE}"})
    t = v.query_region(SkyCoord(ra, dec, unit="deg"), radius=radius_deg * u.deg, catalog=GAIA_DR3)[0]
    t = t.to_pandas()
    return t[t.Plx / t.e_Plx > MIN_PARALLAX_SNR].reset_index(drop=True)


def box(region: GalacticBox, min_parallax: float | None = None) -> pd.DataFrame:
    """Stars in a galactic box with RUWE < 1.4 and parallax_over_error > 5, in the same column
    names as `cone`. No parallax floor unless given: a region holds clusters at every distance."""
    floor = f"AND parallax > {min_parallax}" if min_parallax is not None else ""
    query = f"""
        SELECT source_id, ra, dec, parallax, parallax_error, pmra, pmdec,
               phot_g_mean_mag, phot_bp_mean_mag - phot_rp_mean_mag AS bp_rp, ruwe
        FROM gaia.dr3lite
        WHERE {_box_where(region)} {floor}"""
    t = _async(GAVO_TAP, query).to_pandas()
    t.columns = COLUMNS
    t["GLON"], t["GLAT"] = galactic(t.RA_ICRS, t.DE_ICRS)
    return t[region.contains(t.GLON, t.GLAT)].reset_index(drop=True)


def distances(region: GalacticBox) -> pd.DataFrame:
    """Bailer-Jones et al. (2021) geometric distances in pc (median, 16th and 84th percentiles) for
    the stars `box` returns: same polygon, same cuts. The exact (l, b) filter comes when merging on
    Source. Geometric, not photogeometric: every star with a parallax has one."""
    query = f"""
        SELECT source_id, r_med_geo, r_lo_geo, r_hi_geo
        FROM gedr3dist.litewithdist
        WHERE {_box_where(region)}"""
    t = _async(GAVO_TAP, query).to_pandas()
    t.columns = ["Source", "r_med_geo", "r_lo_geo", "r_hi_geo"]
    return t


def _box_where(region: GalacticBox) -> str:
    poly = ", ".join(f"{ra:.6f}, {dec:.6f}" for ra, dec in region.icrs_corners())
    return f"""1 = CONTAINS(POINT('ICRS', ra, dec), POLYGON('ICRS', {poly}))
          AND ruwe < {MAX_RUWE} AND parallax > {MIN_PARALLAX_SNR} * parallax_error"""


def _async(url: str, query: str, poll_s: int = 30):
    """Run an async TAP job, polling its phase. pyvo's own wait long-polls, and GAVO drops long
    polls after a while even though the job keeps running."""
    job = pyvo.dal.TAPService(url).submit_job(query, maxrec=20_000_000)
    try:
        job.execution_duration = 4 * 3600
    except pyvo.dal.DALServiceError:
        pass  # the service keeps its own limit
    job.run()
    print(f"TAP job {job.url}, time limit {job.execution_duration} s", flush=True)
    phase = "QUEUED"
    while phase not in ("COMPLETED", "ERROR", "ABORTED"):
        time.sleep(poll_s)
        try:
            phase = job.phase
        except (pyvo.dal.DALServiceError, requests.ConnectionError) as e:
            print(f"poll failed, retrying: {e}", flush=True)
    job.raise_if_error()
    return job.fetch_result().to_table()


def hr23_clusters() -> pd.DataFrame:
    return Vizier(row_limit=-1, timeout=300, columns=["**"]).get_catalogs(HR23_CLUSTERS)[0].to_pandas()


def hr23_members(cluster: str) -> pd.DataFrame:
    """Members of one cluster by its Hunt & Reffert name, e.g. 'Melotte_22' for the Pleiades."""
    v = Vizier(row_limit=-1, timeout=300)
    return v.query_constraints(catalog=HR23_MEMBERS, Name=cluster)[0].to_pandas()


def hr23_members_in(region: GalacticBox) -> pd.DataFrame:
    """Every Hunt & Reffert member star inside the box, whatever cluster it belongs to."""
    v = Vizier(row_limit=-1, timeout=600, columns=["Name", "GaiaDR3", "inrt", "Prob", "GLON", "GLAT"])
    t = v.query_constraints(catalog=HR23_MEMBERS, GLON=f"{region.l_min}..{region.l_max}",
                            GLAT=f"{region.b_min}..{region.b_max}")[0].to_pandas()
    return t[region.contains(t.GLON, t.GLAT)].reset_index(drop=True)
