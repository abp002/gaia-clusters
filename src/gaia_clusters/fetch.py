"""Download Gaia DR3 stars and Hunt & Reffert (2023) members from VizieR.

ESA's Gaia TAP returned 500/408 on these cone queries in October 2026 (archive release 3.15),
so both come from VizieR. A 6° cone takes about 15 minutes.
"""
import astropy.units as u
import pandas as pd
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier

GAIA_DR3 = "I/355/gaiadr3"
HR23_MEMBERS = "J/A+A/673/A114/members"
COLUMNS = ["Source", "RA_ICRS", "DE_ICRS", "Plx", "e_Plx", "pmRA", "pmDE", "Gmag", "BP-RP", "RUWE"]

# Quality cuts: parallax > 1 mas (~1 kpc), RUWE < 1.4, parallax_over_error > 5
MIN_PARALLAX, MAX_RUWE, MIN_PARALLAX_SNR = 1, 1.4, 5


def cone(ra: float, dec: float, radius_deg: float) -> pd.DataFrame:
    v = Vizier(columns=COLUMNS, row_limit=-1, timeout=1200,
               column_filters={"Plx": f">{MIN_PARALLAX}", "RUWE": f"<{MAX_RUWE}"})
    t = v.query_region(SkyCoord(ra, dec, unit="deg"), radius=radius_deg * u.deg, catalog=GAIA_DR3)[0]
    t = t.to_pandas()
    return t[t.Plx / t.e_Plx > MIN_PARALLAX_SNR].reset_index(drop=True)


def hr23_members(cluster: str) -> pd.DataFrame:
    """Members of one cluster by its Hunt & Reffert name, e.g. 'Melotte_22' for the Pleiades."""
    v = Vizier(row_limit=-1, timeout=300)
    return v.query_constraints(catalog=HR23_MEMBERS, Name=cluster)[0].to_pandas()
