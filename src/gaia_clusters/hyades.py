"""The Hyades' tidal tails by the convergent point, as fixed in docs/protocol-hyades.md: stars whose
tangential velocity matches the cluster's space velocity projected at their position."""
import numpy as np

K = 4.740470446  # km/s per (mas/yr) / mas: 1 au/yr

# Hunt & Reffert (2023), Melotte 25: position (deg), parallax (mas), proper motion (mas/yr), radial velocity (km/s)
HYADES = {"ra": 66.70864584, "dec": 16.0807978, "plx": 21.23204937, "pmra": 104.13554468, "pmdec": -28.73230917,
          "rv": 39.26349773, "dist": 47.19065976, "rt": 8.61953181}
DOWNLOAD_DV = 10.0  # km/s: the download keeps stars within this of the cluster's tangential velocity
DV, ERRORS = 3.0, 2.0  # the selection: Δ ≤ √(DV² + (ERRORS σ)²) km/s


def basis(ra, dec):
    """Unit vectors towards the star (r), and of increasing α (p) and δ (q), in ICRS Cartesian; degrees in."""
    a, d = np.radians(ra), np.radians(dec)
    r = np.stack([np.cos(d) * np.cos(a), np.cos(d) * np.sin(a), np.sin(d)], axis=-1)
    p = np.stack([-np.sin(a), np.cos(a), np.zeros_like(a)], axis=-1)
    q = np.stack([-np.sin(d) * np.cos(a), -np.sin(d) * np.sin(a), np.cos(d)], axis=-1)
    return r, p, q


def space_velocity(ra, dec, plx, pmra, pmdec, rv) -> np.ndarray:
    """ICRS Cartesian velocity in km/s from position, parallax, proper motion (μ_α*, μ_δ) and radial velocity."""
    r, p, q = basis(np.asarray(ra, float), np.asarray(dec, float))
    return rv * r + K * pmra / plx * p + K * pmdec / plx * q


def residual(ra, dec, plx, e_plx, pmra, pmdec, v) -> tuple[np.ndarray, np.ndarray]:
    """How far each star's tangential velocity is from what velocity `v` would give at its position, in km/s,
    and the error its parallax puts on its tangential velocity."""
    _, p, q = basis(np.asarray(ra, float), np.asarray(dec, float))
    va, vd = K * np.asarray(pmra) / plx, K * np.asarray(pmdec) / plx
    dv = np.hypot(va - p @ v, vd - q @ v)
    return dv, np.hypot(va, vd) * np.asarray(e_plx) / np.asarray(plx)


def selected(dv, sigma, dv_max: float = DV, errors: float = ERRORS) -> np.ndarray:
    return dv <= np.sqrt(dv_max ** 2 + (errors * sigma) ** 2)


def adql_residual_below(v, dv: float = DOWNLOAD_DV) -> str:
    """The ADQL condition of `residual` < dv, without the error term, for gaia.dr3lite's columns."""
    vx, vy, vz = (float(x) for x in v)
    a, d = "RADIANS(ra)", "RADIANS(dec)"
    pa = f"({-vx} * SIN({a}) + {vy} * COS({a}))"
    pd = f"({-vx} * SIN({d}) * COS({a}) + {-vy} * SIN({d}) * SIN({a}) + {vz} * COS({d}))"
    return (f"POWER({K} * pmra / parallax - {pa}, 2) + POWER({K} * pmdec / parallax - {pd}, 2) < {dv ** 2}")


def heliocentric_icrs(ra, dec, plx) -> np.ndarray:
    """ICRS Cartesian position in pc from position (deg) and parallax (mas)."""
    r, _, _ = basis(np.asarray(ra, float), np.asarray(dec, float))
    return r * (1000 / np.asarray(plx, float))[..., None]
