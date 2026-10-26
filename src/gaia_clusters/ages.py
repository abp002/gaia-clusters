"""A cluster's age from PARSEC isochrones fitted to its colour-magnitude diagram, as fixed in
docs/protocol-ages.md: for each age on the grid, extinction and distance modulus by maximum likelihood;
the age with the highest likelihood; the interval from bootstrap resamples of the stars."""
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.spatial import cKDTree

# A_λ / A_V in Gaia EDR3 passbands for a G2V star, Cardelli et al. (1989) + O'Donnell (1994), R_V = 3.1:
# what the CMD 3.9 service quotes with its isochrones
A_G, A_BP, A_RP = 0.83627, 1.08337, 0.63439
ZERO_POINT = 0.017  # mas added to Gaia DR3 parallaxes (Lindegren et al. 2021, median)
PLX_SYSTEMATIC = 0.01  # mas, added in quadrature to the median parallax's error
MAX_AV = 6.0


@dataclass(frozen=True)
class Settings:
    """What calibration fixes (docs/protocol-ages.md)."""
    sigma_colour: float = 0.05  # mag: width of the sequence in BP − RP
    sigma_g: float = 0.15       # mag: width in G
    outliers: float = 0.2       # share of stars from the uniform component
    phases: int = 7             # highest PARSEC label kept: up to the early AGB, no thermal pulses
    step: float = 0.2           # densify the isochrone to this many widths between points
    faintest: float = 18.0      # apparent G: fainter stars are left out, and the model is cut there
    brightest: float = 5.0      # apparent G: likewise at the bright end, where Gaia saturates
    neighbours: int = 32        # isochrone points summed for each star


class Grid:
    """The isochrones by age, each as dense points in (BP − RP, M_G) with the share of stars the IMF puts
    at each (from PARSEC's int_IMF), and a KD-tree of the points scaled by the widths."""

    def __init__(self, table: pd.DataFrame, settings: Settings = Settings()):
        self.settings = s = settings
        self.ages, self.trees, self.points, self.weights = [], [], [], []
        for age, iso in table[table.label <= s.phases].groupby(table.logAge.round(2)):
            xy = np.column_stack([iso.G_BPmag - iso.G_RPmag, iso.Gmag])
            dense, w = densify(xy / [s.sigma_colour, s.sigma_g], np.diff(iso.int_IMF.to_numpy()), s.step)
            self.ages.append(float(age))
            self.points.append(dense * [s.sigma_colour, s.sigma_g])
            self.weights.append(w / w.sum())
            self.trees.append(cKDTree(dense))
        self.ages = np.array(self.ages)


def densify(xy: np.ndarray, weight: np.ndarray, step: float) -> tuple[np.ndarray, np.ndarray]:
    """Points along a polyline no further than `step` apart, each segment's `weight` shared among its points."""
    seg = np.diff(xy, axis=0)
    n = np.maximum(1, np.ceil(np.hypot(*seg.T) / step).astype(int))
    t = np.concatenate([(np.arange(k) + 0.5) / k for k in n])
    start = np.repeat(np.arange(len(seg)), n)
    return xy[start] + seg[start] * t[:, None], np.repeat(np.clip(weight, 0, None) / n, n)


@dataclass
class Fit:
    age: float           # log10(age / yr), best on the grid
    av: float
    dm: float            # distance modulus
    lo: float = np.nan   # bootstrap 16th and 84th percentiles of the age
    hi: float = np.nan
    n: int = 0
    reason: str | None = None  # why no age is given, if not

    @property
    def given(self) -> bool:
        return self.reason is None


def prior(parallax: np.ndarray) -> tuple[float, float]:
    """The distance modulus from the stars' median parallax (mas), corrected for the zero point, and its
    uncertainty: the median's error and the systematics."""
    plx = np.median(parallax) + ZERO_POINT
    err = np.hypot(1.2533 * np.std(parallax) / np.sqrt(len(parallax)), PLX_SYSTEMATIC)
    return 5 * np.log10(1000 / plx) - 5, 5 / np.log(10) * err / plx


class Likelihood:
    """−2 ln L of a cluster's stars for one isochrone, as a function of (A_V, distance modulus). A star's
    density is a mixture: the isochrone's points weighted by the IMF, each blurred by the widths, normalised
    over the points that fall between the brightest and faintest magnitudes kept; and a uniform background.
    So an isochrone is penalised for predicting stars where there are none, such as a turnoff brighter
    than the brightest star."""

    def __init__(self, grid: Grid, colour: np.ndarray, g: np.ndarray, dm0: float, sdm: float):
        s = grid.settings
        self.grid, self.colour, self.g, self.dm0, self.sdm = grid, colour, g, dm0, sdm
        area = max(np.ptp(colour), 0.5) * (s.faintest - s.brightest)
        self.background = s.outliers / area  # the uniform component's density in the diagram
        self.norm = (1 - s.outliers) / (2 * np.pi * s.sigma_colour * s.sigma_g)

    def __call__(self, k: int, av: float, dm: float, idx=slice(None)) -> float:
        if not 0 <= av <= MAX_AV:
            return 1e12
        s, grid = self.grid.settings, self.grid
        shown = grid.points[k][:, 1] + dm + A_G * av
        inside = grid.weights[k][(shown >= s.brightest) & (shown <= s.faintest)].sum()
        if inside <= 0:
            return 1e12
        x = (self.colour[idx] - (A_BP - A_RP) * av) / s.sigma_colour
        y = (self.g[idx] - dm - A_G * av) / s.sigma_g
        d, j = grid.trees[k].query(np.column_stack([x, y]), k=s.neighbours, distance_upper_bound=6)
        w = np.where(np.isfinite(d), grid.weights[k][np.minimum(j, len(grid.weights[k]) - 1)], 0.0)
        density = self.norm * (w * np.exp(-0.5 * np.where(np.isfinite(d), d, 0) ** 2)).sum(axis=1) / inside
        return -2 * np.log(density + self.background).sum() + ((dm - self.dm0) / self.sdm) ** 2

    def best(self, k: int, starts, idx=slice(None)) -> tuple[float, float, float]:
        """The lowest −2 ln L over (A_V, dm) from each start, by Nelder–Mead."""
        out = (np.inf, 0.0, self.dm0)
        for av, dm in starts:
            r = minimize(lambda p: self(k, p[0], p[1], idx), [av, dm], method="Nelder-Mead",
                         options={"xatol": 1e-3, "fatol": 1e-2, "initial_simplex": [[av, dm], [av + 0.3, dm], [av, dm + 0.3]]})
            if r.fun < out[0]:
                out = (r.fun, *r.x)
        return out


def fit(grid: Grid, colour, g, parallax, boots: int = 50, window: float = 0.5, max_width: float = 0.5,
        seed: int = 0) -> Fit:
    """The age, extinction and distance modulus of one cluster's stars (BP − RP, G, parallax in mas).
    Bootstrap fits search ages within `window` dex of the best, starting from its A_V and dm there."""
    s = grid.settings
    keep = np.isfinite(colour) & np.isfinite(g) & (g <= s.faintest) & (g >= s.brightest)
    colour, g, parallax = np.asarray(colour)[keep], np.asarray(g)[keep], np.asarray(parallax)[keep]
    dm0, sdm = prior(parallax)
    like = Likelihood(grid, colour, g, dm0, sdm)
    starts = [(av, dm0) for av in (0.2, 1.0, 2.5, 4.5)]
    profile = [like.best(k, starts) for k in range(len(grid.ages))]
    k = int(np.argmin([p[0] for p in profile]))
    out = Fit(grid.ages[k], profile[k][1], profile[k][2], n=len(colour))

    near = np.flatnonzero(np.abs(grid.ages - out.age) <= window + 1e-9)
    rng, ages = np.random.default_rng(seed), []
    for _ in range(boots):
        idx = rng.integers(0, len(colour), len(colour))
        ages.append(grid.ages[near[np.argmin([like.best(j, [(out.av, out.dm)], idx)[0] for j in near])]])
    out.lo, out.hi = np.percentile(ages, [16, 84]) if boots else (np.nan, np.nan)
    if k in (0, len(grid.ages) - 1):
        out.reason = "at the edge of the grid"
    elif boots and out.hi - out.lo > max_width:
        out.reason = f"interval wider than {max_width} dex"
    return out
