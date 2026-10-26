from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from gaia_clusters.ages import A_BP, A_G, A_RP, ZERO_POINT, Grid, densify, fit, prior

ISOCHRONES = pd.read_csv(Path(__file__).parent / "data" / "parsec_solar_7.9_9.3.csv.gz")


def cluster(age, av, dm, n=400, outliers=40, seed=0):
    """Stars drawn from the isochrone by its IMF, blurred, reddened and shifted, plus field stars; parallaxes
    that agree with `dm`."""
    rng = np.random.default_rng(seed)
    iso = ISOCHRONES[ISOCHRONES.logAge.round(2) == age]
    w = np.clip(np.diff(iso.int_IMF.to_numpy()), 0, None)
    j = rng.choice(len(w), n, p=w / w.sum())
    colour = (iso.G_BPmag - iso.G_RPmag).to_numpy()[j] + (A_BP - A_RP) * av + rng.normal(0, 0.03, n)
    g = iso.Gmag.to_numpy()[j] + dm + A_G * av + rng.normal(0, 0.05, n)
    bright = (g > 6) & (g < 17.5)
    colour, g = np.r_[colour[bright], rng.uniform(0, 3, outliers)], np.r_[g[bright], rng.uniform(8, 17.5, outliers)]
    plx = 1000 / 10 ** (dm / 5 + 1) - ZERO_POINT + rng.normal(0, 0.05, len(g))
    return colour, g, plx


@pytest.mark.parametrize("age,av", [(8.3, 0.4), (8.9, 1.2)])
def test_a_synthetic_cluster_gives_back_its_age_and_extinction(age, av):
    grid = Grid(ISOCHRONES)
    f = fit(grid, *cluster(age, av, dm=10.0), boots=10)
    assert f.given and f.age == pytest.approx(age, abs=0.051)
    assert f.av == pytest.approx(av, abs=0.1) and f.dm == pytest.approx(10.0, abs=0.15)
    assert f.lo <= f.age <= f.hi


def test_no_age_at_the_edge_of_the_grid():
    grid = Grid(ISOCHRONES[ISOCHRONES.logAge.round(2) >= 8.5])  # the cluster is younger than the grid
    assert fit(grid, *cluster(8.0, 0.3, dm=9.0), boots=0).reason == "at the edge of the grid"


def test_the_prior_is_the_distance_of_the_corrected_median_parallax():
    dm, sigma = prior(np.full(100, 1.0 - ZERO_POINT))
    assert dm == pytest.approx(10.0) and sigma == pytest.approx(5 / np.log(10) * 0.01, rel=1e-6)


def test_densified_points_keep_the_weight_and_the_path():
    xy, w = densify(np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 3.0]]), np.array([0.25, 0.75]), step=0.5)
    assert w.sum() == pytest.approx(1.0) and len(xy) == 2 + 6
    assert (np.diff(xy[:2, 0]) > 0).all() and (xy[2:, 0] == 1.0).all()
