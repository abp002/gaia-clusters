import numpy as np
import pytest

from gaia_clusters.hyades import HYADES, K, adql_residual_below, basis, heliocentric_icrs, residual, selected, space_velocity

V = space_velocity(*(HYADES[k] for k in ("ra", "dec", "plx", "pmra", "pmdec", "rv")))


def test_the_hyades_velocity_is_the_published_one():
    assert V == pytest.approx([-5.735, 45.478, 4.712], abs=1e-3)  # Gaia Collaboration (2018): (−6.1, 45.6, 5.6)


def moving_with(v, ra, dec, plx, extra=(0.0, 0.0)):
    """The proper motion (mas/yr) of a star at (ra, dec, plx) moving with `v`, plus `extra` km/s along α and δ."""
    _, p, q = basis(np.asarray(ra, float), np.asarray(dec, float))
    return (p @ v + extra[0]) * plx / K, (q @ v + extra[1]) * plx / K


def test_a_star_moving_with_the_cluster_anywhere_on_the_sky_has_no_residual():
    rng = np.random.default_rng(0)
    ra, dec, plx = rng.uniform(0, 360, 200), rng.uniform(-89, 89, 200), rng.uniform(2, 30, 200)
    pmra, pmdec = moving_with(V, ra, dec, plx)
    dv, _ = residual(ra, dec, plx, 0.1, pmra, pmdec, V)
    assert dv == pytest.approx(0, abs=1e-9)
    pmra, pmdec = moving_with(V, ra, dec, plx, extra=(3.0, -4.0))
    assert residual(ra, dec, plx, 0.1, pmra, pmdec, V)[0] == pytest.approx(5.0)


def test_the_selection_widens_with_the_parallax_error():
    assert selected(np.array([2.9, 3.1, 3.1]), np.array([0.0, 0.0, 1.0])).tolist() == [True, False, True]


def test_the_adql_condition_computes_the_same_residual():
    """Evaluates the ADQL text as Python, with its SQL functions, against residual()."""
    rng = np.random.default_rng(1)
    ra, dec, parallax = rng.uniform(0, 360, 50), rng.uniform(-80, 80, 50), rng.uniform(2, 30, 50)
    pmra, pmdec = moving_with(V, ra, dec, parallax, extra=(rng.normal(0, 8, 50), rng.normal(0, 8, 50)))
    env = {"RADIANS": np.radians, "SIN": np.sin, "COS": np.cos, "POWER": np.power, "ra": ra, "dec": dec,
           "parallax": parallax, "pmra": pmra, "pmdec": pmdec}
    inside = eval(adql_residual_below(V, 10.0), env)
    assert (inside == (residual(ra, dec, parallax, 0.1, pmra, pmdec, V)[0] < 10.0)).all()


def test_positions_are_at_one_over_the_parallax():
    xyz = heliocentric_icrs([0, 90, 0], [0, 0, 90], [10, 2, 4])
    assert xyz == pytest.approx(np.array([[100, 0, 0], [0, 500, 0], [0, 0, 250]]), abs=1e-9)
