import numpy as np
import pytest

from gaia_clusters.cluster import tangent_plane
from gaia_clusters.region import GalacticBox, galactic


def test_centre_projects_to_origin():
    x, y = tangent_plane(56.75, 24.12, 56.75, 24.12)
    assert (x, y) == pytest.approx((0, 0), abs=1e-12)


def test_one_degree_is_one_degree_on_both_axes_away_from_the_equator():
    # At dec 60°, one degree of RA is only half a degree on the sky: the projection must undo that
    x, _ = tangent_plane(60 + 2, 60, 60, 60)
    _, y = tangent_plane(60, 60 + 1, 60, 60)
    assert x == pytest.approx(1.0, rel=2e-3)
    assert y == pytest.approx(1.0, rel=1e-3)


def test_east_is_positive_x_and_north_is_positive_y():
    x, y = tangent_plane(np.array([10.5, 10]), np.array([0, 0.5]), 10, 0)
    assert x[0] > 0 and y[1] > 0


def test_box_corners_map_back_to_the_box():
    box = GalacticBox(170, 180, -5, 5)
    l, b = galactic(*zip(*box.icrs_corners()))
    assert np.allclose(l, [170, 180, 180, 170], atol=1e-6)
    assert np.allclose(b, [-5, -5, 5, 5], atol=1e-6)


def test_box_wraps_longitude_and_excludes_upper_edges():
    box = GalacticBox(170, 180, -5, 5)
    assert box.contains([170, 179.99, 180, 175, 175, 530], [0, 0, 0, -5, 5, 0]).tolist() == \
        [True, True, False, True, False, True]
