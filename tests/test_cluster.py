import numpy as np
import pytest

from gaia_clusters.cluster import best_match, tangent_plane


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


def test_recovery_counts_only_reference_members_present_in_the_input():
    ids = [1, 2, 3, 4, 5]
    labels = [0, 0, 0, -1, 1]
    m = best_match(ids, labels, reference_ids=[1, 2, 4, 99])  # 99 failed our cuts
    assert (m.label, m.size, m.hits, m.reference) == (0, 3, 2, 3)
    assert m.recovery == pytest.approx(2 / 3)
    assert m.purity == pytest.approx(2 / 3)


def test_picks_the_cluster_with_most_members_not_the_biggest():
    ids = list(range(10))
    labels = [0] * 7 + [1] * 3  # the field is cluster 0, the real cluster is 1
    m = best_match(ids, labels, reference_ids=[0, 7, 8, 9])
    assert m.label == 1 and m.hits == 3


def test_noise_never_matches():
    assert best_match([1, 2], [-1, -1], reference_ids=[1, 2]) is None
