import pytest

from gaia_clusters.predict import predict_share


def test_each_bin_counts_as_often_as_the_new_region_has_clusters_in_it():
    # The reference recovered 90 % of its small clusters and 10 % of its big ones
    k_ref, n_ref = [90, 10], [100, 100]
    assert predict_share([50, 0], k_ref, n_ref)[0] == pytest.approx(0.9, abs=0.01)
    assert predict_share([0, 50], k_ref, n_ref)[0] == pytest.approx(0.1, abs=0.01)
    assert predict_share([30, 10], k_ref, n_ref)[0] == pytest.approx((30 * 0.9 + 10 * 0.1) / 40, abs=0.01)


def test_a_smaller_reference_sample_widens_the_range():
    _, lo_small, hi_small = predict_share([100], [8], [10])
    _, lo_big, hi_big = predict_share([100], [800], [1000])
    assert hi_small - lo_small > 1.5 * (hi_big - lo_big)


def test_a_perfect_reference_bin_still_allows_misses():
    mean, lo, _ = predict_share([100], [10], [10])
    assert mean < 1 and lo < 1
