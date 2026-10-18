import pandas as pd
import pytest

from gaia_clusters.score import how_missed, score, shown_run, unmatched


def members(**clusters):
    return pd.DataFrame([(name, sid) for name, ids in clusters.items() for sid in ids], columns=["Name", "GaiaDR3"])


def row(df, name):
    return df.set_index("Name").loc[name]


def test_recovery_and_purity_count_only_members_in_the_sample():
    ids = list(range(20))
    labels = [0] * 12 + [-1] * 8
    ref = members(A=list(range(10)) + [99, 98])  # 98 and 99 failed our cuts
    r = row(score(ids, labels, ref, ["A"]), "A")
    assert (r.n, r.hits, r["size"]) == (10, 10, 12)
    assert r.recovery == pytest.approx(1.0)
    assert r.purity == pytest.approx(10 / 12)
    assert r.recovered


def test_picks_the_cluster_with_most_members_not_the_biggest():
    ids = list(range(60))
    labels = [0] * 40 + [1] * 20  # the field came out as cluster 0
    ref = members(A=[0, 1] + list(range(42, 60)))
    r = row(score(ids, labels, ref, ["A"]), "A")
    assert r.label == 1 and r.hits == 18 and r.recovered


def test_a_split_cluster_is_not_recovered():
    ids = list(range(30))
    labels = [0] * 15 + [1] * 15
    r = row(score(ids, labels, members(A=ids), ["A"]), "A")
    assert r.recovery == pytest.approx(0.5) and not r.recovered  # exactly half is not more than half


def test_two_reference_clusters_merged_into_one_recover_neither():
    ids = list(range(30))
    labels = [0] * 30
    s = score(ids, labels, members(A=ids[:15], B=ids[15:]), ["A", "B"]).set_index("Name")
    assert s.recovery.tolist() == [1.0, 1.0]
    assert not s.recovered.any()


def test_noise_never_matches():
    ids = list(range(20))
    r = row(score(ids, [-1] * 20, members(A=ids), ["A"]), "A")
    assert r.label == -1 and r.hits == 0 and not r.recovered


def test_clusters_with_too_few_members_in_the_sample_are_not_recoverable():
    ids = list(range(12))
    labels = [0] * 12
    s = score(ids, labels, members(A=ids[:9], B=[500, 501]), ["A", "B"]).set_index("Name")
    assert not s.recoverable.any() and not s.recovered.any()
    assert s.n.tolist() == [9, 0]


def test_unmatched_flags_clusters_without_a_majority_owner():
    ids = list(range(40))
    labels = [0] * 10 + [1] * 10 + [2] * 10 + [-1] * 10
    ref = members(A=ids[:6], B=ids[10:15], Outside=ids[30:40])  # B holds exactly half of cluster 1
    assert unmatched(ids, labels, ref).tolist() == [1, 2]  # cluster 2 has no reference member at all


def runs(**by_size):
    """A score_region row from {size: (recovery, purity)}; a run recovers when both exceed one half."""
    r = {}
    for k, (rec, pur) in by_size.items():
        r[f"recovery_{k}"], r[f"purity_{k}"], r[f"recovered_{k}"] = rec, pur, rec > 0.5 and pur > 0.5
    r["recovered"] = any(r[f"recovered_{k}"] for k in by_size)
    return pd.Series(r)


SIZES = [10, 20, 40, 80]


def test_a_recovered_cluster_is_shown_with_a_run_that_recovers_it():
    # min_cluster_size 10 has the best product, 0.5, but its purity is not above one half
    r = runs(**{"10": (1.0, 0.5), "20": (0.6, 0.6), "40": (0.55, 0.8), "80": (0.0, 0.0)})
    assert shown_run(r, SIZES) == 40


def test_a_missed_cluster_is_shown_with_the_run_that_came_closest():
    r = runs(**{"10": (0.9, 0.3), "20": (0.45, 0.45), "40": (0.2, 1.0), "80": (0.0, 0.0)})
    assert shown_run(r, SIZES) == 20
    assert how_missed(r, SIZES) == "not found"
    assert how_missed(runs(**{"10": (0.8, 0.4), "20": (0, 0), "40": (0, 0), "80": (0, 0)}), SIZES).startswith("inside")
    assert how_missed(runs(**{"10": (0.4, 0.8), "20": (0, 0), "40": (0, 0), "80": (0, 0)}), SIZES).startswith("only partly")


def test_ties_go_to_the_smallest_min_cluster_size():
    same = (0.8, 0.9)
    assert shown_run(runs(**{"10": same, "20": same, "40": same, "80": same}), SIZES) == 10
