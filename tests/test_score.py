import pandas as pd
import pytest

from gaia_clusters.score import score, unmatched


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
