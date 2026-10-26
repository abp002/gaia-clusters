import json

import numpy as np
import pandas as pd
import pytest

from gaia_clusters.score import score
from gaia_clusters.region import galactic
from gaia_clusters.viewer import (PM_STEP, build, cmd_stars, drift_reference, galactic_pm, near, rates, heliocentric, level_bounds, owners, quantize,
                                   quantize_pm, view, write)

SIZES = [10, 20, 40, 80]


def test_heliocentric_axes_point_to_the_centre_the_rotation_and_the_pole():
    xyz = heliocentric([0, 90, 180, 33], [0, 0, 0, 90], 1000)
    assert xyz == pytest.approx(np.array([[1000, 0, 0], [0, 1000, 0], [-1000, 0, 0], [0, 0, 1000]]), abs=1e-9)


def test_galactic_proper_motion_is_how_far_the_star_moves_in_l_and_b():
    rng = np.random.default_rng(1)
    ra, dec = rng.uniform(0, 360, 50), rng.uniform(-80, 80, 50)
    pmra, pmdec = rng.normal(0, 20, 50), rng.normal(0, 20, 50)
    pm_l, pm_b = galactic_pm(ra, dec, pmra, pmdec)
    years, mas = 100, 1 / 3.6e6  # degrees per mas
    l0, b0 = galactic(ra, dec)
    l1, b1 = galactic(ra + pmra * years * mas / np.cos(np.radians(dec)), dec + pmdec * years * mas)
    dl = ((l1 - l0 + 180) % 360 - 180) * np.cos(np.radians(b0))
    # a finite step on a curved sky: good to about 0.1 %, while a wrong sign or axis is off by the whole motion
    assert dl / mas / years == pytest.approx(pm_l, rel=2e-3, abs=1e-3)
    assert (b1 - b0) / mas / years == pytest.approx(pm_b, rel=2e-3, abs=1e-3)


def test_quantized_proper_motions_come_back_within_half_a_step_and_clip():
    pm = quantize_pm([1.234, -500.0, np.nan], [-0.005, 2.0, 1.0]) * PM_STEP
    assert pm[0] == pytest.approx([1.23, -0.0], abs=PM_STEP / 2 + 1e-12)
    assert pm[1, 0] == pytest.approx(-327.67) and pm[2, 0] == 0


def test_the_reference_is_the_field_median_by_distance_and_thin_bins_take_a_nearer_one():
    rng = np.random.default_rng(2)
    d = np.r_[rng.uniform(0, 1000, 4000), rng.uniform(1500, 1750, 20)]  # nothing in 1000–1500, too few past 1500
    pm_l, pm_b = -5 - d / 1000 + rng.normal(0, 0.05, len(d)), np.full(len(d), 0.5)
    ref = drift_reference(d, pm_l, pm_b, far=1800, step=250, least=100)
    assert ref.shape == (8, 2)
    assert ref[:4, 0] == pytest.approx(-5 - np.array([125, 375, 625, 875]) / 1000, abs=0.03)
    assert (ref[4:] == ref[3]).all() and (ref[:, 1] == 0.5).all()


def test_the_field_near_a_cluster_shares_its_patch_and_distance_and_leaves_out_its_own():
    l = np.array([100.0, 100.5, 101.0, 100.2, 100.2, 100.1, 100.0])
    b = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.0])
    d = np.array([1000, 1000, 1000, 700, 1300, 1000, 1000.0])
    # 3 tidal radii of 4 pc at 1 kpc: 0.69°. Out: 1° away, 30 % nearer, 30 % farther, and star 6 (its own)
    assert near(l, b, d, (100.0, 0.0), 1000, 4.0, exclude=np.array([6])).tolist() == [0, 1, 5]
    assert near(l, b, d, (100.0, 0.0), 1000, 0.1, exclude=np.array([], int)).tolist() == [0, 5, 6]  # 0.3° at least
    many = near(np.full(1000, 100.0), np.zeros(1000), np.full(1000, 1000.0), (100.0, 0.0), 1000, 4.0, np.array([], int), most=50)
    assert len(many) == 50 and (np.diff(many) > 0).all()


def test_rates_count_only_recoverable_clusters_in_every_bin_in_order():
    table = pd.DataFrame({"recoverable": [True, True, True, False, True], "recovered": [True, False, True, False, False],
                          "members": ["10–29", "10–29", "≥ 100", "10–29", "≥ 100"]})
    assert rates(table, "members", ["10–29", "30–99", "≥ 100"]) == [
        {"bin": "10–29", "recovered": 1, "total": 2}, {"bin": "30–99", "recovered": 0, "total": 0},
        {"bin": "≥ 100", "recovered": 1, "total": 2}]


def test_each_level_doubles_and_the_last_ends_at_the_field():
    assert level_bounds(300_000) == [65_536, 131_072, 262_144, 300_000]
    assert level_bounds(1_000) == [1_000]


def test_quantized_positions_come_back_within_half_a_step():
    rng = np.random.default_rng(0)
    lo, hi = np.array([-8000.0, -1.0, -600.0]), np.array([-10.0, 1300.0, 600.0])
    xyz = lo + rng.random((10_000, 3)) * (hi - lo)
    back = lo + quantize(xyz, lo, hi) * (hi - lo) / 65535
    assert np.abs(back - xyz).max() <= (hi - lo).max() / 65535 / 2 + 1e-9


def test_a_star_in_nested_groups_goes_to_the_smallest():
    big, small = np.arange(0, 50), np.arange(10, 20)
    for groups in ([big, small], [small, big]):
        owner = owners(groups, 60)
        i_big, i_small = (0, 1) if groups[0] is big else (1, 0)
        assert (owner[10:20] == i_small).all() and (owner[:10] == i_big).all() and (owner[50:] == -1).all()


def test_the_diagram_categories_give_back_recovery_and_purity():
    ids = list(range(40))
    labels = [0] * 25 + [-1] * 15
    members = pd.DataFrame({"Name": "A", "GaiaDR3": list(range(5, 35))})  # 20 in the group, 10 outside it
    s = score(ids, labels, members, ["A"]).iloc[0]
    _, cat = cmd_stars(np.arange(25), np.arange(5, 35))
    both, ours, theirs = (cat == 0).sum(), (cat == 1).sum(), (cat == 2).sum()
    assert both / (both + ours) == pytest.approx(s.purity)
    assert both / (both + theirs) == pytest.approx(s.recovery)


def region():
    """A small region: cluster A recovered by every run, cluster B inside a bigger group at
    min_cluster_size 10 and missed, a group of 12 stars that nobody owns, and the field."""
    rng = np.random.default_rng(3)
    n = 600
    stars = pd.DataFrame({"Source": np.arange(1, n + 1), "GLON": rng.uniform(170, 180, n), "GLAT": rng.uniform(-5, 5, n),
                          "Gmag": rng.uniform(10, 19, n), "BP-RP": rng.uniform(0, 3, n), "r_med_geo": rng.uniform(100, 4000, n)})
    a, b, extra, free = np.arange(0, 40), np.arange(40, 60), np.arange(60, 85), np.arange(85, 97)
    stars.loc[a, ["GLON", "GLAT"]] = np.column_stack([rng.normal(175, 0.1, 40), rng.normal(1, 0.1, 40)])
    stars.loc[a[0], "BP-RP"] = np.nan    # no colour: left out of the diagram
    stars.loc[n - 1, "r_med_geo"] = np.nan  # no distance: left out of the field
    stars["pm_l"], stars["pm_b"] = stars.r_med_geo / 100, -stars.GLAT  # so a star's motion can be told from its place
    labels = pd.DataFrame(-1, index=stars.index, columns=[str(k) for k in SIZES])
    labels.loc[a] = 0
    labels.loc[np.r_[b, extra], "10"] = 1
    labels.loc[free, "10"] = 2
    members = pd.DataFrame({"Name": ["A"] * 40 + ["B"] * 20, "GaiaDR3": stars.Source.to_numpy()[np.r_[a, b]]})
    table = pd.DataFrame({"GLON": [175.0, 177.0], "GLAT": [1.0, -2.0], "dist50": [1500.0, 2500.0],
                          "reliable": [True, False], "rtpc": [8.0, 12.5], "pm_l": [-3.5, 1.25], "pm_b": [0.5, -2.0],
                          "logAge16": [7.9, 8.5], "logAge50": [8.0, 8.7], "logAge84": [8.1, 8.9]},
                         index=pd.Index(["A", "B"], name="Name"))
    return stars, labels, members, rescore(table, stars, labels, members)


def rescore(table, stars, labels, members):
    """score_region's columns for every run in `labels`, and recovered by any of them."""
    for k in labels.columns:
        s = score(stars.Source, labels[k], members, table.index).set_index("Name")
        table["n"], table["recoverable"] = s.n, s.recoverable
        table[f"recovered_{k}"], table[f"label_{k}"] = s.recovered, s.label
        table[f"recovery_{k}"], table[f"purity_{k}"] = s.recovery, s.purity
    table["recovered"] = table[[f"recovered_{k}" for k in labels.columns]].any(axis=1)
    return table


def test_the_written_files_draw_every_star_once_and_in_the_right_place(tmp_path):
    stars, labels, members, table = region()
    meta = write(tmp_path, build(stars, labels, members, table, {}, unverified_run=10), {"region": "test"})
    a, b = meta["clusters"]
    assert (a["recovered"], b["recovered"], b["how"]) == (True, False, "inside a bigger group (purity ≤ 0.5)")
    assert (a["radius"], b["radius"]) == (8.0, 12.5)

    raw = (tmp_path / "members.bin").read_bytes()
    m = meta["members"]["count"]
    xyz = np.frombuffer(raw, "<f4", 3 * m).reshape(-1, 3)
    cluster = np.frombuffer(raw, "<u2", m, offset=12 * m)
    assert m == 40 and (cluster == 0).all() and a["members"] == [0, 40] and b["members"] is None
    assert np.linalg.norm(xyz, axis=1) == pytest.approx(1500, rel=1e-6)  # at the cluster's distance

    lo, hi = np.array(meta["quant"]["lo"]), np.array(meta["quant"]["hi"])
    field = np.vstack([np.frombuffer((tmp_path / lv["file"]).read_bytes(), "<u2").reshape(-1, 3) for lv in meta["field"]["levels"]])
    assert len(field) == meta["field"]["count"] == len(stars) - 40 - 1
    own = np.sort(np.linalg.norm(lo + field * (hi - lo) / 65535, axis=1))
    expected = np.sort(stars.r_med_geo.to_numpy()[40:-1])  # everyone but A, and the star without a distance
    assert own == pytest.approx(expected, abs=(hi - lo).max() / 65535)

    raw, k = (tmp_path / "cmd.bin").read_bytes(), meta["cmd"]["count"]
    cat = np.frombuffer(raw, "u1", k, offset=4 * k)
    counts = lambda c: np.bincount(cat[c["cmd"][0]:sum(c["cmd"])], minlength=3).tolist()
    assert counts(a) == [39, 0, 0] and counts(b) == [20, 25, 0]

    assert meta["unverified"]["count"] == 2 and sorted(meta["unverified"]["n"]) == [12, 45]


def test_every_star_keeps_its_own_proper_motion(tmp_path):
    stars, labels, members, table = region()
    meta = write(tmp_path, build(stars, labels, members, table, {}, unverified_run=10), {"region": "test"})
    step = meta["pm_step"]
    assert [c["pm"] for c in meta["clusters"]] == [[-3.5, 0.5], [1.25, -2.0]]
    ref = np.array(meta["reference"]["pm"])
    assert len(ref) * meta["reference"]["step"] >= 4000 and np.isfinite(ref).all()

    lo, hi = np.array(meta["quant"]["lo"]), np.array(meta["quant"]["hi"])
    for lv in meta["field"]["levels"]:  # the field is shuffled: each level's motions must follow its positions
        xyz = lo + np.frombuffer((tmp_path / lv["file"]).read_bytes(), "<u2").reshape(-1, 3) * (hi - lo) / 65535
        pm = np.frombuffer((tmp_path / lv["pm"]).read_bytes(), "<i2").reshape(-1, 2) * step
        assert len(pm) == lv["count"]
        assert pm[:, 0] == pytest.approx(np.linalg.norm(xyz, axis=1) / 100, abs=step / 2 + (hi - lo).max() / 65535 / 100)
        assert pm[:, 1] == pytest.approx(-np.degrees(np.arcsin(xyz[:, 2] / np.linalg.norm(xyz, axis=1))), abs=0.01)

    pm = np.frombuffer((tmp_path / meta["members"]["pm"]).read_bytes(), "<i2").reshape(-1, 2) * step
    assert np.sort(pm[:, 1]) == pytest.approx(np.sort(-stars.GLAT.to_numpy()[:40]), abs=step / 2)

    a, b = meta["clusters"]  # the diagram's stars carry their own motion: A's are its 39 stars with a colour
    cmd = np.frombuffer((tmp_path / meta["cmd"]["pm"]).read_bytes(), "<i2").reshape(-1, 2) * step
    assert len(cmd) == meta["cmd"]["count"]
    assert np.sort(cmd[a["cmd"][0]:sum(a["cmd"]), 1]) == pytest.approx(np.sort(-stars.GLAT.to_numpy()[1:40]), abs=step / 2)
    around = np.frombuffer((tmp_path / meta["near"]["file"]).read_bytes(), "<i2").reshape(-1, 2) * step
    assert len(around) == meta["near"]["count"] == a["near"][1] + b["near"][1] and b["near"][0] == a["near"][1]
    assert (np.abs(around[:a["near"][1], 0] * 100 - 1500) < 1500 * 0.25 + 1).all()  # pm_l is distance / 100

    u = np.array(meta["unverified"]["pm"]).reshape(-1, 2)
    free = stars.iloc[85:97]  # the group of 12 that nobody owns
    assert [free.pm_l.median(), free.pm_b.median()] == pytest.approx(u[np.array(meta["unverified"]["n"]) == 12][0], abs=1e-3)
    assert json.loads((tmp_path / "meta.json").read_text())["region"] == "test"


def test_a_single_run_shows_what_that_run_recovers_and_the_score_any_run(tmp_path):
    stars, labels, members, table = region()
    labels.loc[np.arange(0, 40), "80"] = -1  # A, recovered by the other runs, is lost at 80
    table = rescore(table, stars, labels, members)
    runs = {k: view(stars, labels, members, table, run=k) for k in SIZES}
    meta = write(tmp_path, build(stars, labels, members, table, {}, unverified_run=10), {"region": "test"}, runs, SIZES)

    assert [c["recovered"] for c in meta["clusters"]] == [True, False]  # the score: any run
    assert [(s["run"], s["recovered"], s["groups"], s["unverified"]) for s in meta["runs"]["sweep"]] == [
        (10, 1, 3, 2), (20, 1, 1, 0), (40, 1, 1, 0), (80, 0, 0, 0)]
    assert meta["runs"]["protocol"] == SIZES

    lost = json.loads((tmp_path / "run-80.json").read_text())
    a, b = lost["clusters"]
    assert (a["recovered"], a["how"], a["members"], lost["members"]["count"]) == (False, "not found", None, 0)
    assert b["how"] == "not found"
    assert (tmp_path / lost["members"]["file"]).read_bytes() == b""

    fine = json.loads((tmp_path / "run-10.json").read_text())
    a, b = fine["clusters"]
    assert (a["recovered"], a["members"], b["how"]) == (True, [0, 40], "inside a bigger group (purity ≤ 0.5)")
    assert sorted(fine["unverified"]["n"]) == [12, 45]
    k = fine["cmd"]["count"]
    cat = np.frombuffer((tmp_path / fine["cmd"]["file"]).read_bytes(), "u1", k, offset=4 * k)
    assert np.bincount(cat[b["cmd"][0]:sum(b["cmd"])], minlength=3).tolist() == [20, 25, 0]

    write(tmp_path, build(stars, labels, members, table, {}, unverified_run=10), {"region": "test"})
    assert not list(tmp_path.glob("run-*"))  # a write without runs leaves none from before
