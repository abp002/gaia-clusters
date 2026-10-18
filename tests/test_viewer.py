import json

import numpy as np
import pandas as pd
import pytest

from gaia_clusters.score import score
from gaia_clusters.viewer import build, cmd_stars, heliocentric, level_bounds, owners, quantize, write

SIZES = [10, 20, 40, 80]


def test_heliocentric_axes_point_to_the_centre_the_rotation_and_the_pole():
    xyz = heliocentric([0, 90, 180, 33], [0, 0, 0, 90], 1000)
    assert xyz == pytest.approx(np.array([[1000, 0, 0], [0, 1000, 0], [-1000, 0, 0], [0, 0, 1000]]), abs=1e-9)


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
    labels = pd.DataFrame(-1, index=stars.index, columns=[str(k) for k in SIZES])
    labels.loc[a] = 0
    labels.loc[np.r_[b, extra], "10"] = 1
    labels.loc[free, "10"] = 2
    members = pd.DataFrame({"Name": ["A"] * 40 + ["B"] * 20, "GaiaDR3": stars.Source.to_numpy()[np.r_[a, b]]})
    table = pd.DataFrame({"GLON": [175.0, 177.0], "GLAT": [1.0, -2.0], "dist50": [1500.0, 2500.0],
                          "reliable": [True, False]}, index=pd.Index(["A", "B"], name="Name"))
    for k in SIZES:
        s = score(stars.Source, labels[str(k)], members, table.index).set_index("Name")
        table["n"], table["recoverable"] = s.n, s.recoverable
        table[f"recovered_{k}"], table[f"label_{k}"] = s.recovered, s.label
        table[f"recovery_{k}"], table[f"purity_{k}"] = s.recovery, s.purity
    table["recovered"] = table[[f"recovered_{k}" for k in SIZES]].any(axis=1)
    return stars, labels, members, table


def test_the_written_files_draw_every_star_once_and_in_the_right_place(tmp_path):
    stars, labels, members, table = region()
    meta = write(tmp_path, build(stars, labels, members, table, {}, unverified_run=10), {"region": "test"})
    a, b = meta["clusters"]
    assert (a["recovered"], b["recovered"], b["how"]) == (True, False, "inside a bigger group (purity ≤ 0.5)")

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
    assert json.loads((tmp_path / "meta.json").read_text())["region"] == "test"
