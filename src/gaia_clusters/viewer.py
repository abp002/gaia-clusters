"""What the 3D viewer draws, and the compact form it is sent in. scripts/export_viewer.py writes the files."""
import json
from pathlib import Path

import astropy.units as u
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord

from gaia_clusters.score import how_missed, shown_run, unmatched

FIRST_LEVEL = 1 << 16  # field stars in the coarsest level; each further level doubles them
PM_STEP = 0.01         # mas/yr per step of the 16-bit proper motions; ±327 mas/yr, past which they clip
REFERENCE_STEP = 250   # pc: the distance bins of the field's typical proper motion
REFERENCE_LEAST = 100  # stars a bin needs for its own median; a thinner one takes its nearer neighbour's
NEAR_MOST = 400        # field stars at most behind a cluster's proper-motion diagram


def heliocentric(l, b, d) -> np.ndarray:
    """Galactic Cartesian coordinates, in the unit of `d`, with the Sun at the origin: x towards the
    Galactic centre, y towards l = 90°, z towards the north Galactic pole. One row per star."""
    l, b, d = np.radians(l), np.radians(b), np.asarray(d, dtype=float)
    return np.column_stack([d * np.cos(b) * np.cos(l), d * np.cos(b) * np.sin(l), d * np.sin(b)])


def galactic_pm(ra, dec, pmra, pmdec) -> tuple[np.ndarray, np.ndarray]:
    """Proper motion along Galactic longitude (μ_l cos b) and latitude, in mas/yr, from Gaia's ICRS
    μ_α* (already times cos δ) and μ_δ."""
    c = SkyCoord(ra=np.asarray(ra) * u.deg, dec=np.asarray(dec) * u.deg, pm_ra_cosdec=np.asarray(pmra) * u.mas / u.yr,
                 pm_dec=np.asarray(pmdec) * u.mas / u.yr).galactic
    return c.pm_l_cosb.to_value(u.mas / u.yr), c.pm_b.to_value(u.mas / u.yr)


def drift_reference(d, pm_l, pm_b, far: float, step: float = REFERENCE_STEP,
                    least: int = REFERENCE_LEAST) -> np.ndarray:
    """The field's median proper motion (μ_l cos b, μ_b) in bins of distance `step` wide, from 0 to `far`.
    The viewer moves stars relative to it: the Sun's own motion and the Galaxy's rotation shift every star
    at a distance alike, by several times the spread between them. A bin with fewer than `least` stars
    takes the median of the nearest bin closer to the Sun that has them (the first such bin if none);
    with no such bin at all, every bin takes the whole field's median."""
    d, pm = np.asarray(d, float), np.column_stack([pm_l, pm_b])
    edges = np.arange(0, far + step, step)
    out = np.full((len(edges) - 1, 2), np.nan)
    for i, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
        inside = (d >= lo) & (d < hi)
        if inside.sum() >= least:
            out[i] = np.median(pm[inside], axis=0)
    have = np.flatnonzero(~np.isnan(out[:, 0]))
    if not len(have):  # too few stars anywhere: the whole field's median
        return np.tile(np.median(pm, axis=0), (len(out), 1))
    take = np.maximum.accumulate(np.where(~np.isnan(out[:, 0]), np.arange(len(out)), have[0]))
    return out[take]


def near(l, b, d, centre: tuple[float, float], dist: float, radius_pc: float, exclude: np.ndarray,
         most: int = NEAR_MOST, rng: np.random.Generator | None = None) -> np.ndarray:
    """The stars that share a cluster's patch of sky and distance, less `exclude` (its own): the field its
    proper-motion diagram is drawn against. The patch reaches three tidal radii from the centre, kept
    between 0.3° and 3°; the distance, a quarter of the cluster's either way, since a single star's
    distance is that uncertain a few kiloparsecs out. At most `most`, picked at random, in index order."""
    theta = np.clip(np.degrees(3 * radius_pc / dist), 0.3, 3.0)
    l, b, (l0, b0) = np.radians(l), np.radians(b), np.radians(centre)
    sep = 2 * np.degrees(np.arcsin(np.sqrt(np.sin((b - b0) / 2) ** 2 + np.cos(b) * np.cos(b0) * np.sin((l - l0) / 2) ** 2)))
    pick = np.setdiff1d(np.flatnonzero((sep < theta) & (d > dist / 1.25) & (d < dist * 1.25)), exclude)
    if len(pick) > most:
        pick = np.sort((rng or np.random.default_rng(0)).choice(pick, most, replace=False))
    return pick


def quantize_pm(pm_l, pm_b) -> np.ndarray:
    """Proper motions as pairs of 16-bit steps of PM_STEP, one row per star; the viewer multiplies back."""
    steps = np.round(np.column_stack([pm_l, pm_b]) / PM_STEP)
    return np.clip(np.nan_to_num(steps), -32767, 32767).astype(np.int16)


def level_bounds(n: int, first: int = FIRST_LEVEL) -> list[int]:
    """Where each level of detail ends in a field stored in random order. Level k is the first
    first·2^k stars, so every level is a uniform random sample that contains the previous one."""
    bounds = [min(first, n)]
    while bounds[-1] < n:
        bounds.append(min(2 * bounds[-1], n))
    return bounds


def quantize(xyz: np.ndarray, lo: np.ndarray, hi: np.ndarray) -> np.ndarray:
    """Positions as 16-bit steps from `lo` to `hi`; the viewer decodes lo + q · (hi − lo) / 65535."""
    return np.round((xyz - lo) / (hi - lo) * 65535).astype(np.uint16)


def owners(groups: list[np.ndarray], n: int) -> np.ndarray:
    """For each of `n` stars, the index of the group it is drawn in, or -1. Groups shown from different
    runs can nest: a small cluster's group at min_cluster_size 10 sits inside a neighbour's group at
    20. A star goes to the smallest group that holds it, so the small cluster keeps its stars."""
    owner = np.full(n, -1)
    size = np.full(n, np.iinfo(np.int64).max)
    for i, idx in enumerate(groups):
        take = idx[len(idx) < size[idx]]
        owner[take], size[take] = i, len(idx)
    return owner


def cmd_stars(group: np.ndarray, reference: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """The stars of a cluster's colour-magnitude diagram, as indices and a category: 0 in our group and
    in Hunt & Reffert, 1 in our group only, 2 a Hunt & Reffert member in our sample that the group
    misses. Purity is the share of 0 among 0 and 1; recovery, the share of 0 among 0 and 2."""
    both = np.intersect1d(group, reference)
    ours, theirs = np.setdiff1d(group, reference), np.setdiff1d(reference, group)
    return np.concatenate([both, ours, theirs]), np.repeat([0, 1, 2], [len(both), len(ours), len(theirs)])


def rates(table: pd.DataFrame, column: str, order: list[str]) -> list[dict]:
    """Per bin of `column` (score_region's distance or members), in `order`: how many recoverable clusters
    it holds and how many of them were recovered by at least one run."""
    rec = table[table.recoverable]
    return [{"bin": b, "recovered": int(rec.recovered[rec[column] == b].sum()), "total": int((rec[column] == b).sum())}
            for b in order]


def view(stars: pd.DataFrame, labels: pd.DataFrame, members: pd.DataFrame, scored: pd.DataFrame,
         run: int | None = None, unverified_run: int | None = None) -> dict:
    """The clusters as one view of the region draws them. With `run` None, each scored cluster comes from
    the run it is shown with (shown_run) and counts as recovered if any run recovers it: the score. With a
    min_cluster_size, every cluster comes from that run alone. `labels` holds the runs, one column per
    min_cluster_size; `scored` needs recovered_<k>, label_<k>, recovery_<k> and purity_<k> for each of them.

    A recovered cluster's group is drawn at the cluster's catalogue distance; its stars are `member_index`.
    The groups no Hunt & Reffert cluster owns come from `unverified_run` (`run` if not given)."""
    sizes = [int(c) for c in labels.columns]
    pm = lambda idx: quantize_pm(stars.pm_l.to_numpy()[idx], stars.pm_b.to_numpy()[idx])
    by_name = members.groupby("Name").GaiaDR3.apply(np.asarray)
    clusters, groups, cmd = [], [], []
    for name, r in scored.iterrows():
        k = shown_run(r, sizes) if run is None else run
        recovered = bool(r.recovered if run is None else r[f"recovered_{k}"])
        group = np.flatnonzero(labels[str(k)].to_numpy() == r[f"label_{k}"]) if r[f"label_{k}"] >= 0 else np.array([], int)
        reference = np.flatnonzero(stars.Source.isin(by_name.get(name, [])).to_numpy())
        idx, cat = cmd_stars(group, reference)
        keep = np.isfinite(stars["BP-RP"].to_numpy()[idx])  # a star without BP − RP has no place in the diagram
        cmd.append((idx[keep], cat[keep]))
        if recovered:
            groups.append(group)
        if run is not None:  # how_missed judges by any run's verdict; here only this run's counts
            r = r.copy()
            r["recovered"] = recovered
        clusters.append({"recovered": recovered, "run": k, "size": len(group), "recovery": float(r[f"recovery_{k}"]),
                         "purity": float(r[f"purity_{k}"]), "runs": [s for s in sizes if r[f"recovered_{s}"]],
                         "how": None if recovered else how_missed(r, sizes if run is None else [k])})

    owner = owners(groups, len(stars))
    recovered = [i for i, c in enumerate(clusters) if c["recovered"]]
    drawn = np.flatnonzero(owner >= 0)
    drawn = drawn[np.argsort(owner[drawn], kind="stable")]
    member_cluster = np.asarray(recovered, dtype=int)[owner[drawn]]
    member_xyz = heliocentric(stars.GLON.to_numpy()[drawn], stars.GLAT.to_numpy()[drawn],
                              scored.dist50.to_numpy()[member_cluster])
    for i, c in enumerate(clusters):
        hit = np.flatnonzero(member_cluster == i)
        c["members"] = [int(hit[0]), len(hit)] if len(hit) else None
    offset = np.cumsum([0] + [len(i) for i, _ in cmd])
    for c, o, (i, _) in zip(clusters, offset, cmd):
        c["cmd"] = [int(o), len(i)]
    cmd_idx = np.concatenate([i for i, _ in cmd])

    lab = labels[str(unverified_run or run)].to_numpy()
    in_free = np.flatnonzero(np.isin(lab, unmatched(stars.Source, lab, members)))
    g = pd.DataFrame({"label": lab[in_free], "l": stars.GLON.to_numpy()[in_free], "b": stars.GLAT.to_numpy()[in_free],
                      "d": stars.r_med_geo.to_numpy()[in_free], "pm_l": stars.pm_l.to_numpy()[in_free],
                      "pm_b": stars.pm_b.to_numpy()[in_free]}).groupby("label")
    med = g.median()
    return {"clusters": clusters, "member_xyz": member_xyz, "member_cluster": member_cluster, "member_index": drawn,
            "member_pm": pm(drawn), "cmd_index": cmd_idx, "cmd_bp_rp": stars["BP-RP"].to_numpy()[cmd_idx],
            "cmd_g": stars.Gmag.to_numpy()[cmd_idx], "cmd_cat": np.concatenate([c for _, c in cmd]), "cmd_pm": pm(cmd_idx),
            "unverified_xyz": heliocentric(med.l, med.b, med.d), "unverified_n": g.size().to_numpy(),
            "unverified_pm": med[["pm_l", "pm_b"]].to_numpy(), "groups": int(lab.max() + 1) if run else None}


def build(stars: pd.DataFrame, labels: pd.DataFrame, members: pd.DataFrame, table: pd.DataFrame,
          named: dict[str, str], unverified_run: int, seed: int = 0) -> dict:
    """Everything the viewer draws for one region, in parsecs: the score's view (view() with no run) and
    what every view shares.

    `stars` needs Source, GLON, GLAT, Gmag, BP-RP, r_med_geo (Bailer-Jones) and pm_l, pm_b (galactic_pm);
    `labels` the leaf runs, one column per min_cluster_size, row for row with `stars`; `members` the Hunt &
    Reffert members in the region; `table` score_region's table, indexed by name, with Hunt & Reffert's
    tidal radius (rtpc), age (logAge16, logAge50, logAge84) and the cluster's proper motion as pm_l, pm_b.
    The members of the score's view are left out of the field; every other star is field, at its own distance.
    """
    scored = table[table.recoverable]
    built = view(stars, labels, members, scored, unverified_run=unverified_run)
    for c, (name, r) in zip(built["clusters"], scored.iterrows()):
        c |= {"name": name, "label": named.get(name, name.replace("_", " ")), "reliable": bool(r.reliable),
              "dist": round(float(r.dist50)), "n": int(r.n), "radius": float(r.rtpc),
              "xyz": heliocentric(r.GLON, r.GLAT, r.dist50)[0], "pm": [float(r.pm_l), float(r.pm_b)],
              "age": [round(float(r[f"logAge{q}"]), 2) for q in (16, 50, 84)]}
    clusters = built["clusters"]

    in_group = np.zeros(len(stars), bool)
    in_group[built["member_index"]] = True
    field = np.flatnonzero(~in_group & np.isfinite(stars.r_med_geo.to_numpy()))
    field = field[np.random.default_rng(seed).permutation(len(field))]
    field_xyz = heliocentric(stars.GLON.to_numpy()[field], stars.GLAT.to_numpy()[field], stars.r_med_geo.to_numpy()[field])

    every = np.vstack([field_xyz, built["member_xyz"], np.array([c["xyz"] for c in clusters]).reshape(-1, 3),
                       built["unverified_xyz"]])
    pm = lambda idx: quantize_pm(stars.pm_l.to_numpy()[idx], stars.pm_b.to_numpy()[idx])
    rng, around = np.random.default_rng(seed), []
    for c in clusters:
        own = built["cmd_index"][c["cmd"][0]:sum(c["cmd"])]
        around.append(near(stars.GLON.to_numpy(), stars.GLAT.to_numpy(), stars.r_med_geo.to_numpy(),
                           (table.GLON[c["name"]], table.GLAT[c["name"]]), c["dist"], c["radius"], own, rng=rng))
    start = np.cumsum([0] + [len(a) for a in around])
    for c, o, a in zip(clusters, start, around):
        c["near"] = [int(o), len(a)]
    far = max(np.nanmax(stars.r_med_geo), max(np.linalg.norm(c["xyz"]) for c in clusters))
    reference = drift_reference(stars.r_med_geo.to_numpy()[field], stars.pm_l.to_numpy()[field],
                                stars.pm_b.to_numpy()[field], far)
    return built | {"field_xyz": field_xyz, "field_pm": pm(field), "near_pm": pm(np.concatenate(around).astype(int)),
                    "reference": reference, "lo": np.floor(every.min(axis=0)) - 1, "hi": np.ceil(every.max(axis=0)) + 1,
                    "stars": len(stars), "reach": float(np.nanpercentile(stars.r_med_geo, 99))}


def write_view(out: Path, v: dict, prefix: str = "") -> dict:
    """One view's files (view()): its members (members.bin), its colour-magnitude points (cmd.bin) and their
    proper motions, named after `prefix`. Returns how meta.json points at them, and its unverified groups."""
    name = lambda f: f"{prefix}{f}"
    (out / name("members.bin")).write_bytes(v["member_xyz"].astype("<f4").tobytes() + v["member_cluster"].astype("<u2").tobytes())
    (out / name("members-pm.bin")).write_bytes(v["member_pm"].astype("<i2").tobytes())
    (out / name("cmd-pm.bin")).write_bytes(v["cmd_pm"].astype("<i2").tobytes())  # row for row with cmd.bin
    (out / name("cmd.bin")).write_bytes(np.round(v["cmd_bp_rp"] * 1000).astype("<i2").tobytes()
                                        + np.round(v["cmd_g"] * 1000).astype("<u2").tobytes()
                                        + v["cmd_cat"].astype("u1").tobytes())
    return {"members": {"file": name("members.bin"), "pm": name("members-pm.bin"), "count": len(v["member_xyz"])},
            "cmd": {"file": name("cmd.bin"), "pm": name("cmd-pm.bin"), "count": len(v["cmd_cat"])},
            "unverified": {"count": len(v["unverified_n"]), "n": v["unverified_n"].tolist(),
                           "xyz": np.round(v["unverified_xyz"], 1).ravel().tolist(),
                           "pm": np.round(v["unverified_pm"], 3).ravel().tolist()}}


def write(out: Path, built: dict, summary: dict, runs: dict[int, dict] | None = None, protocol=()) -> dict:
    """The files one region's page loads: meta.json, the field in nested levels (field-<k>.bin), and the
    score's view (write_view). Little-endian, as browsers read. Proper motions go in files of their own
    (field-<k>-pm.bin, members-pm.bin): the page loads them only to move the stars in time.
    `runs` are the single-run views (view() with a run) by min_cluster_size, for the slider: each in
    run-<k>.json and run-<k>-*.bin; `protocol` the sizes the score uses.
    The field is quantized; members are not, since the page comes close enough to see the steps."""
    out.mkdir(parents=True, exist_ok=True)
    for old in [*out.glob("field-*.bin"), *out.glob("run-*")]:  # also the -pm files
        old.unlink()
    lo, hi = built["lo"], built["hi"]
    field = quantize(built["field_xyz"], lo, hi)
    levels, start = [], 0
    for k, end in enumerate(level_bounds(len(field))):
        (out / f"field-{k}.bin").write_bytes(field[start:end].astype("<u2").tobytes())
        (out / f"field-{k}-pm.bin").write_bytes(built["field_pm"][start:end].astype("<i2").tobytes())
        levels.append({"file": f"field-{k}.bin", "pm": f"field-{k}-pm.bin", "count": end - start})
        start = end
    (out / "near-pm.bin").write_bytes(built["near_pm"].astype("<i2").tobytes())
    rounded = lambda c: c | {"recovery": round(c["recovery"], 3), "purity": round(c["purity"], 3)}
    clusters = [rounded(c) | {"xyz": [round(float(v), 1) for v in c["xyz"]], "radius": round(c["radius"], 1),
                              "pm": [round(v, 3) for v in c["pm"]]} for c in built["clusters"]]
    sweep = []
    for k, v in (runs or {}).items():
        files = write_view(out, v, f"run-{k}-")
        recovered = sum(c["recovered"] for c in v["clusters"])
        per = [{key: c[key] for key in ("recovered", "size", "recovery", "purity", "how", "members", "cmd")}
               for c in map(rounded, v["clusters"])]
        (out / f"run-{k}.json").write_text(json.dumps(files | {"run": k, "clusters": per}, ensure_ascii=False,
                                                      separators=(",", ":")))
        sweep.append({"run": k, "file": f"run-{k}.json", "recovered": recovered, "groups": v["groups"],
                      "unverified": files["unverified"]["count"]})
    meta = summary | {
        "stars": built["stars"], "reach": round(built["reach"]),
        "quant": {"lo": lo.tolist(), "hi": hi.tolist()},
        "field": {"count": len(field), "levels": levels},
        "pm_step": PM_STEP,
        "reference": {"step": REFERENCE_STEP, "pm": np.round(built["reference"], 3).tolist()},
        "near": {"file": "near-pm.bin", "count": len(built["near_pm"])},
        "clusters": clusters} | write_view(out, built) | {"runs": {"protocol": list(protocol), "sweep": sweep}}
    (out / "meta.json").write_text(json.dumps(meta, ensure_ascii=False, separators=(",", ":")))
    return meta


def isochrone_curves(table: pd.DataFrame, phases: int, step: float = 0.03) -> dict:
    """PARSEC isochrones for the cards' diagrams: per age, (BP − RP, M_G) points at least `step` mag apart
    along the curve, rounded to the millimagnitude. The page reddens and shifts them by a fit's A_V and dm."""
    ages, curves = [], []
    for age, iso in table[table.label <= phases].groupby(table.logAge.round(2)):
        xy = np.column_stack([iso.G_BPmag - iso.G_RPmag, iso.Gmag])
        keep, last = [0], xy[0]
        for j in range(1, len(xy)):
            if np.hypot(*(xy[j] - last)) >= step or j == len(xy) - 1:
                keep.append(j)
                last = xy[j]
        ages.append(float(age))
        curves.append(np.round(xy[keep], 3).ravel().tolist())
    return {"ages": ages, "curves": curves}
