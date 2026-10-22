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


def build(stars: pd.DataFrame, labels: pd.DataFrame, members: pd.DataFrame, table: pd.DataFrame,
          named: dict[str, str], unverified_run: int, seed: int = 0) -> dict:
    """Everything the viewer draws for one region, in parsecs.

    `stars` needs Source, GLON, GLAT, Gmag, BP-RP, r_med_geo (Bailer-Jones) and pm_l, pm_b (galactic_pm);
    `labels` the leaf runs, one column per min_cluster_size, row for row with `stars`; `members` the Hunt &
    Reffert members in the region; `table` score_region's table, indexed by name, with Hunt & Reffert's
    tidal radius (rtpc) and the cluster's proper motion as pm_l, pm_b. A recovered cluster's group, from the run
    it is shown with, is drawn at the cluster's catalogue distance: a member's own distance is uncertain
    by far more than the cluster is big. Every other star is field, at its own distance.
    """
    sizes = [int(c) for c in labels.columns]
    scored = table[table.recoverable]
    by_name = members.groupby("Name").GaiaDR3.apply(np.asarray)
    clusters, groups, cmd = [], [], []
    for name, r in scored.iterrows():
        k = shown_run(r, sizes)
        group = np.flatnonzero(labels[str(k)].to_numpy() == r[f"label_{k}"]) if r[f"label_{k}"] >= 0 else np.array([], int)
        reference = np.flatnonzero(stars.Source.isin(by_name.get(name, [])).to_numpy())
        idx, cat = cmd_stars(group, reference)
        keep = np.isfinite(stars["BP-RP"].to_numpy()[idx])  # a star without BP − RP has no place in the diagram
        cmd.append((idx[keep], cat[keep]))
        if r.recovered:
            groups.append(group)
        clusters.append({
            "name": name, "label": named.get(name, name.replace("_", " ")), "recovered": bool(r.recovered),
            "reliable": bool(r.reliable), "dist": round(float(r.dist50)), "n": int(r.n), "run": k,
            "size": len(group), "recovery": float(r[f"recovery_{k}"]), "purity": float(r[f"purity_{k}"]),
            "runs": [s for s in sizes if r[f"recovered_{s}"]], "how": None if r.recovered else how_missed(r, sizes),
            "radius": float(r.rtpc), "xyz": heliocentric(r.GLON, r.GLAT, r.dist50)[0], "pm": [float(r.pm_l), float(r.pm_b)]})

    owner = owners(groups, len(stars))
    recovered = [i for i, c in enumerate(clusters) if c["recovered"]]
    drawn = np.flatnonzero(owner >= 0)
    drawn = drawn[np.argsort(owner[drawn], kind="stable")]
    member_cluster = np.asarray(recovered)[owner[drawn]]
    member_xyz = heliocentric(stars.GLON.to_numpy()[drawn], stars.GLAT.to_numpy()[drawn],
                              scored.dist50.to_numpy()[member_cluster])
    for i, c in enumerate(clusters):
        hit = np.flatnonzero(member_cluster == i)
        c["members"] = [int(hit[0]), len(hit)] if len(hit) else None

    field = np.flatnonzero((owner < 0) & np.isfinite(stars.r_med_geo.to_numpy()))
    field = field[np.random.default_rng(seed).permutation(len(field))]
    field_xyz = heliocentric(stars.GLON.to_numpy()[field], stars.GLAT.to_numpy()[field], stars.r_med_geo.to_numpy()[field])

    lab = labels[str(unverified_run)].to_numpy()
    free = unmatched(stars.Source, lab, members)
    in_free = np.flatnonzero(np.isin(lab, free))
    g = pd.DataFrame({"label": lab[in_free], "l": stars.GLON.to_numpy()[in_free], "b": stars.GLAT.to_numpy()[in_free],
                      "d": stars.r_med_geo.to_numpy()[in_free], "pm_l": stars.pm_l.to_numpy()[in_free],
                      "pm_b": stars.pm_b.to_numpy()[in_free]}).groupby("label")
    med = g.median()
    unverified_xyz = heliocentric(med.l, med.b, med.d)

    offset = np.cumsum([0] + [len(i) for i, _ in cmd])
    for c, o, (i, _) in zip(clusters, offset, cmd):
        c["cmd"] = [int(o), len(i)]
    cmd_idx = np.concatenate([i for i, _ in cmd])
    every = np.vstack([field_xyz, member_xyz, np.array([c["xyz"] for c in clusters]).reshape(-1, 3), unverified_xyz])
    pm = lambda idx: quantize_pm(stars.pm_l.to_numpy()[idx], stars.pm_b.to_numpy()[idx])
    far = max(np.nanmax(stars.r_med_geo), max(np.linalg.norm(c["xyz"]) for c in clusters))
    reference = drift_reference(stars.r_med_geo.to_numpy()[field], stars.pm_l.to_numpy()[field],
                                stars.pm_b.to_numpy()[field], far)
    return {"clusters": clusters, "field_xyz": field_xyz, "member_xyz": member_xyz, "member_cluster": member_cluster,
            "field_pm": pm(field), "member_pm": pm(drawn), "unverified_pm": med[["pm_l", "pm_b"]].to_numpy(),
            "reference": reference, "pm_centre": np.median(np.column_stack([stars.pm_l.to_numpy()[field], stars.pm_b.to_numpy()[field]]), axis=0),
            "member_index": drawn,
            "unverified_xyz": unverified_xyz, "unverified_n": g.size().to_numpy(),
            "cmd_bp_rp": stars["BP-RP"].to_numpy()[cmd_idx], "cmd_g": stars.Gmag.to_numpy()[cmd_idx],
            "cmd_cat": np.concatenate([c for _, c in cmd]), "lo": np.floor(every.min(axis=0)) - 1,
            "hi": np.ceil(every.max(axis=0)) + 1, "stars": len(stars),
            "reach": float(np.nanpercentile(stars.r_med_geo, 99))}


def write(out: Path, built: dict, summary: dict) -> dict:
    """The files one region's page loads: meta.json, the field in nested levels (field-<k>.bin), the
    members (members.bin) and the colour-magnitude points (cmd.bin). Little-endian, as browsers read.
    Proper motions go in files of their own (field-<k>-pm.bin, members-pm.bin): the page loads them only
    to move the stars in time.
    The field is quantized; members are not, since the page comes close enough to see the steps."""
    out.mkdir(parents=True, exist_ok=True)
    for old in out.glob("field-*.bin"):  # also the -pm files
        old.unlink()
    lo, hi = built["lo"], built["hi"]
    field = quantize(built["field_xyz"], lo, hi)
    levels, start = [], 0
    for k, end in enumerate(level_bounds(len(field))):
        (out / f"field-{k}.bin").write_bytes(field[start:end].astype("<u2").tobytes())
        (out / f"field-{k}-pm.bin").write_bytes(built["field_pm"][start:end].astype("<i2").tobytes())
        levels.append({"file": f"field-{k}.bin", "pm": f"field-{k}-pm.bin", "count": end - start})
        start = end
    (out / "members.bin").write_bytes(built["member_xyz"].astype("<f4").tobytes()
                                      + built["member_cluster"].astype("<u2").tobytes())
    (out / "members-pm.bin").write_bytes(built["member_pm"].astype("<i2").tobytes())
    (out / "cmd.bin").write_bytes(np.round(built["cmd_bp_rp"] * 1000).astype("<i2").tobytes()
                                  + np.round(built["cmd_g"] * 1000).astype("<u2").tobytes()
                                  + built["cmd_cat"].astype("u1").tobytes())
    clusters = [c | {"xyz": [round(float(v), 1) for v in c["xyz"]], "recovery": round(c["recovery"], 3),
                     "purity": round(c["purity"], 3), "radius": round(c["radius"], 1),
                     "pm": [round(v, 3) for v in c["pm"]]} for c in built["clusters"]]
    meta = summary | {
        "stars": built["stars"], "reach": round(built["reach"]),
        "quant": {"lo": lo.tolist(), "hi": hi.tolist()},
        "field": {"count": len(field), "levels": levels},
        "members": {"file": "members.bin", "pm": "members-pm.bin", "count": len(built["member_xyz"])},
        "pm_step": PM_STEP,
        "reference": {"step": REFERENCE_STEP, "pm": np.round(built["reference"], 3).tolist()},
        "pm_centre": np.round(built["pm_centre"], 3).tolist(),  # the field's median motion: the middle of the motion view
        "cmd": {"file": "cmd.bin", "count": len(built["cmd_cat"])},
        "clusters": clusters,
        "unverified": {"count": len(built["unverified_n"]), "n": built["unverified_n"].tolist(),
                       "xyz": np.round(built["unverified_xyz"], 1).ravel().tolist(),
                       "pm": np.round(built["unverified_pm"], 3).ravel().tolist()}}
    (out / "meta.json").write_text(json.dumps(meta, ensure_ascii=False, separators=(",", ":")))
    return meta
