"""Clustering in (tangent-plane position, parallax, proper motion), the five axes Hunt & Reffert use."""
import tempfile

import hdbscan
import numpy as np
import pandas as pd
from joblib import Memory
from sklearn.preprocessing import RobustScaler


def tangent_plane(lon, lat, lon0: float, lat0: float) -> tuple[np.ndarray, np.ndarray]:
    """Gnomonic projection around (lon0, lat0), all in degrees, so one degree means the same on both axes.
    Works for any frame: (ra, dec) or (l, b)."""
    lon, lat = np.radians(lon), np.radians(lat)
    lon0, lat0 = np.radians(lon0), np.radians(lat0)
    cosc = np.sin(lat0) * np.sin(lat) + np.cos(lat0) * np.cos(lat) * np.cos(lon - lon0)
    x = np.cos(lat) * np.sin(lon - lon0) / cosc
    y = (np.cos(lat0) * np.sin(lat) - np.sin(lat0) * np.cos(lat) * np.cos(lon - lon0)) / cosc
    return np.degrees(x), np.degrees(y)


def features(stars: pd.DataFrame, lon0: float, lat0: float, lon="RA_ICRS", lat="DE_ICRS") -> np.ndarray:
    """Every axis rescaled to the same median and interquartile range, as Hunt & Reffert do."""
    x, y = tangent_plane(stars[lon], stars[lat], lon0, lat0)
    return RobustScaler().fit_transform(np.column_stack([x, y, stars.Plx, stars.pmRA, stars.pmDE]))


def hdbscan_runs(X: np.ndarray, min_cluster_sizes, min_samples: int = 10, method: str = "eom",
                 cache_dir: str | None = None) -> dict[int, np.ndarray]:
    """Labels for each min_cluster_size. The minimum spanning tree depends only on min_samples, so it
    is built once and cached; each extra size or selection method only re-cuts the tree. With
    `cache_dir` the tree outlives the call."""
    with tempfile.TemporaryDirectory() as tmp:
        memory = Memory(cache_dir or tmp, verbose=0)
        return {size: hdbscan.HDBSCAN(min_cluster_size=size, min_samples=min_samples, memory=memory,
                                      cluster_selection_method=method, core_dist_n_jobs=-1).fit_predict(X)
                for size in min_cluster_sizes}
