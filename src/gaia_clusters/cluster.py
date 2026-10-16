"""Clustering in (tangent-plane position, parallax, proper motion) and scoring against a reference."""
from dataclasses import dataclass

import numpy as np
import pandas as pd
from sklearn.cluster import HDBSCAN
from sklearn.preprocessing import RobustScaler


def tangent_plane(ra, dec, ra0: float, dec0: float) -> tuple[np.ndarray, np.ndarray]:
    """Gnomonic projection around (ra0, dec0), all in degrees, so one degree means the same on both axes."""
    ra, dec = np.radians(ra), np.radians(dec)
    ra0, dec0 = np.radians(ra0), np.radians(dec0)
    cosc = np.sin(dec0) * np.sin(dec) + np.cos(dec0) * np.cos(dec) * np.cos(ra - ra0)
    x = np.cos(dec) * np.sin(ra - ra0) / cosc
    y = (np.cos(dec0) * np.sin(dec) - np.sin(dec0) * np.cos(dec) * np.cos(ra - ra0)) / cosc
    return np.degrees(x), np.degrees(y)


def hdbscan_labels(stars: pd.DataFrame, ra0: float, dec0: float,
                   min_cluster_size: int = 40, min_samples: int = 10) -> np.ndarray:
    x, y = tangent_plane(stars.RA_ICRS, stars.DE_ICRS, ra0, dec0)
    features = np.column_stack([x, y, stars.Plx, stars.pmRA, stars.pmDE])
    X = RobustScaler().fit_transform(features)
    return HDBSCAN(min_cluster_size=min_cluster_size, min_samples=min_samples, copy=True).fit_predict(X)


@dataclass
class Match:
    label: int
    size: int
    hits: int
    reference: int  # reference members present in the input (they passed our cuts)

    @property
    def recovery(self) -> float:
        return self.hits / self.reference if self.reference else 0.0

    @property
    def purity(self) -> float:
        return self.hits / self.size if self.size else 0.0


def best_match(source_ids, labels, reference_ids) -> Match | None:
    """The cluster holding the most reference members; noise (label -1) never counts."""
    df = pd.DataFrame({"label": labels, "ref": pd.Series(source_ids).isin(set(reference_ids)).to_numpy()})
    reference = int(df.ref.sum())
    groups = df[df.label >= 0].groupby("label").agg(size=("ref", "size"), hits=("ref", "sum"))
    if groups.empty or groups.hits.max() == 0:
        return None
    label = int(groups.hits.idxmax())
    return Match(label, int(groups.loc[label, "size"]), int(groups.loc[label, "hits"]), reference)
