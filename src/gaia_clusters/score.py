"""Which Hunt & Reffert clusters a clustering run recovers. The rule is fixed in docs/protocol-auriga.md."""
import numpy as np
import pandas as pd

MIN_MEMBERS = 10  # the smallest min_cluster_size: with fewer members in the sample nothing can find it


def score(source_ids, labels, members: pd.DataFrame, reference) -> pd.DataFrame:
    """One row per reference cluster name.

    `members` holds Hunt & Reffert members (columns Name, GaiaDR3) of any cluster. Only members that
    are in our sample count. A cluster is recovered when one of ours holds more than half of them
    (recovery > 0.5) and more than half of that cluster's stars are its members (purity > 0.5).
    Both are strict, so the match is one-to-one. Noise (label -1) never matches.
    """
    ours = pd.DataFrame({"GaiaDR3": np.asarray(source_ids), "label": np.asarray(labels)})
    sizes = ours.label[ours.label >= 0].value_counts()
    m = members.loc[members.Name.isin(reference), ["Name", "GaiaDR3"]].merge(ours, on="GaiaDR3")
    best = (m[m.label >= 0].groupby(["Name", "label"]).size().rename("hits").reset_index()
            .sort_values(["hits", "label"], ascending=[False, True]).drop_duplicates("Name").set_index("Name"))

    out = pd.DataFrame(index=pd.Index(list(reference), name="Name"))
    out["n"] = m.groupby("Name").size().reindex(out.index, fill_value=0)
    out["label"] = best.label.reindex(out.index).fillna(-1).astype(int)
    out["hits"] = best.hits.reindex(out.index).fillna(0).astype(int)
    out["size"] = out.label.map(sizes).fillna(0).astype(int)
    out["recovery"] = (out.hits / out.n.where(out.n > 0)).fillna(0.0)
    out["purity"] = (out.hits / out["size"].where(out["size"] > 0)).fillna(0.0)
    out["recoverable"] = out.n >= MIN_MEMBERS
    out["recovered"] = out.recoverable & (out.recovery > 0.5) & (out.purity > 0.5)
    return out.reset_index()


def unmatched(source_ids, labels, members: pd.DataFrame) -> np.ndarray:
    """Labels of our clusters in which no Hunt & Reffert cluster, wherever its centre, holds more than half the stars."""
    ours = pd.DataFrame({"GaiaDR3": np.asarray(source_ids), "label": np.asarray(labels)})
    ours = ours[ours.label >= 0]
    sizes = ours.label.value_counts()
    top = members[["Name", "GaiaDR3"]].merge(ours, on="GaiaDR3").groupby(["label", "Name"]).size().groupby(level=0).max()
    share = top.reindex(sizes.index, fill_value=0) / sizes
    return np.sort(share.index[share <= 0.5].to_numpy())
