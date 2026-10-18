"""What one region's recovery rates predict for another, so a replication can fail by a stated margin."""
import numpy as np


def predict_share(n_new, k_ref, n_ref, draws: int = 100_000, seed: int = 0) -> tuple[float, float, float]:
    """Share of a new region's clusters that should be recovered if every bin keeps the reference rate.

    Bins are matched by position: the new region has `n_new[i]` clusters in bin i, of which the
    reference region recovered `k_ref[i]` of `n_ref[i]`. Each reference rate is uncertain (Jeffreys
    prior, Beta(k + 1/2, n - k + 1/2)) and the new region's count is binomial on top of it.
    Returns the mean share and the central 95 % range of the predictive distribution.
    """
    n_new, k_ref, n_ref = (np.asarray(a) for a in (n_new, k_ref, n_ref))
    rng = np.random.default_rng(seed)
    p = rng.beta(k_ref + 0.5, n_ref - k_ref + 0.5, size=(draws, len(n_new)))
    share = rng.binomial(n_new, p).sum(axis=1) / n_new.sum()
    lo, hi = np.percentile(share, [2.5, 97.5])
    return float(share.mean()), float(lo), float(hi)
