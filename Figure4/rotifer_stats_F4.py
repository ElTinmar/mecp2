"""Figure 4A-D: cluster-based permutation tests (Maris & Oostenveld 2007) and bootstrap CIs of the mean.

Random draws are made in the same order as in the original analysis notebook, so results are identical.
"""
import warnings

import numpy as np
from scipy.stats import t as t_dist, ttest_ind

N_RESAMPLES = 1000
ALPHA = 0.05
SEED = 0
COMPARISONS = {
    '2 lux: wild type vs mecp2-mutant': ('WT_2lux', 'MeCP2_2lux'),
    '30 lux: wild type vs mecp2-mutant': ('WT_30lux', 'MeCP2_30lux'),
    'wild type: 2 vs 30 lux': ('WT_2lux', 'WT_30lux'),
    'mecp2-mutant: 2 vs 30 lux': ('MeCP2_2lux', 'MeCP2_30lux'),
}
Cluster = tuple[int, int, float]  # first and last time index, p-value


def rotifer_stats_F4(data: dict[str, np.ndarray]) -> tuple[dict[str, np.ndarray], dict[str, list[Cluster]]]:
    """95% bootstrap CI of the mean per condition, and clusters (with p-values) per comparison."""
    rng = np.random.default_rng(SEED)
    T = len(next(iter(data.values())))
    threshold = {k: t_dist.ppf(1 - ALPHA / 2, data[a].shape[1] + data[b].shape[1] - 2)
                 for k, (a, b) in COMPARISONS.items()}
    boots = {k: np.empty((N_RESAMPLES, T)) for k in data}
    null = {k: np.empty(N_RESAMPLES) for k in COMPARISONS}
    for i in range(N_RESAMPLES):
        for k, d in data.items():
            boots[k][i] = d[:, rng.integers(0, d.shape[1], size=d.shape[1])].mean(axis=1)
        for k, (a, b) in COMPARISONS.items():
            pooled = np.hstack([data[a], data[b]])
            idx = rng.permutation(pooled.shape[1])
            n = data[a].shape[1]
            t_perm = welch_t(pooled[:, idx[:n]], pooled[:, idx[n:]])
            null[k][i] = max((mass(t_perm, c) for c in clusters(t_perm, threshold[k])), default=0)

    found = {}
    for k, (a, b) in COMPARISONS.items():
        t_obs = welch_t(data[a], data[b])
        found[k] = [(c[0], c[-1], (np.sum(null[k] >= mass(t_obs, c)) + 1) / (N_RESAMPLES + 1))
                    for c in clusters(t_obs, threshold[k])]
    ci = {k: np.quantile(b, [ALPHA / 2, 1 - ALPHA / 2], axis=0) for k, b in boots.items()}
    return ci, found


def welch_t(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Welch's t at each time point (time x chamber arrays)."""
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)  # undefined at t=0, where every count is 1
        return ttest_ind(a, b, axis=1, equal_var=False).statistic


def clusters(t_vals: np.ndarray, threshold: float) -> list[np.ndarray]:
    """Runs of consecutive time points with |t| above threshold."""
    edges = np.diff(np.r_[0, (np.abs(t_vals) > threshold).astype(int), 0])
    return [np.arange(s, e) for s, e in zip(np.flatnonzero(edges == 1), np.flatnonzero(edges == -1))]


def mass(t_vals: np.ndarray, cluster: np.ndarray) -> float:
    return np.abs(t_vals[cluster]).sum()
