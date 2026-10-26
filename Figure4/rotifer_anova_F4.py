"""Figure 4A-D: genotype x illuminance ANOVA on the rotifers eaten after 10 min (one value per chamber)."""
import numpy as np
import pandas as pd
from scipy.stats import f as f_dist, mannwhitneyu
from statsmodels.stats.multitest import multipletests

from rotifer_stats_F4 import COMPARISONS


def rotifer_anova_F4(eaten: dict[str, np.ndarray]) -> None:
    y = np.concatenate(list(eaten.values()))
    wt = np.concatenate([np.full(len(e), 1 if k.startswith('WT') else -1) for k, e in eaten.items()])
    dim = np.concatenate([np.full(len(e), 1 if k.endswith('_2lux') else -1) for k, e in eaten.items()])
    X = np.column_stack([np.ones_like(y), wt, dim, wt * dim])  # effects coding: Type III sums of squares

    def sse(cols: list[int]) -> float:
        beta = np.linalg.lstsq(X[:, cols], y, rcond=None)[0]
        return np.sum((y - X[:, cols] @ beta) ** 2)

    error, df = sse([0, 1, 2, 3]), len(y) - 4
    rows = []
    for term, col in [('genotype', 1), ('illuminance', 2), ('genotype x illuminance', 3)]:
        ss = sse([c for c in range(4) if c != col]) - error
        F = ss / (error / df)
        rows.append((term, ss, 1, F, f_dist.sf(F, 1, df), ss / (ss + error)))
    rows.append(('Error', error, df, np.nan, np.nan, np.nan))
    print('Two-way ANOVA, rotifers eaten at 10 min (%), genotype x illuminance, between chambers (Type III SS):')
    print(pd.DataFrame(rows, columns=['term', 'SumSq', 'DF', 'F', 'pValue', 'pEta2'])
          .to_string(index=False, float_format='{:.4g}'.format, na_rep=''))

    tests = [mannwhitneyu(eaten[a], eaten[b], alternative='two-sided') for a, b in COMPARISONS.values()]
    p = [t.pvalue for t in tests]
    print('Post-hoc (Mann-Whitney U, Holm corrected over the 4 comparisons):')
    for k, t, p_holm in zip(COMPARISONS, tests, multipletests(p, method='holm')[1]):
        print(f'    {k}: U={t.statistic:.1f}, p={t.pvalue:.4g}, Holm p={p_holm:.4g}')
    print()
