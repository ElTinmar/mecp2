"""Mixed-design (split-plot) ANOVA, as MATLAB fitrm/ranova."""
import numpy as np
import pandas as pd
from scipy.linalg import null_space
from scipy.stats import f as f_dist


def mixed_anova(groups: list[np.ndarray]) -> tuple[pd.DataFrame, float]:
    """Genotype between subjects (one array per group), stimulus within subjects (columns).

    Effects coding; within effects are Greenhouse-Geisser corrected.
    Effect sizes: partial (pEta2) and generalized (gEta2) eta squared.
    Returns the ANOVA table and the Greenhouse-Geisser epsilon.
    """
    Y = np.vstack(groups)
    n, k = Y.shape
    X = np.column_stack([np.ones(n), np.repeat([1.0, -1.0], [len(g) for g in groups])])
    XtX_inv = np.linalg.inv(X.T @ X)
    df_err = n - X.shape[1]

    def sscp(Z: np.ndarray, L: np.ndarray) -> tuple[float, np.ndarray]:
        B = XtX_inv @ X.T @ Z
        LB = L @ B
        R = Z - X @ B
        return np.trace(np.outer(LB, LB)) / (L @ XtX_inv @ L), R.T @ R

    intercept, genotype = np.array([1.0, 0.0]), np.array([0.0, 1.0])
    ss_g, E_between = sscp(Y.sum(axis=1, keepdims=True) / np.sqrt(k), genotype)
    Z = Y @ null_space(np.ones((1, k)))
    ss_stim, E = sscp(Z, intercept)
    ss_int, _ = sscp(Z, genotype)
    eps = np.trace(E) ** 2 / ((k - 1) * np.trace(E @ E))

    rows = []
    for term, ss, df, ss_e, df_e, e in [
        ('genotype', ss_g, 1, np.trace(E_between), df_err, 1.0),
        ('(Intercept):stim', ss_stim, k - 1, np.trace(E), (k - 1) * df_err, eps),
        ('genotype:stim', ss_int, k - 1, np.trace(E), (k - 1) * df_err, eps),
    ]:
        F = (ss / df) / (ss_e / df_e)
        rows.append({'term': term, 'SumSq': ss, 'DF': df, 'MeanSq': ss / df, 'F': F,
                     'pValue': f_dist.sf(F, df, df_e), 'pValueGG': f_dist.sf(F, e * df, e * df_e),
                     'pEta2': ss / (ss + ss_e), 'gEta2': ss / (ss + np.trace(E_between) + np.trace(E))})
    rows.insert(1, {'term': 'Error', 'SumSq': np.trace(E_between), 'DF': df_err, 'MeanSq': np.trace(E_between) / df_err})
    rows.append({'term': 'Error(stim)', 'SumSq': np.trace(E), 'DF': (k - 1) * df_err, 'MeanSq': np.trace(E) / ((k - 1) * df_err)})
    return pd.DataFrame(rows).set_index('term'), eps
