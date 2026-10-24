"""Figure 4E: ipsi/contralateral J-turn frequency vs UV intensity, mixed-design ANOVA and post-hoc tests."""
import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from scipy.stats import mannwhitneyu, sem
from statsmodels.stats.multitest import multipletests

from load_data_F4 import SIDES, UV_INTENSITIES
from mixed_anova import mixed_anova
from style_F4 import COLORS, FISH_MEAN, label


def panel_E_F4(freq: dict[str, dict[str, np.ndarray]], ax: Axes) -> None:
    mecp2, wt = freq['mecp2-mutant']['ipsi'], freq['wild type']['ipsi']

    table, eps = mixed_anova([mecp2, wt])
    print(f'J-turn frequency, ipsilateral (n = {len(mecp2)} mutant, {len(wt)} wild type)')
    print(f'Mixed-design ANOVA, Greenhouse-Geisser epsilon: {eps:.4f}')
    print(table.to_string(float_format='{:.4g}'.format, na_rep=''))

    stats = posthoc(mecp2, wt)
    print('\nPost-hoc tests (Mann-Whitney U, Holm corrected):')
    print(stats.to_string(index=False, formatters={
        'Mecp2_Mean_Hz': '{:.3f}'.format, 'WT_Mean_Hz': '{:.3f}'.format,
        'Raw_p': '{:.4e}'.format, 'Holm_Corrected_p': '{:.4e}'.format}))
    print()

    for side, ls in zip(SIDES, ['-', '--']):
        for g, sides in freq.items():
            mu, err = sides[side].mean(axis=0), sem(sides[side], axis=0)
            ax.plot(UV_INTENSITIES, mu, color=COLORS[g], ls=ls, lw=2, label=f'{label(g)} ({side})')
            ax.fill_between(UV_INTENSITIES, mu - err, mu + err, color=COLORS[g], alpha=0.5, lw=0)
    edges = log_midpoints(UV_INTENSITIES)
    change = np.diff(np.r_[0, (stats.Holm_Corrected_p.to_numpy() < 0.05).astype(int), 0])
    for start, end in zip(np.flatnonzero(change == 1), np.flatnonzero(change == -1)):
        ax.axvspan(edges[start], edges[end], color=(0.3, 0.3, 0.3), alpha=0.3, lw=0, zorder=0)
    ax.set(xscale='log', xlabel='UV intensity', ylabel=f'<J-turn frequency (Hz)>{FISH_MEAN}', ylim=(0, None))
    ticks = [0.01, 0.03, 0.1, 0.3, 1.0]
    ax.set_xticks(ticks, [str(t) for t in ticks])
    ax.minorticks_off()
    ax.legend(frameon=False, loc='upper left')


def posthoc(mecp2: np.ndarray, wt: np.ndarray) -> pd.DataFrame:
    """Mann-Whitney U per intensity across fish, Holm corrected."""
    tests = [mannwhitneyu(m, w, alternative='two-sided') for m, w in zip(mecp2.T, wt.T)]
    p = np.array([t.pvalue for t in tests])
    return pd.DataFrame({
        'UV_Intensity': UV_INTENSITIES,
        'Mecp2_Mean_Hz': mecp2.mean(axis=0),
        'WT_Mean_Hz': wt.mean(axis=0),
        'U_Statistic': [t.statistic for t in tests],
        'Raw_p': p,
        'Holm_Corrected_p': multipletests(p, method='holm')[1],
    })


def log_midpoints(x: np.ndarray) -> np.ndarray:
    """Edges halfway (in log space) between consecutive values."""
    log_x = np.log10(x)
    half = np.diff(log_x) / 2
    return 10 ** np.concatenate([[log_x[0] - half[0]], log_x[:-1] + half, [log_x[-1] + half[-1]]])
