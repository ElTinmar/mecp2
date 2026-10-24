"""Supplementary Figure 2C: bootstrapped sigmoid (ipsi) and linear (contra) fits, EC50 comparison."""
import numpy as np
from matplotlib import legend_handler
from matplotlib.axes import Axes

from bootstrap_fits_F4 import INTENSITY_STOP
from load_data_F4 import UV_INTENSITIES
from style_F4 import COLORS, FISH_MEAN, label


def SupPanel_C_F4(freq: dict[str, dict[str, np.ndarray]], boot: dict[str, np.ndarray], ax: Axes) -> None:
    diff = boot['mecp2-mutant ipsi param'] - boot['wild type ipsi param']
    print('Bootstrapped sigmoid (ipsi) / linear (contra) fits:')
    for g in freq:
        print(f'{g} EC50: {median_ci(boot[f"{g} ipsi param"])}')
    print(f'EC50 difference (mutant - wild type): {median_ci(diff)}')
    print(f'Empirical p-value: {2 * min(np.mean(diff > 0), np.mean(diff < 0)):.4f}')
    for g in freq:
        print(f'{g} contra slope: {median_ci(boot[f"{g} contra param"])}')
    for metric in ['r2', 'rmse']:
        for side in ['ipsi', 'contra']:
            for g in freq:
                print(f'{g} {side} {metric}: {median_ci(boot[f"{g} {side} {metric}"])}')
    print()

    x, x_smooth = boot['x'], boot['x_smooth']
    log_all = np.log10(UV_INTENSITIES)
    handles, labels = [], []
    for side, marker, lw, ls in [('ipsi', 'o', 2, '-'), ('contra', '^', 1.2, '--')]:
        for g, sides in freq.items():
            low, med, high = np.percentile(boot[f'{g} {side} curve'], [2.5, 50, 97.5], axis=0)
            line, = ax.plot(x_smooth, med, color=COLORS[g], lw=lw, ls=ls, zorder=1)
            ax.fill_between(x_smooth, low, high, color=COLORS[g], alpha=0.15, edgecolor=None, zorder=0)
            mean = sides[side].mean(axis=0)
            ax.scatter(log_all[INTENSITY_STOP:], mean[INTENSITY_STOP:], marker=marker, s=15,
                       color=COLORS[g], edgecolor='k', zorder=5, alpha=0.4)
            scat = ax.scatter(x, mean[:INTENSITY_STOP], marker=marker, s=20, color=COLORS[g], edgecolor='k', zorder=5)
            handles.append((line, scat))
            labels.append(f'{label(g)} ({side})')
            if side == 'ipsi':
                lo, mid, hi = np.log10(np.percentile(boot[f'{g} ipsi param'], [2.5, 50, 97.5]))
                ax.errorbar(mid, (med.min() + med.max()) / 2, xerr=[[mid - lo], [hi - mid]], fmt='none',
                            color=COLORS[g], capsize=3, elinewidth=1, capthick=1, zorder=5)
    ticks = [0.01, 0.03, 0.1, 0.3, 1.0]
    ax.set_xticks(np.log10(ticks), [str(t) for t in ticks])
    ax.set(xlim=(np.log10(0.008), np.log10(1.2)), xlabel='UV intensity', ylabel=f'<J-turn frequency (Hz)>{FISH_MEAN}')
    ax.legend(handles, labels, loc='lower center', bbox_to_anchor=(0.5, 1.0), ncol=2, frameon=False,
              handler_map={tuple: legend_handler.HandlerTuple(ndivide=None)})


def median_ci(values: np.ndarray) -> str:
    return f'{np.median(values):.4f} [95% CI: {np.percentile(values, 2.5):.4f}, {np.percentile(values, 97.5):.4f}]'
