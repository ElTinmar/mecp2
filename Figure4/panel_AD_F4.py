"""Figure 4A-D: rotifers remaining over time, by genotype (A: 30 lux, B: 2 lux) and illuminance (C, D)."""
import numpy as np
from matplotlib.axes import Axes

from load_rotifers_F4 import FPS
from rotifer_stats_F4 import ALPHA, COMPARISONS, Cluster
from style_F4 import COLORS, label

STYLE = {  # color, line style, legend label
    'WT_2lux': (COLORS['wild type'], '-', '2 lux'),
    'WT_30lux': ('#2EC4F1', '--', '30 lux'),
    'MeCP2_2lux': (COLORS['mecp2-mutant'], '-', '2 lux'),
    'MeCP2_30lux': ('#EA8D3B', '--', '30 lux'),
}
PANELS = [  # comparison, title, legend by genotype (else by illuminance)
    ('30 lux: wild type vs mecp2-mutant', '30 lux', True),
    ('2 lux: wild type vs mecp2-mutant', '2 lux', True),
    ('mecp2-mutant: 2 vs 30 lux', label('mecp2-mutant', bold=True), False),
    ('wild type: 2 vs 30 lux', 'wild type', False),
]
CHAMBER_MEAN = r'$_{\mathrm{chamber}}$'


def panel_AD_F4(data: dict[str, np.ndarray], ci: dict[str, np.ndarray],
                found: dict[str, list[Cluster]], axs: list[Axes]) -> None:
    eaten = {k: 100 * (1 - d[-1]) for k, d in data.items()}
    print('Rotifers eaten after 10 min (%, mean ± SD across chambers, SD with ddof=0 as in the original analysis)')
    for k, e in eaten.items():
        print(f'    {k}: {e.mean():.2f} ± {e.std():.2f} (n = {e.size} chambers)')
    print('\nCluster-based permutation tests (Welch t, two-sided, p = (r+1)/(N+1)):')
    for k, cs in found.items():
        print(f'    {k}: ' + ('; '.join(f'{minutes(s):.2f}-{minutes(e):.2f} min, p={p:.4f}' for s, e, p in cs) or 'no cluster'))
    print()

    minute = np.arange(len(next(iter(data.values())))) / FPS / 60
    lo = min(c[0].min() for c in ci.values())
    for ax, (comparison, title, by_genotype) in zip(axs, PANELS):
        for k in COMPARISONS[comparison]:
            color, ls, lux = STYLE[k]
            name = label('wild type' if k.startswith('WT') else 'mecp2-mutant') if by_genotype else lux
            ax.plot(minute, data[k].mean(axis=1), color=color, ls=ls, lw=2, label=name)
            ax.fill_between(minute, *ci[k], color=color, alpha=0.3, lw=0)
        for s, e, p in found[comparison]:
            if p < ALPHA:
                ax.axvspan(minute[s], minute[e], color=(0.3, 0.3, 0.3), alpha=0.3, lw=0, zorder=0)
        ax.text(0.98, 0.95, title, transform=ax.transAxes, ha='right', va='top', fontsize=10, fontweight='bold')
        ax.legend(frameon=False, loc='lower left')
        ax.set(xlim=(0, 10), ylim=(lo - 0.02, 1.02), xticks=range(11))
    for ax in axs[::2]:
        ax.set_ylabel(f'<Rotifers (norm.)>{CHAMBER_MEAN}')
    for ax in axs[2:]:
        ax.set_xlabel('Time (min)')


def minutes(index: int) -> float:
    return index / FPS / 60
