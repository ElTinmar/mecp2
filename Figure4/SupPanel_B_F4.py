"""Supplementary Figure 2B: J-turn frequency heatmaps (fish x UV intensity)."""
import numpy as np
import seaborn as sns
from matplotlib.axes import Axes
import matplotlib.pyplot as plt

from load_data_F4 import SIDES, UV_INTENSITIES
from style_F4 import COLORS, label


def SupPanel_B_F4(freq: dict[str, dict[str, np.ndarray]], ax: np.ndarray, cax: Axes) -> None:
    """ax: 2 x 2 axes (rows: groups, columns: ipsi / contra), cax: colorbar axes."""
    vmax = np.percentile(np.concatenate([a.ravel() for sides in freq.values() for a in sides.values()]), 99.5)
    for row, g in enumerate(['wild type', 'mecp2-mutant']):
        sides = freq[g]
        for col, side in enumerate(SIDES):
            a = ax[row, col]
            data = sides[side]
            sns.heatmap(data, ax=a, cmap='gist_yarg', vmin=0, vmax=vmax, cbar=False,
                        xticklabels=[str(i) if i in (0.01, 0.1, 1.0) else '' for i in UV_INTENSITIES],
                        yticklabels=False)
            for spine in a.spines.values():
                spine.set(visible=True, color=COLORS[g], linewidth=1.2)
            if row == 0:
                a.set_title(f'{side}lateral', fontweight='bold')
                a.set_xticklabels([])
            else:
                a.set_xlabel('UV intensity')
                a.tick_params(axis='x', labelrotation=0)
        ax[row, 0].set_ylabel(f'{label(g, bold=True)}\n(n = {len(sides["ipsi"])} larvae)', color=COLORS[g], fontweight='bold')
    sm = plt.cm.ScalarMappable(cmap='gist_yarg', norm=plt.Normalize(0, vmax))
    cb = plt.colorbar(sm, cax=cax)
    cb.set_label('J-turn frequency (Hz)', size=8)
    cb.outline.set(visible=True, edgecolor='black', linewidth=1.0)
