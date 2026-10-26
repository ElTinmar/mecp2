"""Figure 4 layout."""
import numpy as np
from matplotlib.figure import Figure

from panel_AD_F4 import panel_AD_F4
from panel_E_F4 import panel_E_F4
from rotifer_stats_F4 import Cluster
from style_F4 import a4_figure, axes, panel


def plot_F4(rotifers: dict[str, np.ndarray], ci: dict[str, np.ndarray], found: dict[str, list[Cluster]],
            freq: dict[str, dict[str, np.ndarray]]) -> Figure:
    f = a4_figure()
    axs = [axes(f, panel(f, pos, letter), (0.2, 0.18, 0.75, 0.72))
           for letter, pos in zip('ABCD', [(0, 0.75, 0.5, 0.25), (0.5, 0.75, 0.5, 0.25),
                                           (0, 0.5, 0.5, 0.25), (0.5, 0.5, 0.5, 0.25)])]
    panel_AD_F4(rotifers, ci, found, axs)
    pE = panel(f, (0, 0.27, 1, 0.22), 'E')
    panel_E_F4(freq, axes(f, pE, (0.08, 0.15, 0.9, 0.78)))
    return f
