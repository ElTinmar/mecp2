"""Supplementary Figure 2 (Figure 4 data)."""
import numpy as np
from matplotlib.figure import Figure

from style_F4 import a4_figure, axes, panel
from SupPanel_A_F4 import SupPanel_A_F4
from SupPanel_B_F4 import SupPanel_B_F4
from SupPanel_C_F4 import SupPanel_C_F4


def plot_SupF4(freq: dict[str, dict[str, np.ndarray]], boot: dict[str, np.ndarray]) -> Figure:
    f = a4_figure()
    pA = panel(f, (0, 0.72, 1, 0.26), 'A')
    SupPanel_A_F4(axes(f, pA, (0, 0, 1, 0.95)))
    pB = panel(f, (0, 0.34, 0.5, 0.36), 'B')
    ax_B = np.array([[axes(f, pB, (0.18 + 0.37 * c, 0.52 - 0.44 * r, 0.31, 0.42)) for c in range(2)]
                     for r in range(2)])
    SupPanel_B_F4(freq, ax_B, axes(f, pB, (0.92, 0.3, 0.03, 0.4)))
    pC = panel(f, (0.58, 0.38, 0.42, 0.32), 'C')
    SupPanel_C_F4(freq, boot, axes(f, pC, (0.2, 0.15, 0.78, 0.68)))
    return f
