"""Figure 4 layout. Panels A-D are placeholders for now."""
import numpy as np
from matplotlib.figure import Figure

from panel_E_F4 import panel_E_F4
from style_F4 import a4_figure, axes, panel


def plot_F4(freq: dict[str, dict[str, np.ndarray]]) -> Figure:
    f = a4_figure()
    for letter, pos in zip('ABCD', [(0, 0.75, 0.5, 0.25), (0.5, 0.75, 0.5, 0.25),
                                    (0, 0.5, 0.5, 0.25), (0.5, 0.5, 0.5, 0.25)]):
        panel(f, pos, letter)
    pE = panel(f, (0, 0.27, 1, 0.22), 'E')
    panel_E_F4(freq, axes(f, pE, (0.08, 0.15, 0.9, 0.78)))
    return f
