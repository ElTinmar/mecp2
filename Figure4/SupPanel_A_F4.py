"""Supplementary Figure 2A: 3D render of the virtual-reality setup."""
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.axes import Axes

RENDER = Path(__file__).parent / '3d_render' / 'setup_vr_pc.png'


def SupPanel_A_F4(ax: Axes) -> None:
    ax.imshow(plt.imread(RENDER), interpolation='none')
    ax.axis('off')
