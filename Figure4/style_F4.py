"""Shared colors, style and A4 layout (as myfigure/mypanel/myaxes/myaxes_properties in MATLAB)."""
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure

COLORS = {'mecp2-mutant': '#D95319', 'wild type': '#0072BD'}
CM = 1 / 2.54
PAGE = (21, 29.7)
MARGIN = 1.5

plt.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['Arial'],
    'font.size': 8, 'axes.labelsize': 10, 'axes.titlesize': 10, 'legend.fontsize': 8,
    'xtick.direction': 'in', 'ytick.direction': 'in', 'xtick.top': True, 'ytick.right': True,
    'axes.linewidth': 1, 'svg.fonttype': 'none', 'pdf.fonttype': 42,
    'mathtext.fontset': 'custom', 'mathtext.rm': 'Arial', 'mathtext.it': 'Arial:italic',
    'mathtext.bfit': 'Arial:italic:bold',
})


FISH_MEAN = r'$_{\mathrm{fish}}$'


def label(group: str, bold: bool = False) -> str:
    """Group name for display, gene name in italic."""
    return group.replace('mecp2', r'$\mathbfit{mecp2}$' if bold else r'$\mathit{mecp2}$')


def a4_figure() -> Figure:
    return plt.figure(figsize=(PAGE[0] * CM, PAGE[1] * CM), facecolor='w')


def panel(fig: Figure, pos: tuple[float, float, float, float], letter: str) -> tuple[float, float, float, float]:
    """Panel at pos (x, y, w, h, normalized within the page margins), labelled with letter.

    Returns the panel rectangle in figure coordinates.
    """
    mx, my = MARGIN / PAGE[0], MARGIN / PAGE[1]
    rect = (mx + pos[0] * (1 - 2 * mx), my + pos[1] * (1 - 2 * my), pos[2] * (1 - 2 * mx), pos[3] * (1 - 2 * my))
    fig.text(rect[0], rect[1] + rect[3], letter, fontsize=12, fontweight='bold', va='top')
    return rect


def axes(fig: Figure, panel_rect: tuple[float, float, float, float], pos: tuple[float, float, float, float]) -> Axes:
    """Axes at pos (normalized within panel_rect)."""
    x, y, w, h = panel_rect
    return fig.add_axes((x + pos[0] * w, y + pos[1] * h, pos[2] * w, pos[3] * h))


def save_figure(fig: Figure, basename: str) -> None:
    fig.savefig(f'{basename}.pdf')
    fig.savefig(f'{basename}.svg')
