"""Figure 4 and Supplementary Figure 2B-C: rotifer consumption (A-D) and J-turns vs UV intensity (E),
mecp2-mutant vs wild type.

usage: python Figure4/F4.py <datadir>   (datadir contains Behavior_prey/ and MECP2_VR_Behavior/)
Writes fig4.npz, fig4.pdf/svg, figS2_F4.pdf/svg and table_F4.txt (statistics only;
progress messages go to stderr).
"""
import sys
from contextlib import redirect_stdout
from pathlib import Path
from typing import TextIO

import numpy as np

from bootstrap_fits_F4 import bootstrap_fits_F4
from load_data_F4 import load_data_F4
from load_rotifers_F4 import load_rotifers_F4
from plot_F4 import plot_F4
from rotifer_stats_F4 import rotifer_stats_F4
from plot_SupF4 import plot_SupF4
from style_F4 import save_figure


class Tee:
    """Write to several streams at once."""
    def __init__(self, *streams: TextIO) -> None:
        self.streams = streams

    def write(self, text: str) -> None:
        for s in self.streams:
            s.write(text)

    def flush(self) -> None:
        for s in self.streams:
            s.flush()


def F4(datadir: Path) -> None:
    freq = load_data_F4(datadir / 'MECP2_VR_Behavior')
    boot = bootstrap_fits_F4(freq)

    with open('table_F4.txt', 'w') as table, redirect_stdout(Tee(sys.stdout, table)):
        rotifers = load_rotifers_F4(datadir / 'Behavior_prey' / 'processed_videos')
        ci, found = rotifer_stats_F4(rotifers)
        save_figure(plot_F4(rotifers, ci, found, freq), 'fig4')
        save_figure(plot_SupF4(freq, boot), 'figS2_F4')

    np.savez('fig4.npz',
             **{f'rotifers {k}': d for k, d in rotifers.items()},
             **{f'{g} {side}': a for g, sides in freq.items() for side, a in sides.items()},
             **{f'boot {key}': val for key, val in boot.items()})


if __name__ == '__main__':
    F4(Path(sys.argv[1]))
