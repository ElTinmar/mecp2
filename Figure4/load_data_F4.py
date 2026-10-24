"""Load J-turn frequencies per fish from the bout tables."""
from enum import IntEnum
from pathlib import Path

import numpy as np
import pandas as pd

bouts_category_name_short = ["AS", "S1", "S2", "SCS", "LCS", "BS", "JT", "HAT", "RT", "SAT", "O", "LLC", "SLC"]

class Stim(IntEnum):
    PREY_CAPTURE = 6

GROUPS = {'mecp2-mutant': 'mecp2/bouts.csv', 'wild type': 'nacre/bouts.csv'}
SIDES = ['ipsi', 'contra']
JTURN = bouts_category_name_short.index('JT')
PROB_THRESHOLD = 0.5
TRIAL_DURATION_S = 7
PREY_ARC_DEG = 20
UV_INTENSITIES = np.array([0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1,
                           0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.75, 1.0])


def load_data_F4(datadir: Path) -> dict[str, dict[str, np.ndarray]]:
    """freq[group][side]: J-turn frequency (Hz), fish x UV intensity."""
    return {g: load_jt_frequency(datadir / path) for g, path in GROUPS.items()}


def load_jt_frequency(bout_file: Path) -> dict[str, np.ndarray]:
    df = pd.read_csv(bout_file)
    qc_file = bout_file.with_name('qc.csv')
    blacklisted = set(pd.read_csv(qc_file).file) if qc_file.exists() else set()

    prey = df[df.stim == Stim.PREY_CAPTURE]
    fish = [f for f in prey.file.unique() if f not in blacklisted]
    colors = [f'[0.0, 0.0, {i}, 1.0]' for i in UV_INTENSITIES]
    missing = set(colors) - set(prey.foreground_color)
    if missing:
        raise ValueError(f'{bout_file}: no prey-capture trials for {sorted(missing)}')

    jt = prey[
        (prey.category == JTURN) &
        (prey.proba > PROB_THRESHOLD) &
        (prey.trial_time <= TRIAL_DURATION_S) &
        (prey.prey_arc_start_deg.abs() == PREY_ARC_DEG) &
        prey.file.isin(fish)
    ]
    prey_side = np.sign(jt.prey_arc_start_deg)
    duration = 2 * TRIAL_DURATION_S  # two prey sides per laterality
    freq = {}
    for side, mask in zip(SIDES, [jt.sign == prey_side, jt.sign == -prey_side]):
        counts = pd.crosstab(jt.file[mask], jt.foreground_color[mask])
        freq[side] = counts.reindex(index=fish, columns=colors, fill_value=0).to_numpy() / duration
    return freq
