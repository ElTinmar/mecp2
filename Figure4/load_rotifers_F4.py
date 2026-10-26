"""Figure 4A-D: rotifer counts per chamber, denoised (non-increasing isotonic fit) and normalized."""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import isotonic_regression

FPS = 30
START = 5        # s trimmed at the start: closing the room door, binarization frame buffer
ADAPT = 5 * 60   # s of adaptation before the analysed window
END = 15         # s trimmed at the end: edge effect of the isotonic fit (leaves 10 min)
CONDITIONS = ['WT_2lux', 'WT_30lux', 'MeCP2_2lux', 'MeCP2_30lux']


def load_rotifers_F4(datadir: Path) -> dict[str, np.ndarray]:
    """Rotifers remaining (fraction of the count at the end of adaptation), time x chamber, per condition."""
    counts = pd.read_csv(datadir / 'prey_counts.csv').iloc[START * FPS:]
    unknown = [c for c in counts if not any(k in c for k in CONDITIONS)]
    if unknown:
        raise ValueError(f'columns without a condition ({", ".join(CONDITIONS)}): {unknown}')

    denoised = np.column_stack([isotonic_regression(counts[c].to_numpy(float), increasing=False).x
                                for c in counts])[ADAPT * FPS:-END * FPS]
    print(f'Rotifers per chamber after adaptation: {denoised[0].mean():.1f} ± {denoised[0].std(ddof=1):.1f} (mean ± SD)')

    normalized = denoised / denoised[0]
    return {k: normalized[:, [k in c for c in counts.columns]] for k in CONDITIONS}
