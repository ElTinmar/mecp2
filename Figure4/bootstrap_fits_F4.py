"""Bootstrap (over fish) sigmoid fits of ipsilateral and linear fits of contralateral J-turn curves."""
import sys

import numpy as np
from scipy.optimize import curve_fit

from load_data_F4 import SIDES, UV_INTENSITIES

INTENSITY_STOP = 15  # intensities used for the fits
N_BOOTSTRAP = 10_000
SEED = 0


def sigmoid(x: np.ndarray, bottom: float, top: float, log_ec50: float, hill_slope: float) -> np.ndarray:
    return bottom + (top - bottom) / (1 + 10 ** ((log_ec50 - x) * hill_slope))


def linear(x: np.ndarray, slope: float, intercept: float) -> np.ndarray:
    return slope * x + intercept


def bootstrap_fits_F4(freq: dict[str, dict[str, np.ndarray]]) -> dict[str, np.ndarray]:
    """Per group and side: fit curve on x_smooth, r2, rmse and param (EC50 for ipsi, slope for contra)."""
    x = np.log10(UV_INTENSITIES[:INTENSITY_STOP])
    x_smooth = np.linspace(x.min(), x.max(), 200)
    rng = np.random.default_rng(SEED)
    m_ipsi = freq['mecp2-mutant']['ipsi']
    p0 = {'ipsi': [m_ipsi.min(), m_ipsi.max(), -1.0, 1.0],
          'contra': [0.0, freq['mecp2-mutant']['contra'].mean()]}
    bounds = {'ipsi': ([0, 0, x.min() - 1, 0.01], [10, 10, x.max() + 1, 5]),
              'contra': ([-2, 0], [2, 10])}
    model = {'ipsi': sigmoid, 'contra': linear}

    out = {}
    n_failed = 0
    for _ in range(N_BOOTSTRAP):
        sample = {}
        try:
            for g, sides in freq.items():
                idx = rng.integers(len(sides['ipsi']), size=len(sides['ipsi']))
                for side in SIDES:
                    y = sides[side][idx, :INTENSITY_STOP].mean(axis=0)
                    popt, _ = curve_fit(model[side], x, y, p0=p0[side], bounds=bounds[side], maxfev=5000)
                    fit = model[side](x, *popt)
                    ss_tot = np.sum((y - y.mean()) ** 2)
                    sample[f'{g} {side} r2'] = 1 - np.sum((y - fit) ** 2) / ss_tot if ss_tot else 0.0
                    sample[f'{g} {side} rmse'] = np.sqrt(np.mean((y - fit) ** 2))
                    sample[f'{g} {side} curve'] = model[side](x_smooth, *popt)
                    sample[f'{g} {side} param'] = 10 ** popt[2] if side == 'ipsi' else popt[0]
        except RuntimeError:
            n_failed += 1
            continue
        for key, val in sample.items():
            out.setdefault(key, []).append(val)
    print(f'{n_failed} of {N_BOOTSTRAP} bootstrap samples failed to fit', file=sys.stderr)
    return {'x': x, 'x_smooth': x_smooth} | {key: np.array(val) for key, val in out.items()}
