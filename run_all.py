"""Generate all figures: MATLAB (Figures 1-3, S1) and Python (Figure 4), in parallel.

usage: python run_all.py <datadir> [--figures F1 F3] [--jobs 2] [--out results] [--matlab /path/to/matlab]

Each figure runs in its own process; its output goes to <out>/logs/<figure>.log, and
its figures, tables and movies to <out>. Run from the conda environment of environment.yml.
"""
import argparse
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

REPO = Path(__file__).resolve().parent
FIGURES = ['F1', 'F2', 'F3', 'F4', 'FS1']


def command(figure: str, datadir: Path, matlab: str) -> list[str]:
    if figure == 'F4':
        return [sys.executable, str(REPO / 'Figure4' / 'F4.py'), str(datadir)]
    quote = lambda p: "'" + str(p).replace("'", "''") + "'"
    return [matlab, '-batch', f"addpath(genpath({quote(REPO)})); {figure}({quote(datadir)})"]


def run(figure: str, datadir: Path, out: Path, matlab: str) -> tuple[str, int, float]:
    start = time.time()
    with open(out / 'logs' / f'{figure}.log', 'w') as log:
        code = subprocess.run(command(figure, datadir, matlab), cwd=out, stdout=log, stderr=subprocess.STDOUT).returncode
    print(f'{figure}: {"done" if code == 0 else f"FAILED (exit {code})"} in {(time.time() - start) / 60:.1f} min', flush=True)
    return figure, code, time.time() - start


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('datadir', type=Path, help='folder containing Spontaneous/, Stimulation_*/, MECP2_VR_Behavior/ and Immunostaining/')
    parser.add_argument('--figures', nargs='+', choices=FIGURES, default=FIGURES)
    parser.add_argument('--jobs', type=int, default=len(FIGURES),
                        help='figures run at once; lower it if memory or MATLAB licenses are limited')
    parser.add_argument('--out', type=Path, default=Path('results'))
    parser.add_argument('--matlab', default='matlab', help='MATLAB executable (R2019a+ for -batch)')
    args = parser.parse_args()

    datadir = args.datadir.resolve()
    if not datadir.is_dir():
        sys.exit(f'data folder not found: {datadir}')
    out = args.out.resolve()
    (out / 'logs').mkdir(parents=True, exist_ok=True)

    print(f'Running {", ".join(args.figures)} ({args.jobs} at a time); logs in {out / "logs"}', flush=True)
    with ThreadPoolExecutor(args.jobs) as pool:
        results = list(pool.map(lambda f: run(f, datadir, out, args.matlab), args.figures))
    failed = [f for f, code, _ in results if code]
    sys.exit(f'Failed: {", ".join(failed)} (see logs)' if failed else 0)


if __name__ == '__main__':
    main()
